#include "LazyIndex.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Util.h"

#include <carquet/include/carquet/carquet.h>

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {
const int64_t ROWS_PER_PAGE = 1024;
const int64_t ROWS_PER_ROW_GROUP = 16 * ROWS_PER_PAGE;
const char* const META_ENTRIES = "mmseqs.index.entries";
const char* const META_DATA_SIZE = "mmseqs.index.dataSize";
const char* const META_LAST_KEY = "mmseqs.index.lastKey";
const char* const META_MAX_SEQ_LEN = "mmseqs.index.maxSeqLen";
const char* const META_SORTED_BY_OFFSET = "mmseqs.index.sortedByOffset";
const char* const META_TEXT_INDEX_BYTES = "mmseqs.index.textIndexBytes";

// binary file: header, then one record per entry sorted by key
const char BINARY_MAGIC[8] = {'M', 'M', 'S', 'I', 'D', 'X', 'B', '1'};
struct BinaryRecord {
    size_t offset;
    DBKeyType id;
    unsigned int length;
};
struct BinaryHeader {
    char magic[8];
    uint64_t recordSize;
    uint64_t entries;
    uint64_t dataSize;
    uint64_t lastKey;
    uint64_t maxSeqLen;
    uint64_t sortedByOffset;
    uint64_t textIndexBytes;
};

bool lazyIndexDisabled = false;
size_t nextSerial = 1;

// size of a file, SIZE_MAX if it does not exist
size_t fileBytes(const std::string& file) {
    struct stat sb;
    if (stat(file.c_str(), &sb) != 0) {
        return SIZE_MAX;
    }
    return static_cast<size_t>(sb.st_size);
}

void carquetFail(const std::string& message, const carquet_error_t* err) {
    Debug(Debug::ERROR) << message << (err != NULL && err->message[0] != '\0' ? std::string(": ") + err->message : "") << "\n";
    EXIT(EXIT_FAILURE);
}

// ---------------------------------------------------------------------------------------------------------------
// custom decoder: thrift compact page header and DELTA_BINARY_PACKED values, read from the mapped file

struct Cursor {
    const uint8_t* p;
    const uint8_t* end;
    bool ok;
};

inline uint64_t readVarint(Cursor& c) {
    uint64_t v = 0;
    for (unsigned int shift = 0; shift < 64; shift += 7) {
        if (c.p >= c.end) {
            c.ok = false;
            return 0;
        }
        const uint8_t b = *c.p++;
        v |= static_cast<uint64_t>(b & 0x7F) << shift;
        if ((b & 0x80) == 0) {
            return v;
        }
    }
    c.ok = false;
    return 0;
}

inline int64_t zigzag(uint64_t v) {
    return static_cast<int64_t>(v >> 1) ^ -static_cast<int64_t>(v & 1);
}

inline void skipBytes(Cursor& c, uint64_t n) {
    if (static_cast<uint64_t>(c.end - c.p) < n) {
        c.ok = false;
        c.p = c.end;
        return;
    }
    c.p += n;
}

void thriftSkipStruct(Cursor& c, int depth);

// thrift compact types: 1/2 bool, 3 byte, 4 i16, 5 i32, 6 i64, 7 double, 8 binary, 9 list, 10 set, 11 map, 12 struct
void thriftSkip(Cursor& c, int type, bool inContainer, int depth) {
    if (depth > 16) {
        c.ok = false;
        return;
    }
    switch (type) {
        case 1:
        case 2:
            // a bool field has its value in the field header, a bool in a list is one byte
            if (inContainer) skipBytes(c, 1);
            break;
        case 3:
            skipBytes(c, 1);
            break;
        case 4:
        case 5:
        case 6:
            readVarint(c);
            break;
        case 7:
            skipBytes(c, 8);
            break;
        case 8:
            skipBytes(c, readVarint(c));
            break;
        case 9:
        case 10: {
            if (c.p >= c.end) { c.ok = false; return; }
            const uint8_t b = *c.p++;
            uint64_t n = b >> 4;
            if (n == 15) n = readVarint(c);
            for (uint64_t i = 0; i < n && c.ok; i++) thriftSkip(c, b & 0x0F, true, depth + 1);
            break;
        }
        case 11: {
            const uint64_t n = readVarint(c);
            if (n > 0) {
                if (c.p >= c.end) { c.ok = false; return; }
                const uint8_t kv = *c.p++;
                for (uint64_t i = 0; i < n && c.ok; i++) {
                    thriftSkip(c, kv >> 4, true, depth + 1);
                    thriftSkip(c, kv & 0x0F, true, depth + 1);
                }
            }
            break;
        }
        case 12:
            thriftSkipStruct(c, depth + 1);
            break;
        default:
            c.ok = false;
    }
}

// next field of a struct: false at the stop byte
inline bool thriftField(Cursor& c, int& id, int& type) {
    if (c.p >= c.end) {
        c.ok = false;
        return false;
    }
    const uint8_t b = *c.p++;
    if (b == 0) {
        return false;
    }
    type = b & 0x0F;
    if ((b >> 4) == 0) {
        id = static_cast<int>(zigzag(readVarint(c)));
    } else {
        id += b >> 4;
    }
    return c.ok;
}

void thriftSkipStruct(Cursor& c, int depth) {
    int id = 0, type = 0;
    while (thriftField(c, id, type)) {
        thriftSkip(c, type, false, depth);
    }
}

struct PageHeader {
    int type;               // 0 DATA_PAGE, 3 DATA_PAGE_V2
    int64_t uncompressedSize;
    int64_t compressedSize;
    int64_t numValues;
    int encoding;           // 5 DELTA_BINARY_PACKED
    int64_t levelBytes;     // v2: repetition + definition level bytes before the values
    PageHeader() : type(-1), uncompressedSize(0), compressedSize(0), numValues(0), encoding(-1), levelBytes(0) {}
};

// PageHeader: 1 type, 2 uncompressed_page_size, 3 compressed_page_size, 5 data_page_header {1 num_values, 2 encoding},
// 8 data_page_header_v2 {1 num_values, 4 encoding, 5 definition_levels_byte_length, 6 repetition_levels_byte_length}
bool readPageHeader(Cursor& c, PageHeader& h) {
    int id = 0, type = 0;
    while (thriftField(c, id, type)) {
        if (id == 1 && type == 5) {
            h.type = static_cast<int>(zigzag(readVarint(c)));
        } else if (id == 2 && type == 5) {
            h.uncompressedSize = zigzag(readVarint(c));
        } else if (id == 3 && type == 5) {
            h.compressedSize = zigzag(readVarint(c));
        } else if ((id == 5 || id == 8) && type == 12) {
            const bool v2 = (id == 8);
            int sid = 0, stype = 0;
            while (thriftField(c, sid, stype)) {
                if (sid == 1 && stype == 5) {
                    h.numValues = zigzag(readVarint(c));
                } else if (sid == (v2 ? 4 : 2) && stype == 5) {
                    h.encoding = static_cast<int>(zigzag(readVarint(c)));
                } else if (v2 && (sid == 5 || sid == 6) && stype == 5) {
                    h.levelBytes += zigzag(readVarint(c));
                } else {
                    thriftSkip(c, stype, false, 1);
                }
            }
        } else {
            thriftSkip(c, type, false, 0);
        }
    }
    return c.ok;
}

// w bits at bit position bitPos (least significant bit first); fileEnd bounds the fast unaligned read
inline uint64_t readBits(const uint8_t* base, const uint8_t* fileEnd, uint64_t bitPos, unsigned int w) {
    const uint8_t* q = base + (bitPos >> 3);
    const unsigned int shift = static_cast<unsigned int>(bitPos & 7);
    uint64_t v;
#if defined(__BYTE_ORDER__) && __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
    if (q + 9 <= fileEnd) {
        memcpy(&v, q, sizeof(v));
        v >>= shift;
        if (shift + w > 64) {
            v |= static_cast<uint64_t>(q[8]) << (64 - shift);
        }
    } else
#endif
    {
        v = 0;
        const unsigned int bytes = (shift + w + 7) / 8;
        unsigned int have = 0;
        for (unsigned int i = 0; i < bytes && q + i < fileEnd; i++) {
            const uint64_t byte = q[i];
            if (i == 0) {
                v = byte >> shift;
                have = 8 - shift;
            } else {
                if (have < 64) v |= byte << have;
                have += 8;
            }
        }
    }
    return (w == 64) ? v : (v & ((static_cast<uint64_t>(1) << w) - 1));
}

// DELTA_BINARY_PACKED: block size, miniblocks per block, value count, first value; per block: min delta, one bit
// width per miniblock, the miniblocks (deltas - min delta, bit packed). n values are expected
template <typename O>
bool decodeDelta(Cursor& c, const uint8_t* fileEnd, O* out, size_t n) {
    const uint64_t blockSize = readVarint(c);
    const uint64_t miniblocks = readVarint(c);
    const uint64_t total = readVarint(c);
    uint64_t prev = static_cast<uint64_t>(zigzag(readVarint(c)));
    if (c.ok == false || miniblocks == 0 || blockSize == 0 || blockSize % miniblocks != 0 || total != n || n == 0) {
        return false;
    }
    const uint64_t valuesPerMiniblock = blockSize / miniblocks;
    out[0] = static_cast<O>(prev);
    size_t produced = 1;
    while (produced < n) {
        const uint64_t minDelta = static_cast<uint64_t>(zigzag(readVarint(c)));
        const uint8_t* widths = c.p;
        skipBytes(c, miniblocks);
        if (c.ok == false) {
            return false;
        }
        for (uint64_t m = 0; m < miniblocks && produced < n; m++) {
            const unsigned int w = widths[m];
            const size_t count = static_cast<size_t>(std::min<uint64_t>(valuesPerMiniblock, n - produced));
            if (w > 64 || static_cast<uint64_t>(c.end - c.p) < (count * w + 7) / 8) {
                return false;
            }
            if (w == 0) {
                for (size_t i = 0; i < count; i++) {
                    prev += minDelta;
                    out[produced + i] = static_cast<O>(prev);
                }
            } else {
                for (size_t i = 0; i < count; i++) {
                    prev += minDelta + readBits(c.p, fileEnd, i * w, w);
                    out[produced + i] = static_cast<O>(prev);
                }
            }
            produced += count;
            // a miniblock is stored in full, the last one of the page may be cut after its values
            c.p += std::min<uint64_t>(valuesPerMiniblock * w / 8, static_cast<uint64_t>(c.end - c.p));
        }
    }
    return true;
}

// one data page of a REQUIRED column (no levels) at p; n values
template <typename O>
bool decodePage(const uint8_t* p, const uint8_t* fileEnd, O* out, size_t n) {
    Cursor c = {p, fileEnd, true};
    PageHeader h;
    if (readPageHeader(c, h) == false || (h.type != 0 && h.type != 3) || h.encoding != 5
        || h.numValues != static_cast<int64_t>(n) || h.compressedSize != h.uncompressedSize || h.compressedSize < 0
        || static_cast<uint64_t>(fileEnd - c.p) < static_cast<uint64_t>(h.compressedSize)) {
        return false;
    }
    Cursor values = {c.p, c.p + h.compressedSize, true};
    skipBytes(values, static_cast<uint64_t>(h.levelBytes));
    return values.ok && decodeDelta(values, fileEnd, out, n);
}
}

// ---------------------------------------------------------------------------------------------------------------

// decoded page of the three columns; one per thread and LazyIndex
struct LazyIndex::Page {
    size_t serial;
    size_t pageIdx;
    size_t rows;
    std::vector<int64_t> key;
    std::vector<int64_t> offset;
    std::vector<int32_t> length;
    Page() : serial(0), pageIdx(SIZE_MAX), rows(0) {}
};

LazyIndex::LazyIndex() : backend(BACKEND_NONE), serial(0), reader(NULL), rowsPerPage(0), map(NULL), mapSize(0) {
    memset(&meta, 0, sizeof(meta));
#pragma omp critical(lazyIndexSerial)
    serial = nextSerial++;
}

LazyIndex::~LazyIndex() {
    if (reader != NULL) {
        carquet_reader_close(reader);
    }
    if (map != NULL) {
        munmap(map, mapSize);
    }
}

LazyIndex::Backend LazyIndex::backendFromEnv() {
    const char* env = getenv("MMSEQS_INDEX_BACKEND");
    if (env == NULL) {
        return BACKEND_NONE;
    }
    if (strcmp(env, "parquet") == 0) return BACKEND_CARQUET;
    if (strcmp(env, "parquet-custom") == 0) return BACKEND_CUSTOM;
    if (strcmp(env, "binary") == 0) return BACKEND_BINARY;
    return BACKEND_NONE;
}

const char* LazyIndex::backendName(Backend backend) {
    switch (backend) {
        case BACKEND_CARQUET: return "parquet (carquet)";
        case BACKEND_CUSTOM: return "parquet (custom decoder)";
        case BACKEND_BINARY: return "binary";
        default: return "text";
    }
}

std::string LazyIndex::fileName(const std::string& indexFile, Backend backend) {
    return indexFile + (backend == BACKEND_BINARY ? ".bin" : ".parquet");
}

bool LazyIndex::isDisabled() {
    return lazyIndexDisabled;
}

size_t LazyIndex::memoryBytes() const {
    return firstKey.size() * sizeof(int64_t) + firstRow.size() * sizeof(size_t) + rowGroup.size() * sizeof(int32_t)
           + rowGroupFirstRow.size() * sizeof(size_t) + pageOffset.size() * sizeof(uint64_t);
}

void LazyIndex::writeParquet(const std::string& file, const Entry* entries, const Meta& meta) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* schema = carquet_schema_create(&err);
    if (schema == NULL
        || carquet_schema_add_column(schema, "key", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK
        || carquet_schema_add_column(schema, "offset", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK
        || carquet_schema_add_column(schema, "length", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
        carquetFail("Cannot create the schema of " + file, &err);
    }
    carquet_writer_options_t opts;
    carquet_writer_options_init(&opts);
    // DELTA_BINARY_PACKED is the compression: a codec on top saves < 1% and slows down every page read
    opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
    // row groups and pages are cut by rows: small pages make a lookup cheap, small row groups keep carquet's cost
    // per column reader low (it grows with the pages of the row group)
    opts.row_group_size = static_cast<int64_t>(1) << 40;
    opts.page_size = 64 * 1024 * 1024;
    opts.max_rows_per_page = ROWS_PER_PAGE;
    opts.write_statistics = true;
    opts.write_page_index = true;
    opts.write_bloom_filters = false;
    opts.write_crc = false;
    opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    opts.created_by = "MMseqs parquet index";
    const std::string tmp = file + ".tmp";
    carquet_writer_t* writer = carquet_writer_create(tmp.c_str(), schema, &opts, &err);
    carquet_schema_free(schema);
    if (writer == NULL) {
        carquetFail("Cannot write " + tmp, &err);
    }
    bool ok = true;
    // page statistics only for the key: the first key of every page is the page directory of the reader
    ok = ok && carquet_writer_set_column_statistics(writer, 1, false) == CARQUET_OK;
    ok = ok && carquet_writer_set_column_statistics(writer, 2, false) == CARQUET_OK;
    for (int32_t c = 0; c < 3; c++) {
        ok = ok && carquet_writer_set_column_encoding(writer, c, CARQUET_ENCODING_DELTA_BINARY_PACKED) == CARQUET_OK;
    }
    ok = ok && carquet_writer_add_metadata(writer, META_ENTRIES, SSTR(meta.entries).c_str()) == CARQUET_OK;
    ok = ok && carquet_writer_add_metadata(writer, META_DATA_SIZE, SSTR(meta.dataSize).c_str()) == CARQUET_OK;
    ok = ok && carquet_writer_add_metadata(writer, META_LAST_KEY, SSTR(meta.lastKey).c_str()) == CARQUET_OK;
    ok = ok && carquet_writer_add_metadata(writer, META_MAX_SEQ_LEN, SSTR(meta.maxSeqLen).c_str()) == CARQUET_OK;
    ok = ok && carquet_writer_add_metadata(writer, META_SORTED_BY_OFFSET, meta.sortedByOffset ? "1" : "0") == CARQUET_OK;
    ok = ok && carquet_writer_add_metadata(writer, META_TEXT_INDEX_BYTES, SSTR(meta.textIndexBytes).c_str()) == CARQUET_OK;
    if (ok == false) {
        carquetFail("Cannot set the writer options of " + tmp, NULL);
    }
    std::vector<int64_t> key(ROWS_PER_ROW_GROUP), offset(ROWS_PER_ROW_GROUP);
    std::vector<int32_t> length(ROWS_PER_ROW_GROUP);
    for (size_t start = 0; start < meta.entries; start += ROWS_PER_ROW_GROUP) {
        const int64_t n = static_cast<int64_t>(std::min<size_t>(ROWS_PER_ROW_GROUP, meta.entries - start));
        for (int64_t i = 0; i < n; i++) {
            key[i] = static_cast<int64_t>(entries[start + i].id);
            offset[i] = static_cast<int64_t>(entries[start + i].offset);
            length[i] = static_cast<int32_t>(entries[start + i].length);
        }
        if ((start > 0 && carquet_writer_new_row_group(writer) != CARQUET_OK)
            || carquet_writer_write_batch(writer, 0, key.data(), n, NULL, NULL) != CARQUET_OK
            || carquet_writer_write_batch(writer, 1, offset.data(), n, NULL, NULL) != CARQUET_OK
            || carquet_writer_write_batch(writer, 2, length.data(), n, NULL, NULL) != CARQUET_OK) {
            carquetFail("Cannot write the entries of " + tmp, NULL);
        }
    }
    if (carquet_writer_close(writer) != CARQUET_OK) {
        carquetFail("Cannot close " + tmp, NULL);
    }
    FileUtil::move(tmp.c_str(), file.c_str());
}

void LazyIndex::writeBinary(const std::string& file, const Entry* entries, const Meta& meta) {
    const std::string tmp = file + ".tmp";
    FILE* fp = FileUtil::openFileOrDie(tmp.c_str(), "wb", false);
    BinaryHeader header;
    memset(&header, 0, sizeof(header));
    memcpy(header.magic, BINARY_MAGIC, sizeof(BINARY_MAGIC));
    header.recordSize = sizeof(BinaryRecord);
    header.entries = meta.entries;
    header.dataSize = meta.dataSize;
    header.lastKey = meta.lastKey;
    header.maxSeqLen = meta.maxSeqLen;
    header.sortedByOffset = meta.sortedByOffset;
    header.textIndexBytes = meta.textIndexBytes;
    bool ok = fwrite(&header, sizeof(header), 1, fp) == 1;
    const size_t BATCH = 1024 * 1024;
    std::vector<BinaryRecord> records(std::min(BATCH, meta.entries));
    for (size_t start = 0; start < meta.entries && ok; start += BATCH) {
        const size_t n = std::min(BATCH, meta.entries - start);
        for (size_t i = 0; i < n; i++) {
            records[i].offset = entries[start + i].offset;
            records[i].id = entries[start + i].id;
            records[i].length = entries[start + i].length;
        }
        ok = fwrite(records.data(), sizeof(BinaryRecord), n, fp) == n;
    }
    if (fclose(fp) != 0 || ok == false) {
        Debug(Debug::ERROR) << "Cannot write " << tmp << "\n";
        EXIT(EXIT_FAILURE);
    }
    FileUtil::move(tmp.c_str(), file.c_str());
}

void LazyIndex::writeFromTextIndex(const std::string& indexFile, const std::string& parquetFile,
                                   const std::string& binaryFile, int threads) {
    const bool parquet = parquetFile.empty() == false;
    const bool binary = binaryFile.empty() == false;
    if (parquet == false && binary == false) {
        return;
    }
    const std::string suffix = ".index";
    std::string dataFile = indexFile;
    if (dataFile.size() > suffix.size() && dataFile.compare(dataFile.size() - suffix.size(), suffix.size(), suffix) == 0) {
        dataFile.resize(dataFile.size() - suffix.size());
    }
    // the text index is the source: no lazy index while it is read
    lazyIndexDisabled = true;
    DBReader<DBKeyType> reader(dataFile.c_str(), indexFile.c_str(), threads, DBReader<DBKeyType>::USE_INDEX);
    reader.open(DBReader<DBKeyType>::NOSORT);
    lazyIndexDisabled = false;
    Meta meta;
    memset(&meta, 0, sizeof(meta));
    meta.entries = reader.getSize();
    if (meta.entries > 0) {
        // DBReader<DBKeyType>::Index and Entry have the same members
        const DBReader<DBKeyType>::Index* index = reader.getIndex();
        std::vector<Entry> sorted;
        const Entry* entries = reinterpret_cast<const Entry*>(index);
        meta.sortedByOffset = true;
        bool sortedByKey = true;
        for (size_t i = 0; i < meta.entries; i++) {
            meta.dataSize += index[i].length;
            meta.maxSeqLen = std::max(meta.maxSeqLen, index[i].length);
            meta.lastKey = std::max(meta.lastKey, index[i].id);
            meta.sortedByOffset = meta.sortedByOffset && (i == 0 || index[i].offset >= index[i - 1].offset);
            sortedByKey = sortedByKey && (i == 0 || index[i].id >= index[i - 1].id);
        }
        if (sortedByKey == false) {
            Debug(Debug::ERROR) << "Index " << indexFile << " is not sorted by key\n";
            EXIT(EXIT_FAILURE);
        }
        meta.textIndexBytes = fileBytes(indexFile);
        if (parquet) {
            writeParquet(parquetFile, entries, meta);
        }
        if (binary) {
            writeBinary(binaryFile, entries, meta);
        }
    }
    reader.close();
}

void LazyIndex::writeAfterDbWrite(const std::string& indexFile, int threads) {
    const char* env = getenv("MMSEQS_WRITE_LAZY_INDEX");
    if (env == NULL || env[0] == '\0') {
        return;
    }
    const std::string what(env);
    writeFromTextIndex(indexFile, what.find("parquet") != std::string::npos ? fileName(indexFile, BACKEND_CARQUET) : "",
                       what.find("binary") != std::string::npos ? fileName(indexFile, BACKEND_BINARY) : "", threads);
}

LazyIndex* LazyIndex::open(const std::string& indexFile, Backend backend) {
    if (backend == BACKEND_NONE) {
        return NULL;
    }
    const std::string file = fileName(indexFile, backend);
    if (fileBytes(file) == SIZE_MAX) {
        return NULL;
    }
    LazyIndex* idx = new LazyIndex();
    idx->backend = backend;
    const bool ok = (backend == BACKEND_BINARY) ? idx->openBinary(file) : idx->openParquet(file, backend == BACKEND_CUSTOM);
    // made from another text index (the DB was written again without the lazy index)
    const size_t textBytes = fileBytes(indexFile);
    const bool stale = ok && textBytes != SIZE_MAX && textBytes != idx->meta.textIndexBytes;
    if (ok == false || stale || idx->meta.entries == 0) {
        Debug(Debug::WARNING) << "Not using " << file << (stale ? ": it does not fit " + indexFile : ": cannot be read as index")
                              << ", reading the text index\n";
        delete idx;
        return NULL;
    }
    return idx;
}

bool LazyIndex::openBinary(const std::string& file) {
    const size_t bytes = fileBytes(file);
    BinaryHeader header;
    if (bytes < sizeof(header)) {
        return false;
    }
    const int fd = ::open(file.c_str(), O_RDONLY);
    if (fd < 0) {
        return false;
    }
    void* m = mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
    ::close(fd);
    if (m == MAP_FAILED) {
        return false;
    }
    map = static_cast<char*>(m);
    mapSize = bytes;
    memcpy(&header, map, sizeof(header));
    if (memcmp(header.magic, BINARY_MAGIC, sizeof(BINARY_MAGIC)) != 0 || header.recordSize != sizeof(BinaryRecord)
        || bytes != sizeof(header) + header.entries * sizeof(BinaryRecord)) {
        return false;
    }
    meta.entries = header.entries;
    meta.dataSize = header.dataSize;
    meta.lastKey = static_cast<DBKeyType>(header.lastKey);
    meta.maxSeqLen = static_cast<unsigned int>(header.maxSeqLen);
    meta.sortedByOffset = header.sortedByOffset != 0;
    meta.textIndexBytes = header.textIndexBytes;
    return true;
}

bool LazyIndex::openParquet(const std::string& file, bool custom) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t options;
    carquet_reader_options_init(&options);
    options.use_mmap = true;
    reader = carquet_reader_open(file.c_str(), &options, &err);
    if (reader == NULL) {
        return false;
    }
    const carquet_schema_t* schema = carquet_reader_schema(reader);
    const char* entries = carquet_reader_find_metadata(reader, META_ENTRIES);
    const char* dataSize = carquet_reader_find_metadata(reader, META_DATA_SIZE);
    const char* lastKey = carquet_reader_find_metadata(reader, META_LAST_KEY);
    const char* maxSeqLen = carquet_reader_find_metadata(reader, META_MAX_SEQ_LEN);
    const char* sortedByOffset = carquet_reader_find_metadata(reader, META_SORTED_BY_OFFSET);
    const char* textIndexBytes = carquet_reader_find_metadata(reader, META_TEXT_INDEX_BYTES);
    if (schema == NULL || carquet_schema_num_columns(schema) != 3
        || carquet_schema_column_type(schema, 0) != CARQUET_PHYSICAL_INT64
        || carquet_schema_column_type(schema, 1) != CARQUET_PHYSICAL_INT64
        || carquet_schema_column_type(schema, 2) != CARQUET_PHYSICAL_INT32
        || entries == NULL || dataSize == NULL || lastKey == NULL || maxSeqLen == NULL || sortedByOffset == NULL
        || textIndexBytes == NULL) {
        return false;
    }
    meta.entries = strtoull(entries, NULL, 10);
    meta.dataSize = strtoull(dataSize, NULL, 10);
    meta.lastKey = static_cast<DBKeyType>(strtoull(lastKey, NULL, 10));
    meta.maxSeqLen = static_cast<unsigned int>(strtoull(maxSeqLen, NULL, 10));
    meta.sortedByOffset = sortedByOffset[0] == '1';
    meta.textIndexBytes = strtoull(textIndexBytes, NULL, 10);
    if (static_cast<int64_t>(meta.entries) != carquet_reader_num_rows(reader)) {
        return false;
    }

    // page directory: first key and first row of every page of the key column; the custom decoder also needs the
    // file offset of the page in each column, and the pages of the three columns have to hold the same rows
    const int32_t rowGroups = carquet_reader_num_row_groups(reader);
    size_t groupFirstRow = 0;
    bool uniform = true;
    for (int32_t g = 0; g < rowGroups; g++) {
        carquet_row_group_metadata_t groupMeta;
        if (carquet_reader_row_group_metadata(reader, g, &groupMeta) != CARQUET_OK) {
            return false;
        }
        carquet_column_index_t* columnIndex = carquet_reader_get_column_index(reader, g, 0, &err);
        carquet_offset_index_t* offsetIndex[3] = {NULL, NULL, NULL};
        bool ok = columnIndex != NULL;
        for (int32_t c = 0; c < (custom ? 3 : 1) && ok; c++) {
            offsetIndex[c] = carquet_reader_get_offset_index(reader, g, c, &err);
            ok = offsetIndex[c] != NULL;
            if (custom && ok) {
                carquet_column_chunk_metadata_t chunk;
                ok = carquet_reader_column_chunk_metadata(reader, g, c, &chunk) == CARQUET_OK
                     && chunk.codec == CARQUET_COMPRESSION_UNCOMPRESSED && chunk.has_dictionary_page == false;
            }
        }
        const int32_t pages = ok ? carquet_column_index_num_pages(columnIndex) : 0;
        for (int32_t c = 0; c < (custom ? 3 : 1) && ok; c++) {
            ok = carquet_offset_index_num_pages(offsetIndex[c]) == pages;
        }
        for (int32_t p = 0; p < pages && ok; p++) {
            carquet_page_stats_t stats;
            carquet_page_location_t location;
            ok = carquet_column_index_get_page_stats(columnIndex, p, &stats) == CARQUET_OK
                 && stats.min_value != NULL && stats.min_value_size == static_cast<int32_t>(sizeof(int64_t))
                 && carquet_offset_index_get_page_location(offsetIndex[0], p, &location) == CARQUET_OK;
            if (ok == false) {
                break;
            }
            int64_t first;
            memcpy(&first, stats.min_value, sizeof(first));
            const size_t row = groupFirstRow + static_cast<size_t>(location.first_row_index);
            uniform = uniform && row == firstKey.size() * static_cast<size_t>(ROWS_PER_PAGE);
            firstKey.push_back(first);
            firstRow.push_back(row);
            if (custom) {
                pageOffset.push_back(static_cast<uint64_t>(location.offset));
                for (int32_t c = 1; c < 3 && ok; c++) {
                    carquet_page_location_t other;
                    ok = carquet_offset_index_get_page_location(offsetIndex[c], p, &other) == CARQUET_OK
                         && other.first_row_index == location.first_row_index;
                    pageOffset.push_back(static_cast<uint64_t>(other.offset));
                }
            } else {
                rowGroup.push_back(g);
            }
        }
        carquet_column_index_free(columnIndex);
        for (int32_t c = 0; c < 3; c++) {
            carquet_offset_index_free(offsetIndex[c]);
        }
        if (ok == false) {
            return false;
        }
        rowGroupFirstRow.push_back(groupFirstRow);
        groupFirstRow += static_cast<size_t>(groupMeta.num_rows);
    }
    if (groupFirstRow != meta.entries || firstKey.empty()) {
        return false;
    }
    if (uniform) {
        // the usual layout: the page of a row is row / rowsPerPage, no first row per page needed
        rowsPerPage = ROWS_PER_PAGE;
        std::vector<size_t>().swap(firstRow);
    }
    if (custom) {
        // the pages are decoded from an own mapping; carquet was only needed for the footer and the page index
        std::vector<size_t>().swap(rowGroupFirstRow);
        carquet_reader_close(reader);
        reader = NULL;
        const size_t bytes = fileBytes(file);
        const int fd = ::open(file.c_str(), O_RDONLY);
        if (fd < 0) {
            return false;
        }
        void* m = mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        ::close(fd);
        if (m == MAP_FAILED) {
            return false;
        }
        map = static_cast<char*>(m);
        mapSize = bytes;
        for (size_t i = 0; i < pageOffset.size(); i++) {
            if (pageOffset[i] >= mapSize) {
                return false;
            }
        }
    }
    return true;
}

size_t LazyIndex::pageOfRow(size_t row) const {
    if (rowsPerPage != 0) {
        return row / rowsPerPage;
    }
    return static_cast<size_t>(std::upper_bound(firstRow.begin(), firstRow.end(), row) - firstRow.begin()) - 1;
}

size_t LazyIndex::rowsOfPage(size_t pageIdx) const {
    const size_t first = (rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx];
    const size_t next = (pageIdx + 1 == firstKey.size()) ? meta.entries
                        : ((rowsPerPage != 0) ? (pageIdx + 1) * rowsPerPage : firstRow[pageIdx + 1]);
    return next - first;
}

void LazyIndex::loadPage(size_t pageIdx, Page& p) {
    const size_t rows = rowsOfPage(pageIdx);
    p.key.resize(rows);
    p.offset.resize(rows);
    p.length.resize(rows);
    bool ok = true;
    if (backend == BACKEND_CUSTOM) {
        const uint8_t* base = reinterpret_cast<const uint8_t*>(map);
        const uint8_t* end = base + mapSize;
        ok = decodePage(base + pageOffset[3 * pageIdx], end, p.key.data(), rows)
             && decodePage(base + pageOffset[3 * pageIdx + 1], end, p.offset.data(), rows)
             && decodePage(base + pageOffset[3 * pageIdx + 2], end, p.length.data(), rows);
    } else {
        const int32_t g = rowGroup[pageIdx];
        const size_t first = (rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx];
        const int64_t skip = static_cast<int64_t>(first - rowGroupFirstRow[g]);
        void* out[3] = {p.key.data(), p.offset.data(), p.length.data()};
        for (int32_t c = 0; c < 3 && ok; c++) {
            carquet_error_t err = CARQUET_ERROR_INIT;
            carquet_column_reader_t* column = carquet_reader_get_column(reader, g, c, &err);
            ok = column != NULL && (skip == 0 || carquet_column_skip(column, skip) == skip)
                 && carquet_column_read_batch(column, out[c], static_cast<int64_t>(rows), NULL, NULL) == static_cast<int64_t>(rows);
            carquet_column_reader_free(column);
        }
    }
    if (ok == false) {
        Debug(Debug::ERROR) << "Cannot read page " << pageIdx << " of the " << backendName(backend) << " index\n";
        EXIT(EXIT_FAILURE);
    }
    p.serial = serial;
    p.pageIdx = pageIdx;
    p.rows = rows;
}

// the decoded page of this thread; a thread keeps one page per LazyIndex it reads (query DB, target DB, ...)
LazyIndex::Page& LazyIndex::page(size_t pageIdx) {
    static thread_local std::vector<Page> pages;
    static thread_local size_t evict = 0;
    Page* p = NULL;
    for (size_t i = 0; i < pages.size(); i++) {
        if (pages[i].serial == serial) {
            p = &pages[i];
            break;
        }
    }
    if (p == NULL) {
        const size_t MAX_PAGES = 16;
        if (pages.size() < MAX_PAGES) {
            pages.push_back(Page());
            p = &pages.back();
        } else {
            p = &pages[evict++ % MAX_PAGES];
        }
        p->serial = serial;
        p->pageIdx = SIZE_MAX;
    }
    if (p->pageIdx != pageIdx) {
        loadPage(pageIdx, *p);
    }
    return *p;
}

size_t LazyIndex::getId(DBKeyType key) {
    if (backend == BACKEND_BINARY) {
        const BinaryRecord* records = reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader));
        size_t lo = 0, hi = meta.entries;
        while (lo < hi) {
            const size_t mid = lo + (hi - lo) / 2;
            if (records[mid].id < key) lo = mid + 1; else hi = mid;
        }
        return (lo < meta.entries && records[lo].id == key) ? lo : SIZE_MAX;
    }
    const int64_t k = static_cast<int64_t>(key);
    // first page that starts with a key >= k: the first entry with this key is in the page before or starts this page
    const size_t q = static_cast<size_t>(std::lower_bound(firstKey.begin(), firstKey.end(), k) - firstKey.begin());
    if (q > 0) {
        const Page& p = page(q - 1);
        const std::vector<int64_t>::const_iterator it = std::lower_bound(p.key.begin(), p.key.end(), k);
        if (it != p.key.end() && *it == k) {
            const size_t first = (rowsPerPage != 0) ? (q - 1) * rowsPerPage : firstRow[q - 1];
            return first + static_cast<size_t>(it - p.key.begin());
        }
    }
    if (q < firstKey.size() && firstKey[q] == k) {
        return (rowsPerPage != 0) ? q * rowsPerPage : firstRow[q];
    }
    return SIZE_MAX;
}

LazyIndex::Entry LazyIndex::get(size_t id) {
    Entry e;
    if (backend == BACKEND_BINARY) {
        const BinaryRecord* records = reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader));
        e.id = records[id].id;
        e.offset = records[id].offset;
        e.length = records[id].length;
        return e;
    }
    const size_t pageIdx = pageOfRow(id);
    const Page& p = page(pageIdx);
    const size_t i = id - ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
    e.id = static_cast<DBKeyType>(p.key[i]);
    e.offset = static_cast<size_t>(p.offset[i]);
    e.length = static_cast<unsigned int>(p.length[i]);
    return e;
}

void LazyIndex::readAll(Entry* out) {
    if (backend == BACKEND_BINARY) {
        const BinaryRecord* records = reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader));
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < meta.entries; i++) {
            out[i].id = records[i].id;
            out[i].offset = records[i].offset;
            out[i].length = records[i].length;
        }
        return;
    }
#pragma omp parallel
    {
        Page p;
#pragma omp for schedule(static)
        for (size_t pageIdx = 0; pageIdx < firstKey.size(); pageIdx++) {
            loadPage(pageIdx, p);
            Entry* o = out + ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
            for (size_t i = 0; i < p.rows; i++) {
                o[i].id = static_cast<DBKeyType>(p.key[i]);
                o[i].offset = static_cast<size_t>(p.offset[i]);
                o[i].length = static_cast<unsigned int>(p.length[i]);
            }
        }
    }
}
