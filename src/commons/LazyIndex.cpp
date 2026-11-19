#include "LazyIndex.h"
#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "itoa.h"
#include "Util.h"

#include <carquet/include/carquet/carquet.h>

#include <algorithm>
#include <functional>
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
// 256 rows: a random lookup costs as much as in the text index (1024: 1.7x), the file is 1.4x larger (2.3 B/entry)
const int64_t ROWS_PER_PAGE = 256;
const int64_t ROWS_PER_ROW_GROUP = 16384;
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

// rows per page of the files that are written; MMSEQS_LAZY_INDEX_PAGE_ROWS (a divisor of the rows per row group) is
// for benchmarks: smaller pages make a random lookup cheaper and the page directory in memory larger
int64_t rowsPerPageForWrite() {
    const char* env = getenv("MMSEQS_LAZY_INDEX_PAGE_ROWS");
    const int64_t rows = (env != NULL) ? strtoll(env, NULL, 10) : 0;
    return (rows >= 16 && rows <= ROWS_PER_ROW_GROUP && ROWS_PER_ROW_GROUP % rows == 0) ? rows : ROWS_PER_PAGE;
}

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
// width per miniblock, the miniblocks (deltas - min delta, bit packed). The page has n values; decoding stops after
// the miniblock that reaches want values or (stopAt != NULL, ascending values) a value >= *stopAt.
// decoded: number of values written to out
template <typename O>
bool decodeDelta(Cursor& c, const uint8_t* fileEnd, O* out, size_t n, size_t want, const int64_t* stopAt, size_t& decoded) {
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
    decoded = 1;
    while (produced < n && produced < want && (stopAt == NULL || static_cast<int64_t>(prev) < *stopAt)) {
        const uint64_t minDelta = static_cast<uint64_t>(zigzag(readVarint(c)));
        const uint8_t* widths = c.p;
        skipBytes(c, miniblocks);
        if (c.ok == false) {
            return false;
        }
        for (uint64_t m = 0; m < miniblocks && produced < n && produced < want
                                 && (stopAt == NULL || static_cast<int64_t>(prev) < *stopAt); m++) {
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
    decoded = produced;
    return true;
}

// one data page of a REQUIRED column (no levels) at p; n values, see decodeDelta
template <typename O>
bool decodePage(const uint8_t* p, const uint8_t* fileEnd, O* out, size_t n, size_t want, const int64_t* stopAt, size_t& decoded) {
    Cursor c = {p, fileEnd, true};
    PageHeader h;
    if (readPageHeader(c, h) == false || (h.type != 0 && h.type != 3) || h.encoding != 5
        || h.numValues != static_cast<int64_t>(n) || h.compressedSize != h.uncompressedSize || h.compressedSize < 0
        || static_cast<uint64_t>(fileEnd - c.p) < static_cast<uint64_t>(h.compressedSize)) {
        return false;
    }
    Cursor values = {c.p, c.p + h.compressedSize, true};
    skipBytes(values, static_cast<uint64_t>(h.levelBytes));
    return values.ok && decodeDelta(values, fileEnd, out, n, want, stopAt, decoded);
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
    // values decoded so far per column (key, offset, length): a lookup decodes only the column and the rows it needs
    size_t decoded[3];
    Page() : serial(0), pageIdx(SIZE_MAX), rows(0) { decoded[0] = decoded[1] = decoded[2] = 0; }
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

// ---------------------------------------------------------------------------------------------------------------
// writing: entries are added in key order and written in batches, nothing but one batch is kept in memory

namespace {
struct MappedFile {
    char* data;
    size_t size;
    MappedFile() : data(NULL), size(0) {}
    // an empty or missing file stays unmapped
    void open(const std::string& file) {
        const size_t bytes = fileBytes(file);
        if (bytes == SIZE_MAX || bytes == 0) {
            return;
        }
        const int fd = ::open(file.c_str(), O_RDONLY);
        void* p = (fd < 0) ? MAP_FAILED : mmap(NULL, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        if (fd >= 0) {
            ::close(fd);
        }
        if (p == MAP_FAILED) {
            Debug(Debug::ERROR) << "Cannot map " << file << "\n";
            EXIT(EXIT_FAILURE);
        }
        madvise(p, bytes, MADV_SEQUENTIAL);
        data = static_cast<char*>(p);
        size = bytes;
    }
    void close() {
        if (data != NULL) {
            munmap(data, size);
            data = NULL;
            size = 0;
        }
    }
};

// one line "key \t offset \t length \n" of a text index; false at the end of the file
inline bool parseIndexLine(const char*& pos, const char* end, LazyIndex::Entry& e) {
    while (pos < end && (*pos == '\n' || *pos == '\0')) {
        pos++;
    }
    if (pos >= end) {
        return false;
    }
    uint64_t v[3] = {0, 0, 0};
    for (int c = 0; c < 3; c++) {
        while (pos < end && *pos >= '0' && *pos <= '9') {
            v[c] = v[c] * 10 + static_cast<uint64_t>(*pos - '0');
            pos++;
        }
        if (c < 2 && pos < end && *pos == '\t') {
            pos++;
        }
    }
    while (pos < end && *pos != '\n') {
        pos++;
    }
    e.id = static_cast<DBKeyType>(v[0]);
    e.offset = static_cast<size_t>(v[1]);
    e.length = static_cast<unsigned int>(v[2]);
    return true;
}

inline char* formatIndexLine(char* buffer, DBKeyType key, size_t offset, unsigned int length) {
    char* p = Itoa::u64toa_sse2(static_cast<uint64_t>(key), buffer);
    *(p - 1) = '\t';
    p = Itoa::u64toa_sse2(static_cast<uint64_t>(offset), p);
    *(p - 1) = '\t';
    p = Itoa::u32toa_sse2(length, p);
    *(p - 1) = '\n';
    return p;
}

bool entryByKey(const LazyIndex::Entry& a, const LazyIndex::Entry& b) {
    return a.id < b.id;
}

// one index part as a stream of entries in key order: the text part itself if it is sorted, otherwise a file with
// its sorted entries
struct IndexRun {
    MappedFile file;
    std::string runFile;
    bool binary;
    size_t shift;
    const char* pos;
    const char* end;
    // part without text index: its parquet/binary file, read row by row
    LazyIndex* lazy;
    size_t row;
    LazyIndex::Entry cur;
    IndexRun() : binary(false), shift(0), pos(NULL), end(NULL), lazy(NULL), row(0) {}
    void rewind() {
        row = 0;
        pos = file.data;
        end = pos + file.size;
    }
    bool next() {
        if (lazy != NULL) {
            if (row >= lazy->getMeta().entries) {
                return false;
            }
            cur = lazy->get(row++);
        } else if (binary) {
            if (pos >= end) {
                return false;
            }
            memcpy(&cur, pos, sizeof(cur));
            pos += sizeof(cur);
        } else if (parseIndexLine(pos, end, cur) == false) {
            return false;
        }
        cur.offset += shift;
        return true;
    }
};
}

struct LazyIndex::Writer::Impl {
    std::string parquetFile, binaryFile;
    carquet_writer_t* writer;
    FILE* binary;
    Meta meta;
    size_t prevOffset;
    DBKeyType prevKey;
    bool rowGroupOpen;
    std::vector<int64_t> key, offset;
    std::vector<int32_t> length;
    std::vector<BinaryRecord> records;
    size_t fill;
    Impl() : writer(NULL), binary(NULL), prevOffset(0), prevKey(0), rowGroupOpen(false), fill(0) {
        memset(&meta, 0, sizeof(meta));
        meta.sortedByOffset = true;
    }
    void flush() {
        if (fill == 0) {
            return;
        }
        const int64_t n = static_cast<int64_t>(fill);
        if (writer != NULL) {
            if ((rowGroupOpen && carquet_writer_new_row_group(writer) != CARQUET_OK)
                || carquet_writer_write_batch(writer, 0, key.data(), n, NULL, NULL) != CARQUET_OK
                || carquet_writer_write_batch(writer, 1, offset.data(), n, NULL, NULL) != CARQUET_OK
                || carquet_writer_write_batch(writer, 2, length.data(), n, NULL, NULL) != CARQUET_OK) {
                carquetFail("Cannot write the entries of " + parquetFile, NULL);
            }
            rowGroupOpen = true;
        }
        if (binary != NULL) {
            for (size_t i = 0; i < fill; i++) {
                records[i].offset = static_cast<size_t>(offset[i]);
                records[i].id = static_cast<DBKeyType>(key[i]);
                records[i].length = static_cast<unsigned int>(length[i]);
            }
            if (fwrite(records.data(), sizeof(BinaryRecord), fill, binary) != fill) {
                Debug(Debug::ERROR) << "Cannot write " << binaryFile << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
        fill = 0;
    }
};

LazyIndex::Writer::Writer(const std::string& parquetFile, const std::string& binaryFile) : impl(new Impl) {
    impl->parquetFile = parquetFile;
    impl->binaryFile = binaryFile;
    impl->key.resize(ROWS_PER_ROW_GROUP);
    impl->offset.resize(ROWS_PER_ROW_GROUP);
    impl->length.resize(ROWS_PER_ROW_GROUP);
    if (parquetFile.empty() == false) {
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_schema_t* schema = carquet_schema_create(&err);
        if (schema == NULL
            || carquet_schema_add_column(schema, "key", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK
            || carquet_schema_add_column(schema, "offset", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK
            || carquet_schema_add_column(schema, "length", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0) != CARQUET_OK) {
            carquetFail("Cannot create the schema of " + parquetFile, &err);
        }
        carquet_writer_options_t opts;
        carquet_writer_options_init(&opts);
        // DELTA_BINARY_PACKED is the compression: a codec on top saves < 1% and slows down every page read
        opts.compression = CARQUET_COMPRESSION_UNCOMPRESSED;
        // row groups and pages are cut by rows: small pages make a lookup cheap, small row groups keep carquet's cost
        // per column reader low (it grows with the pages of the row group)
        opts.row_group_size = static_cast<int64_t>(1) << 40;
        opts.page_size = 64 * 1024 * 1024;
        opts.max_rows_per_page = rowsPerPageForWrite();
        opts.write_statistics = true;
        opts.write_page_index = true;
        opts.write_bloom_filters = false;
        opts.write_crc = false;
        opts.dictionary_encoding = CARQUET_ENCODING_PLAIN;
        opts.created_by = "MMseqs parquet index";
        const std::string tmp = parquetFile + ".tmp";
        impl->writer = carquet_writer_create(tmp.c_str(), schema, &opts, &err);
        carquet_schema_free(schema);
        if (impl->writer == NULL) {
            carquetFail("Cannot write " + tmp, &err);
        }
        bool ok = true;
        // page statistics only for the key: the first key of every page is the page directory of the reader
        ok = ok && carquet_writer_set_column_statistics(impl->writer, 1, false) == CARQUET_OK;
        ok = ok && carquet_writer_set_column_statistics(impl->writer, 2, false) == CARQUET_OK;
        for (int32_t c = 0; c < 3; c++) {
            ok = ok && carquet_writer_set_column_encoding(impl->writer, c, CARQUET_ENCODING_DELTA_BINARY_PACKED) == CARQUET_OK;
        }
        if (ok == false) {
            carquetFail("Cannot set the writer options of " + tmp, NULL);
        }
    }
    if (binaryFile.empty() == false) {
        impl->records.resize(ROWS_PER_ROW_GROUP);
        impl->binary = FileUtil::openFileOrDie((binaryFile + ".tmp").c_str(), "wb", false);
        // the header is written again by close, when the values are known
        BinaryHeader header;
        memset(&header, 0, sizeof(header));
        if (fwrite(&header, sizeof(header), 1, impl->binary) != 1) {
            Debug(Debug::ERROR) << "Cannot write " << binaryFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }
}

LazyIndex::Writer::~Writer() {
    delete impl;
}

void LazyIndex::Writer::add(DBKeyType key, size_t offset, unsigned int length) {
    Impl& w = *impl;
    if (w.meta.entries > 0 && key < w.prevKey) {
        Debug(Debug::ERROR) << "Entries of the index " << (w.parquetFile.empty() ? w.binaryFile : w.parquetFile)
                            << " are not added in key order\n";
        EXIT(EXIT_FAILURE);
    }
    w.meta.sortedByOffset = w.meta.sortedByOffset && offset >= w.prevOffset;
    w.prevOffset = offset;
    w.prevKey = key;
    w.meta.entries++;
    w.meta.dataSize += length;
    w.meta.maxSeqLen = std::max(w.meta.maxSeqLen, length);
    w.meta.lastKey = std::max(w.meta.lastKey, key);
    w.key[w.fill] = static_cast<int64_t>(key);
    w.offset[w.fill] = static_cast<int64_t>(offset);
    w.length[w.fill] = static_cast<int32_t>(length);
    if (++w.fill == static_cast<size_t>(ROWS_PER_ROW_GROUP)) {
        w.flush();
    }
}

const LazyIndex::Meta& LazyIndex::Writer::getMeta() const {
    return impl->meta;
}

// an index without entries gets no file (as a missing file, the text index is read)
void LazyIndex::Writer::close(size_t textIndexBytes) {
    Impl& w = *impl;
    w.flush();
    w.meta.textIndexBytes = textIndexBytes;
    const Meta& meta = w.meta;
    if (w.writer != NULL) {
        const std::string tmp = w.parquetFile + ".tmp";
        bool ok = true;
        ok = ok && carquet_writer_add_metadata(w.writer, META_ENTRIES, SSTR(meta.entries).c_str()) == CARQUET_OK;
        ok = ok && carquet_writer_add_metadata(w.writer, META_DATA_SIZE, SSTR(meta.dataSize).c_str()) == CARQUET_OK;
        ok = ok && carquet_writer_add_metadata(w.writer, META_LAST_KEY, SSTR(meta.lastKey).c_str()) == CARQUET_OK;
        ok = ok && carquet_writer_add_metadata(w.writer, META_MAX_SEQ_LEN, SSTR(meta.maxSeqLen).c_str()) == CARQUET_OK;
        ok = ok && carquet_writer_add_metadata(w.writer, META_SORTED_BY_OFFSET, meta.sortedByOffset ? "1" : "0") == CARQUET_OK;
        ok = ok && carquet_writer_add_metadata(w.writer, META_TEXT_INDEX_BYTES, SSTR(meta.textIndexBytes).c_str()) == CARQUET_OK;
        if (ok == false || carquet_writer_close(w.writer) != CARQUET_OK) {
            carquetFail("Cannot close " + tmp, NULL);
        }
        w.writer = NULL;
        if (meta.entries > 0) {
            FileUtil::move(tmp.c_str(), w.parquetFile.c_str());
        } else {
            FileUtil::remove(tmp.c_str());
        }
    }
    if (w.binary != NULL) {
        const std::string tmp = w.binaryFile + ".tmp";
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
        const bool ok = fseek(w.binary, 0, SEEK_SET) == 0 && fwrite(&header, sizeof(header), 1, w.binary) == 1;
        if (fclose(w.binary) != 0 || ok == false) {
            Debug(Debug::ERROR) << "Cannot write " << tmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        w.binary = NULL;
        if (meta.entries > 0) {
            FileUtil::move(tmp.c_str(), w.binaryFile.c_str());
        } else {
            FileUtil::remove(tmp.c_str());
        }
    }
}

void LazyIndex::writeParquet(const std::string& file, const Entry* entries, const Meta& meta) {
    Writer writer(file, "");
    for (size_t i = 0; i < meta.entries; i++) {
        writer.add(entries[i].id, entries[i].offset, entries[i].length);
    }
    writer.close(meta.textIndexBytes);
}

void LazyIndex::writeBinary(const std::string& file, const Entry* entries, const Meta& meta) {
    Writer writer("", file);
    for (size_t i = 0; i < meta.entries; i++) {
        writer.add(entries[i].id, entries[i].offset, entries[i].length);
    }
    writer.close(meta.textIndexBytes);
}

// Index parts (text indexes, e.g. one per writer thread) -> one index, without loading it:
//   1. each part becomes a run of entries in key order: the part itself if it is sorted (the usual case), otherwise
//      its entries are sorted in memory (one part at a time) and written to a file next to it
//   2. the runs are merged by key; shifts[i] is added to the offsets of part i
// textOut (may be empty: no text index is written):
//   sortText == true:  the merged entries, as DBWriter::sortIndex writes them
//   sortText == false: the parts one after the other, as DBWriter::mergeIndex writes them
// parquetFile/binaryFile (may be empty) always get the entries in key order.
void LazyIndex::mergeTextIndexes(const std::vector<std::string>& parts, const std::vector<size_t>& shifts,
                                 const std::string& textOut, bool sortText,
                                 const std::string& parquetFile, const std::string& binaryFile,
                                 bool removeParts, size_t textBytesWithoutText) {
    const size_t n = parts.size();
    std::vector<IndexRun> runs(n);
    std::vector<char> sorted(n, 1);
    for (size_t i = 0; i < n; i++) {
        runs[i].shift = shifts[i];
        const Backend backend = backendWithoutText(parts[i]);
        if (backend != BACKEND_NONE) {
            runs[i].lazy = open(parts[i], backend);
        } else {
            runs[i].file.open(parts[i]);
        }
    }
#pragma omp parallel for schedule(dynamic, 1)
    for (size_t i = 0; i < n; i++) {
        const char* pos = runs[i].file.data;
        const char* end = pos + runs[i].file.size;
        Entry e;
        DBKeyType prev = 0;
        while (parseIndexLine(pos, end, e)) {
            if (e.id < prev) {
                sorted[i] = 0;
                break;
            }
            prev = e.id;
        }
    }
    for (size_t i = 0; i < n; i++) {
        if (sorted[i]) {
            continue;
        }
        std::vector<Entry> entries;
        const char* pos = runs[i].file.data;
        const char* end = pos + runs[i].file.size;
        Entry e;
        memset(&e, 0, sizeof(e));
        while (parseIndexLine(pos, end, e)) {
            entries.push_back(e);
        }
        // equal keys keep the order of the part
        std::stable_sort(entries.begin(), entries.end(), entryByKey);
        runs[i].runFile = parts[i] + ".run";
        FILE* fp = FileUtil::openFileOrDie(runs[i].runFile.c_str(), "wb", false);
        const bool ok = fwrite(entries.data(), sizeof(Entry), entries.size(), fp) == entries.size();
        if (fclose(fp) != 0 || ok == false) {
            Debug(Debug::ERROR) << "Cannot write " << runs[i].runFile << "\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const std::string textTmp = textOut + ".merge.tmp";
    FILE* text = textOut.empty() ? NULL : FileUtil::openFileOrDie(textTmp.c_str(), "wb", false);
    char line[128];
    bool textOk = true;
    if (text != NULL && sortText == false) {
        for (size_t i = 0; i < n && textOk; i++) {
            if (i == 0 && runs[i].lazy == NULL) {
                // DBWriter::mergeIndex appends to the first part
                textOk = runs[i].file.size == 0 || fwrite(runs[i].file.data, 1, runs[i].file.size, text) == runs[i].file.size;
                continue;
            }
            // the entries in the order of the part
            runs[i].rewind();
            while (textOk && runs[i].next()) {
                const size_t len = static_cast<size_t>(formatIndexLine(line, runs[i].cur.id, runs[i].cur.offset, runs[i].cur.length) - line);
                textOk = fwrite(line, 1, len, text) == len;
            }
        }
    }

    std::vector<MappedFile> runFiles(n);
    // smallest key first, equal keys in the order of the parts
    std::vector<std::pair<DBKeyType, size_t> > heap;
    const std::greater<std::pair<DBKeyType, size_t> > cmp = std::greater<std::pair<DBKeyType, size_t> >();
    for (size_t i = 0; i < n; i++) {
        runs[i].rewind();
        if (sorted[i] == 0) {
            runs[i].file.close();
            runFiles[i].open(runs[i].runFile);
            runs[i].binary = true;
            runs[i].pos = runFiles[i].data;
            runs[i].end = runs[i].pos + runFiles[i].size;
        }
        if (runs[i].next()) {
            heap.push_back(std::make_pair(runs[i].cur.id, i));
        }
    }
    std::make_heap(heap.begin(), heap.end(), cmp);
    const bool needMerge = parquetFile.empty() == false || binaryFile.empty() == false || (text != NULL && sortText);
    Writer writer(parquetFile, binaryFile);
    while (needMerge && heap.empty() == false && textOk) {
        std::pop_heap(heap.begin(), heap.end(), cmp);
        IndexRun& run = runs[heap.back().second];
        // the run stays on top as long as its next key is not larger than the smallest key of the other runs
        const bool others = heap.size() > 1;
        const std::pair<DBKeyType, size_t> limit = others ? heap.front() : std::make_pair(DBKeyType(0), size_t(0));
        bool more = true;
        while (more) {
            writer.add(run.cur.id, run.cur.offset, run.cur.length);
            if (text != NULL && sortText) {
                const size_t len = static_cast<size_t>(formatIndexLine(line, run.cur.id, run.cur.offset, run.cur.length) - line);
                textOk = textOk && fwrite(line, 1, len, text) == len;
            }
            more = run.next();
            if (more && others && std::make_pair(run.cur.id, heap.back().second) > limit) {
                break;
            }
        }
        if (more) {
            heap.back().first = run.cur.id;
            std::push_heap(heap.begin(), heap.end(), cmp);
        } else {
            heap.pop_back();
        }
    }
    size_t textBytes = textBytesWithoutText;
    if (text != NULL) {
        if (fclose(text) != 0 || textOk == false) {
            Debug(Debug::ERROR) << "Cannot write " << textTmp << "\n";
            EXIT(EXIT_FAILURE);
        }
        textBytes = fileBytes(textTmp);
    }
    for (size_t i = 0; i < n; i++) {
        runs[i].file.close();
        runFiles[i].close();
        delete runs[i].lazy;
        runs[i].lazy = NULL;
        if (runs[i].runFile.empty() == false) {
            FileUtil::remove(runs[i].runFile.c_str());
        }
    }
    // the parts go before the new files appear: a part can have the name of the result
    if (removeParts) {
        for (size_t i = 0; i < n; i++) {
            const std::string files[3] = {parts[i], fileName(parts[i], BACKEND_CARQUET), fileName(parts[i], BACKEND_BINARY)};
            for (size_t f = 0; f < 3; f++) {
                if (FileUtil::fileExists(files[f].c_str())) {
                    FileUtil::remove(files[f].c_str());
                }
            }
        }
    }
    writer.close(textBytes);
    if (text != NULL) {
        FileUtil::move(textTmp.c_str(), textOut.c_str());
    }
}

bool LazyIndex::writeRequested(std::string* parquetFile, std::string* binaryFile, const std::string& indexFile, bool* noText) {
    const char* env = getenv("MMSEQS_WRITE_LAZY_INDEX");
    if (env == NULL || env[0] == '\0') {
        return false;
    }
    const std::string what(env);
    *parquetFile = what.find("parquet") != std::string::npos ? fileName(indexFile, BACKEND_CARQUET) : "";
    *binaryFile = what.find("binary") != std::string::npos ? fileName(indexFile, BACKEND_BINARY) : "";
    if (noText != NULL) {
        *noText = what.find("notext") != std::string::npos;
    }
    return parquetFile->empty() == false || binaryFile->empty() == false;
}

bool LazyIndex::mergeIndexParts(const char** indexFilenames, unsigned int fileCount, const std::vector<size_t>& dataSizes,
                                const char* outIndex, bool sortText) {
    std::string parquetFile, binaryFile;
    bool noText = false;
    if (writeRequested(&parquetFile, &binaryFile, outIndex, &noText) == false) {
        return false;
    }
    std::vector<std::string> parts;
    std::vector<size_t> shifts;
    size_t shift = 0;
    for (unsigned int i = 0; i < fileCount; i++) {
        parts.push_back(indexFilenames[i]);
        shifts.push_back(shift);
        shift += (i < dataSizes.size()) ? dataSizes[i] : 0;
    }
    if (noText) {
        // a text index of an earlier DB with this name would be taken for the current one
        if (FileUtil::fileExists(outIndex)) {
            FileUtil::remove(outIndex);
        }
        mergeTextIndexes(parts, shifts, "", sortText, parquetFile, binaryFile, true, 0);
    } else {
        mergeTextIndexes(parts, shifts, outIndex, sortText, parquetFile, binaryFile, true, 0);
    }
    return true;
}

LazyIndex::Backend LazyIndex::backendWithoutText(const std::string& indexFile) {
    if (fileBytes(indexFile) != SIZE_MAX) {
        return BACKEND_NONE;
    }
    if (fileBytes(fileName(indexFile, BACKEND_CUSTOM)) != SIZE_MAX) {
        const Backend env = backendFromEnv();
        return (env == BACKEND_CARQUET) ? BACKEND_CARQUET : BACKEND_CUSTOM;
    }
    if (fileBytes(fileName(indexFile, BACKEND_BINARY)) != SIZE_MAX) {
        return BACKEND_BINARY;
    }
    return BACKEND_NONE;
}

void LazyIndex::removeIndex(const std::string& indexFile) {
    const std::string files[3] = {indexFile, fileName(indexFile, BACKEND_CARQUET), fileName(indexFile, BACKEND_BINARY)};
    for (size_t f = 0; f < 3; f++) {
        if (FileUtil::fileExists(files[f].c_str())) {
            FileUtil::remove(files[f].c_str());
        }
    }
}

size_t LazyIndex::countEntries(const std::string& indexFile) {
    const Backend backend = backendWithoutText(indexFile);
    if (backend == BACKEND_NONE) {
        return FileUtil::countLines(indexFile.c_str());
    }
    LazyIndex* idx = open(indexFile, backend);
    const size_t entries = (idx != NULL) ? idx->meta.entries : 0;
    delete idx;
    return entries;
}

// the text index is not changed; it does not have to be sorted
void LazyIndex::writeFromTextIndex(const std::string& indexFile, const std::string& parquetFile,
                                   const std::string& binaryFile, int) {
    if (parquetFile.empty() && binaryFile.empty()) {
        return;
    }
    mergeTextIndexes(std::vector<std::string>(1, indexFile), std::vector<size_t>(1, 0), "", false, parquetFile, binaryFile,
                     false, fileBytes(indexFile));
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
    size_t uniformRows = 0;
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
            if (firstKey.size() == 1) {
                uniformRows = row;
            }
            uniform = uniform && row == firstKey.size() * uniformRows;
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
        rowsPerPage = (firstKey.size() > 1) ? uniformRows : std::max<size_t>(meta.entries, 1);
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

// the page without decoded values
void LazyIndex::setPage(size_t pageIdx, Page& p) {
    const size_t rows = rowsOfPage(pageIdx);
    if (p.key.size() < rows) {
        p.key.resize(rows);
        p.offset.resize(rows);
        p.length.resize(rows);
    }
    p.serial = serial;
    p.pageIdx = pageIdx;
    p.rows = rows;
    p.decoded[0] = p.decoded[1] = p.decoded[2] = 0;
}

// at least the first want values of a column (0 key, 1 offset, 2 length) of the page; for the key column instead
// (stopAt != NULL) up to the first key >= *stopAt
void LazyIndex::decodeColumn(Page& p, int column, size_t want, const int64_t* stopAt) {
    want = std::min(want, p.rows);
    if (stopAt == NULL) {
        if (p.decoded[column] >= want) {
            return;
        }
        // a second request for more rows of the same page: all of them
        if (p.decoded[column] > 0) {
            want = p.rows;
        }
    } else if (p.decoded[column] == p.rows || (p.decoded[column] > 0 && p.key[p.decoded[column] - 1] >= *stopAt)) {
        return;
    }
    const size_t pageIdx = p.pageIdx;
    void* out = (column == 0) ? static_cast<void*>(p.key.data())
                : (column == 1 ? static_cast<void*>(p.offset.data()) : static_cast<void*>(p.length.data()));
    bool ok = true;
    size_t decoded = 0;
    if (backend == BACKEND_CUSTOM) {
        const uint8_t* base = reinterpret_cast<const uint8_t*>(map);
        const uint8_t* end = base + mapSize;
        const uint8_t* at = base + pageOffset[3 * pageIdx + column];
        ok = (column == 2) ? decodePage(at, end, static_cast<int32_t*>(out), p.rows, want, stopAt, decoded)
                           : decodePage(at, end, static_cast<int64_t*>(out), p.rows, want, stopAt, decoded);
    } else {
        // carquet reads whole batches: all rows of the page for a key search
        if (stopAt != NULL) {
            want = p.rows;
        }
        const int32_t g = rowGroup[pageIdx];
        const size_t first = (rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx];
        const int64_t skip = static_cast<int64_t>(first - rowGroupFirstRow[g]);
        carquet_error_t err = CARQUET_ERROR_INIT;
        carquet_column_reader_t* reader_ = carquet_reader_get_column(reader, g, column, &err);
        ok = reader_ != NULL && (skip == 0 || carquet_column_skip(reader_, skip) == skip)
             && carquet_column_read_batch(reader_, out, static_cast<int64_t>(want), NULL, NULL) == static_cast<int64_t>(want);
        carquet_column_reader_free(reader_);
        decoded = want;
    }
    if (ok == false) {
        Debug(Debug::ERROR) << "Cannot read page " << pageIdx << " of the " << backendName(backend) << " index\n";
        EXIT(EXIT_FAILURE);
    }
    p.decoded[column] = decoded;
}

// all values of the page
void LazyIndex::loadPage(size_t pageIdx, Page& p) {
    setPage(pageIdx, p);
    for (int c = 0; c < 3; c++) {
        decodeColumn(p, c, p.rows, NULL);
    }
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
        setPage(pageIdx, *p);
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
        Page& p = page(q - 1);
        // the keys of the page up to the first one that is not smaller
        decodeColumn(p, 0, p.rows, &k);
        const std::vector<int64_t>::const_iterator last = p.key.cbegin() + p.decoded[0];
        const std::vector<int64_t>::const_iterator it = std::lower_bound(p.key.cbegin(), last, k);
        if (it != last && *it == k) {
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
    Page& p = page(pageIdx);
    const size_t i = id - ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
    for (int c = 0; c < 3; c++) {
        decodeColumn(p, c, i + 1, NULL);
    }
    e.id = static_cast<DBKeyType>(p.key[i]);
    e.offset = static_cast<size_t>(p.offset[i]);
    e.length = static_cast<unsigned int>(p.length[i]);
    return e;
}

DBKeyType LazyIndex::getKey(size_t id) {
    if (backend == BACKEND_BINARY) {
        return reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader))[id].id;
    }
    const size_t pageIdx = pageOfRow(id);
    Page& p = page(pageIdx);
    const size_t i = id - ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
    decodeColumn(p, 0, i + 1, NULL);
    return static_cast<DBKeyType>(p.key[i]);
}

size_t LazyIndex::getOffset(size_t id) {
    if (backend == BACKEND_BINARY) {
        return reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader))[id].offset;
    }
    const size_t pageIdx = pageOfRow(id);
    Page& p = page(pageIdx);
    const size_t i = id - ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
    decodeColumn(p, 1, i + 1, NULL);
    return static_cast<size_t>(p.offset[i]);
}

unsigned int LazyIndex::getLength(size_t id) {
    if (backend == BACKEND_BINARY) {
        return reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader))[id].length;
    }
    const size_t pageIdx = pageOfRow(id);
    Page& p = page(pageIdx);
    const size_t i = id - ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
    decodeColumn(p, 2, i + 1, NULL);
    return static_cast<unsigned int>(p.length[i]);
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

void LazyIndex::readLengths(unsigned int* out) {
    if (backend == BACKEND_BINARY) {
        const BinaryRecord* records = reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader));
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < meta.entries; i++) {
            out[i] = records[i].length;
        }
        return;
    }
#pragma omp parallel
    {
        Page p;
#pragma omp for schedule(static)
        for (size_t pageIdx = 0; pageIdx < firstKey.size(); pageIdx++) {
            loadPage(pageIdx, p);
            unsigned int* o = out + ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
            for (size_t i = 0; i < p.rows; i++) {
                o[i] = static_cast<unsigned int>(p.length[i]);
            }
        }
    }
}

void LazyIndex::readOffsets(size_t* out) {
    if (backend == BACKEND_BINARY) {
        const BinaryRecord* records = reinterpret_cast<const BinaryRecord*>(map + sizeof(BinaryHeader));
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < meta.entries; i++) {
            out[i] = records[i].offset;
        }
        return;
    }
#pragma omp parallel
    {
        Page p;
#pragma omp for schedule(static)
        for (size_t pageIdx = 0; pageIdx < firstKey.size(); pageIdx++) {
            loadPage(pageIdx, p);
            size_t* o = out + ((rowsPerPage != 0) ? pageIdx * rowsPerPage : firstRow[pageIdx]);
            for (size_t i = 0; i < p.rows; i++) {
                o[i] = static_cast<size_t>(p.offset[i]);
            }
        }
    }
}
