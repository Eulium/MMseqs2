#ifndef MMSEQS_LAZYINDEX_H
#define MMSEQS_LAZYINDEX_H

// An MMseqs2 .index (key, offset, length) that is not loaded into memory.
//
// DBReader normally parses the text index into an array of 24 bytes per entry. A LazyIndex keeps the index on disk in
// a second file next to the text index and reads only what is asked for:
//   <db>.index.parquet  parquet file (written with carquet): columns key (int64), offset (int64), length (int32),
//                       sorted by key, DELTA_BINARY_PACKED, no codec, 256 rows per page, 16384 rows per row group,
//                       page index. Only one first key per page is kept in memory (8 to 32 bytes per 256 entries).
//                       A lookup decodes one page (the three columns) into a per thread cache.
//                       Two readers for the same file: BACKEND_CARQUET (carquet column readers) and BACKEND_CUSTOM
//                       (own decoder of the page header and the DELTA_BINARY_PACKED values on the mapped file).
//   <db>.index.bin      binary file: header + one 16 byte record per entry, sorted by key, memory mapped, binary search.
//                       Reference for the benchmark: no decoding, but 3-6x the size of the parquet file.
// The entry order is the key order, so the id of an entry (row) is the id DBReader uses in the modes NOSORT/SORT_BY_ID.
//
// Prototype switches (environment variables, no module parameters yet):
//   MMSEQS_INDEX_BACKEND=parquet|parquet-custom|binary   DBReader uses the file if it exists and fits the text index
//   MMSEQS_WRITE_LAZY_INDEX=parquet,binary               DBWriter::close also writes these files
//   MMSEQS_WRITE_LAZY_INDEX=parquet,notext               DBWriter::close writes the parquet file INSTEAD of the text
//                                                        index; DBReader then reads the parquet file in every mode
//   MMSEQS_LAZY_INDEX_TEXT_MEMORY=1                     test mode: readers register the memory of a text index
//                                                        (same memory budget and k-mer parts as with the text index)
//   MMSEQS_INDEX_DEBUG=1                                 DBReader reports each index it opens with a backend

#include "IndexTypes.h"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

struct carquet_reader;

class LazyIndex {
public:
    // same layout as DBReader<DBKeyType>::Index
    struct Entry {
        DBKeyType id;
        size_t offset;
        unsigned int length;
    };

    // what DBReader::readIndex computes while parsing the text index
    struct Meta {
        size_t entries;
        size_t dataSize;            // sum of the entry lengths
        DBKeyType lastKey;
        unsigned int maxSeqLen;     // largest entry length
        bool sortedByOffset;        // offsets increase with the key
        size_t textIndexBytes;      // size of the text index the file was made from (stale check)
    };

    enum Backend {
        BACKEND_NONE = 0,
        BACKEND_CARQUET,
        BACKEND_CUSTOM,
        BACKEND_BINARY
    };

    // MMSEQS_INDEX_BACKEND
    static Backend backendFromEnv();
    static const char* backendName(Backend backend);
    static std::string fileName(const std::string& indexFile, Backend backend);

    // writes the entries it is given in key order, one row group at a time
    class Writer {
    public:
        // empty name: that file is not written
        Writer(const std::string& parquetFile, const std::string& binaryFile);
        ~Writer();
        void add(DBKeyType key, size_t offset, unsigned int length);
        // textIndexBytes: size of the text index with the same entries (stale check of the reader)
        void close(size_t textIndexBytes);
        const Meta& getMeta() const;
    private:
        Writer(const Writer&);
        Writer& operator=(const Writer&);
        struct Impl;
        Impl* impl;
    };

    // entries sorted by key; meta.entries entries
    static void writeParquet(const std::string& file, const Entry* entries, const Meta& meta);
    static void writeBinary(const std::string& file, const Entry* entries, const Meta& meta);
    // merges index parts by key without loading them (see LazyIndex.cpp); shifts[i] is added to the offsets of part i.
    // A part is a text index or, if that does not exist, its parquet/binary file.
    // Writes the text index textOut (empty: none) and/or the parquet/binary file. removeParts: the parts are deleted.
    // textBytesWithoutText: stored as size of the text index if none is written (0: there is no text index)
    static void mergeTextIndexes(const std::vector<std::string>& parts, const std::vector<size_t>& shifts,
                                 const std::string& textOut, bool sortText,
                                 const std::string& parquetFile, const std::string& binaryFile,
                                 bool removeParts, size_t textBytesWithoutText);
    // MMSEQS_WRITE_LAZY_INDEX (parquet, binary or both, comma separated; with notext the text index is not written):
    // the files DBReader looks for next to indexFile
    static bool writeRequested(std::string* parquetFile, std::string* binaryFile, const std::string& indexFile,
                               bool* noText = NULL);
    // for DBWriter::mergeResults, instead of mergeIndex + sortIndex: if MMSEQS_WRITE_LAZY_INDEX is set, the index parts
    // of a DB (dataSizes: size of the data of each part) become outIndex and its parquet/binary file. false: not set
    static bool mergeIndexParts(const char** indexFilenames, unsigned int fileCount, const std::vector<size_t>& dataSizes,
                                const char* outIndex, bool sortText);
    // backend of a DB without text index (written with notext): BACKEND_NONE if the text index exists or no file is there
    static Backend backendWithoutText(const std::string& indexFile);
    // removes the index of a DB that is used as plain file (flat output): the text index and/or its parquet/binary file
    static void removeIndex(const std::string& indexFile);
    // entries of the index (text lines, or from the parquet/binary file if there is no text index)
    static size_t countEntries(const std::string& indexFile);
    // writes an existing text index as parquet and/or binary file (empty name: not written)
    static void writeFromTextIndex(const std::string& indexFile, const std::string& parquetFile,
                                   const std::string& binaryFile, int threads);
    // true while the text index must be read: DBReader must not use a lazy index then
    static bool isDisabled();

    // NULL if the file does not exist, does not fit the text index or cannot be read by the backend
    static LazyIndex* open(const std::string& indexFile, Backend backend);
    ~LazyIndex();

    const Meta& getMeta() const { return meta; }
    Backend getBackend() const { return backend; }
    // bytes kept in memory for the page directory
    size_t memoryBytes() const;

    // id (row in key order) of the first entry with this key, SIZE_MAX if there is none; thread safe
    size_t getId(DBKeyType key);
    // entry id (0 .. entries - 1); thread safe
    Entry get(size_t id);
    // one member of the entry: only that column of the page is decoded, up to the row
    DBKeyType getKey(size_t id);
    size_t getOffset(size_t id);
    unsigned int getLength(size_t id);
    // all entries in key order
    void readAll(Entry* out);
    // one column of all entries in key order (for the orders DBReader builds from it)
    void readLengths(unsigned int* out);
    void readOffsets(size_t* out);

private:
    LazyIndex();
    LazyIndex(const LazyIndex&);
    LazyIndex& operator=(const LazyIndex&);

    struct Page;
    Page& page(size_t pageIdx);
    void loadPage(size_t pageIdx, Page& p);
    void setPage(size_t pageIdx, Page& p);
    void decodeColumn(Page& p, int column, size_t want, const int64_t* stopAt);
    size_t pageOfRow(size_t row) const;
    size_t rowsOfPage(size_t pageIdx) const;
    bool openParquet(const std::string& file, bool custom);
    bool openBinary(const std::string& file);

    Backend backend;
    Meta meta;
    size_t serial;                      // unique per instance, the per thread page cache is found by it

    // parquet: one entry per page (pages hold the same rows in the three columns)
    carquet_reader* reader;
    size_t rowsPerPage;                 // rows of every page but the last of a row group; 0 = not uniform
    std::vector<int64_t> firstKey;      // first key of each page
    std::vector<size_t> firstRow;       // first row of each page, only if rowsPerPage == 0
    std::vector<int32_t> rowGroup;      // carquet backend: row group of each page
    std::vector<size_t> rowGroupFirstRow;
    std::vector<uint64_t> pageOffset;   // custom backend: 3 file offsets per page (key, offset, length)

    // mapped file (custom backend: the parquet file, binary backend: the record file)
    char* map;
    size_t mapSize;
};

#endif
