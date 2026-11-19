// Compares every access mode of DBReader on two copies of the same DB:
//   reference: the DB with its text index only
//   test:      the DB as the lazy index reads it (parquet/binary file, with or without text index; see LazyIndex.h)
// For each mode the two readers must give the same key, lengths, offset and data for every local id and the same
// local id for every key (and for keys that are not in the DB).
// usage: test_lazyindex <reference DB> <test DB> [threads]
// MMSEQS_INDEX_BACKEND selects the backend of the test DB (the reference DB has no lazy index file).
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include "DBReader.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LazyIndex.h"

#ifdef OPENMP
#include <omp.h>
#endif

const char* binary_name = "test_lazyindex";

static const char* MODE_NAME[] = {"NOSORT", "SORT_BY_LENGTH", "LINEAR_ACCCESS", "SORT_BY_ID", "SORT_BY_LINE", "SHUFFLE",
                                  "HARDNOSORT", "SORT_BY_ID_OFFSET", "SORT_BY_OFFSET", "SORT_BY_WEIGHTS"};

// number of differences between the two readers (0: the same)
static size_t compare(DBReader<DBKeyType>& ref, DBReader<DBKeyType>& tst, bool lookups, int threads) {
    if (ref.getSize() != tst.getSize()) {
        std::cout << "    size " << ref.getSize() << " != " << tst.getSize() << "\n";
        return 1;
    }
    size_t diff = 0;
    if (ref.getMaxSeqLen() != tst.getMaxSeqLen() || ref.getLastKey() != tst.getLastKey()
        || ref.getAminoAcidDBSize() != tst.getAminoAcidDBSize() || ref.isSortedByOffset() != tst.isSortedByOffset()) {
        std::cout << "    maxSeqLen/lastKey/dataSize/sortedByOffset differ\n";
        diff++;
    }
    const size_t n = ref.getSize();
#pragma omp parallel for schedule(dynamic, 10000) reduction(+:diff) num_threads(threads)
    for (size_t id = 0; id < n; id++) {
        int thread = 0;
#ifdef OPENMP
        thread = omp_get_thread_num();
#endif
        const DBKeyType key = ref.getDbKey(id);
        bool same = key == tst.getDbKey(id) && ref.getSeqLen(id) == tst.getSeqLen(id)
                    && ref.getEntryLen(id) == tst.getEntryLen(id) && ref.getOffset(id) == tst.getOffset(id);
        if (same) {
            const char* a = ref.getData(id, thread);
            const char* b = tst.getData(id, thread);
            const size_t len = ref.getEntryLen(id);
            same = memcmp(a, b, len > 0 ? len - 1 : 0) == 0;
        }
        if (same && lookups) {
            // the id of a key (the first entry with this key), and of a key next to it that may not exist
            same = ref.getId(key) == tst.getId(key) && ref.getId(key + 1) == tst.getId(key + 1);
        }
        diff += same ? 0 : 1;
    }
    return diff;
}

int main(int argc, const char** argv) {
    if (argc < 3) {
        std::cerr << "usage: test_lazyindex <reference DB> <test DB> [threads]\n";
        return EXIT_FAILURE;
    }
    const std::string refDb = argv[1], tstDb = argv[2];
    const int threads = argc > 3 ? atoi(argv[3]) : 1;
#ifdef OPENMP
    omp_set_num_threads(threads);
#endif
    const int dataMode = DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA;
    size_t failed = 0;
    for (int mode = 0; mode <= 9; mode++) {
        DBReader<DBKeyType> ref(refDb.c_str(), (refDb + ".index").c_str(), threads, dataMode);
        DBReader<DBKeyType> tst(tstDb.c_str(), (tstDb + ".index").c_str(), threads, dataMode);
        if (mode == DBReader<DBKeyType>::SORT_BY_WEIGHTS) {
            // as Clustering does: open sorted by key, then order by weights
            ref.open(DBReader<DBKeyType>::SORT_BY_ID);
            tst.open(DBReader<DBKeyType>::SORT_BY_ID);
            std::vector<float> weights(ref.getSize());
            for (size_t i = 0; i < weights.size(); i++) {
                weights[i] = static_cast<float>((i * 2654435761u) % 1000);
            }
            ref.sortIndex(weights.data());
            tst.sortIndex(weights.data());
        } else {
            ref.open(mode);
            tst.open(mode);
        }
        // SORT_BY_OFFSET and HARDNOSORT keep no key order: getId is not defined there
        const bool lookups = mode != DBReader<DBKeyType>::SORT_BY_OFFSET && mode != DBReader<DBKeyType>::HARDNOSORT;
        const size_t diff = compare(ref, tst, lookups, threads);
        std::cout << (diff == 0 ? "  OK    " : "  FAIL  ") << MODE_NAME[mode] << ": " << ref.getSize() << " entries, "
                  << diff << " differences\n";
        failed += diff == 0 ? 0 : 1;
        // the index as array (materialized from the lazy index)
        if (mode == DBReader<DBKeyType>::NOSORT) {
            const DBReader<DBKeyType>::Index* a = ref.getIndex();
            const DBReader<DBKeyType>::Index* b = tst.getIndex();
            size_t d = 0;
            for (size_t i = 0; i < ref.getSize(); i++) {
                d += (a[i].id == b[i].id && a[i].offset == b[i].offset && a[i].length == b[i].length) ? 0 : 1;
            }
            d += compare(ref, tst, true, threads);
            std::cout << (d == 0 ? "  OK    " : "  FAIL  ") << "getIndex (array): " << d << " differences\n";
            failed += d == 0 ? 0 : 1;
        }
        ref.close();
        tst.close();
    }
    std::cout << (failed == 0 ? "== all modes the same\n" : "== FAILED modes: ") << (failed == 0 ? "" : SSTR(failed) + "\n");
    return failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}
