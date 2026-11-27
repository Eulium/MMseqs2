// Checksum of the content of a DB that does not depend on the order of the entries in the index or in the data file:
// the sum over all entries of a hash of (key, data). Two runs of a module that wrote the same entries give the same sum.
// usage: test_dbchecksum <DB> [<DB> ...]   (the index is read as the DBReader finds it: text or lazy index)
// With DBCHECKSUM_ENTRIES set, one line per entry (key, hash, length) is written to stderr, to find the entries
// that differ between two DBs.
#include <cstdint>
#include <cstdlib>
#include <cstdio>
#include <string>

#include "DBReader.h"
#include "Debug.h"

const char* binary_name = "test_dbchecksum";

int main(int argc, const char** argv) {
    for (int a = 1; a < argc; a++) {
        const std::string db = argv[a];
        DBReader<DBKeyType> reader(db.c_str(), (db + ".index").c_str(), 1,
                                   DBReader<DBKeyType>::USE_INDEX | DBReader<DBKeyType>::USE_DATA);
        reader.open(DBReader<DBKeyType>::NOSORT);
        uint64_t sum = 0, bytes = 0;
        const bool entries = getenv("DBCHECKSUM_ENTRIES") != NULL;
        for (size_t id = 0; id < reader.getSize(); id++) {
            const unsigned char* data = reinterpret_cast<const unsigned char*>(reader.getData(id, 0));
            const size_t len = reader.getEntryLen(id);
            uint64_t h = 1469598103934665603ULL ^ (static_cast<uint64_t>(reader.getDbKey(id)) * 0x9E3779B97F4A7C15ULL);
            for (size_t i = 0; i + 1 < len; i++) {
                h = (h ^ data[i]) * 1099511628211ULL;
            }
            sum += h;
            if (entries) {
                fprintf(stderr, "%llu\t%016llx\t%zu\n", static_cast<unsigned long long>(reader.getDbKey(id)),
                        static_cast<unsigned long long>(h), len);
            }
            bytes += len;
        }
        printf("%016llx\t%zu\t%llu\t%s\n", static_cast<unsigned long long>(sum), reader.getSize(),
               static_cast<unsigned long long>(bytes), db.c_str());
        reader.close();
    }
    return 0;
}
