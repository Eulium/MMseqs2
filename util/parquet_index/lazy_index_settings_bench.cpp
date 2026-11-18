// throwaway benchmark: an MMseqs2 .index (key, offset, length) as parquet, different write settings
// size, load of the page directory, full scan, random point lookups; baselines: array in RAM, binary file via mmap
#include <carquet/carquet.h>
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <string>
#include <vector>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }
static void die(const char* m, const carquet_error_t* e = NULL) { fprintf(stderr, "ERROR %s %s\n", m, e ? e->message : ""); exit(1); }

struct Idx { std::vector<int64_t> key, off; std::vector<int32_t> len; };
struct Cfg { const char* name; carquet_compression_t codec; int level; int enc; int rowsPerPage; int64_t rowsPerRG; };
// enc: 0 all PLAIN, 1 all DELTA_BINARY_PACKED, 2 key+offset delta, length plain, 3 all BYTE_STREAM_SPLIT

static void gen(Idx& d, size_t n, bool seqLike) {
    std::mt19937_64 rng(42);
    d.key.resize(n); d.off.resize(n); d.len.resize(n);
    if (seqLike) {   // sequence DB: consecutive keys, offsets increase with the key
        int64_t o = 0;
        for (size_t i = 0; i < n; i++) {
            d.key[i] = i; d.len[i] = 30 + (int32_t)(rng() % 700); d.off[i] = o; o += d.len[i];
        }
    } else {         // result DB: sparse keys, written by 16 threads in chunks, so offsets jump between thread files
        const int T = 16; const size_t chunk = 1000;
        std::vector<int> thr((n + chunk - 1) / chunk);
        std::vector<int64_t> tsize(T, 0), tbase(T, 0), tpos(T, 0);
        int64_t k = 0;
        for (size_t i = 0; i < n; i++) { k += 1 + rng() % 3; d.key[i] = k; d.len[i] = 1 + (int32_t)(rng() % 20000); }
        for (size_t c = 0; c < thr.size(); c++) thr[c] = rng() % T;
        for (size_t i = 0; i < n; i++) tsize[thr[i / chunk]] += d.len[i];
        for (int t = 1; t < T; t++) tbase[t] = tbase[t - 1] + tsize[t - 1];
        for (size_t i = 0; i < n; i++) { int t = thr[i / chunk]; d.off[i] = tbase[t] + tpos[t]; tpos[t] += d.len[i]; }
    }
}

static void writePq(const Idx& d, const std::string& file, const Cfg& c) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_schema_t* s = carquet_schema_create(&err);
    carquet_schema_add_column(s, "key", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "offset", CARQUET_PHYSICAL_INT64, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_schema_add_column(s, "length", CARQUET_PHYSICAL_INT32, NULL, CARQUET_REPETITION_REQUIRED, 0, 0);
    carquet_writer_options_t o; carquet_writer_options_init(&o);
    o.compression = c.codec; o.compression_level = c.level;
    o.row_group_size = (int64_t)1 << 40;   // row groups are cut by hand
    o.page_size = 64 * 1024 * 1024; o.max_rows_per_page = c.rowsPerPage;
    o.write_statistics = true; o.write_page_index = true; o.write_bloom_filters = false; o.write_crc = false;
    o.dictionary_encoding = CARQUET_ENCODING_PLAIN;
    carquet_writer_t* w = carquet_writer_create(file.c_str(), s, &o, &err);
    if (!w) die("writer", &err);
    carquet_schema_free(s);
    carquet_writer_set_column_statistics(w, 1, false); carquet_writer_set_column_statistics(w, 2, false);
    carquet_encoding_t e[3] = {CARQUET_ENCODING_PLAIN, CARQUET_ENCODING_PLAIN, CARQUET_ENCODING_PLAIN};
    if (c.enc == 1) e[0] = e[1] = e[2] = CARQUET_ENCODING_DELTA_BINARY_PACKED;
    if (c.enc == 2) e[0] = e[1] = CARQUET_ENCODING_DELTA_BINARY_PACKED;
    if (c.enc == 3) e[0] = e[1] = e[2] = CARQUET_ENCODING_BYTE_STREAM_SPLIT;
    for (int i = 0; i < 3; i++) if (carquet_writer_set_column_encoding(w, i, e[i]) != CARQUET_OK) die("encoding");
    const size_t n = d.key.size();
    for (size_t p = 0; p < n; p += c.rowsPerRG) {
        const int64_t m = std::min<size_t>(c.rowsPerRG, n - p);
        if (p > 0 && carquet_writer_new_row_group(w) != CARQUET_OK) die("new row group");
        if (carquet_writer_write_batch(w, 0, d.key.data() + p, m, NULL, NULL) != CARQUET_OK) die("write key");
        if (carquet_writer_write_batch(w, 1, d.off.data() + p, m, NULL, NULL) != CARQUET_OK) die("write off");
        if (carquet_writer_write_batch(w, 2, d.len.data() + p, m, NULL, NULL) != CARQUET_OK) die("write len");
    }
    if (carquet_writer_close(w) != CARQUET_OK) die("close");
}

struct PageDir { std::vector<int64_t> firstKey, firstRow, rows; };   // one entry per page of the key column

static carquet_reader_t* openPq(const std::string& file, bool mm) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    carquet_reader_options_t ro; carquet_reader_options_init(&ro); ro.use_mmap = mm;
    carquet_reader_t* r = carquet_reader_open(file.c_str(), &ro, &err);
    if (!r) die("open", &err);
    return r;
}
static void loadDir(carquet_reader_t* r, PageDir& pd) {
    carquet_error_t err = CARQUET_ERROR_INIT;
    int64_t rgFirst = 0;
    for (int32_t g = 0; g < carquet_reader_num_row_groups(r); g++) {
        carquet_row_group_metadata_t md; carquet_reader_row_group_metadata(r, g, &md);
        carquet_column_index_t* ci = carquet_reader_get_column_index(r, g, 0, &err);
        carquet_offset_index_t* oi = carquet_reader_get_offset_index(r, g, 0, &err);
        if (!ci || !oi) die("page index", &err);
        const int32_t np = carquet_column_index_num_pages(ci);
        for (int32_t p = 0; p < np; p++) {
            carquet_page_stats_t st; carquet_page_location_t loc;
            if (carquet_column_index_get_page_stats(ci, p, &st) != CARQUET_OK) die("page stats");
            if (carquet_offset_index_get_page_location(oi, p, &loc) != CARQUET_OK) die("page loc");
            int64_t k; memcpy(&k, st.min_value, sizeof(k));
            if (!pd.rows.empty() && pd.firstRow.back() >= rgFirst) pd.rows.back() = rgFirst + loc.first_row_index - pd.firstRow.back();
            pd.firstKey.push_back(k); pd.firstRow.push_back(rgFirst + loc.first_row_index); pd.rows.push_back(0);
        }
        rgFirst += md.num_rows;
        pd.rows.back() = rgFirst - pd.firstRow.back();
        carquet_column_index_free(ci); carquet_offset_index_free(oi);
    }
}

int main(int argc, char** argv) {
    const size_t n = argc > 1 ? strtoull(argv[1], NULL, 10) : 30000000;
    const bool seqLike = argc > 2 ? atoi(argv[2]) != 0 : true;
    const std::string dir = argc > 3 ? argv[3] : ".";
    const size_t Q = argc > 4 ? strtoull(argv[4], NULL, 10) : 200000;
    Idx d; gen(d, n, seqLike);
    std::mt19937_64 rng(7);
    std::vector<size_t> qrow(Q);
    for (size_t i = 0; i < Q; i++) qrow[i] = rng() % n;
    printf("# %zu entries, %s index, %zu random lookups\n", n, seqLike ? "sequence-DB-like" : "result-DB-like", Q);

    // baseline 1: text .index as MMseqs2 writes it, parsed into the 24 byte DBReader::Index array (single thread)
    {
        const std::string f = dir + "/text.index";
        FILE* fp = fopen(f.c_str(), "w");
        for (size_t i = 0; i < n; i++) fprintf(fp, "%lld\t%lld\t%d\n", (long long)d.key[i], (long long)d.off[i], d.len[i]);
        fclose(fp);
        struct stat sb; stat(f.c_str(), &sb);
        struct Index { unsigned int id; size_t offset; unsigned int length; };
        double t = now();
        int fd = open(f.c_str(), O_RDONLY);
        char* m = (char*)mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
        std::vector<Index> ix(n);
        const char* p = m;
        for (size_t i = 0; i < n; i++) {
            uint64_t v = 0; while (*p != '\t') v = v * 10 + (*p++ - '0'); ix[i].id = v; p++;
            v = 0; while (*p != '\t') v = v * 10 + (*p++ - '0'); ix[i].offset = v; p++;
            v = 0; while (*p != '\n') v = v * 10 + (*p++ - '0'); ix[i].length = v; p++;
        }
        const double tl = now() - t;
        munmap(m, sb.st_size); close(fd);
        t = now(); uint64_t sum = 0;
        for (size_t i = 0; i < Q; i++) {
            Index val; val.id = (unsigned int)d.key[qrow[i]];
            const Index* e = std::lower_bound(ix.data(), ix.data() + n, val, [](const Index& a, const Index& b) { return a.id < b.id; });
            sum += e->offset + e->length;
        }
        const double tq = now() - t;
        printf("%-34s size %8.1f MB  ram %7.1f MB  load %6.2f s  scan      -      lookup %8.3f us  (chk %llu)\n",
               "text .index -> Index[] in RAM", sb.st_size / 1048576.0, n * sizeof(Index) / 1048576.0, tl, tq / Q * 1e6, (unsigned long long)sum);
    }
    // baseline 2: binary records (key u32, length u32, offset u64 = 16 bytes), sorted by key, mmap + binary search
    {
        struct Rec { uint32_t key, len; uint64_t off; };
        const std::string f = dir + "/binary.idx";
        { std::vector<Rec> r(n); for (size_t i = 0; i < n; i++) { r[i].key = d.key[i]; r[i].len = d.len[i]; r[i].off = d.off[i]; }
          FILE* fp = fopen(f.c_str(), "wb"); fwrite(r.data(), sizeof(Rec), n, fp); fclose(fp); }
        double t = now();
        int fd = open(f.c_str(), O_RDONLY);
        const Rec* r = (const Rec*)mmap(NULL, n * sizeof(Rec), PROT_READ, MAP_PRIVATE, fd, 0);
        const double tl = now() - t;
        t = now(); uint64_t sum = 0;
        for (size_t i = 0; i < n; i++) sum += r[i].off + r[i].len;
        const double ts = now() - t;
        t = now(); sum = 0;
        for (size_t i = 0; i < Q; i++) {
            Rec val; val.key = (uint32_t)d.key[qrow[i]];
            const Rec* e = std::lower_bound(r, r + n, val, [](const Rec& a, const Rec& b) { return a.key < b.key; });
            sum += e->off + e->len;
        }
        const double tq = now() - t;
        printf("%-34s size %8.1f MB  ram %7.1f MB  load %6.2f s  scan %6.2f s  lookup %8.3f us  (chk %llu)\n",
               "binary 16 B records, mmap", n * sizeof(Rec) / 1048576.0, 0.0, tl, ts, tq / Q * 1e6, (unsigned long long)sum);
        munmap((void*)r, n * sizeof(Rec)); close(fd);
    }

    const carquet_compression_t U = CARQUET_COMPRESSION_UNCOMPRESSED, Z = CARQUET_COMPRESSION_ZSTD,
                                L = CARQUET_COMPRESSION_LZ4_RAW, S = CARQUET_COMPRESSION_SNAPPY;
    const int64_t RG = 4 * 1024 * 1024;
    const Cfg cfgs[] = {
        {"delta none  p4096 rg4M", U, 0, 1, 4096, RG},       {"delta zstd3 p4096 rg4M", Z, 3, 1, 4096, RG},
        {"delta none  p4096 rg256k", U, 0, 1, 4096, 262144}, {"delta zstd3 p4096 rg256k", Z, 3, 1, 4096, 262144},
        {"delta none  p4096 rg64k", U, 0, 1, 4096, 65536},   {"delta zstd3 p4096 rg64k", Z, 3, 1, 4096, 65536},
        {"delta lz4   p4096 rg64k", L, 0, 1, 4096, 65536},
        {"delta none  p1024 rg64k", U, 0, 1, 1024, 65536},   {"delta zstd3 p1024 rg64k", Z, 3, 1, 1024, 65536},
        {"delta none  p1024 rg16k", U, 0, 1, 1024, 16384},   {"delta zstd3 p1024 rg16k", Z, 3, 1, 1024, 16384},
        {"delta none  p256  rg16k", U, 0, 1, 256, 16384},    {"delta zstd3 p256  rg16k", Z, 3, 1, 256, 16384},
        {"plain none  p1024 rg16k", U, 0, 0, 1024, 16384},   {"plain zstd3 p1024 rg16k", Z, 3, 0, 1024, 16384},
        {"plain lz4   p1024 rg16k", L, 0, 0, 1024, 16384},
        {"dlt-ko zstd3 p1024 rg16k", Z, 3, 2, 1024, 16384},  {"bss   zstd3 p1024 rg16k", Z, 3, 3, 1024, 16384},
        {"delta none  p4096 rg4096", U, 0, 1, 4096, 4096},   {"delta zstd3 p4096 rg4096", Z, 3, 1, 4096, 4096},
    };
    for (const Cfg& c : cfgs) {
        const std::string f = dir + "/idx.parquet";
        double t = now(); writePq(d, f, c); const double tw = now() - t;
        struct stat sb; stat(f.c_str(), &sb);
        for (int mm = 1; mm <= 1; mm++) {
            t = now();
            carquet_reader_t* r = openPq(f, mm);
            PageDir pd; loadDir(r, pd);
            const double tl = now() - t;
            carquet_error_t err = CARQUET_ERROR_INIT;
            carquet_batch_reader_config_t bc; carquet_batch_reader_config_init(&bc);
            bc.batch_size = 65536; bc.num_threads = 1; bc.use_mmap = mm;
            carquet_batch_reader_t* b = carquet_batch_reader_create(r, &bc, &err);
            if (!b) die("batch", &err);
            // full scan, one thread
            t = now(); uint64_t sum = 0; int64_t seen = 0; carquet_row_batch_t* batch = NULL;
            while (carquet_batch_reader_next(b, &batch) == CARQUET_OK && batch != NULL) {
                const void* k; const void* o; const void* l; const uint8_t* nb; int64_t m;
                carquet_row_batch_column(batch, 0, &k, &nb, &m); carquet_row_batch_column(batch, 1, &o, &nb, &m);
                carquet_row_batch_column(batch, 2, &l, &nb, &m);
                for (int64_t i = 0; i < m; i++) sum += ((const int64_t*)o)[i] + ((const int32_t*)l)[i];
                seen += m; carquet_row_batch_free(batch); batch = NULL;
            }
            const double ts = now() - t;
            if ((size_t)seen != n) die("scan rows");
            carquet_batch_reader_free(b);
            bc.batch_size = c.rowsPerPage;
            b = carquet_batch_reader_create(r, &bc, &err);
            // random point lookups: page from the directory, decode that page, binary search in it
            t = now(); sum = 0; size_t bad = 0;
            for (size_t i = 0; i < Q; i++) {
                const int64_t key = d.key[qrow[i]];
                const size_t p = std::upper_bound(pd.firstKey.begin(), pd.firstKey.end(), key) - pd.firstKey.begin() - 1;
                if (carquet_batch_reader_set_row_range(b, pd.firstRow[p], pd.rows[p]) != CARQUET_OK) die("row range");
                bool found = false;
                while (carquet_batch_reader_next(b, &batch) == CARQUET_OK && batch != NULL) {
                    const void* k; const void* o; const void* l; const uint8_t* nb; int64_t m;
                    carquet_row_batch_column(batch, 0, &k, &nb, &m); carquet_row_batch_column(batch, 1, &o, &nb, &m);
                    carquet_row_batch_column(batch, 2, &l, &nb, &m);
                    const int64_t* kk = (const int64_t*)k;
                    const int64_t* e = std::lower_bound(kk, kk + m, key);
                    if (e < kk + m && *e == key) { found = true; sum += ((const int64_t*)o)[e - kk] + ((const int32_t*)l)[e - kk]; }
                    carquet_row_batch_free(batch); batch = NULL;
                }
                bad += !found;
            }
            const double tq = now() - t;
            // method B: column readers of the row group; key column: skip to the page, read the page, binary search;
            // offset and length: skip to the row, read one value
            std::vector<int64_t> rgFirst(1, 0);
            for (int32_t g = 0; g < carquet_reader_num_row_groups(r); g++) {
                carquet_row_group_metadata_t md; carquet_reader_row_group_metadata(r, g, &md); rgFirst.push_back(rgFirst.back() + md.num_rows);
            }
            std::vector<int64_t> kbuf(c.rowsPerPage);
            t = now(); uint64_t sumB = 0; size_t badB = 0;
            for (size_t i = 0; i < Q; i++) {
                const int64_t key = d.key[qrow[i]];
                const size_t p = std::upper_bound(pd.firstKey.begin(), pd.firstKey.end(), key) - pd.firstKey.begin() - 1;
                const int32_t g = std::upper_bound(rgFirst.begin(), rgFirst.end(), pd.firstRow[p]) - rgFirst.begin() - 1;
                const int64_t inRg = pd.firstRow[p] - rgFirst[g];
                carquet_column_reader_t* ck = carquet_reader_get_column(r, g, 0, &err);
                if (!ck) die("col", &err);
                if (inRg > 0 && carquet_column_skip(ck, inRg) != inRg) die("skip key");
                const int64_t m = carquet_column_read_batch(ck, kbuf.data(), pd.rows[p], NULL, NULL);
                carquet_column_reader_free(ck);
                const int64_t* e = std::lower_bound(kbuf.data(), kbuf.data() + m, key);
                if (e == kbuf.data() + m || *e != key) { badB++; continue; }
                const int64_t row = inRg + (e - kbuf.data());
                int64_t o = 0; int32_t l = 0;
                carquet_column_reader_t* co = carquet_reader_get_column(r, g, 1, &err);
                carquet_column_reader_t* cl = carquet_reader_get_column(r, g, 2, &err);
                if (row > 0) { carquet_column_skip(co, row); carquet_column_skip(cl, row); }
                if (carquet_column_read_batch(co, &o, 1, NULL, NULL) != 1 || carquet_column_read_batch(cl, &l, 1, NULL, NULL) != 1) die("read 1");
                carquet_column_reader_free(co); carquet_column_reader_free(cl);
                sumB += o + l;
            }
            const double tqB = now() - t;
            printf("%-34s size %8.1f MB  load %6.3f s  scan %6.2f s  lookup A %8.2f us  B %8.2f us  write %5.1f s  pages %zu rg %d%s%s\n",
                   c.name, sb.st_size / 1048576.0, tl, ts, tq / Q * 1e6, tqB / Q * 1e6, tw,
                   pd.firstKey.size(), carquet_reader_num_row_groups(r), bad ? "  A MISSES!" : "", (badB || sumB != sum) ? "  B WRONG!" : "");
            fflush(stdout);
            carquet_batch_reader_free(b); carquet_reader_close(r);
        }
    }
    return 0;
}
