# Session summary: parquet index for MMseqs2 (branch parquet_index, written 2026-11-18)

State of the work on reading the MMseqs2 .index from a parquet file instead of loading it into memory, so that it
can be continued (also by a new Claude Code session, on the laptop or on spark004). Everything below was true when
written; verify paths and running jobs before relying on them.

Files in this folder (`util/parquet_index/`):
- `SESSION_SUMMARY.md`: this file
- `lazy_index_test.sh`, `lazy_index_bench.sh`, `lazy_index_cluster_bench.sh`: correctness test and benchmarks
- `lazy_index_spark.sh`, `lazy_index_spark2.sh`: drivers of the benchmark runs on spark004
- `peak_anon.py`: peak RssAnon of a process tree, used by the benchmarks (has to sit next to them)
- `lazy_index_settings_bench.cpp`: throwaway benchmark of the parquet write settings (links libcarquet.a, zstd, lz4)
- `NOTES_parquet_writer.md`: results and ideas of the earlier parquet writer / reader work
The same scripts are in `examples/parquet_benchmark/` of the working copy (ignored by git), together with the other
benchmark scripts named in section 4.

## 0. How to continue on spark004

- Copy of the repo with all uncommitted changes: `~/mmseqs_session/MMseqs2` on spark004 (user `user2`).
  Start the session there: `cd ~/mmseqs_session/MMseqs2 && claude`, first message e.g.
  "Read util/parquet_index/SESSION_SUMMARY.md and continue" (that copy was made before the branch existed and has
  the summary as `examples/parquet_benchmark/SESSION_HANDOVER.md`).
- The copy has no build directory. Build it inside the repo (so it is within the session's working directory):
  `mkdir -p build && cd build && cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_SHELLCHECK=0 .. && make -j20 mmseqs`
  (about 3-5 minutes on 20 cores; not while a timed benchmark is running). The binary used by the running
  benchmarks is `~/mmseqs_lazy_index/build2/src/mmseqs`, built from identical source.
- Not copied (too large, laptop only): `build/`, `examples/parquet_benchmark/{data,out,raw,logs,stress_557114}`,
  `examples/parquet_tests`, `examples/tmp`. Test data on spark004 is listed in section 5.
- The benchmark working folders are outside that directory (`~/mmseqs_lazy_index`, see section 5).
- The persistent memory of the laptop session is not copied. Its only entry: "TODO: cA3M parquet column as BLOB —
  deferred fix, cA3M written as UTF-8 string is unreadable".

## 2. Repository state

- Repo: MMseqs2 fork (origin Eulium/MMseqs2, upstream soedinglab). Branch `parquet_index`, created from `parquet`
  and not pushed; `parquet` has base MMseqs2 plus the parquet writer and `convert2parquet`. carquet 0.8.0 is a subtree in
  `lib/carquet`.
- Base commit: `59feac4d1 small changes that improve performance on subseted databases` (committed by the user;
  contains the subset index fix and the three parsing changes in ParquetDBWriter).
- The parquet index feature (section 3) on top of it:
  - modified: `src/commons/DBReader.h`, `DBReader.cpp`, `DBReaderSortIndex.cpp`, `DBWriter.cpp`,
    `src/commons/CMakeLists.txt`, `src/util/convert2parquet.cpp`
  - new: `src/commons/LazyIndex.h`, `LazyIndex.cpp`
  - untracked and paused (not built, ignore): `src/commons/ParquetDBReader.{h,cpp}`
- `examples/parquet_benchmark/` is in `.gitignore`: scripts and results there are not pushed.
- Build: C++14 (`-std=c++1y`), CMake. Laptop build dir `build/Integrated_Native` (not copied to spark004).

## 3. Current task: index on disk instead of in memory ("lazy index")

### Request
MMseqs2 loads the whole `.index` into RAM (24 bytes per entry: `DBReader<DBKeyType>::Index`). New feature: after a DB
is written, also write the index as a compressed parquet file; when it exists, modules look entries up in that file
instead of loading the index. The user asked for: best write settings, alternatives to carquet, a test module, a
benchmark, the right place in the DB handler, the implementation, stability/speed tests, other parquet readers.
User decisions: build BOTH a carquet reader and a custom decoder; include a binary mmap index as reference; benchmark
on spark004 with a medium-sized DB (UniRef50); also test the clustering workflows (largest inputs).

### Findings before implementation (20M synthetic entries, laptop, warm cache)
- Text index 443-508 MB on disk, 458 MB in RAM; binary 16 B records 305 MB; parquet with DELTA_BINARY_PACKED 57-89 MB.
- The encoding does the compression: zstd/lz4/snappy on top of delta saves < 1 % and slows lookups.
- carquet's cost per lookup grows with the pages per row group -> small row groups (16384 rows) and pages (1024 rows).
- Random lookup: 0.6 us in RAM; 15-19 us via carquet with these settings; plain+zstd 48-81 us.
- Alternatives (web survey): no small C/C++ library beats carquet; Rust `parquet` crate is the mature option but
  needs a Rust toolchain + C shim; DuckDB is 50-60 MB and needs SQL per lookup; nanoparquet/miniparquet unsuitable.
- Settings benchmark source: `lazy_index_settings_bench.cpp` in this folder.

### Implementation
- `LazyIndex` (`src/commons/LazyIndex.{h,cpp}`), three backends for a file next to the text index:
  - `<db>.index.parquet` via carquet column readers (BACKEND_CARQUET)
  - the same file via an own decoder: thrift compact page header + DELTA_BINARY_PACKED on the mmapped file
    (BACKEND_CUSTOM; carquet is only used for the footer and page index at open)
  - `<db>.index.bin`: header + 16 byte records, mmap, binary search (BACKEND_BINARY)
- Parquet layout: columns key int64, offset int64, length int32; sorted by key; delta encoding, no codec; 1024
  rows/page; 16384 rows/row group; page index; statistics on key only; key-value metadata `mmseqs.index.*`
  (entries, dataSize, lastKey, maxSeqLen, sortedByOffset, textIndexBytes).
- In RAM: the first key of each page (8-32 bytes per 1024 entries). Each thread caches one decoded page per LazyIndex
  (thread_local, found by an instance serial).
- `DBReader::open` -> `openLazyIndex()`: used for NOSORT, SORT_BY_ID, and LINEAR_ACCCESS if offsets increase with the
  key. Other modes read the text index as before. `getIndex()`, `sortIndex(weights)`, `getDataOffset`,
  `findNextOffsetid`, `serialize`, `decomposeDomainByAminoAcid` call `materializeIndex()` (read everything into the
  array, close the lazy index). `DBReader<std::string>` never uses it.
- Stale check: the file stores the size of the text index it was made from; mismatch -> warning, text index is read.
- Writing: `DBWriter::close` -> `LazyIndex::writeAfterDbWrite` when `MMSEQS_WRITE_LAZY_INDEX` is set;
  `mmseqs convert2parquet <db>.index <db>.index.parquet` for existing DBs (with `MMSEQS_WRITE_LAZY_INDEX=binary` also
  `<db>.index.bin`).
- Switches are environment variables (prototype, no module parameters yet):
  - `MMSEQS_INDEX_BACKEND=parquet|parquet-custom|binary`
  - `MMSEQS_WRITE_LAZY_INDEX=parquet,binary`
  - `MMSEQS_INDEX_DEBUG=1`: DBReader prints "Index backend ..." for every index it opens (also for text, with
    entries, bytes and access mode). MMseqs2 prints warnings to stdout; use `-v 2` or higher to see them.

### Verified
- `lazy_index_test.sh`: 0 failed checks on the laptop (Swiss-Prot, 574,627 entries) and on spark004 (1.29M entries,
  sparse keys, ARM): createsubdb, view, convert2fasta, search + convertalis byte-identical for all backends; also
  with a lazy index for every DB the search writes; 1 and all threads; stale files rejected.
- Clustering (laptop): linclust (Swiss-Prot, 337,226 clusters) and cluster (20,000 sequences, 8,420 clusters) give
  identical clusters for all backends, input-only and all-DBs. Not yet run on spark004.
- DuckDB reads the parquet index and returns exactly the rows of the text index.
- Size, Swiss-Prot: text 11.6 MB, parquet 1.8 MB, binary 9.2 MB; 7 KB in RAM instead of 13 MB.
- Speed, laptop, one round only: random lookup ~12 us carquet, ~6 us custom, 0.4 us binary, 0.5 us text; lookups
  in key order cost the same for all.

### Known limits / open
- Clustering only partly benefits: `align2clust` opens the sequence DB with SORT_BY_LENGTH (mode 1) and result DBs
  are opened with LINEAR_ACCCESS (mode 2) with offsets not in key order -> text index. Supporting these needs the
  lazy index combined with the local2id/id2local permutation (8 B/entry instead of 32). Not built.
- Custom decoder always decodes all three columns of a page; decoding only up to the requested row is a possible
  speedup. Not built.
- Rust `parquet` crate not tried (step 7 covers DuckDB CLI and R arrow only).
- Cold-cache behaviour not measured (no root to drop caches). 1e9-entry scale not tested (61,000 row groups: carquet
  footer parsing cost unknown).
- Module parameters instead of environment variables; default on/off: not decided.
- In `DBReader::open` the variable `isSortedById` is shadowed in the text branch (upstream behaviour, open() always
  returns false); left as it is.

## 4. Benchmark scripts (this folder; all of them in examples/parquet_benchmark/ of the working copy)

| Script | Purpose |
|---|---|
| `lazy_index_test.sh <mmseqs> <seqDB with _h> <workdir> [queries]` | correctness of all backends incl. linclust/cluster |
| `lazy_index_bench.sh <mmseqs> <seqDB> <workdir> [rounds] [threads] [queries]` | lookup benchmark: view, createsubdb with 1e3/1e5/1e7/all keys (random and sorted), convertalis (20 and 1 threads), DuckDB/R arrow on the same file |
| `lazy_index_cluster_bench.sh <mmseqs> <seqDB> <workdir> <sequences> [rounds] [threads] [workflow]` | linclust/cluster per backend, counts lazy vs text index opens |
| `lazy_index_spark.sh <base>` | driver 1 on spark004: wait for download, createdb, lookup benchmark |
| `lazy_index_spark2.sh <base> <small DB>` | driver 2: waits for driver 1, then correctness test, linclust 20M, cluster 500k |
| `peak_anon.py <out> <cmd...>` | peak RssAnon of a process tree (allocated memory without mapped files); must sit next to the benchmark scripts |
| `run_parquet_benchmark.sh`, `run_stress_test.sh`, `run_stress_test.sbatch` | TSV vs parquet conversion/query benchmark and stress test; measure and plot RssAnon |
| `subset_index_bench.sh`, `parse_checks_test.sh`, `carquet_regression_test.sh`, `anon_test_spark.sh` | earlier tests (subset index fix, parsing changes, carquet update, RssAnon) |
| `make_stacked_db.sh`, `setup_conda_env.sh`, `tsv_merge_proto_bench.sh` | helpers; the createtsv merge prototype was never finished |
| `TODO.md` | results and ideas of the earlier writer/reader work (section 6b: writer benchmarks) |

Both lazy index benchmarks record `other_cores`: CPU cores used by OTHER processes during each run (spark004 is
shared). Workloads above 2 cores are marked as disturbed.

## 5. spark004 (Brev Connect node, aarch64, 20 cores, 121 GB RAM, shared user user2)

- Access from the laptop: `ssh spark004`. `scp` and `brev copy` hang; upload with `ssh spark004 'cat > file' < file`.
  Start long jobs with `setsid nohup ... &`. Do not `pkill -f <pattern>` with a pattern that occurs in your own ssh
  command line (it kills the ssh shell) - kill by pid.
- The machine is shared: user3 runs batch jobs; other Claude sessions run as user2 (`~/mmseqs_c2p_opt` and
  `~/MMseqs2` belong to another session - do not touch).
- conda: `source ~/miniforge3/etc/profile.d/conda.sh && conda activate rbench` (R 4.4.3, arrow, data.table, dplyr,
  duckdb CLI 1.5.1; no r-duckdb on aarch64). Rust in `~/.cargo`.
- Folders of this work:
  - `~/mmseqs_session/MMseqs2`: repo copy for a new session (this handover)
  - `~/mmseqs_lazy_index/`: `MMseqs2` (source used for the builds), `build` (first build), `build2` (current
    binary, with the text-path debug line), `data/` (UniRef50 download `parts/`, then `uniref50` DB), scripts,
    `test1` (correctness run, passed), later `bench_uniref50`, `test2`, `cluster_linclust`, `cluster_cluster`
  - older: `~/mmseqs_subset_index_test`, `~/mmseqs_anon_test`, `~/mmseqs_carquet080_test`,
    `~/mmseqs_parquet_writer_bench`, `~/mmseqs_parquet_reader_tests` (test data such as the stacked result DB
    `mmseqs_parquet_writer_bench/stress_conv/stack/res8` and the stress DBs under
    `mmseqs_carquet080_test/stress_full/scratch/stress/db/`). The user decides when to delete them.

### Jobs running at handover (check before starting anything)
- UniRef50 download from the Expasy mirror in 16 parts (`~/mmseqs_lazy_index/data/parts`, total 9,182,717,121
  bytes); short parts were resumed, the last missing ranges as sub-chunks `part_NN.r0` .. `.r7` (they sort
  right after `part_NN`, so `cat part_*` gives the file; the driver checks the total size and `gzip -t`). The EBI mirror was too slow (30-80 KB/s);
  single connections drop to ~80 KB/s after a burst, parallel ranges reach 6-12 MB/s.
- `lazy_index_spark.sh` (log `~/mmseqs_lazy_index/spark_run.log`): waits for the full size (max 3 h), `gzip -t`,
  `mmseqs createdb` -> `data/uniref50`, then `lazy_index_bench.sh` (3 rounds, 20 threads, 2000 queries) into
  `bench_uniref50/` (`lazy_index.csv`, `index_files.csv`, `other_readers.csv`).
- `lazy_index_spark2.sh` (log `spark_run2.log`): waits for driver 1, then `test2`, `cluster_linclust` (linclust of
  20M sequences, 1 round per variant), `cluster_cluster` (cluster of 500k sequences).
- Check: `pgrep -u user2 -fa lazy_index_spark; tail ~/mmseqs_lazy_index/spark_run*.log`.
- No timed measurement on spark004 for the lazy index existed yet at handover.

## 6. Earlier results of this project (short; details in NOTES_parquet_writer.md)

- Parquet writer (`ParquetDBWriter`, `convert2parquet`): at 20 GiB uncompressed 1.35x and zstd 1.23x faster than
  createtsv; parquet zstd 1.30x faster and 20 % smaller than `createtsv --db-output --compressed`.
- carquet 0.8.0 update: regression test 0 failures; files 0-2.4 % larger (max_rows_per_page fix).
- Subset index fix (in HEAD): convert2parquet on a subset index (createsubdb --subdb-mode 1) 3.7-4.1x faster, output
  byte-identical to the compact DB; before, 17 of 20 output files were empty.
- Three parsing changes (in HEAD): one compare per byte in the TSV parser, parseInteger without per-digit branch,
  projection always filled; identical output, no measurable speed change.
- Peak RSS of GNU time counts memory mapped DB pages; RssAnon (peak_anon.py) is the memory that must fit into RAM.
  Example 2 GiB conversion: createtsv 766 MB anon / 3,430 MB RSS; convert2parquet zstd 816 / 2,860; uncompressed
  2,263 / 4,283.
- `--parquet-schema` never converts types: non-generic DBs must match the layout type; generic DBs parse strictly.

## 7. Open items and unanswered questions

1. Report the UniRef50 lookup and clustering numbers when the spark004 runs finish; repeat disturbed runs.
2. Lazy index for SORT_BY_LENGTH / unsorted LINEAR_ACCCESS (needed for a real memory saving in clustering).
3. `availableMemory()` in ParquetDBWriter ignores cgroup/slurm limits; proposal: `--split-memory-limit` + walk the
   cgroup hierarchy (`/proc/self/cgroup`), min with MemAvailable and RLIMIT_AS. Not applied.
4. Proposals not applied: `writeCell(const char*, size_t)` overload (no temporary std::string); MADV_RANDOM or
   pread for subset indexes; faster R `best_hit` (setorder + unique); collinear/parallel merge for createtsv.
5. Unanswered: which spark004 test folders to delete.
