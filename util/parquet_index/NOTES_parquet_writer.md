---
name: todo-ca3m-blob-column
description: "Parquet project status and open todos (writer, convert2parquet, benchmarks, stress test, ParquetDBReader, spark004 setup, cA3M BLOB, Linclust step 6)"
metadata:
  node_type: memory
  type: project
  originSessionId: 30d6be04-d7b1-465c-9663-7e39898f5323
  modified: 2026-11-16T02:50:42.512Z
---

# Parquet for MMseqs2: status, results, ideas, todos (as of 2026-11-16)

Repo: /home/me/Documents/SteineggerLab/parquet_writer/MMseqs2, branch `parquet`. Nothing of the work below is committed.
Uncommitted: src/commons/ParquetDBReader.{h,cpp} (new), ParquetDBWriter.{h,cpp}, CMakeLists.txt, src/util/filterdb.cpp,
src/util/createsubdb.cpp, examples/parquet_benchmark/*.sh. (`git status` also shows examples/compare_columns.py deleted
and util/regression modified; not from this work.)

## 1. Writer / convert2parquet (done)
- `mmseqs convert2parquet <DB> <out>`: one file, or `<out>/part-NNNNN.parquet` (~2 GB raw DB per file).
  Sort modes auto/linear/keyorder/spillsort (auto: keyorder if raw DB < MemAvailable/2, else spillsort).
- Options: --parquet-compression (0 none, 1 zstd default, 2 snappy, 3 gzip, 4 LZ4_RAW), --parquet-compression-level
  (default 3), --parquet-row-group-size (128M), --parquet-schema (required for generic DBs), --lookup-db.
- Sidecars: lookup `_lookup.pq` / `_<name>.lookup.pq` (key, accession, file_number); NEW: empty entries
  `_empty.pq` / `_<name>.empty.pq` (keys of entries without rows, needed by a reader).
- NEW: key-value metadata `mmseqs.dbtype`; NEW: max 20000 rows per page (cheap row-range reads).
- `--lookup-db` is not inferred from a result DB: `<resultDB>.lookup` does not exist, pass the query seq DB.
- convertalis --format-mode 5 uses the same machinery with its own file/column naming.
- Writer build note: CARQUET_ENABLE_BROTLI OFF. Upstream workflow scripts (easyproteomesearch/cluster.sh) fail
  shellcheck -> build with -DHAVE_SHELLCHECK=0 (build/Integrated_Native).

## 2. Query benchmarks TSV vs parquet (done)
- Scripts: examples/parquet_benchmark/run_parquet_benchmark.sh (+ setup_conda_env.sh). Variants: R data.table
  (fread all / needed columns), R arrow (load all / lazy), R duckdb, DuckDB CLI; queries A (best hit per query),
  B (hits of 100 queries, now data.table join), C (cluster size distribution, needs --cluster-db).
  Each variant runs each query in its own process (per-query peak memory). OOM/failed are recorded as results.
- Main finding: most of the speedup is the engine (DuckDB vs data.table); parquet adds 1.5-4x and ~40% less memory,
  most for selective queries and loading.
- Stress test 1-30 GiB TrEMBL search subsets (job 555745, super011, 3h11): at 30 GiB (634 M lines)
  A: fread 449 s vs R duckdb parquet 11.7 s (38x); B: 299 s vs 0.94 s (317x); memory 103 GB vs 25 GB;
  conversion createtsv 21.9 s vs convert2parquet 7.9 s. Per-size folders (size_*) were not copied locally.

## 3. Stress test scripts (ready, NOT yet run in the new version)
- run_stress_test.sh + run_stress_test.sbatch (`sbatch run_stress_test.sbatch /path/to/mmseqs`): copies TrEMBL to
  $SCRATCH, sizes 1 5 10 20 30 40 50 GiB, linclust per size (--min-seq-id 0.5 -c 0.8) for query C,
  convert2parquet with zstd AND uncompressed (queries use only zstd), results copied to --out after each size.
- Plots: stress_plot (row 1 time, row 2 peak memory; conversion, A, B, C), stress_filesize (TSV vs parquet zstd vs
  parquet uncompressed, alignments + clusters).
- Expect OOM for fread/arrow-full at 40-50 GiB (~170/145 GB). Estimated runtime 10-12 h.
- Open: enable `#SBATCH --exclusive`?; confirm linclust parameters; add GNU `time` to setup_conda_env.sh;
  copy updated run_parquet_benchmark.sh + run_stress_test.sh to the server.

## 4. ParquetDBReader (paused by the user on 2026-11-16)
- Design: (a) compatibility getData(id) rebuilds the exact DB text (alignment, cluster, prefilter, AA/NUC);
  (b) typed entryRows() + numericValue(); EntryReader wraps DBReader or parquet (magic PAR1 / part dir).
  Used in filterdb (text path + typed numeric-comparison path, MMSEQS_PARQUET_COMPAT=1 forces text) and
  createsubdb (parquet input). Output byte-identical to DB input in all tests (local + spark004, incl. empty entries,
  seqId boundary 0.845).
- Benchmark on spark004 (4.46 GiB alignment DB, 100.7 M lines, 20 threads, median of 3, noisy shared machine):

  | case | DB | pq baseline | fix1b | fix2 |
  |---|---|---|---|---|
  | filterdb none, typed, warm | 1.1-1.4 s | 15.2 s | 2.2 s | 3.0 s |
  | filterdb none, typed, cold | 2.0-2.5 s | 16.6 s | 4.6 s | 2.1 s |
  | filterdb selective, typed, cold | 2.3-3.2 s | 17.7 s | 6.2 s | 3.0 s |
  | filterdb permissive, text, cold | 2.4-3.6 s | 19.8 s | 8.3 s | 5.7 s |
  | createsubdb 1 %, warm | 0.06-0.18 s | 4.6 s | 9.2 s | 3.5 s |
  | createsubdb 10 %, warm | 0.3-0.9 s | 14.0 s | 26.9 s | 10.5 s |

  fix1 = OpenMP chunks sized to the decode window + no CRC check on read; fix1b = window rule "forward = full window"
  (fix1 entry-only windows made createsubdb 7-15x slower); fix2 = 20000 rows per page.
- Compressed-input round (DB via `mmseqs compress` = zstd 3 per entry, 1.18 GB in 20 files; parquet zstd 3,
  1.14 GB; raw DB 4.79 GB) finished on spark004 (label `compressed` in ~/mmseqs_parquet_bench/results/results.csv),
  results not yet read.
- Bottlenecks, least invasive first (1-3 done):
  1 overlapping per-thread windows (14x decode -> 1.5x) DONE; 2 CRC verify on read DONE; 3 1 MB pages DONE;
  4 windows not exactly aligned to chunks (1.5x instead of 1.0x): reader hands out row-sized chunks;
  5 sparse key access decodes 65k-row windows (~20x more than needed): prefetch API with the key list;
  6 open decodes the whole key column: write an entry index sidecar (key, part, first row, rows);
  7 copies of decoded values / string allocations: zero-copy views on carquet batch buffers;
  8 text path: snprintf("%.3E") + string appends: own evalue formatter, preallocated buffer;
  9 PLAIN encoding incl. repeated key: DELTA/RLE keys, dictionary for small ints (changes files);
  10 modules write text DBs: ParquetDBWriter as DBWriter output backend (large);
  11 typed record API (result_t rows) for all readers (large rework of DBReader + module loops);
  12 columnar native DB format (largest).
- Next when resumed: read the `compressed` results; implement 4 + 5; rerun quieter (8 threads, 5 runs).

## 5. Compression notes
- `createtsv --compressed 1` only compresses with `--db-output 1` (generic DB); otherwise plain TSV.
- DBWriter compression: zstd level 3 fixed, one frame per entry, entries < 60 bytes stay uncompressed,
  dbtype bit 31 marks compressed. Parquet compresses per column/page (better ratio).

## 6. spark004 (Brev Connect node, user2, shared account, other users active)
- aarch64, 20 cores, 121 GB RAM. Rust installed for user2 in ~/.cargo, ~/.rustup (user approved).
- ~/mmseqs_parquet_bench: MMseqs2 (source copy without .git; build in MMseqs2/build_arm), swissprot_db, data
  (res, resbig, aln5g*, aln5gb, aln5gb_c, aln5gb_p20k.pq, small*), results/results.csv, scripts bench.sh,
  prepare*.sh, correctness.sh, run_fix2.sh, run_compressed.sh, summarize.py, dbdump.py.
- Transfer: `brev copy` fails for Connect nodes, `brev exec` does not forward stdin; `scp`/`ssh spark004` work.
  Jobs started with setsid nohup survive SSH drops. Nothing was deleted there.
- Moving the Claude session to spark004 (tmux + /remote-control) was discussed, not done.

## 6b. Writer (convert2parquet) optimisations, 2026-11-16 (applied, uncommitted)
Test: 2 GiB alignment DB (2.15 GB, 119023 entries, 45.3 M lines, Swiss-Prot result subset), local laptop, 8 threads,
`convert2parquet --parquet-compression 0`; content checked with a DuckDB hash over all columns (identical in every
step), plus cluster, prefilter, alignment-with-backtrace and sequence DBs (identical before/after fix 4/6) and
ParquetDBReader/filterdb on a delta-encoded file (identical output).
- Cost split (single thread, 10 M lines, measured with a standalone carquet benchmark): text split 11 %, number
  parsing 13 %, std::map column lookups 13 %, carquet encode 48 % (bloom filters ~25 % of total, statistics ~10 %,
  CRC not measurable).
- Applied:
  - fix 2: `record_batch` std::map -> std::vector (2 map lookups per cell removed); freeRecordBatch clears it.
  - fix 3: statistics only on key/query/target/representative/member/score/seqId/evalue/bits/fident/pident/
    prefScore/taxid/accession (carquet_writer_set_column_statistics); no measurable gain locally.
  - fix 4: line/field ends via memchr in writeEntry (entry ends at its first \0 or len); ~2 % CPU.
  - fix 6: key column (key/query/representative, INT64) DELTA_BINARY_PACKED: query column 362 MB -> 2.2 MB,
    file 2567 -> 2207 MB (-14 %). Bug found on the way: the existing "all columns PLAIN" loop in openWriter must run
    before the per-column override.
- Results (laptop is noisy, load up to 10.9 on 8 logical cores):
  sequential, 3 runs: baseline 9.62 s / 64.7 s CPU, fix 2 8.02 s / 51.5 s, fix 3 9.17 s / 53.0 s,
  fix 4 7.14 s / 39.6 s (noise), fix 6 6.07 s / 37.2 s; createtsv 9.94 s / 49.2 s.
  interleaved, 5 rounds: createtsv 15.34 s / 67.4 s, fix 3 9.21 s / 53.4 s, fix 4 9.22 s / 52.4 s,
  fix 6 8.57 s / 51.3 s. Baseline/fix 2 binaries were not kept, so only fix 3 -> 6 are compared interleaved.
  - fix 1 (applied later the same day): bloom filters only for the `query` and `target` columns. carquet detail:
    the plain per-column setter (carquet_writer_set_column_bloom_filter) is IGNORED while the global flag is on;
    only carquet_writer_set_column_bloom_filter_options switches to per-column control, so it is used for every
    column (ndv 0 / fpp 0 = defaults). Cluster files (representative/member) get no bloom filter.
- spark004 rerun (2026-11-16, quiet machine, 20 threads, 2 GiB alignment DB cut from resbig: 2.15 GB, 257847
  entries, 45.1 M lines; 3 interleaved rounds; stage binaries built from sources with the fixes undone step by step,
  ~/mmseqs_parquet_writer_bench, results/writer_results.csv):

  | variant | wall median | CPU | peak MB | size MB | vs baseline | vs createtsv |
  |---|---|---|---|---|---|---|
  | createtsv | 1.72 s | 11.9 s | 2798 | 2473 | - | 1.00x |
  | baseline | 1.44 s | 9.2 s | 4289 | 2560 | 1.00x | 1.19x |
  | fix 2 (vector) | 1.29 s | 7.7 s | 4226 | 2560 | 1.12x | 1.33x |
  | fix 3 (+ stats) | 1.28 s | 7.6 s | 4262 | 2560 | 1.12x | 1.34x |
  | fix 4 (+ memchr) | 1.34 s | 8.1 s | 4258 | 2560 | 1.07x | 1.28x |
  | fix 6 (+ delta key) | 1.32 s | 7.9 s | 4036 | 2201 | 1.09x | 1.30x |
  | fix 1 (+ bloom query/target) | 1.20 s | 7.0 s | 4011 | 2176 | 1.20x | 1.43x |

  Runs are ~1-2 s (input in page cache), so fixed costs weigh in; fix 3/4 differences are within noise.
  Only 8 of 20 cores are used (2 GB -> 8 files, see fix 5).
- spark004 rerun with 20 GiB (2026-11-16): resbig stacked 5x with concatdbs (alignments, sequences, headers with
  the same key shift), cut to 20 GiB: 21.5 GB, 2576166 entries, 451.4 M lines; 20 threads, 3 interleaved rounds,
  results/writer_results_20g.csv. createtsv --compressed only works with --db-output 1 (else plain TSV), so the
  compressed TSV variant is `createtsv --db-output 1 --compressed 1` (zstd 3 per entry).

  | variant | mode | wall median | CPU | size MB | vs baseline | vs createtsv |
  |---|---|---|---|---|---|---|
  | createtsv | plain | 24.91 s | 109.7 s | 24726 | - | 1.00x |
  | createtsv --db-output --compressed | zstd 3 | 9.48 s | 170.0 s | 6044 | - | 2.63x |
  | baseline | raw | 8.27 s | 107.2 s | 25547 | 1.00x | 3.01x |
  | fix 2 | raw | 6.27 s | 93.5 s | 25547 | 1.32x | 3.97x |
  | fix 3 | raw | 6.63 s | 85.2 s | 25547 | 1.25x | 3.76x |
  | fix 4 | raw | 8.61 s | 90.6 s | 25547 | 0.96x | 2.89x |
  | fix 6 | raw | 6.38 s | 89.9 s | 21961 | 1.30x | 3.90x |
  | fix 1 | raw | 6.11 s | 77.2 s | 21745 | 1.35x | 4.08x |
  | baseline | zstd 3 | 9.00 s | 144.3 s | 5040 | 1.00x | 2.77x |
  | fix 2 | zstd 3 | 7.65 s | 122.2 s | 5040 | 1.18x | 3.26x |
  | fix 3 | zstd 3 | 7.53 s | 122.2 s | 5040 | 1.20x | 3.31x |
  | fix 4 | zstd 3 | 8.02 s | 129.1 s | 5040 | 1.12x | 3.11x |
  | fix 6 | zstd 3 | 8.26 s | 130.1 s | 5038 | 1.09x | 3.02x |
  | fix 1 | zstd 3 | 7.29 s | 116.2 s | 4822 | 1.23x | 3.42x |

  Findings: all fixes: raw 1.35x (CPU -28 %), zstd 1.23x (CPU -19 %); raw file -15 %, zstd file -4 %.
  Parquet zstd vs createtsv --db-output --compressed: 1.30x faster, 32 % less CPU, 20 % smaller.
  fix 4 (memchr) is SLOWER on spark004/ARM in both modes (raw 6.63 -> 8.61 s, zstd 7.53 -> 8.02 s), likely many
  short memchr calls on ~5-byte fields. REVERTED on 2026-11-16 (byte loop restored with a comment; content of
  seq/cluster/prefilter/alignment+backtrace parquet identical; fixes 1, 2, 3, 6 kept). Peak memory ~21-24 GB is mostly the mmapped input.
- Conversion stress test on spark004 (2026-11-16, writer with fixes 1/2/3/6, fix 4 reverted): Swiss-Prot result
  stacked 8x without copying (make_stacked_db.sh: index + symlinks, 35.7 GiB), sizes 1/5/10/20/30 GiB, 3 runs,
  20 threads, run_stress_test.sh --result-db/--seq-db --conversion-only --shared-output; plotted locally with
  --summary-only (no R on spark004). Results in ~/mmseqs_parquet_writer_bench/stress_conv/out on spark004.
  At 30 GiB (677 M lines): createtsv 43.8 s / 34.5 GiB output; parquet zstd 11.5 s (3.8x faster) / 6.8 GiB
  (5x smaller); parquet uncompressed 15.5 s / 30.5 GiB. Peak memory 27-32 GiB for all (mostly the mmapped
  input); uncompressed parquet needs more at 5-10 GiB (10.5 / 17.6 GiB vs 6.2 / 11.8 GiB zstd). createtsv has a
  large spread at 20-30 GiB (sd 7-9 s, writing). Cluster conversions < 1 s at all sizes.
- New script options (2026-11-16): run_parquet_benchmark.sh --conversion-only; run_stress_test.sh --result-db,
  --seq-db, --conversion-only, --shared-output, --summary-only; make_stacked_db.sh.
- carquet 0.7.2 -> 0.8.0 (subtree updated by the user, 2026-11-17): stability/regression test on spark004
  (examples/parquet_benchmark/carquet_regression_test.sh, ~/mmseqs_carquet080_test, report.txt): 0 failed checks.
  Release + LTO build links (archiver gcc-ar-13: the 0.7.3 fix for the cluster link error). Old (0.7.2) vs new
  binary on small alignment with 764 empty entries, 5 GiB alignment (113 M lines, zstd and uncompressed), cluster,
  prefilter, alignment with backtrace, sequences, convertalis --format-mode 5: carquet validate OK, content
  identical (CSV export hash), same rows, metadata and empty-entry sidecar present. Stability: 5 repeats, 1/4/20
  threads, keyorder and spillsort give identical content, no warnings. Performance 5 GiB (3 rounds): zstd 2.03 ->
  1.88 s (1.08x), uncompressed 1.76 -> 1.79 s (0.98x, noise). Files 0-2.4 % larger with 0.8.0 (zstd 5 GiB
  1.315 -> 1.347 GB), consistent with the max_rows_per_page fix (more pages); page counts not recorded.
- Full stress test with carquet 0.8.0 on spark004 (2026-11-17, 1/5/10/20 GiB, conversions + queries A/B/C, conda env
  rbench without R duckdb -> R duckdb variants "failed" as agreed): results and plots copied to
  examples/parquet_benchmark/stress_full_carquet080/. At 20 GiB: conversion createtsv 23.1 s / 23.0 GiB, parquet
  zstd 7.3 s / 4.8 GiB, uncompressed 8.5 s / 20.3 GiB; A: DuckDB CLI parquet 20.3 s vs fread 117.8 s; C identical
  rows in all variants. CAVEAT: A/B row counts differ between TSV and parquet variants because of the stacked test
  data (the 8 copies share header names, the lookup has unique _cN accessions; B ids come from TSV names, so parquet B
  only matches copy 0) -> B timings TSV vs parquet are not comparable here. Fix for future stacked tests: suffix the
  headers per copy too (needs a copy of the header DB).
- Reader-test data on spark004 moved (not deleted) to ~/mmseqs_parquet_reader_tests; writer benchmark in
  ~/mmseqs_parquet_writer_bench (MMseqs2 source with all fixes, build_arm, bin/mmseqs_<stage>, data/aln2g).
- Not applied (user decision): fix 5 more files/workers on many-core machines (2 GB DB -> max 8 workers today):
  REJECTED 2026-11-16, the file layout must not change (no more, smaller files); only parallelism inside a file
  (parse in parallel, one writer per file) would be acceptable if ever revisited. Fix 7 string columns without
  temporary std::string: not applied.
- To do: quieter rerun (e.g. 4 threads, more rounds) for reliable numbers; consider DELTA for `target` and
  dictionary for small integer columns.

## 7. Other open todos / ideas
- cA3M (dbtype 10): msa column must be BLOB (BYTE_ARRAY without STRING logical type), split the shared CA3M/MSA case
  in dbTypeToFields; currently cA3M (and INDEX_DB) are marked unsupported in convert2parquet.
- Linclust paper step 6 (Pfam + GO functional consistency): original evaluation code not public (linclust-analysis,
  uniclust-pipeline lack it). Needs UniProt 2016_03 (sequences + .dat with DR GO / DR Pfam), go.obo 2016,
  clusterings; scripts 1-6 (extract annotations, cluster + convert2parquet, GO graph closure, IC, scores,
  aggregate/plot). Open decisions: Pfam set/multiset/ordered, EXP_F evidence codes, unannotated members and
  singletons, is_a vs part_of. Suggested pilot: Swiss-Prot 2016_03, Pfam score only.
- Possible commit of the parquet reader/writer changes (user decides; I may not run git add/commit/push).

**Why:** the user asked to collect all results and ideas in this todo file (2026-11-16); this memory file is the
only todo file, written outside the working directory with the user's explicit permission.
**How to apply:** use as the starting point when the user resumes any of these topics; verify file paths and code
before relying on details, numbers are from noisy shared-machine runs.
