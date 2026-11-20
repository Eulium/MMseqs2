# Parquet index for MMseqs2: benchmarks, findings, fixes

Working log of the work on branch `parquet_index` (started 2026-11-18, machine spark004: aarch64, 20 cores, 121 GB RAM,
shared with other users). Newest results are added at the end of each section. Numbers marked *exploratory* were
measured while something else of ours was running (download, compile) and are only used to find where time goes.

## 1. What is being built

MMseqs2 loads the `.index` of every database into memory (24 bytes per entry) and loads it again whenever a database
is written (merge of the per-thread parts, sort by key). For very large clustering inputs this no longer fits.
The index is therefore kept in a parquet file `<db>.index.parquet` (key, offset, length; sorted by key; delta encoded)
and read page by page.

State of the code (all uncommitted in `~/MMseqs2` unless noted):

| Part | State |
|---|---|
| Lazy reader for NOSORT, SORT_BY_ID, key-ordered LINEAR_ACCCESS | committed (`92d4ef92`) |
| Lazy reader for SORT_BY_LENGTH, unsorted LINEAR_ACCCESS, SHUFFLE (permutation of 8 B/entry in memory) | written, tested |
| Streaming index merge in `DBWriter` (k-way merge of the per-thread parts, no full load) | written, tested |
| Parquet-only mode: no text `.index` is written (`MMSEQS_WRITE_LAZY_INDEX=parquet,notext`) | written, tested |
| Per-column partial page decoding, configurable page size | written, not yet measured |

Switches (environment variables, prototype):

- `MMSEQS_INDEX_BACKEND=parquet|parquet-custom|binary`: reader backend when a text index exists next to the file
- `MMSEQS_WRITE_LAZY_INDEX=parquet[,binary][,notext]`: what `DBWriter` writes; `notext` leaves the text index out
- `MMSEQS_INDEX_DEBUG=1`: every index open is reported (kind, entries, access mode)
- `MMSEQS_LAZY_INDEX_PAGE_ROWS=<n>`: rows per page of written files (benchmark only, default 1024)

A database without text index is read through its parquet file automatically (custom decoder), in every access mode.

## 2. How the index is used in the clustering workflows (code trace)

`linclust` in this fork runs `kmermatcher` -> `align2clust` -> `mergeclusters` -> `createsubdb` -> second round.

| Step | DB and open mode | Access pattern |
|---|---|---|
| `kmermatcher` | sequence DB, NOSORT | one pass in key order |
| `align2clust` | sequence DB, SORT_BY_LENGTH | random: per cluster the representative, then every k-mer hit by key |
| `align2clust` | k-mer result DB, LINEAR_ACCCESS | random `getId(key)` (loop runs in order of hit-list size) |
| `mergeclusters` | sequence DB, NOSORT (index only) | random `getId(member)` for every member; one pass at the end |
| `mergeclusters` | cluster DBs, LINEAR_ACCCESS | sequential |
| `createsubdb --subdb-mode 1` | sequence DB, NOSORT | `getId` per representative, in key order |
| `DBWriter::close` (every output) | per-thread index parts, merged index | was: load all, shift offsets, load all again, sort |

Other per-sequence memory that stays: `align2clust` 36 B/sequence (assignment, hit-list sizes, final pairs),
`kmermatcher` 3-5 B/key, `mergeclusters` a `std::list` per sequence. Default builds use 32-bit ids (limit 4.29e9).

## 3. Correctness tests

| Test | Result |
|---|---|
| `lazy_index_test.sh` (createsubdb, view, convert2fasta, search + convertalis, linclust, cluster; 3 backends; 1.29M entries) | 0 failed checks, last run with the streaming writer |
| `test_lazyindex` (`src/test/TestLazyIndex.cpp`): all 10 access modes + `getIndex`, every id and key, 20 threads; text reference against parquet-only and against text + parquet/custom/binary; sequence DB and a result DB whose offsets are not in key order | all modes the same (2 DBs x 4 configurations) |
| Old write path against streaming write path: `createsubdb` (random keys, both sub-DB modes), `filterdb` 20 threads | single-thread writers byte-identical; multi-thread writers same content per key (their index bytes also differ between two runs of the old path: thread scheduling) |
| Parquet-only mode: linclust, cluster, search + convertalis, createtsv, convert2fasta, mvdb, rmdb; input DB also parquet-only | same clusters/hits as the text build; 61 index opens, all lazy, none loaded |

## 4. Benchmarks

### 4.1 Lookup benchmark, UniRef50 (41.2M entries), committed build, 3 rounds

All workloads were disturbed by other users (2-7 cores). Index on disk: text 979 MB, parquet 64 MB, binary 660 MB.

| Workload | text | parquet (carquet) | parquet (custom) | binary |
|---|---|---|---|---|
| view, 1 key | 0.74 s | 0.01 s | 0.01 s | 0.00 s |
| createsubdb, 100k random keys | 1.45 s | 0.71 s | 0.20 s | 0.04 s |
| createsubdb, 10M random keys | 8.57 s | 69.18 s | 20.53 s | 5.97 s |
| createsubdb, 10M sorted keys | 2.58 s | 1.35 s | 1.22 s | 1.20 s |
| createsubdb, all keys sorted | 6.02 s | 5.20 s | 4.84 s | 4.83 s |
| convertalis, 20 threads | 1.73 s | 0.27 s | 0.26 s | 0.22 s |
| allocated memory, most workloads | 944 MB | 0-8 MB | 8 MB | 0 MB |

Per random key: text 0.70 us, carquet 6.90 us, custom 2.04 us, binary 0.58 us. Per sorted key all 0.12-0.13 us.

### 4.2 linclust, 20M UniRef50 sequences, 20 threads, one round each (whole workflow)

| Build | Variant | Wall | Allocated | Large text-index loads | Other load |
|---|---|---|---|---|---|
| committed | text | 283.6 s | 45,944 MB | 15 | 1.2 cores |
| committed | parquet-custom, all DBs | 329.8 s | 45,493 MB | 13 | 1.9 |
| + permutation modes | text | 272.4 s | 45,944 MB | 15 | 0.7 |
| + permutation modes | parquet-custom, all DBs | 318.9 s | 45,491 MB | 9 | 1.7 |
| + streaming writer | text | 329.8 s | 45,944 MB | 15 | 2.3 (download + compiles running) |
| + streaming writer | parquet-custom, all DBs | 310.5 s | 45,493 MB | 0 | 1.8 (same) |

Same clusters in every run. The peak memory is not set by the index at this size (457 MB of 46 GB); it does not move.
These single rounds with changing other load cannot rank the variants; section 4.4 repeats them properly.

### 4.3 cluster (cascaded), 500k sequences, committed build, one round

text 237.9 s, parquet 246.2 s, parquet-custom 246.1 s, binary 220.9 s, parquet-custom all DBs 200.9 s; 2.28-2.29 GB
allocated in all; same clusters. Differences follow the other load (0.4-1.5 cores), not the backend.

### 4.4 Per-module benchmark (`lazy_index_steps_bench.sh`)

Runs the workflow through `mmseqs_steps.py`, which records wall, CPU and peak RssAnon of every module call.
Variants: `text` (MMseqs2 as it is), `both` (text + parquet written, read through parquet), `notext` (parquet only).

#### 4.4.1 Exploratory profile, 20M UniRef50, parquet-only build before any reader fix (1 round, download running)

Whole workflow: text 277.5 s, both 282.3 s, notext 278.5 s; other load 0.8-0.9 cores; same clusters.
Index opens: text variant 72 text loads (15 of them of at least 10M entries); both/notext 12 lazy opens, 0 loads.

| Module | text: wall / CPU / peak anon | notext: wall / CPU / peak anon | Change |
|---|---|---|---|
| kmermatcher (1) | 41.9 s / 685 s / 19,149 MB | 41.8 s / 682 s / 18,696 MB | level |
| align2clust (1) | 79.1 s / 1,309 s / 2,834 MB | 79.9 s / 1,392 s / 1,926 MB | CPU +6 %, memory -908 MB |
| createsubdb | 2.4 s / 2 s / 458 MB | 4.0 s / 4 s / 5 MB | +1.6 s, memory -453 MB |
| kmermatcher (2) | 59.9 s / 956 s / 45,945 MB | 61.4 s / 961 s / 45,490 MB | level |
| align2clust (2) | 82.5 s / 1,371 s / 4,213 MB | 84.7 s / 1,485 s / 2,424 MB | CPU +8 %, memory -1,789 MB |
| mergeclusters | 11.6 s / 103 s / 2,269 MB | 6.5 s / 100 s / 1,257 MB | -5.1 s, memory -1,012 MB |

The peak of the whole workflow is `kmermatcher` (2) with 45.9 GB: its k-mer buffers, sized from the available memory.
That is why the whole-workflow peak never moved in section 4.2. Per module the index saving is visible.

#### 4.4.2 Timed run, 20M UniRef50, 20 threads, 3 rounds, nothing else of ours running (2026-11-19 01:41-02:38)

Variants: `text`; `notext_v1` (parquet only, reader before the lookup fix of 5.1, 1,024-row pages);
`notext` (parquet only, partial decoding, 256-row pages).

| Run | text | notext_v1 | notext |
|---|---|---|---|
| Round 1: wall s (other processes, cores) | 269.9 (0.3) | 303.5 (1.5) | 290.8 (1.3) |
| Round 2 | 304.6 (1.4) | 299.4 (1.3) | 286.0 (1.1) |
| Round 3 | 304.6 (1.4) | 302.0 (1.0) | 286.9 (1.3) |
| CPU s (median) | 4,623 | 4,823 | 4,586 |
| Peak RssAnon | 45,951 MB | 45,499 MB | 45,501 MB |
| Index opens: lazy / text loads (large) | 0 / 72 (15) | 12 / 0 | 12 / 0 |

**This run is not a clean comparison.** After the first run (text, 269.9 s) another user started a single-threaded
job that ran until the end. With 20 threads on 20 cores one busy outside core costs about 12 % (text: 269.9 s
alone, 304.6 s twice with the outside job). Only the eight runs with the outside job compare with each other:
text 304.6 s, notext_v1 about 302 s, notext about 287 s. The run is repeated with fewer threads than cores (4.4.3).

Per module (median of 3 rounds; text median is from the two disturbed rounds):

| Module | text: wall / CPU / peak anon | notext_v1 | notext |
|---|---|---|---|
| kmermatcher (1) | 45.8 s / 705 s / 19,148 MB | 46.5 s / 705 s / 18,694 MB | 41.4 s / 678 s / 18,699 MB |
| align2clust (1) | 83.7 s / 1,348 s / 2,837 MB | 84.4 s / 1,436 s / 1,926 MB | 80.9 s / 1,364 s / 1,930 MB |
| createsubdb | 2.7 s / 3 s / 458 MB | 4.5 s / 4 s / 5 MB | 4.8 s / 5 s / 8 MB |
| kmermatcher (2) | 69.2 s / 1,027 s / 45,944 MB | 67.8 s / 1,029 s / 45,491 MB | 68.5 s / 1,013 s / 45,494 MB |
| align2clust (2) | 88.2 s / 1,420 s / 4,197 MB | 90.2 s / 1,544 s / 2,401 MB | 84.0 s / 1,438 s / 2,414 MB |
| mergeclusters | 13.4 s / 108 s / 2,280 MB | 7.3 s / 101 s / 1,257 MB | 7.4 s / 101 s / 1,258 MB |

What holds despite the disturbance (CPU time and memory are much less sensitive to it than wall time):

- The lookup fix removes the extra CPU of `align2clust`: +6.5 % / +8.7 % CPU before the fix, +1.2 % / +1.3 % after.
- Memory per module is the same before and after the fix: `align2clust` -0.9 / -1.8 GB, `mergeclusters` -1.0 GB,
  `createsubdb` -450 MB against text. The workflow peak (`kmermatcher` 2) is 450 MB lower: one index array of 20M.
- `mergeclusters` is 6 s faster (no index loads), `createsubdb` 2 s slower (5.2).

#### 4.4.3 Clean timed run, 20M UniRef50, 16 of 20 threads, 3 rounds (2026-11-19 04:41-05:37)

Four cores left free for other users' processes; outside load 0.1-0.6 cores, no swapping, at least 70 GB of memory
available throughout. All nine runs give the same clusters as text.

| | text | notext_v1 (before lookup fix) | notext (current) |
|---|---|---|---|
| Wall, median (min-max) | 300.4 s (300.1-300.7) | 302.8 s (300.9-305.0) | 294.6 s (294.3-298.3) |
| Against text | - | +0.8 % | -1.9 % |
| CPU | 4,122 s | 4,262 s | 4,128 s |
| Peak RssAnon (workflow) | 45,942 MB | 45,491 MB | 45,490 MB |
| Index opens: lazy / text loads (large) | 0 / 64 (15) | 12 / 0 | 12 / 0 |

| Module | text: wall / CPU / peak anon | notext_v1 | notext |
|---|---|---|---|
| kmermatcher (1) | 46.7 s / 642 s / 19,138 MB | 47.4 s / 641 s / 18,683 MB | 47.6 s / 641 s / 18,686 MB |
| align2clust (1) | 84.5 s / 1,211 s / 2,788 MB | 85.6 s / 1,281 s / 1,876 MB | 82.2 s / 1,222 s / 1,881 MB |
| createsubdb | 2.4 s / 2 s / 458 MB | 4.0 s / 4 s / 5 MB | 4.5 s / 4 s / 8 MB |
| kmermatcher (2) | 66.4 s / 890 s / 45,934 MB | 68.2 s / 893 s / 45,484 MB | 68.3 s / 894 s / 45,483 MB |
| align2clust (2) | 88.7 s / 1,270 s / 4,084 MB | 92.0 s / 1,372 s / 2,233 MB | 86.0 s / 1,280 s / 2,265 MB |
| mergeclusters | 11.5 s / 108 s / 2,269 MB | 5.1 s / 69 s / 1,257 MB | 6.4 s / 93 s / 1,258 MB |

Reading:

- Whole workflow: parquet only is 2 % faster than text at this size (6 s of 300 s), and was 1 % slower before
  the lookup fix.
- `align2clust`: 2-3 s faster than text after the fix (CPU +1 %; it no longer loads three index arrays), 1-3 s
  slower before it (CPU +6 % / +8 %). Memory -0.9 GB and -1.8 GB.
- `mergeclusters`: 5 s faster, -1.0 GB.
- `kmermatcher`: 1-2 s slower (2-3 %) with the same CPU time. Not analysed yet; candidates are the result index
  merge (text part sorted in memory, then written as parquet) and reading keys and lengths page by page.
- `createsubdb`: 2 s slower; fixed after this run (5.2).

**One run gave a different clustering:** `notext_v1` round 3 has 19,957,551 clusters instead of 19,957,552. The other
eight runs are identical to text. Cause: a known race in `kmermatcher`, not the index (6.1).

### 4.5 UniRef100 (release 2026_04): 227,086,028 sequences

DB: 94.6 GB sequences, 28.0 GB headers; text index 5.76 GB (5.2 GB as array per load), parquet index 0.89 GB;
`createdb` 17 min, 33.6 GB peak RSS.

#### 4.5.1 All 227M sequences, text index, 16 threads, `--split-memory-limit 80G`: stopped after 6 h 22 min

| Step | Result |
|---|---|
| kmermatcher (1) | 1,420 s wall, 18,700 s CPU, 76.4 GB peak RssAnon; result 92 GB in 16 files |
| align2clust (1) | stopped by hand after 5 h 58 min, no output written yet |

`align2clust` was reading from disk at about 5 GB/s the whole time (100 TB in total), at 2-4 of 16 cores busy, in
I/O wait. 7 GB of the process were swapped out. The memory guard did not fire: page cache counts as available
memory.

Correction (2026-11-20): the first explanation here was that the sequence DB did not fit into the page cache. That
was not the cause. The 120M run (4.5.2) showed the same reading with 40 GB of memory free; the cause is a wrong
read hint in `align2clust` (5.4). Whether all 227M sequences run in useful time with the fix has not been tested:
the sequence DB (94.6 GB) and the k-mer result (92 GB) are both read in random order and together do not fit into
121 GB, so some disk reading will remain.

- The build for the `createsubdb` fix (5.2) and its test ran during the first minutes of this run.

#### 4.5.2 120M UniRef100 sequences, 16 threads, `--split-memory-limit 70G`

Subset: 49.9 GB sequences, text index 3.0 GB (2.7 GB as array per load), parquet index 0.5 GB. k-mer result of
round one: 40.1 GB in 16 files.

Before the read hint fix (5.4), text index, one round (2026-11-19 12:06 - 2026-11-20 01:18): 47,282 s (13 h 8 min),
of which `align2clust` 31,800 s + 13,917 s at about 2 of 16 cores. The parquet-only round was stopped after
25 min of `align2clust` with the same behaviour.

With the fix (binary `mmseqs_v4_advice`), round 1 (2026-11-20 01:57 - 06:34; a job of another user used 1.4 cores):

| Module | text (s) | parquet-only (s) | text RssAnon (MB) | parquet-only RssAnon (MB) |
|---|---|---|---|---|
| kmermatcher | 699.5 | 652.6 | 68,344 | 68,255 |
| align2clust | 3,851.5 | 3,994.6 | 14,726 | 9,691 |
| createsubdb | 11.5 | 13.2 | 2,747 | 53 |
| kmermatcher (2) | 666.5 | 640.8 | 68,180 | 68,241 |
| align2clust (2) | 2,741.0 | 2,780.8 | 18,022 | 9,948 |
| mergeclusters | 50.9 | 27.6 | 11,733 | 7,077 |
| **linclust total** | **8,021** | **8,110 (+1.1 %)** | peak 68,346 | peak 68,263 |

- The fix makes the text run 5.9 times faster (47,282 s to 8,021 s) with the same clustering (72,128,220
  clusters, checksum `2b82388f6a78` before and after).
- Parquet-only against text: the total is the same within what one round on a shared machine can show (+1.1 %).
  `align2clust` is 1.5-3.7 % slower in wall time and needs 2.6-5.5 % more CPU (key lookups through pages, as in
  5.1); `kmermatcher` was faster in this round (the other user's job makes single-module differences of a few
  percent unreliable); `mergeclusters` takes half the time.
- Memory: `align2clust` needs 5.0 GB and 8.1 GB less, `mergeclusters` 4.7 GB less, `createsubdb` 2.7 GB less. The
  peak of the workflow is `kmermatcher` and does not change. The text run swapped out 97.8 GB over the run, the
  parquet-only run 9.2 GB.
- **Open problem:** both parquet-only rounds gave 72,128,221 clusters (checksum `3bcd9dbe2537` both times), text
  72,128,220 (`2b82388f6a78`, before and after the fix). Two identical results are not the `kmermatcher` race
  (6.1): at 120M parquet-only and text differ by one cluster, reproducibly. Not found at 20M. Cause not known yet.

Round 2: parquet-only 8,009 s (same result as round 1); the text round was stopped by hand after about 40 min for
the comparison with `master` (5.4).

## 5. Bottlenecks and fixes

### 5.1 Random lookups decode whole pages (found in `align2clust`)

Finding: with the parquet index `align2clust` needs 6-8 % more CPU. It looks up every k-mer hit by key
(`getId`, `getData`, `getSeqLen`); each of these decoded all three columns of a 1,024-row page.

Fix A, partial decoding (`LazyIndex::decodeColumn`): a lookup decodes only the column it needs and only up to the
row it needs; a key search stops at the first key that is not smaller. `DBReader` asks for one member
(`getKey`/`getOffset`/`getLength`) instead of the whole entry.

Fix B, page size: 256 rows per page instead of 1,024 (row groups stay at 16,384 rows). Files written with either
size are read by both reader versions.

Lookup timing, `createsubdb` with 10M random keys on UniRef50 (41.2M entries), one run each, download running:

| Reader | Rows per page | Time | Parquet index size |
|---|---|---|---|
| text index | - | 8.29 s | (979 MB text) |
| custom decoder before the fix | 1,024 | 19.70 s | 64.5 MB |
| custom decoder, fix A | 1,024 | 13.67 s | 64.5 MB |
| custom decoder, fix A + B | 256 | 8.21 s | 92.8 MB |
| custom decoder, fix A | 128 | 8.10 s | 130.5 MB |
| custom decoder, fix A | 64 | 7.95 s | 200.4 MB |
| carquet reader, fix A | 1,024 | 68.18 s | 64.5 MB |

All seven outputs are identical. 256 rows is the default now: random lookups cost the same as with the text index,
the file grows from 1.6 to 2.3 bytes per entry, the page directory in memory from 8 to 32 bytes per 256 entries
(about 5 MB for 41M entries with the custom decoder). Smaller pages gain almost nothing more.
The carquet reader was not measured with 256-row pages; its cost grows with the pages per row group.

Tests after the fix: access-mode test (old 1,024-row and new 256-row files, custom and carquet reader), the
correctness suite and the parquet-only workflows all pass with identical results.

### 5.2 `createsubdb` with a parquet-only key DB wrote the keys to a temporary text file

Finding: `createsubdb <clustering DB> <sequence DB>` takes the keys from the index of the first DB. Without a text
index the keys were printed to a temporary file and parsed again: 4.0-4.5 s instead of 2.4 s at 20M.

Fix: the keys are read from the lazy index directly (`LazyIndex::getKey`), the text form is only built for the
lookup-by-name mode and for warnings. Key files and DBs with text index go the old way.

| Build | createsubdb, 19.96M keys of 20M (2 runs) | Result checksum |
|---|---|---|
| before | 4.40 s, 4.37 s | 7708520a0a4c304b |
| after | 3.19 s, 3.25 s | 7708520a0a4c304b (same, and same as in the linclust runs) |

A text key file gives byte-identical output with both builds. The remaining 0.8 s against text is one key lookup
per entry in the sequence DB and writing the parquet index. Binary: `~/mmseqs_lazy_index/bin/mmseqs_v3_subdb`
(the UniRef100 run below was started before this fix and uses the build without it).

### 5.3 Benchmark method: threads and memory on a shared machine

- With 20 threads on 20 cores a single outside process slows a run by about 12 % (4.4.2). The benchmarks now use
  16 threads; the spread between rounds dropped to 0.2-1.4 %.
- The benchmark script waits for free memory before a run, stops and repeats a run if memory gets short, and
  records the lowest available memory and the swapped MB of every run (after a queued job was stopped by a memory
  shortage caused by other users' processes).

### 5.4 `align2clust` read its hit lists with the sequential read hint (the cause of the hours of disk reading)

Not an index problem, but it hid every index effect at large scale.

Finding (120M sequences, sampled once a minute during `align2clust`): 6.2 GB/s read from disk, all 16 worker
threads waiting for a page (`folio_wait_bit_common`), 1 core busy, 40 GB of memory free, 1.6 million pages per
second read that had been in the cache before (`workingset_refault_file`), 2,000 major faults per second, so
about 3 MB read per fault. `/proc/<pid>/smaps` showed the flag `sr` (sequential read) on the 16 files of the
k-mer result.

Cause: `align2clust` opens the k-mer result with `LINEAR_ACCCESS`, and `DBReader::open` sets
`POSIX_MADV_SEQUENTIAL` on the data for this mode. But the module reads the hit lists by key, in the order of
the clustering (set cover: first one pass in order of sequence length to count the hits, then by list size), which
is random in the file. With the sequential hint each page fault reads a large window ahead, and the kernel does
not count such pages as used and drops them again, even with memory free. Below about 20M sequences everything
is cached so fast that it does not show.

Fix: `DBReader::setNormalAdvice()` (new) takes the hint back; `align2clust` calls it after opening the k-mer
result and the two filter DBs of round two. The open mode itself is unchanged, so ids and output are the same.

Effect at 120M (4.5.2): `align2clust` 31,800 s to 3,851 s and 13,917 s to 2,741 s, linclust 5.9 times faster,
16 of 16 cores busy, disk reads 0.25-0.5 GB/s and falling while the cache fills; same clustering checksum.

Check against `master` (`d401e78c`, own build), 20M UniRef50, text index, 16 threads, 2 rounds each
(2026-11-20 09:36-10:14, other users' processes used 1.0-2.4 cores, so times are only rough):

| Build | linclust wall (s) | align2clust 1 / 2 (s) | Clusters, checksum |
|---|---|---|---|
| this branch with the hint change | 312.6, 350.9 | 85.2, 92.3 / 92.8, 103.3 | 19,957,552, `8a66e0ba3b62` |
| `master` | 376.1, 344.1 | 104.4, 99.1 / 104.9, 92.8 | 19,957,552, `8a66e0ba3b62` |
| this branch without the hint change | 338.4, 319.0 | 92.1, 87.8 / 96.7, 91.9 | 19,957,552, `8a66e0ba3b62` |

Same output in all six runs. No time difference beyond the noise: at 20M the k-mer result is in the page cache
and the hint has no visible cost. `master` has the same `align2clust` and the same hint code; that it shows the
same disk reading at 120M follows from the code but was not run.

Not tested: `POSIX_MADV_RANDOM` instead of normal (reads only the page asked for; better when the data does not
fit into memory, slower to fill a cold cache); other modules that open a result with `LINEAR_ACCCESS` and read
it by key (not searched for).

## 6. Known problems and limits

- All of UniRef100 (227M) has not been clustered: the one attempt was stopped (4.5.1) before the read hint fix (5.4).
- Shell workflows that read `.index` with awk break in parquet-only mode: `linclust.sh` (clusthash branch, old
  `linclust` path), `cascaded_clustering.sh` (reassign branch), `blastp.sh` (sensitivity steps, one prefilter branch).
  The default paths tested here do not use them.
- A parquet-only DB cannot be read by upstream MMseqs2, Foldseek or scripts; there is no parquet-to-text converter yet.
- The parquet footer is written last: a crash during the index merge leaves no usable index.
- Callers that need the whole index array (`getIndex`, sorting by weights, split merge in `prefilter`,
  `SORT_BY_LINE`, `SORT_BY_OFFSET`) load it completely from the parquet file (24 B/entry, as before).
- `SORT_BY_LINE` on a parquet-only DB gives key order, not the order in which the entries were written.
- An unsorted per-thread index part is sorted in memory as one piece in the streaming merge.
- Merging more than 16 parquet-only DBs at once (split/MPI results) thrashes the per-thread page cache.
- 14 modules had their "remove the index of a flat output file" call replaced; only `convertalis` and `createtsv`
  were run.
- `createRenumberedDB` and `sortDatafileByIdOrder` still load the index.
- `mvdb`/`rmdb`/copy/link did not handle the parquet file in the committed build: a moved DB silently fell back to
  the text index (fixed in the working tree).

### 6.1 One run with a different clustering: known race in `kmermatcher`, not the index

Observed once: timed run 4.4.2, `notext_v1` round 3 (parquet only, reader *before* the lookup fix, 20 threads, another
user's job on one core). Index sizes in the log of that run against all other runs:

| DB | this run | every other run |
|---|---|---|
| `pref` (first `kmermatcher`) | 20,015,494 entries | 20,000,000 |
| `clu`, `input_rep`, `pref_rep` | 19,959,560 | 19,959,558 |
| `clu_rep` and final clusters | 19,957,551 | 19,957,552 |

The first `kmermatcher` wrote 15,494 entries too many (some keys twice); everything after follows from that.
`kmermatcher` writes its result with one thread: a result block per representative, then one entry for every
sequence not marked as written. Extra entries mean the sorted k-mer list had a representative in more than one
block, so the fault is in or before the parallel k-mer stage, which reads key, length and data of every sequence
through the index. The index merge after it only copies entries and has a single non-empty part here.

Attempts to reproduce it (20M sequences, 20 threads, results compared with an order-independent checksum of every
DB, `src/test/TestDbChecksum.cpp`):

| Test | Runs | Different |
|---|---|---|
| whole linclust, reader before the fix, every tmp DB checked | 8 | 0 |
| `kmermatcher` alone, reader before the fix | 30 | 0 |
| `kmermatcher` alone, current reader | 30 | 0 |
| `kmermatcher` alone, text index (run in parallel with the above: 40 threads on 20 cores) | 60 | 0 |

Resolved (project owner, 2026-11-19): this is a known race condition in this fork's `kmermatcher`, independent of
the index. Single off cases of this kind are to be ignored when results are compared; it is not counted against
the parquet path.

## 7. Options for later

- Module parameters instead of environment variables; decide the default.
- Command that prints a parquet index as text (for the awk uses and for exchange with other tools).
- Hard memory cap for the run sort in the streaming merge (cut unsorted parts into chunks).
- Reduce the remaining per-sequence arrays of `align2clust` and the `std::list` per sequence in `mergeclusters`.
- 64-bit id build (`MMSEQS_INT64_IDS`) for more than 4.29e9 entries: permutation becomes 16 B/entry.
- Datasets beyond UniRef100: Metaclust nr (about 710M sequences, 103 GB gz), BFD (about 2.2B, 292 GB gz).
