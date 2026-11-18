#!/bin/bash
# Benchmark of the lazy index (src/commons/LazyIndex.h): text index in memory (MMseqs2 as it is) vs the index on disk
# as parquet (carquet reader), parquet (custom decoder) and binary records (mmap).
#   1. index files: write time and size
#   2. workloads, each as one mmseqs process per backend (wall time, CPU time, peak RssAnon = allocated memory,
#      peak RSS incl. mapped files), interleaved rounds, warm page cache:
#        view1          view of one entry: opening the DB + one lookup
#        sub_<K>_<o>    createsubdb --subdb-mode 1 (index only) with K keys in random (rnd) or key (srt) order:
#                       getId + getOffset + getEntryLen per key, one thread
#        sub_all_srt    all keys in key order (a full pass over the index)
#        m8_t<threads>  convertalis of a search result: a target lookup per hit from all threads (random access),
#                       plus the header DB (the search itself runs against a 5M entry subset of a large DB)
#   3. other parquet readers on the same file: DuckDB CLI and R arrow (if available), single key and 1000 keys
# usage: lazy_index_bench.sh <mmseqs> <sequence DB with _h> <workdir> [rounds] [threads] [search queries]
# needs peak_anon.py next to this script and GNU time
set -uo pipefail
M="$1"; SRC="$2"; W="$3"; ROUNDS="${4:-3}"; T="${5:-$(nproc)}"; NQ="${6:-2000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
BACKENDS="text parquet parquet-custom binary"
mkdir -p "$W"; cd "$W"
log() { echo "[$(date '+%F %T')] $*"; }
CSV="$W/lazy_index.csv"
echo "workload,backend,round,wall_s,cpu_s,peak_rss_mb,anon_mb,other_cores" > "$CSV"
# CPU time of the whole machine so far (jiffies: user nice system irq softirq steal); the difference over a run minus
# the run's own CPU time is what other processes used: the spark machine is shared, such runs are marked
busy() { awk '/^cpu / { print $2 + $3 + $4 + $7 + $8 + $9 }' /proc/stat; }
HZ=$(getconf CLK_TCK)

# run <workload> <backend> <round> <command...>: one measured process
run() {
    local wl="$1" b="$2" r="$3"; shift 3
    local env=(env)
    [ "$b" != text ] && env=(env MMSEQS_INDEX_BACKEND="$b")
    local b0; b0=$(busy)
    python3 -I "$HERE/peak_anon.py" "$W/anon.txt" "$TIME_BIN" -o "$W/t.txt" -f "%e %U %S %M" "${env[@]}" "$@" > "$W/log_${wl}_$b.txt" 2>&1
    local rc=$? w u s k a b1
    b1=$(busy)
    read -r w u s k < <(tail -1 "$W/t.txt"); read -r a < "$W/anon.txt"
    [ "$rc" -eq 0 ] || { log "FAILED $wl $b (exit $rc, see $W/log_${wl}_$b.txt)"; return; }
    [ "$r" -gt 0 ] && echo "$wl,$b,$r,$w,$(awk -v a="$u" -v b="$s" 'BEGIN { print a + b }'),$((k / 1024)),$((a / 1024)),$(awk -v d="$((b1 - b0))" -v hz="$HZ" -v c="$u" -v c2="$s" -v w="$w" 'BEGIN { o = (d / hz - c - c2) / (w > 0.05 ? w : 0.05); printf "%.2f", (o > 0 ? o : 0) }')" >> "$CSV"
}

log "0. DB $SRC"
for f in "" .index .dbtype .lookup _h _h.index _h.dbtype; do
    [ -e "$SRC$f" ] && [ ! -e "db$f" ] && ln -s "$(readlink -f "$SRC$f")" "db$f"
done
for f in "$SRC".[0-9]* "${SRC}_h".[0-9]*; do   # split data files
    [ -e "$f" ] && [ ! -e "db${f#"$SRC"}" ] && ln -s "$(readlink -f "$f")" "db${f#"$SRC"}"
done
N=$(wc -l < db.index)
log "   $N entries, text index $(stat -Lc %s db.index) B, as array in memory: $((N * 24 / 1048576)) MB"

log "1. index files"
{
    echo "file,write_s,bytes"
    for d in db db_h; do
        "$TIME_BIN" -o t.txt -f "%e" "$M" convert2parquet "$d.index" "$d.index.parquet" --threads "$T" -v 1
        echo "$d.index.parquet,$(cat t.txt),$(stat -c %s "$d.index.parquet")"
        MMSEQS_WRITE_LAZY_INDEX=binary "$TIME_BIN" -o t.txt -f "%e" "$M" convert2parquet "$d.index" "$d.index.parquet" --threads "$T" -v 1
        echo "$d.index.bin (+ parquet again),$(cat t.txt),$(stat -c %s "$d.index.bin")"
        echo "$d.index (text),,$(stat -Lc %s "$d.index")"
    done
} > index_files.csv
cat index_files.csv

log "2. key lists and search result"
KS=""
for k in 1000 100000 10000000; do
    [ "$k" -lt "$N" ] || continue
    KS="$KS $k"
    [ -e "keys_${k}_rnd" ] || cut -f1 db.index | shuf -n "$k" --random-source=<(yes) > "keys_${k}_rnd"
    [ -e "keys_${k}_srt" ] || sort -n "keys_${k}_rnd" > "keys_${k}_srt"
done
[ -e keys_all_srt ] || cut -f1 db.index | sort -n > keys_all_srt
ONE=$(sed -n "$((N / 3))p" db.index | cut -f1)
if [ ! -e res.dbtype ]; then
    cut -f1 db.index | shuf -n "$NQ" --random-source=<(yes 7) > query.keys
    "$M" createsubdb query.keys db q --subdb-mode 0 -v 1
    "$M" createsubdb query.keys db_h q_h --subdb-mode 0 -v 1
    # the search only makes the hits for convertalis: against a subset of the DB if it is large (createsubdb keeps
    # the keys, so the hits are entries spread over the whole DB)
    TGT=db
    if [ "$N" -gt 6000000 ]; then
        { cat query.keys; cut -f1 db.index | shuf -n 5000000 --random-source=<(yes 3); } | sort -nu > tsub.keys
        "$M" createsubdb tsub.keys db tsub --subdb-mode 0 -v 1
        TGT=tsub
    fi
    "$TIME_BIN" -o t.txt -f "%e s, peak %M kB" "$M" search q "$TGT" res tmp -s 1 --max-seqs 300 --threads "$T" -v 1 > search.log 2>&1
    log "   search of $NQ queries against $TGT: $(cat t.txt), $(cat res res.[0-9]* 2> /dev/null | tr -d '\0' | wc -l) hits"
fi

log "3. workloads ($ROUNDS rounds, round 0 = warm-up)"
for r in $(seq 0 "$ROUNDS"); do
    for b in $BACKENDS; do
        run view1 "$b" "$r" "$M" view db --id-list "$ONE" -v 1
        for k in $KS; do
            for o in rnd srt; do
                run "sub_${k}_$o" "$b" "$r" "$M" createsubdb "keys_${k}_$o" db "out_sub_$b" --subdb-mode 1 -v 1
            done
        done
        run sub_all_srt "$b" "$r" "$M" createsubdb keys_all_srt db "out_sub_$b" --subdb-mode 1 -v 1
        run "m8_t$T" "$b" "$r" "$M" convertalis q db res "out_$b.m8" --threads "$T" -v 1
        run m8_t1 "$b" "$r" "$M" convertalis q db res "out_${b}_t1.m8" --threads 1 -v 1
    done
    log "   round $r done"
done
# same output with every backend
for b in $BACKENDS; do
    printf "   %-15s createsubdb index %s  convertalis %s (%s lines)\n" "$b" "$(md5sum < "out_sub_$b.index" | cut -c1-12)" \
        "$(sort "out_$b.m8" | md5sum | cut -c1-12)" "$(wc -l < "out_$b.m8")"
done

log "4. other parquet readers on db.index.parquet (wall time of one process, median of 3)"
med3() { for i in 1 2 3; do "$TIME_BIN" -o t.txt -f "%e" "$@" > /dev/null 2>&1; cat t.txt; done | sort -n | sed -n 2p; }
{
    echo "reader,query,wall_s"
    if command -v duckdb > /dev/null; then
        echo "duckdb_cli,1 key,$(med3 duckdb -c "select * from 'db.index.parquet' where key = $ONE")"
        echo "duckdb_cli,1000 keys (join),$(med3 duckdb -c "select count(*), sum(p.\"offset\") from 'db.index.parquet' p join read_csv('keys_1000_rnd', header = false, columns = {'key': 'BIGINT'}) k using (key)")"
        echo "duckdb_cli,full scan,$(med3 duckdb -c "select count(*), sum(length) from 'db.index.parquet'")"
        echo "duckdb_cli,rows read correctly,$(duckdb -csv -noheader -c "select count(*) from 'db.index.parquet'") of $N"
    fi
    if command -v Rscript > /dev/null && Rscript -e 'library(arrow)' > /dev/null 2>&1; then
        echo "r_arrow,1 key,$(med3 Rscript -e "library(arrow); library(dplyr); print(open_dataset('db.index.parquet') |> filter(key == $ONE) |> collect())")"
        echo "r_arrow,full scan,$(med3 Rscript -e "library(arrow); library(dplyr); print(open_dataset('db.index.parquet') |> summarise(n = n(), s = sum(length)) |> collect())")"
        echo "r_arrow,start-up only,$(med3 Rscript -e "library(arrow); library(dplyr)")"
    fi
    echo "mmseqs parquet (carquet),1 key (view),$(med3 env MMSEQS_INDEX_BACKEND=parquet "$M" view db --id-list "$ONE" -v 1)"
    echo "mmseqs parquet (custom),1 key (view),$(med3 env MMSEQS_INDEX_BACKEND=parquet-custom "$M" view db --id-list "$ONE" -v 1)"
} > other_readers.csv
cat other_readers.csv

log "5. summary (median over rounds; x = time relative to the text index)"
python3 -I - "$CSV" $BACKENDS <<'PY'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1]))); backends = sys.argv[2:]
wls = list(dict.fromkeys(r["workload"] for r in rows))
med = lambda wl, b, k: (lambda x: st.median(x) if x else float("nan"))([float(r[k]) for r in rows if r["workload"] == wl and r["backend"] == b])
other = lambda wl: max([float(r["other_cores"]) for r in rows if r["workload"] == wl] or [0])
print(f'{"workload":18}' + "".join(f'{b:>34}' for b in backends) + "   other load")
print(f'{"":18}' + "".join(f'{"wall s (x)  anon MB  rss MB":>34}' for b in backends))
for wl in wls:
    base = med(wl, "text", "wall_s")
    print(f'{wl:18}' + "".join(f'{med(wl, b, "wall_s"):9.2f} ({med(wl, b, "wall_s") / base:5.2f}) {med(wl, b, "anon_mb"):8.0f} {med(wl, b, "peak_rss_mb"):7.0f}' for b in backends)
          + f'   {other(wl):5.1f} cores max' + ("  <- disturbed" if other(wl) > 2 else ""))
print("other load: CPU cores used by other processes of the machine during a run (largest of all runs of the workload)")
# cost per lookup: slope between the smallest and the largest random key list
ks = sorted({int(w.split("_")[1]) for w in wls if w.startswith("sub_") and w.split("_")[1].isdigit()})
if len(ks) >= 2:
    print("\ncreatesubdb, time per key (largest minus smallest key list):")
    for o in ("rnd", "srt"):
        print(f'  {o}: ' + "  ".join(f'{b} {(med(f"sub_{ks[-1]}_{o}", b, "wall_s") - med(f"sub_{ks[0]}_{o}", b, "wall_s")) / (ks[-1] - ks[0]) * 1e6:.3f} us' for b in backends))
PY
log "done"
