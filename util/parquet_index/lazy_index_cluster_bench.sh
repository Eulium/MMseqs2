#!/bin/bash
# Lazy index (src/commons/LazyIndex.h) in the clustering workflows, which read the largest DBs:
# linclust of a random subset of a sequence DB with the text index (MMseqs2 as it is) and with each backend.
#   variants: text | parquet | parquet-custom | binary  (lazy index only for the input DB)
#             parquet-custom+all                        (MMSEQS_WRITE_LAZY_INDEX: also for every DB the workflow writes)
#   measured per run: wall time, CPU time, peak RssAnon (allocated memory), peak RSS; from MMSEQS_INDEX_DEBUG: how many
#   index opens used the backend and how many read a text index of at least half the input size (and its size in RAM)
#   check: same clusters (sorted representative/member pairs) as with the text index
# usage: lazy_index_cluster_bench.sh <mmseqs> <sequence DB with _h> <workdir> <sequences> [rounds] [threads] [workflow]
set -uo pipefail
M="$1"; SRC="$2"; W="$3"; NSEQ="$4"; ROUNDS="${5:-1}"; T="${6:-$(nproc)}"; WF="${7:-linclust}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
mkdir -p "$W"; cd "$W"
log() { echo "[$(date '+%F %T')] $*"; }
CSV="$W/cluster_$WF.csv"
echo "variant,round,wall_s,cpu_s,peak_rss_mb,anon_mb,lazy_opens,big_text_opens,largest_text_index_mb,clusters,result,other_cores" > "$CSV"
# CPU time of the whole machine so far (jiffies: user nice system irq softirq steal); the difference over a run minus
# the run's own CPU time is what other processes used: the spark machine is shared, such runs are marked
busy() { awk '/^cpu / { print $2 + $3 + $4 + $7 + $8 + $9 }' /proc/stat; }
HZ=$(getconf CLK_TCK)

if [ ! -e cdb.dbtype ]; then
    log "subset of $NSEQ sequences of $SRC"
    cut -f1 "$SRC.index" | shuf -n "$NSEQ" --random-source=<(yes 5) | sort -n > cdb.keys
    "$M" createsubdb cdb.keys "$SRC" cdb --subdb-mode 0 -v 1
    "$M" createsubdb cdb.keys "${SRC}_h" cdb_h --subdb-mode 0 -v 1
    for d in cdb cdb_h; do MMSEQS_WRITE_LAZY_INDEX=binary "$M" convert2parquet "$d.index" "$d.index.parquet" -v 1; done
fi
N=$(wc -l < cdb.index)
log "$WF of $N sequences, $T threads: text index $(stat -c %s cdb.index) B ($((N * 24 / 1048576)) MB as array), parquet $(stat -c %s cdb.index.parquet) B, binary $(stat -c %s cdb.index.bin) B"

for r in $(seq 1 "$ROUNDS"); do
    for v in text parquet parquet-custom binary parquet-custom+all; do
        env=(env MMSEQS_INDEX_DEBUG=1)
        [ "$v" != text ] && env+=(MMSEQS_INDEX_BACKEND="${v%+all}")
        [ "$v" = parquet-custom+all ] && env+=(MMSEQS_WRITE_LAZY_INDEX=parquet)
        out="clu_${v//+/_}_$r"
        b0=$(busy)
        python3 -I "$HERE/peak_anon.py" anon.txt "$TIME_BIN" -o t.txt -f "%e %U %S %M" "${env[@]}" \
            "$M" "$WF" cdb "$out" "tmp_${v//+/_}_$r" --threads "$T" -v 2 > "log_${WF}_${v//+/_}.txt" 2>&1
        rc=$?
        b1=$(busy)
        read -r w u s k < <(tail -1 t.txt); read -r a < anon.txt
        other=$(awk -v d="$((b1 - b0))" -v hz="$HZ" -v c="$u" -v c2="$s" -v w="$w" 'BEGIN { o = (d / hz - c - c2) / (w > 0.05 ? w : 0.05); printf "%.2f", (o > 0 ? o : 0) }')
        L="log_${WF}_${v//+/_}.txt"
        lazy=$(grep "Index backend" "$L" | grep -vc "Index backend text")
        # text indexes of at least half the input: "... (<entries> entries, <bytes> bytes in memory, access mode <m>)"
        read -r big largest < <(grep "Index backend text" "$L" | sed 's/.*(\([0-9]*\) entries, \([0-9]*\) bytes.*/\1 \2/' \
            | awk -v n="$N" '$1 * 2 >= n { c++; if ($2 > m) m = $2 } END { printf "%d %.0f\n", c, m / 1048576 }')
        res="failed (exit $rc)"; clusters=NA
        if [ "$rc" -eq 0 ]; then
            "$M" createtsv cdb cdb "$out" "$out.tsv" --threads "$T" -v 1 > /dev/null 2>&1
            res=$(sort -S 4G --parallel="$T" "$out.tsv" | md5sum | cut -c1-12); clusters=$(cut -f1 "$out.tsv" | sort -S 4G -u | wc -l)
        fi
        echo "$v,$r,$w,$(awk -v a="$u" -v b="$s" 'BEGIN { print a + b }'),$((k / 1024)),$((a / 1024)),$lazy,$big,$largest,$clusters,$res,$other" >> "$CSV"
        log "   round $r $v: $w s, anon $((a / 1024)) MB, $clusters clusters, $res, other processes used $other cores"
    done
done
log "access modes of the text index opens of at least half the input (variant parquet-custom+all): 0 NOSORT, 1 SORT_BY_LENGTH, 2 LINEAR_ACCCESS, 3 SORT_BY_ID"
grep "Index backend text" "log_${WF}_parquet-custom_all.txt" | sed 's/.* for \(.*\) (\([0-9]*\) entries.*access mode \([0-9]*\)).*/\2 \3 \1/' \
    | awk -v n="$N" '$1 * 2 >= n { sub(/.*\//, "", $3); print "   mode " $2 ": " $3 " (" $1 " entries)" }' | sort | uniq -c
python3 -I - "$CSV" <<'PY'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1])))
ref = next((r["result"] for r in rows if r["variant"] == "text"), None)
base = st.median(float(r["wall_s"]) for r in rows if r["variant"] == "text")
print(f'{"variant":20} {"wall s":>9} {"x":>6} {"cpu s":>9} {"anon MB":>8} {"rss MB":>8} {"lazy opens":>10} {"big text opens":>14} {"largest MB":>10} {"other cores":>11}  clusters')
for v in dict.fromkeys(r["variant"] for r in rows):
    x = [r for r in rows if r["variant"] == v]
    m = lambda k: st.median(float(r[k]) for r in x)
    same = "same clusters" if all(r["result"] == ref for r in x) else "DIFFERENT: " + x[0]["result"]
    print(f'{v:20} {m("wall_s"):9.1f} {m("wall_s") / base:6.2f} {m("cpu_s"):9.0f} {m("anon_mb"):8.0f} {m("peak_rss_mb"):8.0f} '
          f'{m("lazy_opens"):10.0f} {m("big_text_opens"):14.0f} {m("largest_text_index_mb"):10.0f} {max(float(r["other_cores"]) for r in x):11.1f}  {x[0]["clusters"]} {same}')
PY
log "done"
