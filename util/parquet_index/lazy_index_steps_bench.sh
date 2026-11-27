#!/bin/bash
# Clustering workflow with the text index (MMseqs2 as it is) against the parquet index (src/commons/LazyIndex.h),
# with the time and memory of every module of the workflow (mmseqs_steps.py).
#   variants (label=binary:kind[:rows per index page] runs a kind with another build):
#             text    text index for all DBs (no switch set)
#             both    text index and parquet file written for every DB, read through the parquet file
#             notext  parquet file only: MMSEQS_WRITE_LAZY_INDEX=parquet,notext, input DB without text index
#             both_textmem, notext_textmem  the same with MMSEQS_LAZY_INDEX_TEXT_MEMORY=1: kmermatcher plans its
#                     memory as with the text index, for the same clusters as text when it works in parts
#   per run:  wall and CPU time, peak RssAnon of the process tree, CPU used by other processes of the machine,
#             index opens by kind (MMSEQS_INDEX_DEBUG), per module: wall, CPU, peak RssAnon
#   check:    same clusters (sorted representative/member pairs) as the first text run
# The order of the variants alternates between the rounds. Nothing else should run on the machine.
# Memory of a shared machine: a run starts when NEED_MB of memory are available (default 0: no wait) and is stopped
# and repeated later if less than MIN_AVAIL_MB (default 6000) are left while it runs, before the machine swaps or
# the OOM killer picks a process. Per run the lowest available memory and the MB swapped out are recorded.
# Runs that are in the result file already are skipped: a stopped benchmark continues when it is started again.
# usage: lazy_index_steps_bench.sh <mmseqs> <sequence DB with _h> <workdir> <sequences, 0 = all> [rounds] [threads]
#                                  [workflow] [variants, e.g. "text notext"] [extra workflow arguments]
set -uo pipefail
M="$1"; SRC="$2"; W="$3"; NSEQ="$4"; ROUNDS="${5:-1}"; T="${6:-$(nproc)}"; WF="${7:-linclust}"; VARIANTS="${8:-text both notext}"
EXTRA="${9:-}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TIME_BIN="${TIME_BIN:-/usr/bin/time}"
mkdir -p "$W"; cd "$W"
log() { echo "[$(date '+%F %T')] $*"; }
# plain read/write copy: cp uses the NFS 4.2 server-side copy, which hung a job for 23 h (nfs42_proc_copy)
copy() { cat -- "$1" > "$2"; }
CSV="$W/runs_$WF.csv"; STEPS="$W/steps_$WF.csv"
[ -e "$CSV" ] || echo "variant,round,wall_s,cpu_s,peak_rss_mb,anon_mb,lazy_opens,text_opens,loaded_opens,big_text_opens,clusters,result,other_cores,exit,min_avail_mb,swapout_mb" > "$CSV"
[ -e "$STEPS" ] || echo "variant,round,start,module,wall_s,user_s,sys_s,peak_rss_kb,anon_kb,exit" > "$STEPS"
busy() { awk '/^cpu / { print $2 + $3 + $4 + $7 + $8 + $9 }' /proc/stat; }
HZ=$(getconf CLK_TCK)
NEED_MB="${NEED_MB:-0}"; MIN_AVAIL_MB="${MIN_AVAIL_MB:-6000}"
avail() { awk '/^MemAvailable:/ { print int($2 / 1024) }' /proc/meminfo; }
swapout() { awk '/^pswpout / { print $2 }' /proc/vmstat; }
# wait until NEED_MB are available at three checks in a row
wait_mem() {
    local ok=0 told=0
    [ "$NEED_MB" -le 0 ] && return
    while [ "$ok" -lt 3 ]; do
        if [ "$(avail)" -ge "$NEED_MB" ]; then ok=$((ok + 1)); else
            ok=0; [ "$told" -eq 0 ] && log "   waiting for $NEED_MB MB of available memory ($(avail) MB now)"; told=1
        fi
        [ "$ok" -lt 3 ] && sleep 20
    done
}
# watches the run (process group $1): lowest available memory to minavail.txt, stops the run if memory gets short
guard() {
    local min a; min=$(avail); echo "$min" > minavail.txt
    while kill -0 "$1" 2> /dev/null; do
        a=$(avail)
        [ "$a" -lt "$min" ] && { min=$a; echo "$min" > minavail.txt; }
        if [ "$a" -lt "$MIN_AVAIL_MB" ]; then
            echo "$a" > killed.txt
            kill -TERM -- "-$1" 2> /dev/null; sleep 5; kill -KILL -- "-$1" 2> /dev/null
            return
        fi
        sleep 1
    done
}
RUN=""
trap '[ -n "$RUN" ] && kill -KILL -- "-$RUN" 2> /dev/null; exit 143' TERM INT

# input: cdb (text index + parquet file) and cdbp (the same data, parquet file only)
if [ ! -e cdb.dbtype ]; then
    if [ "$NSEQ" -gt 0 ]; then
        log "subset of $NSEQ sequences of $SRC"
        cut -f1 "$SRC.index" | shuf -n "$NSEQ" --random-source=<(yes 5) | sort -n > cdb.keys
        "$M" createsubdb cdb.keys "$SRC" cdb --subdb-mode 0 -v 1
        "$M" createsubdb cdb.keys "${SRC}_h" cdb_h --subdb-mode 0 -v 1
    else
        for s in "" _h; do
            for f in "$SRC$s" "$SRC$s".[0-9]*; do [ -e "$f" ] && ln -s "$(readlink -f "$f")" "cdb$s${f#"$SRC$s"}"; done
            copy "$SRC$s.index" "cdb$s.index"; copy "$SRC$s.dbtype" "cdb$s.dbtype"
        done
    fi
    for s in "" _h; do
        "$M" convert2parquet "cdb$s.index" "cdb$s.index.parquet" -v 1
        for f in "cdb$s" "cdb$s".[0-9]*; do [ -e "$f" ] && ln -s "$(readlink -f "$f")" "cdbp$s${f#"cdb$s"}"; done
        copy "cdb$s.dbtype" "cdbp$s.dbtype"; copy "cdb$s.index.parquet" "cdbp$s.index.parquet"
    done
fi
N=$(wc -l < cdb.index)
log "$WF of $N sequences, $T threads, variants: $VARIANTS; text index $(stat -c %s cdb.index) B ($((N * 24 / 1048576)) MB as array), parquet $(stat -c %s cdb.index.parquet) B"
log "load before the benchmark: $(cut -d' ' -f1-3 /proc/loadavg)"

for r in $(seq 1 "$ROUNDS"); do
    order="$VARIANTS"
    [ $((r % 2)) -eq 0 ] && order=$(echo "$VARIANTS" | tr ' ' '\n' | tac | tr '\n' ' ')
    for v in $order; do
        # a variant is a kind (text, both, notext) or label=binary:kind[:rows per index page] for another build
        bin="$M"; kind="$v"; rows=""
        if [ "${v#*=}" != "$v" ]; then
            spec="${v#*=}"; v="${v%%=*}"
            IFS=: read -r bin kind rows <<< "$spec"
        fi
        env=(env MMSEQS_INDEX_DEBUG=1 MMSEQS_REAL="$bin" MMSEQS_STEPS="$W/steps.tmp"); in=cdb
        grep -q "^$v,$r," "$CSV" && { log "   round $r $v: done before"; continue; }
        # <kind>_textmem: the lazy index registers the memory of a text index (same k-mer parts as the text run)
        if [ "${kind%_textmem}" != "$kind" ]; then kind="${kind%_textmem}"; env+=(MMSEQS_LAZY_INDEX_TEXT_MEMORY=1); fi
        [ "$kind" = both ] && env+=(MMSEQS_INDEX_BACKEND=parquet-custom MMSEQS_WRITE_LAZY_INDEX=parquet)
        [ "$kind" = notext ] && { env+=(MMSEQS_WRITE_LAZY_INDEX=parquet,notext); in=cdbp; }
        if [ -n "$rows" ]; then
            # input index with this page size, and the same for every index the workflow writes
            env+=(MMSEQS_LAZY_INDEX_PAGE_ROWS="$rows"); in="cdbp$rows"
            if [ ! -e "$in.dbtype" ]; then
                for s in "" _h; do
                    for f in "cdb$s" "cdb$s".[0-9]*; do [ -e "$f" ] && ln -s "$(readlink -f "$f")" "$in$s${f#"cdb$s"}"; done
                    copy "cdb$s.dbtype" "$in$s.dbtype"
                    MMSEQS_LAZY_INDEX_PAGE_ROWS="$rows" "$M" convert2parquet "cdb$s.index" "$in$s.index.parquet" -v 1
                done
            fi
        fi
        out="clu_${v}_$r"; L="log_${WF}_${v}_$r.txt"
        while :; do
            wait_mem
            rm -f steps.tmp killed.txt; sync
            b0=$(busy); s0=$(swapout)
            # own process group, so that the guard can stop the whole run
            # shellcheck disable=SC2086
            setsid python3 -I "$HERE/peak_anon.py" anon.txt "$TIME_BIN" -o t.txt -f "%e %U %S %M" "${env[@]}" \
                "$HERE/mmseqs_steps.py" "$WF" "$in" "$out" "tmp_${v}_$r" --threads "$T" -v 2 $EXTRA > "$L" 2>&1 &
            RUN=$!
            guard "$RUN" & g=$!
            wait "$RUN"; rc=$?
            wait "$g" 2> /dev/null; RUN=""
            [ -e killed.txt ] || break
            log "   round $r $v stopped: $(cat killed.txt) MB of memory available (limit $MIN_AVAIL_MB MB), run again later"
            mv "$L" "$L.stopped"; rm -rf "tmp_${v}_$r"; "$M" rmdb "$out" -v 1 > /dev/null 2>&1; rm -f "$out".* "$out"
            [ "$NEED_MB" -lt $((MIN_AVAIL_MB * 4)) ] && NEED_MB=$((MIN_AVAIL_MB * 4))
        done
        b1=$(busy); s1=$(swapout); minav=$(cat minavail.txt); swp=$(((s1 - s0) * $(getconf PAGESIZE) / 1048576))
        read -r w u s k < <(tail -1 t.txt); read -r a < anon.txt
        other=$(awk -v d="$((b1 - b0))" -v hz="$HZ" -v c="$u" -v c2="$s" -v w="$w" 'BEGIN { o = (d / hz - c - c2) / (w > 0.05 ? w : 0.05); printf "%.2f", (o > 0 ? o : 0) }')
        lazy=$(grep "Index backend" "$L" | grep -vc "Index backend \(text\|parquet-loaded\)")
        text=$(grep -c "Index backend text" "$L"); loaded=$(grep -c "Index backend parquet-loaded" "$L")
        big=$(grep "Index backend \(text\|parquet-loaded\)" "$L" | sed 's/.*(\([0-9]*\) entries.*/\1/' | awk -v n="$N" '$1 * 2 >= n { c++ } END { print c + 0 }')
        res="failed"; clusters=NA
        if [ "$rc" -eq 0 ]; then
            "$M" createtsv cdb cdb "$out" "$out.tsv" --threads "$T" -v 1 > /dev/null 2>&1
            res=$(sort -S 8G --parallel="$T" "$out.tsv" | md5sum | cut -c1-12); clusters=$(cut -f1 "$out.tsv" | sort -S 8G -u | wc -l)
            rm -f "$out.tsv"
        fi
        awk -v v="$v" -v r="$r" '{ print v "," r "," $0 }' steps.tmp >> "$STEPS"
        echo "$v,$r,$w,$(awk -v a="$u" -v b="$s" 'BEGIN { print a + b }'),$((k / 1024)),$((a / 1024)),$lazy,$text,$loaded,$big,$clusters,$res,$other,$rc,$minav,$swp" >> "$CSV"
        log "   round $r $v: $w s, anon $((a / 1024)) MB, $clusters clusters, $res, index opens lazy/text/loaded $lazy/$text/$loaded (large not lazy: $big), other processes $other cores, lowest available memory $minav MB, swapped out $swp MB"
        # the results are compared by their checksum; the DBs are not needed again
        rm -rf "tmp_${v}_$r"; "$M" rmdb "$out" -v 1 > /dev/null 2>&1
    done
done
python3 -I - "$CSV" "$STEPS" <<'PY'
import csv, sys, statistics as st
rows = list(csv.DictReader(open(sys.argv[1]))); steps = list(csv.DictReader(open(sys.argv[2])))
variants = list(dict.fromkeys(r["variant"] for r in rows))
ref = next((r["result"] for r in rows if r["variant"] == "text"), rows[0]["result"])
base = st.median(float(r["wall_s"]) for r in rows if r["variant"] == variants[0])
print(f'\n{"variant":8} {"runs":>4} {"wall s":>9} {"min":>8} {"max":>8} {"x":>6} {"cpu s":>9} {"anon MB":>8} {"lazy":>5} {"text":>5} {"loaded":>6} {"large not lazy":>14} {"other cores":>11}  result')
for v in variants:
    x = [r for r in rows if r["variant"] == v]
    m = lambda k: st.median(float(r[k]) for r in x)
    w = [float(r["wall_s"]) for r in x]
    same = "same clusters" if all(r["result"] == ref for r in x) else "DIFFERENT: " + ",".join(sorted(set(r["result"] for r in x)))
    print(f'{v:8} {len(x):4d} {m("wall_s"):9.1f} {min(w):8.1f} {max(w):8.1f} {m("wall_s") / base:6.2f} {m("cpu_s"):9.0f} {m("anon_mb"):8.0f} '
          f'{m("lazy_opens"):5.0f} {m("text_opens"):5.0f} {m("loaded_opens"):6.0f} {m("big_text_opens"):14.0f} {max(float(r["other_cores"]) for r in x):11.1f}  {x[0]["clusters"]} {same}')
# modules in the order of the first run; the n-th call of a module in a run is one step
key = lambda s, i: (s["module"], i)
table = {}; order = []
for v in variants:
    for r in sorted(set(s["round"] for s in steps if s["variant"] == v)):
        seen = {}
        for s in sorted((s for s in steps if s["variant"] == v and s["round"] == r), key=lambda s: float(s["start"])):
            i = seen.get(s["module"], 0); seen[s["module"]] = i + 1
            k = key(s, i)
            if k not in order: order.append(k)
            table.setdefault((k, v), []).append(s)
print("\nper module, median over the rounds: wall s | CPU s | peak RssAnon MB")
print(f'{"module":28}' + "".join(f'{v:>30}' for v in variants) + "   wall vs " + variants[0])
for k in order:
    line = f'{k[0] + ("#" + str(k[1] + 1) if k[1] else ""):28}'
    walls = []
    for v in variants:
        x = table.get((k, v), [])
        if not x:
            line += f'{"-":>30}'; walls.append(None); continue
        w = st.median(float(s["wall_s"]) for s in x); c = st.median(float(s["user_s"]) + float(s["sys_s"]) for s in x)
        a = st.median(float(s["anon_kb"]) for s in x) / 1024
        line += f'{w:12.1f} {c:9.0f} {a:7.0f}'; walls.append(w)
    if walls[0] and all(w is not None for w in walls):
        line += "   " + " ".join(f'{w - walls[0]:+7.1f}' for w in walls[1:])
    print(line)
PY
log "done"
