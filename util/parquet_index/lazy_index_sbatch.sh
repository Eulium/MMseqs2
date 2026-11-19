#!/bin/bash
#SBATCH --job-name=mmseqs_lazy_index
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00
#SBATCH --exclusive
#SBATCH --output=lazy_index_%j.log
#
# Builds MMseqs2 (branch parquet_index) and runs all correctness tests and benchmarks of the parquet index
# (src/commons/LazyIndex.h). Results: one CSV with every measured value and a summary text file.
#
#   sbatch [sbatch options] lazy_index_sbatch.sh <MMseqs2 repo> <database> <tmp dir> [output dir]
#   bash lazy_index_sbatch.sh <MMseqs2 repo> <database> <tmp dir> [output dir]          (without slurm)
#
#   <MMseqs2 repo>  checkout of the branch; the helper scripts are taken from <repo>/util/parquet_index
#                   (sbatch runs a copy of this script, so they cannot be found next to it)
#   <database>      MMseqs2 sequence DB (made with createdb, with its _h DB) or a FASTA file (.gz allowed);
#                   a FASTA file is converted with createdb in <tmp dir>
#   <tmp dir>       build directory, subsets, workflow temporary files: local scratch with space for about
#                   four times the largest subset plus the build (1 GB)
#   [output dir]    default <tmp dir>/results: lazy_index_results.csv, lazy_index_summary.txt, logs/, csv/
#
# Build requirements: cmake >= 3.15, gcc/g++ with OpenMP, make, zlib; a Rust toolchain (cargo) for block-aligner;
# network access at configure time unless lz4 is installed (carquet downloads lz4, see lib/carquet/CMakeLists.txt).
# Load them before the build if the cluster uses modules, e.g. in LAZY_INDEX_SETUP:
#   LAZY_INDEX_SETUP="module load gcc cmake rust" sbatch lazy_index_sbatch.sh ...
# The build is:
#   cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=1 -DHAVE_MPI=0 -DCMAKE_INSTALL_PREFIX=. <repo>
#   make -j <threads> mmseqs test_lazyindex test_dbchecksum
# Run requirements: python3 (>= 3.9), GNU time (/usr/bin/time or TIME_BIN), Linux /proc.
#
# Settings (environment variables):
#   STEPS          steps to run, default "build test lookup cluster steps"
#   THREADS        default SLURM_CPUS_PER_TASK, else all cores. On a node shared with other jobs use fewer threads
#                  than cores: one busy outside core slowed a 20 thread run on 20 cores by 12 %.
#   ROUNDS         rounds of every benchmark, default 3
#   TEST_SEQS      sequences of the subset for the correctness tests, default 200000
#   LOOKUP_SEQS    sequences for the lookup benchmark, 0 = whole database (default)
#   CLUSTER_SEQS   sizes of the clustering benchmarks, space separated, 0 = whole database; default "20000000".
#                  Sizes above the database size use the whole database.
#   STEP_VARIANTS  variants of the per-module benchmark, default "text both notext"
#   WORKFLOWS      workflows of the per-module benchmark, default "linclust"
#   EXTRA_ARGS     extra arguments for the workflows of the per-module benchmark, e.g. "--split-memory-limit 80G"
#   BUILD_DIR      build directory, default <tmp dir>/build (an existing build is updated)
#   NEED_MB, MIN_AVAIL_MB   memory wait/guard of lazy_index_steps_bench.sh, default 0 = off
#   CLEAN          1 (default): remove the work directories in <tmp dir> at the end (not the build, not a DB made
#                  from FASTA); 0: keep them
#
# Steps:
#   build    cmake + make
#   test     lazy_index_test.sh (same output with every reader, on the TEST_SEQS subset), test_lazyindex (all
#            access modes of DBReader, text index against carquet/custom/binary readers and against a DB with a
#            parquet index only), cluster workflow text against parquet only (same clusters)
#   lookup   lazy_index_bench.sh: index file sizes, open + lookups, createsubdb with random and sorted keys,
#            convertalis from all threads; readers text, parquet (carquet), parquet (custom decoder), binary
#   cluster  lazy_index_cluster_bench.sh: linclust end to end with the four readers (first size of CLUSTER_SEQS)
#   steps    lazy_index_steps_bench.sh: workflow per module (kmermatcher, align2clust, ...) for STEP_VARIANTS,
#            every size of CLUSTER_SEQS
# Exit code: 0 if every step ran and every check passed.
set -uo pipefail

if [ $# -lt 3 ]; then
    sed -n '11,22p' "$0" | sed 's/^# \{0,1\}//'
    exit 2
fi
REPO="$(readlink -f "$1")"; DB_ARG="$(readlink -f "$2")"; TMP="$3"; OUT="${4:-$3/results}"
mkdir -p "$TMP" "$OUT" || exit 2
TMP="$(readlink -f "$TMP")"; OUT="$(readlink -f "$OUT")"
U="$REPO/util/parquet_index"
STEPS="${STEPS:-build test lookup cluster steps}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
ROUNDS="${ROUNDS:-3}"
TEST_SEQS="${TEST_SEQS:-200000}"
LOOKUP_SEQS="${LOOKUP_SEQS:-0}"
CLUSTER_SEQS="${CLUSTER_SEQS:-20000000}"
STEP_VARIANTS="${STEP_VARIANTS:-text both notext}"
WORKFLOWS="${WORKFLOWS:-linclust}"
EXTRA_ARGS="${EXTRA_ARGS:-}"
BUILD_DIR="${BUILD_DIR:-$TMP/build}"
CLEAN="${CLEAN:-1}"
export NEED_MB="${NEED_MB:-0}" MIN_AVAIL_MB="${MIN_AVAIL_MB:-0}"
export TIME_BIN="${TIME_BIN:-/usr/bin/time}"
WORK="$TMP/lazy_index_work"
LOGS="$OUT/logs"; CSVS="$OUT/csv"
mkdir -p "$WORK" "$LOGS" "$CSVS"
M="$BUILD_DIR/src/mmseqs"
STATUS="$OUT/steps_status.csv"; echo "step,exit,seconds" > "$STATUS"
TESTS="$CSVS/tests.csv"; echo "group,test,result" > "$TESTS"
log() { echo "[$(date '+%F %T')] $*"; }
has() { case " $STEPS " in *" $1 "*) return 0 ;; *) return 1 ;; esac; }
# step <name> <command...>: output to logs/<name>.log, exit code and time to steps_status.csv
step() {
    local name="$1"; shift
    local t0=$SECONDS rc
    log "step $name"
    "$@" > "$LOGS/$name.log" 2>&1; rc=$?
    echo "$name,$rc,$((SECONDS - t0))" >> "$STATUS"
    log "step $name: exit $rc, $((SECONDS - t0)) s"
    return "$rc"
}
# OK/FAIL lines of a test log ("  OK    name") to tests.csv
collect_checks() {   # group, log
    awk -v g="$1" '/^  (OK|FAIL)  / { r = $1; sub(/^  (OK|FAIL)  +/, ""); gsub(/"/, "'"'"'"); printf "%s,\"%s\",%s\n", g, $0, r }' "$2" >> "$TESTS"
}

for f in lazy_index_test.sh lazy_index_bench.sh lazy_index_cluster_bench.sh lazy_index_steps_bench.sh mmseqs_steps.py peak_anon.py; do
    [ -e "$U/$f" ] || { log "missing $U/$f"; exit 2; }
done
[ -x "$TIME_BIN" ] || { log "GNU time not found ($TIME_BIN), set TIME_BIN"; exit 2; }
command -v python3 > /dev/null || { log "python3 not found"; exit 2; }
if [ -n "${LAZY_INDEX_SETUP:-}" ]; then eval "$LAZY_INDEX_SETUP"; fi

{
    echo "date:      $(date '+%F %T')"
    echo "host:      $(hostname), $(nproc) cores available, $(awk '/^MemTotal:/ { printf "%.0f GB", $2 / 1048576 }' /proc/meminfo) RAM, $(uname -srm)"
    echo "cpu:       $(lscpu 2> /dev/null | awk -F: '/^Model name/ { gsub(/^ +/, "", $2); print $2; exit }')"
    echo "slurm:     job ${SLURM_JOB_ID:-none}, cpus per task ${SLURM_CPUS_PER_TASK:-none}"
    echo "repo:      $REPO, $(git -C "$REPO" rev-parse --abbrev-ref HEAD 2> /dev/null) $(git -C "$REPO" rev-parse --short HEAD 2> /dev/null), $(git -C "$REPO" status --porcelain 2> /dev/null | grep -vc '^??') modified files"
    echo "database:  $DB_ARG"
    echo "tmp:       $TMP ($(df -Ph "$TMP" | awk 'NR == 2 { print $4 }') free)"
    echo "settings:  steps \"$STEPS\", threads $THREADS, rounds $ROUNDS, test $TEST_SEQS, lookup $LOOKUP_SEQS, cluster \"$CLUSTER_SEQS\", variants \"$STEP_VARIANTS\", workflows \"$WORKFLOWS\", extra \"$EXTRA_ARGS\""
} > "$OUT/environment.txt"
cat "$OUT/environment.txt"

# ---------------------------------------------------------------- build
build() {
    mkdir -p "$BUILD_DIR" && cd "$BUILD_DIR" || return 1
    if [ ! -e CMakeCache.txt ]; then
        cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=1 -DHAVE_MPI=0 -DCMAKE_INSTALL_PREFIX=. "$REPO" || return 1
    fi
    make -j "$THREADS" mmseqs test_lazyindex test_dbchecksum || return 1
    "$M" version
}
if has build; then
    step build build || { log "build failed, see $LOGS/build.log"; }
fi
if [ ! -x "$M" ]; then
    log "no binary $M"
    echo "binary,127,0" >> "$STATUS"
    STEPS=""
fi

# ---------------------------------------------------------------- database and subsets
SRC=""
prepare_db() {
    if [ -e "$DB_ARG.dbtype" ] && [ -e "$DB_ARG.index" ] && [ -e "${DB_ARG}_h.index" ]; then
        SRC="$DB_ARG"
    elif [ -f "$DB_ARG" ]; then
        SRC="$TMP/lazy_index_db/seqdb"
        if [ ! -e "$SRC.dbtype" ]; then
            mkdir -p "$TMP/lazy_index_db"
            log "createdb $DB_ARG"
            "$TIME_BIN" -f "createdb %e s, peak %M kB" "$M" createdb "$DB_ARG" "$SRC" -v 2 > "$LOGS/createdb.log" 2>&1 \
                || { log "createdb failed, see $LOGS/createdb.log"; return 1; }
        fi
    else
        log "$DB_ARG is neither an MMseqs2 sequence DB (with .index, .dbtype and _h) nor a file"
        return 1
    fi
    NALL=$(wc -l < "$SRC.index")
    log "database $SRC: $NALL sequences, text index $(stat -Lc %s "$SRC.index") B"
    echo "sequences: $NALL in $SRC, text index $(stat -Lc %s "$SRC.index") B, data $(du -shL "$SRC" | cut -f1)" >> "$OUT/environment.txt"
}
# subset <sequences> <name>: path of a DB with that many random sequences (the source DB if it is not larger)
subset() {
    local n="$1" name="$2"
    if [ "$n" -le 0 ] || [ "$n" -ge "$NALL" ]; then echo "$SRC"; return; fi
    local d="$WORK/$name"
    if [ ! -e "$d.dbtype" ]; then
        cut -f1 "$SRC.index" | shuf -n "$n" --random-source=<(yes 5) | sort -n > "$d.keys"
        "$M" createsubdb "$d.keys" "$SRC" "$d" --subdb-mode 0 -v 1 >&2
        "$M" createsubdb "$d.keys" "${SRC}_h" "${d}_h" --subdb-mode 0 -v 1 >&2
    fi
    echo "$d"
}
# copy_db <src DB> <dst DB> <index kind>: the data files as links, with index
#   text: text index only; lazy: text index + parquet + binary file; notext: parquet file only
copy_db() {
    local s="$1" d="$2" kind="$3" f
    for f in "$s" "$s".[0-9]*; do [ -e "$f" ] && ln -sf "$(readlink -f "$f")" "$d${f#"$s"}"; done
    cp "$s.dbtype" "$d.dbtype"
    case "$kind" in
        text) cp "$s.index" "$d.index" ;;
        lazy) cp "$s.index" "$d.index"; MMSEQS_WRITE_LAZY_INDEX=binary "$M" convert2parquet "$d.index" "$d.index.parquet" -v 1 ;;
        notext) "$M" convert2parquet "$s.index" "$d.index.parquet" -v 1 ;;
    esac
}
if [ -n "$STEPS" ] && [ "$STEPS" != build ]; then
    prepare_db || { echo "database,1,0" >> "$STATUS"; STEPS=""; }
fi

# ---------------------------------------------------------------- correctness tests
test_modules() {
    "$U/lazy_index_test.sh" "$M" "$TDB" "$WORK/test" 300
}
test_access_modes() {
    local W="$WORK/modes" rc=0 src name b
    mkdir -p "$W"
    # a sequence DB, its header DB and a search result (written by many threads: offsets not in key order)
    for src in "$TDB" "${TDB}_h" "$WORK/test/res_text"; do
        [ -e "$src.index" ] || { echo "  FAIL  $src: DB for the access mode test is missing"; rc=1; continue; }
        name="$(basename "$src")"
        copy_db "$src" "$W/${name}_ref" text
        copy_db "$src" "$W/${name}_lazy" lazy
        copy_db "$src" "$W/${name}_notext" notext
        for b in parquet parquet-custom binary; do
            echo "== $name, text index against $b"
            MMSEQS_INDEX_BACKEND=$b "$BUILD_DIR/src/test/test_lazyindex" "$W/${name}_ref" "$W/${name}_lazy" "$THREADS" 2>&1 \
                | sed "s/^  \(OK\|FAIL\)  /  \1  $name, $b: /"
            [ "${PIPESTATUS[0]}" -eq 0 ] || rc=1
        done
        echo "== $name, text index against a DB with parquet index only"
        "$BUILD_DIR/src/test/test_lazyindex" "$W/${name}_ref" "$W/${name}_notext" "$THREADS" 2>&1 \
            | sed "s/^  \(OK\|FAIL\)  /  \1  $name, parquet only: /"
        [ "${PIPESTATUS[0]}" -eq 0 ] || rc=1
    done
    return "$rc"
}
if has test; then
    TDB="$(subset "$TEST_SEQS" testdb)"
    step test_modules test_modules; collect_checks modules "$LOGS/test_modules.log"
    step test_access_modes test_access_modes; collect_checks access_modes "$LOGS/test_access_modes.log"
    # parquet only against text in the cascaded cluster workflow and in linclust: same clusters
    for wf in cluster linclust; do
        step "test_notext_$wf" "$U/lazy_index_steps_bench.sh" "$M" "$TDB" "$WORK/test_notext_$wf" 0 1 "$THREADS" "$wf" "text notext"
        if [ -e "$WORK/test_notext_$wf/runs_$wf.csv" ]; then
            awk -F, -v wf="$wf" 'NR == 1 { next } $1 == "text" && ref == "" { ref = $12 }
                { r = ($14 == 0 && $12 == ref && ref != "failed") ? "OK" : "FAIL"
                  printf "notext_workflow,\"%s, %s: same clusters as text (%s clusters, %s lazy opens, %s text loads)\",%s\n", wf, $1, $11, $7, $8, r }' \
                "$WORK/test_notext_$wf/runs_$wf.csv" >> "$TESTS"
        else
            echo "notext_workflow,\"$wf: no result\",FAIL" >> "$TESTS"
        fi
    done
fi

# ---------------------------------------------------------------- benchmarks
if has lookup; then
    LDB="$(subset "$LOOKUP_SEQS" lookupdb)"
    step lookup "$U/lazy_index_bench.sh" "$M" "$LDB" "$WORK/lookup" "$ROUNDS" "$THREADS" 2000
    for f in lazy_index.csv index_files.csv other_readers.csv; do
        [ -e "$WORK/lookup/$f" ] && cp "$WORK/lookup/$f" "$CSVS/lookup_$f"
    done
fi
FIRST_SIZE="${CLUSTER_SEQS%% *}"
if has cluster; then
    CDB="$(subset "$FIRST_SIZE" "clusterdb_$FIRST_SIZE")"
    step "cluster_$FIRST_SIZE" "$U/lazy_index_cluster_bench.sh" "$M" "$CDB" "$WORK/cluster_$FIRST_SIZE" "$(wc -l < "$CDB.index")" "$ROUNDS" "$THREADS" linclust
    [ -e "$WORK/cluster_$FIRST_SIZE/cluster_linclust.csv" ] && cp "$WORK/cluster_$FIRST_SIZE/cluster_linclust.csv" "$CSVS/cluster_${FIRST_SIZE}_linclust.csv"
fi
if has steps; then
    for n in $CLUSTER_SEQS; do
        CDB="$(subset "$n" "clusterdb_$n")"
        for wf in $WORKFLOWS; do
            step "steps_${n}_$wf" "$U/lazy_index_steps_bench.sh" "$M" "$CDB" "$WORK/steps_${n}_$wf" 0 "$ROUNDS" "$THREADS" "$wf" "$STEP_VARIANTS" "$EXTRA_ARGS"
            for f in runs steps; do
                [ -e "$WORK/steps_${n}_$wf/${f}_$wf.csv" ] && cp "$WORK/steps_${n}_$wf/${f}_$wf.csv" "$CSVS/steps_${n}_${wf}_$f.csv"
            done
        done
    done
fi

# ---------------------------------------------------------------- one CSV and the summary
python3 -I - "$OUT" <<'PY'
import csv, glob, os, re, statistics as st, sys
out = sys.argv[1]; cdir = os.path.join(out, "csv")
rows = []   # suite, dataset, benchmark, variant, round, module, metric, value
def add(suite, dataset, bench, variant, rnd, module, metric, value):
    rows.append([suite, dataset, bench, variant, rnd, module, metric, value])
def read(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))

for r in read(os.path.join(cdir, "tests.csv")):
    add("test", "", r["group"], "", "", "", r["test"], r["result"])
p = os.path.join(cdir, "lookup_lazy_index.csv")
if os.path.exists(p):
    for r in read(p):
        for m in ("wall_s", "cpu_s", "peak_rss_mb", "anon_mb", "other_cores"):
            add("lookup", "", r["workload"], r["backend"], r["round"], "", m, r[m])
p = os.path.join(cdir, "lookup_index_files.csv")
if os.path.exists(p):
    for r in read(p):
        for m in ("write_s", "bytes"):
            if r[m] != "":
                add("lookup_index_files", "", "index_file", r["file"], "", "", m, r[m])
p = os.path.join(cdir, "lookup_other_readers.csv")
if os.path.exists(p):
    for r in read(p):
        add("lookup_other_readers", "", r["query"], r["reader"], "", "", "wall_s", r["wall_s"])
for p in sorted(glob.glob(os.path.join(cdir, "cluster_*_linclust.csv"))):
    n = re.match(r"cluster_(\d+)_", os.path.basename(p)).group(1)
    for r in read(p):
        for m in r:
            if m not in ("variant", "round"):
                add("cluster", n, "linclust", r["variant"], r["round"], "workflow", m, r[m])
for p in sorted(glob.glob(os.path.join(cdir, "steps_*_runs.csv"))):
    n, wf = re.match(r"steps_(\d+)_(.+)_runs\.csv", os.path.basename(p)).groups()
    for r in read(p):
        for m in r:
            if m not in ("variant", "round"):
                add("steps", n, wf, r["variant"], r["round"], "workflow", m, r[m])
    sp = p.replace("_runs.csv", "_steps.csv")
    if os.path.exists(sp):
        seen = {}
        for r in sorted(read(sp), key=lambda r: (r["variant"], int(r["round"]), float(r["start"]))):
            # the n-th call of a module in a run is one step: kmermatcher, kmermatcher#2
            k = (r["variant"], r["round"], r["module"]); seen[k] = seen.get(k, 0) + 1
            module = r["module"] + ("#" + str(seen[k]) if seen[k] > 1 else "")
            add("steps", n, wf, r["variant"], r["round"], module, "wall_s", r["wall_s"])
            add("steps", n, wf, r["variant"], r["round"], module, "cpu_s", f'{float(r["user_s"]) + float(r["sys_s"]):.2f}')
            add("steps", n, wf, r["variant"], r["round"], module, "peak_rss_mb", f'{int(r["peak_rss_kb"]) / 1024:.0f}')
            add("steps", n, wf, r["variant"], r["round"], module, "anon_mb", f'{int(r["anon_kb"]) / 1024:.0f}')
            add("steps", n, wf, r["variant"], r["round"], module, "exit", r["exit"])
with open(os.path.join(out, "lazy_index_results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["suite", "dataset_sequences", "benchmark", "variant", "round", "module", "metric", "value"])
    w.writerows(rows)

def num(v):
    try:
        return float(v)
    except ValueError:
        return None
L = []
L.append("Parquet index for MMseqs2: tests and benchmarks")
L.append("=" * 47)
L += open(os.path.join(out, "environment.txt")).read().splitlines()
L.append("")
L.append("Steps (exit code, seconds)")
status = read(os.path.join(out, "steps_status.csv"))
for r in status:
    L.append(f'  {r["step"]:28} exit {r["exit"]:>3}  {r["seconds"]:>7} s')
tests = [r for r in rows if r[0] == "test"]
failed = [r for r in tests if r[7] != "OK"]
L.append("")
L.append(f"Correctness: {len(tests) - len(failed)} of {len(tests)} checks passed")
for g in dict.fromkeys(r[2] for r in tests):
    x = [r for r in tests if r[2] == g]
    L.append(f'  {g:18} {sum(r[7] == "OK" for r in x)} OK, {sum(r[7] != "OK" for r in x)} failed')
for r in failed:
    L.append(f"  FAIL  {r[2]}: {r[6]}")
# median over the rounds per (suite, dataset, benchmark, module), one column per variant, ratio to the first variant
def table(suite, metrics, title):
    sel = [r for r in rows if r[0] == suite and r[6] in metrics and num(r[7]) is not None]
    if not sel:
        return
    L.append("")
    L.append(title)
    for ds, bench_prefix in dict.fromkeys((r[1], r[2] if suite != "lookup" else "") for r in sel):
        grp = [r for r in sel if r[1] == ds and (suite == "lookup" or r[2] == bench_prefix)]
        variants = list(dict.fromkeys(r[3] for r in grp))
        if bench_prefix or ds:
            L.append(f'  {(bench_prefix + ", " if bench_prefix else "")}{ds + " sequences" if ds else ""}'.rstrip(", "))
        for metric in metrics:
            keys = list(dict.fromkeys((r[2], r[5]) for r in grp if r[6] == metric))
            if not keys:
                continue
            L.append(f'    {metric + " (median)":34}' + "".join(f"{v:>22}" for v in variants) + f"   x {variants[0]}")
            for b, mod in keys:
                med = []
                for v in variants:
                    x = [num(r[7]) for r in grp if r[6] == metric and r[2] == b and r[5] == mod and r[3] == v]
                    med.append(st.median(x) if x else None)
                name = b if suite == "lookup" else mod
                line = f"    {name:34}" + "".join(f"{m:22.2f}" if m is not None else f'{"-":>22}' for m in med)
                if med[0]:
                    line += "   " + " ".join(f"{m / med[0]:5.2f}" if m is not None else "    -" for m in med[1:])
                L.append(line)
table("lookup", ["wall_s", "anon_mb"], "Lookup benchmark (one process per workload)")
table("cluster", ["wall_s", "cpu_s", "anon_mb", "big_text_opens"], "Clustering workflow, four readers")
table("steps", ["wall_s", "cpu_s", "anon_mb"], "Workflow per module (workflow = whole run; anon_mb = peak allocated memory)")
for suite in ("cluster", "steps"):
    res = [r for r in rows if r[0] == suite and r[6] == "result"]
    for ds, b in dict.fromkeys((r[1], r[2]) for r in res):
        x = [r for r in res if r[1] == ds and r[2] == b]
        kinds = sorted(set(r[7] for r in x))
        L.append("")
        L.append(f'{suite}, {b}, {ds} sequences: ' + ("all runs gave the same clusters" if len(kinds) == 1 and "failed" not in kinds[0]
                 else "DIFFERENT results: " + ", ".join(f'{r[3]} round {r[4]}: {r[7]}' for r in x)))
        other = [num(r[7]) for r in rows if r[0] == suite and r[1] == ds and r[2] == b and r[6] == "other_cores" and num(r[7]) is not None]
        if other:
            L.append(f"  CPU used by other processes during the runs: at most {max(other):.1f} cores" + ("  <- disturbed" if max(other) > 1 else ""))
bad = [r for r in status if r["exit"] != "0"]
L.append("")
L.append("RESULT: " + ("all steps ran, all checks passed" if not bad and not failed else
                       f'{len(bad)} steps with errors ({", ".join(r["step"] for r in bad)}), {len(failed)} failed checks'))
L.append(f"All values: {os.path.join(out, 'lazy_index_results.csv')} (suite, dataset_sequences, benchmark, variant, round, module, metric, value)")
L.append(f"Logs of the steps with their own tables: {os.path.join(out, 'logs')}")
open(os.path.join(out, "lazy_index_summary.txt"), "w").write("\n".join(L) + "\n")
print("\n".join(L))
sys.exit(1 if bad or failed else 0)
PY
RC=$?
if [ "$CLEAN" = 1 ] && [ -d "$WORK" ]; then
    log "removing $WORK"
    rm -rf "$WORK"
fi
log "results: $OUT/lazy_index_results.csv, $OUT/lazy_index_summary.txt"
exit "$RC"
