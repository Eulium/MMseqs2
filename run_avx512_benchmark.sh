#!/usr/bin/env bash
# =============================================================================
# MMseqs2 AVX2 vs AVX512 reproduction + analysis harness  (branch: AI_Agent)
#
# Recreates the benchmark dataset, builds AVX2 / AVX512 / AVX512-without-P3,
# runs runtime benchmarks (search, align=striped-SW, fwbw@blocklen {16,32,64}),
# checks AVX512-vs-AVX2 output accuracy (delta 0.01), and produces perf + llvm-mca
# analyses. Designed to run on a real x86-64 AVX512 Linux machine.
#
# What P3 is: the striped-SW kernels (StripedSmithWaterman.cpp) are compiled at
# 256-bit width with AVX512VL mask features even in the AVX512 build (to dodge the
# 512-bit striped-SW penalties). Toggle with -DMMSEQS_NO_SIMD256. This script
# builds "avx512" (P3 on) and "avx512_noP3" (P3 off) so you can measure P3 directly.
#
# Config is via environment variables (SLURM-friendly):
#   MMSEQS_DIR  path to the MMseqs2 source tree (the repo; may be anywhere)   [required*]
#   SCRATCH     fast work filesystem; WORK defaults to $SCRATCH/mmseqs_avx_bench
#   WORK        override the working dir (builds/data/results) explicitly
#   THREADS/JOBS default to $SLURM_CPUS_PER_TASK when set
#   REPEATS FWBW_REPS SENS QUERY_N MCPU SPROT_URL  (see defaults below)
#  *REPO/MMSEQS_DIR falls back to this script's own dir only if that dir is a repo.
#
# Standalone:  MMSEQS_DIR=~/src/MMseqs2 ./run_avx512_benchmark.sh
# SLURM (sbatch):
#   #!/bin/bash
#   #SBATCH -c 32
#   #SBATCH --constraint=avx512          # request an AVX512 node
#   #SBATCH -t 04:00:00
#   export MMSEQS_DIR=$HOME/src/MMseqs2   # repo, separate from scratch
#   srun bash "$MMSEQS_DIR/run_avx512_benchmark.sh"   # WORK auto = $SCRATCH/mmseqs_avx_bench
#
# Requires: gcc/g++, cmake, make, cargo/rustc (>=1.78 for block-aligner), curl,
#           awk, python3, perf, llvm-mca (clang optional). sudo only (optional) to
#           relax perf_event_paranoid; skipped if unavailable (common on clusters).
# =============================================================================
set -euo pipefail

# --- repo location: explicit env (MMSEQS_DIR or REPO), else this script's dir ---
_self_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="${MMSEQS_DIR:-${REPO:-$_self_dir}}"
# --- working dir: prefer $SCRATCH (HPC), then $WORK, then $HOME ---
WORK="${WORK:-${SCRATCH:+$SCRATCH/mmseqs_avx_bench}}"; WORK="${WORK:-$HOME/mmseqs_avx_bench}"
THREADS="${THREADS:-${SLURM_CPUS_PER_TASK:-16}}"
JOBS="${JOBS:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
REPEATS="${REPEATS:-5}"
FWBW_REPS="${FWBW_REPS:-3}"
SENS="${SENS:-7.5}"
QUERY_N="${QUERY_N:-5000}"
MCPU="${MCPU:-native}"                              # llvm-mca / -march target
SPROT_URL="${SPROT_URL:-https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz}"

# --- validate ---
[ -f "$REPO/CMakeLists.txt" ] && [ -f "$REPO/src/mmseqs.cpp" ] || {
  echo "ERROR: MMSEQS_DIR does not point to an MMseqs2 source tree: '$REPO'"
  echo "       set MMSEQS_DIR=/path/to/MMseqs2 (repo may live anywhere, separate from \$SCRATCH)"; exit 1; }
mkdir -p "$WORK" 2>/dev/null || { echo "ERROR: cannot create WORK='$WORK' (set \$SCRATCH or \$WORK)"; exit 1; }
[ -w "$WORK" ] || { echo "ERROR: WORK not writable: '$WORK'"; exit 1; }
REPO="$(cd "$REPO" && pwd)"; WORK="$(cd "$WORK" && pwd)"   # normalize to absolute
echo "REPO=$REPO"; echo "WORK=$WORK  THREADS=$THREADS JOBS=$JOBS"

DATA="$WORK/data"; RES="$WORK/results"; LOG="$WORK/logs"; ASM="$WORK/asm"
mkdir -p "$DATA" "$RES" "$LOG" "$ASM"
say(){ printf '\n\033[1m== %s ==\033[0m\n' "$*"; }
have(){ command -v "$1" >/dev/null 2>&1; }

# ---- 0. prereqs -------------------------------------------------------------
say "0. Environment"
grep -o 'avx512[a-z]*' /proc/cpuinfo | sort -u | tr '\n' ' '; echo
for t in cmake make g++ cargo python3 curl awk; do have "$t" || echo "  WARNING: '$t' not found"; done
LLVMMCA=""; for c in llvm-mca llvm-mca-18 llvm-mca-17 llvm-mca-16 llvm-mca-15; do have "$c" && { LLVMMCA="$c"; break; }; done
[ -z "$LLVMMCA" ] && echo "  note: llvm-mca not found -> static analysis section skipped"
PERF=""; have perf && PERF=perf || echo "  note: perf not found -> profiling section skipped"
sudo sysctl -w kernel.perf_event_paranoid=1 >/dev/null 2>&1 || true

# ---- 1. builds --------------------------------------------------------------
build(){  # name  <extra cmake args...>
  local name="$1"; shift
  local bdir="$WORK/build_$name"
  [ -x "$bdir/src/mmseqs" ] && { echo "  [$name] already built"; return; }
  say "Build $name ($*)"; mkdir -p "$bdir"
  ( cd "$bdir" && cmake -DCMAKE_BUILD_TYPE=Release -DHAVE_TESTS=0 -DENABLE_WERROR=0 "$@" "$REPO" \
      > "$LOG/cmake_$name.log" 2>&1 && make -j"$JOBS" mmseqs > "$LOG/build_$name.log" 2>&1 ) \
    && echo "  [$name] OK" || { echo "  [$name] BUILD FAILED (see $LOG/build_$name.log)"; tail -15 "$LOG/build_$name.log"; exit 1; }
}
build avx2         -DHAVE_AVX2=1
build avx512       -DHAVE_AVX512=1 -DCMAKE_CXX_FLAGS="-mavx512dq -mavx512vl" -DCMAKE_C_FLAGS="-mavx512dq -mavx512vl"
build avx512_noP3  -DHAVE_AVX512=1 -DCMAKE_CXX_FLAGS="-mavx512dq -mavx512vl -DMMSEQS_NO_SIMD256" -DCMAKE_C_FLAGS="-mavx512dq -mavx512vl"
# P4 (ungapped-scan prefetch, opt-in): avx512 + P3 + prefetch. Its effect shows in `search`
# (the ungapped diagonal scan runs inside prefilter); compare avx512_pf vs avx512 on search.
build avx512_pf    -DHAVE_AVX512=1 -DCMAKE_CXX_FLAGS="-mavx512dq -mavx512vl -DMMSEQS_UNGAPPED_PREFETCH" -DCMAKE_C_FLAGS="-mavx512dq -mavx512vl"
MM_AVX2="$WORK/build_avx2/src/mmseqs"
VARIANTS="avx2 avx512 avx512_noP3 avx512_pf"

# ---- 2. dataset -------------------------------------------------------------
say "2. Dataset (UniProt Swiss-Prot target + $QUERY_N-query sample)"
sprot="$DATA/uniprot_sprot.fasta"
[ -f "$sprot" ] || { curl -fsSL "$SPROT_URL" -o "$DATA/sp.gz" && gunzip -f "$DATA/sp.gz" && mv "$DATA/uniprot_sprot.fasta" "$sprot" 2>/dev/null || mv "$DATA/sp" "$sprot" 2>/dev/null || true; }
[ -f "$sprot" ] || { echo "download failed"; exit 1; }
q="$DATA/query.fasta"
if [ ! -f "$q" ]; then
  total=$(grep -c '^>' "$sprot"); step=$(( total / QUERY_N )); (( step<1 )) && step=1
  awk -v s="$step" 'BEGIN{n=-1}/^>/{n++;k=(n%s==0)}k' "$sprot" > "$q"
fi
echo "  queries=$(grep -c '^>' "$q")  targets=$(grep -c '^>' "$sprot")"
[ -f "$DATA/targetDB.dbtype" ] || "$MM_AVX2" createdb "$sprot" "$DATA/targetDB" -v1 >/dev/null
[ -f "$DATA/queryDB.dbtype"  ] || "$MM_AVX2" createdb "$q"     "$DATA/queryDB"  -v1 >/dev/null
# inputs for align (prefilter result) and fwbw (alignment result subset)
[ -f "$DATA/prefDB.dbtype" ] || "$MM_AVX2" prefilter "$DATA/queryDB" "$DATA/targetDB" "$DATA/prefDB" -s "$SENS" --threads "$THREADS" -v1 >/dev/null
if [ ! -f "$DATA/aln200.dbtype" ]; then
  rm -rf "$WORK/t"; "$MM_AVX2" search "$DATA/queryDB" "$DATA/targetDB" "$DATA/fullaln" "$WORK/t" -s "$SENS" --threads "$THREADS" -v1 >/dev/null; rm -rf "$WORK/t"
  head -200 "$DATA/fullaln.index" | awk '{print $1}' > "$WORK/k200"; "$MM_AVX2" createsubdb "$WORK/k200" "$DATA/fullaln" "$DATA/aln200" >/dev/null 2>&1
fi

# ---- 3. runtime benchmarks --------------------------------------------------
CSV="$RES/timings.csv"; echo "variant,bench,param,rep,wall_s,max_rss_kb" > "$CSV"
timed(){ # variant bench param  cmd...
  local v="$1" b="$2" p="$3" rep="$4"; shift 4
  /usr/bin/time -f "%e %M" -o "$RES/.t" "$@" >"$LOG/${v}_${b}_${p}_${rep}.log" 2>&1
  read -r w m < "$RES/.t"; echo "$v,$b,$p,$rep,$w,$m" >> "$CSV"
  printf "  %-13s %-6s %-4s r%s  %8ss\n" "$v" "$b" "$p" "$rep" "$w"
}
say "3. Runtime benchmarks (${THREADS} threads, ${REPEATS} reps; fwbw ${FWBW_REPS} reps)"
for v in $VARIANTS; do
  MM="$WORK/build_$v/src/mmseqs"
  for r in $(seq 1 "$REPEATS"); do
    rm -rf "$WORK/ts"; timed "$v" search - "$r" "$MM" search "$DATA/queryDB" "$DATA/targetDB" "$WORK/o_s_$v" "$WORK/ts" -s "$SENS" --threads "$THREADS" -v1; rm -rf "$WORK/o_s_$v"* "$WORK/ts"
    timed "$v" align - "$r" "$MM" align "$DATA/queryDB" "$DATA/targetDB" "$DATA/prefDB" "$WORK/o_a_$v" --threads "$THREADS" -v1; rm -f "$WORK/o_a_$v"*
  done
  for bl in 16 32 64; do for r in $(seq 1 "$FWBW_REPS"); do
    timed "$v" fwbw "$bl" "$r" "$MM" fwbw "$DATA/queryDB" "$DATA/targetDB" "$DATA/aln200" "$WORK/o_f_${v}_$bl" --blocklen "$bl" --threads "$THREADS" -v1; rm -f "$WORK/o_f_${v}_$bl"*
  done; done
done

# ---- 4. accuracy: AVX512 vs AVX2 (delta 0.01) -------------------------------
say "4. Accuracy (AVX512 vs AVX2, tolerance 0.01)"
acc(){ # bench cmdA... :: cmdB...   compares two result DBs numerically
  :; }
check_pair(){ # label  outA_prefix outB_prefix
  python3 - "$1" "$2" "$3" <<'PY'
import sys,glob,struct,re
lbl,a,b=sys.argv[1],sys.argv[2],sys.argv[3]
def blob(pfx):
    d=b''
    for f in sorted(glob.glob(pfx+'.*')):
        if re.search(r'\.(index|dbtype|lookup|source)$',f): continue
        d+=open(f,'rb').read()
    return d
da,db=blob(a),blob(b)
if da==db: print(f"  {lbl}: IDENTICAL ({len(da)} bytes)"); sys.exit(0)
# fall back to float32 tolerance
fa=[struct.unpack('<f',da[i:i+4])[0] for i in range(0,len(da)-3,4)]
fb=[struct.unpack('<f',db[i:i+4])[0] for i in range(0,len(db)-3,4)]
import math; n=min(len(fa),len(fb))
md=max((abs(fa[i]-fb[i]) for i in range(n) if math.isfinite(fa[i]) and math.isfinite(fb[i])), default=0.0)
print(f"  {lbl}: maxAbsDiff={md:.6g} -> {'PASS' if md<=0.01 else 'FAIL'}  (sizes {len(da)} vs {len(db)})")
PY
}
rm -rf "$WORK/ts"
"$WORK/build_avx2/src/mmseqs"   search "$DATA/queryDB" "$DATA/targetDB" "$WORK/acc_s_a2" "$WORK/ts" -s "$SENS" --threads "$THREADS" -v1 >/dev/null; rm -rf "$WORK/ts"
"$WORK/build_avx512/src/mmseqs" search "$DATA/queryDB" "$DATA/targetDB" "$WORK/acc_s_a5" "$WORK/ts" -s "$SENS" --threads "$THREADS" -v1 >/dev/null; rm -rf "$WORK/ts"
"$WORK/build_avx2/src/mmseqs"   align "$DATA/queryDB" "$DATA/targetDB" "$DATA/prefDB" "$WORK/acc_a_a2" --threads "$THREADS" -v1 >/dev/null
"$WORK/build_avx512/src/mmseqs" align "$DATA/queryDB" "$DATA/targetDB" "$DATA/prefDB" "$WORK/acc_a_a5" --threads "$THREADS" -v1 >/dev/null
"$WORK/build_avx2/src/mmseqs"   fwbw "$DATA/queryDB" "$DATA/targetDB" "$DATA/aln200" "$WORK/acc_f_a2" --blocklen 64 --threads "$THREADS" -v1 >/dev/null
"$WORK/build_avx512/src/mmseqs" fwbw "$DATA/queryDB" "$DATA/targetDB" "$DATA/aln200" "$WORK/acc_f_a5" --blocklen 64 --threads "$THREADS" -v1 >/dev/null
check_pair "search" "$WORK/acc_s_a2" "$WORK/acc_s_a5"
check_pair "align (striped SW / P3)" "$WORK/acc_a_a2" "$WORK/acc_a_a5"
check_pair "fwbw@64" "$WORK/acc_f_a2" "$WORK/acc_f_a5"
rm -f "$WORK"/acc_*

# ---- 5. perf hotspots + counters --------------------------------------------
if [ -n "$PERF" ]; then
  say "5. perf (hotspots + HW counters if available)"
  # small query subset for a quick, representative single-thread profile
  awk 'BEGIN{n=-1}/^>/{n++}n<500{print}n>=500{exit}' "$q" > "$DATA/profq.fasta"
  [ -f "$DATA/profqDB.dbtype" ] || "$MM_AVX2" createdb "$DATA/profq.fasta" "$DATA/profqDB" -v1 >/dev/null
  for v in avx2 avx512; do
    MM="$WORK/build_$v/src/mmseqs"; rm -rf "$WORK/tp"
    echo "  -- $v: perf stat (search, 1 thread) --"
    $PERF stat -o "$RES/perfstat_$v.txt" -- "$MM" search "$DATA/profqDB" "$DATA/targetDB" "$WORK/o_p_$v" "$WORK/tp" -s "$SENS" --threads 1 -v1 >/dev/null 2>&1 || true
    rm -rf "$WORK/o_p_$v"* "$WORK/tp"
    grep -E 'cycles|instructions|GHz|insn per cycle|elapsed' "$RES/perfstat_$v.txt" 2>/dev/null | sed 's/^/    /' || echo "    (HW counters unavailable on this host)"
    rm -rf "$WORK/tp"
    echo "  -- $v: perf record hotspots (align, 1 thread) --"
    $PERF record -e cpu-clock -F 999 -g -o "$RES/perf_$v.data" -- "$MM" align "$DATA/queryDB" "$DATA/targetDB" "$DATA/prefDB" "$WORK/o_pa_$v" --threads 1 -v1 >/dev/null 2>&1 || true
    $PERF report -i "$RES/perf_$v.data" --stdio -g none --percent-limit 2 2>/dev/null | grep -E '^\s+[0-9]+\.' | head -10 | sed 's/^/    /' || true
    rm -f "$WORK/o_pa_$v"*
  done
fi

# ---- 6. llvm-mca static analysis of the SW inner loop -----------------------
if [ -n "$LLVMMCA" ] && have clang++; then
  say "6. llvm-mca (${MCPU}) — striped-SW inner loop: AVX2 vs AVX512(512-bit) vs AVX512VL(256-bit,P3)"
  # emit assembly of the SW TU per width and let llvm-mca analyze the hottest region.
  # (For a focused number, we also mca a representative lazy-F micro-kernel below.)
  cat > "$ASM/lazyf.cpp" <<'EOF'
#include <immintrin.h>
#if defined(K512)
bool f(__m512i vF,__m512i vH){ return _mm512_cmp_epi16_mask(vF,vH,_MM_CMPINT_NLE)!=0; }
#elif defined(KVL256)
bool f(__m256i vF,__m256i vH){ return _mm256_cmp_epi16_mask(vF,vH,_MM_CMPINT_NLE)!=0; }
#else
bool f(__m256i vF,__m256i vH){ int m=_mm256_movemask_epi8(_mm256_cmpgt_epi16(vF,vH)); return m!=0; }
#endif
EOF
  for tag in AVX2:  K512:-DK512 KVL256:-DKVL256; do
    name="${tag%%:*}"; def="${tag#*:}"
    clang++ -O2 -march="$MCPU" $def -S -D__MM_MALLOC_H "$ASM/lazyf.cpp" -o "$ASM/lazyf_$name.s" 2>/dev/null \
      && { printf "  lazy-F exit [%s]: " "$name"; $LLVMMCA -mcpu="$MCPU" "$ASM/lazyf_$name.s" 2>/dev/null | grep -E 'Block RThroughput|Total uOps'|tr '\n' ' '; echo; }
  done
  echo "  (full per-ISA SW/fwbw assembly emitted under $ASM for deeper inspection)"
  for v in "" "-mavx512f -mavx512bw -mavx512dq -mavx512vl"; do :; done
fi

# ---- 7. report --------------------------------------------------------------
say "7. SUMMARY (median wall seconds; speedup = AVX2 / variant)"
python3 - "$CSV" <<'PY'
import csv,sys,statistics as st
rows=list(csv.DictReader(open(sys.argv[1])))
key=lambda r:(r['bench'],r['param'],r['variant'])
from collections import defaultdict
d=defaultdict(list)
for r in rows: d[(r['bench'],r['param'],r['variant'])].append(float(r['wall_s']))
benches=sorted({(r['bench'],r['param']) for r in rows})
print(f"{'bench':<10}{'param':<6}{'avx2':>9}{'avx512':>9}{'noP3':>9}{'pf(P4)':>9}{'  512/2':>8}{' P3g':>7}{' P4g':>7}")
for b,p in benches:
    m={v:st.median(d[(b,p,v)]) for v in ('avx2','avx512','avx512_noP3','avx512_pf') if (b,p,v) in d}
    if 'avx2' not in m: continue
    a2=m['avx2']; a5=m.get('avx512'); n3=m.get('avx512_noP3'); pf=m.get('avx512_pf')
    sp = f"{a2/a5:.3f}x" if a5 else "-"
    p3 = f"{n3/a5:.3f}x" if (a5 and n3) else "-"   # noP3/P3: >1 => P3 helped (align)
    p4 = f"{a5/pf:.3f}x" if (a5 and pf) else "-"   # P3/(P3+prefetch): >1 => P4 helped (search)
    print(f"{b:<10}{p:<6}{a2:>9.2f}{(a5 or 0):>9.2f}{(n3 or 0):>9.2f}{(pf or 0):>9.2f}{sp:>8}{p3:>7}{p4:>7}")
print("\nLegend: 512/2=AVX2/AVX512 (>1 AVX512 faster). P3g=noP3/P3 (>1 P3 helped, align). P4g=P3/pf (>1 prefetch helped, search).")
PY
echo
echo "Raw: $CSV | perf: $RES/perf* | asm+mca: $ASM | logs: $LOG"
echo "Dataset + results preserved under $WORK"
