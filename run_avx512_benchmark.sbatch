#!/bin/bash
# =============================================================================
# Example SLURM submission for the MMseqs2 AVX2-vs-AVX512 benchmark harness.
# Submit with:   sbatch run_avx512_benchmark.sbatch
# Edit the SBATCH directives, module loads, and MMSEQS_DIR for your cluster.
# =============================================================================
#SBATCH --job-name=mmseqs_avx512
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --constraint=avx512            # request an AVX512-capable node (site-specific tag)
#SBATCH --output=%x-%j.out             # <jobname>-<jobid>.out
#SBATCH --error=%x-%j.err

set -euo pipefail

# --- 1. point at the MMseqs2 source tree (separate from $SCRATCH) ------------
# The repo can live anywhere (home/project space); results go to $SCRATCH.
export MMSEQS_DIR="${MMSEQS_DIR:-$HOME/src/MMseqs2}"

# --- 2. $SCRATCH is the working dir; the harness auto-uses $SCRATCH/mmseqs_avx_bench.
#     Fall back to the job tmp dir if the site does not define $SCRATCH.
export SCRATCH="${SCRATCH:-${SLURM_TMPDIR:-/tmp/$USER}}"

# --- 3. toolchain (cluster-specific; adjust module names or delete) ----------
# module load gcc cmake llvm rust    2>/dev/null || true
# If block-aligner's Rust is too old system-wide, a user rustup works too:
#   curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y; . "$HOME/.cargo/env"

# THREADS/JOBS auto-derive from $SLURM_CPUS_PER_TASK inside the harness.
echo "node=$(hostname)  MMSEQS_DIR=$MMSEQS_DIR  SCRATCH=$SCRATCH  cpus=${SLURM_CPUS_PER_TASK:-?}"

# --- 4. run -------------------------------------------------------------------
srun bash "$MMSEQS_DIR/run_avx512_benchmark.sh"

echo "Done. Results under $SCRATCH/mmseqs_avx_bench (timings.csv, perf*, asm/, logs/)."
