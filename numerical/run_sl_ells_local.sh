#!/usr/bin/env bash
set -euo pipefail

# Local parallel runner for sl_solve_ell_parallel.
# Usage:
#   ./run_sl_ells_local.sh ell_min ell_max nproc prefix output_h5
# Example:
#   ./run_sl_ells_local.sh 0 80 6 sl_part sl_spectrum_merged.h5

ELL_MIN=${1:-0}
ELL_MAX=${2:-5}
NPROC=${3:-4}
PREFIX=${4:-sl_part}
OUTPUT=${5:-sl_spectrum_merged.h5}

FC=${FC:-h5fc}
FFLAGS=${FFLAGS:--O3 -march=native}

mkdir -p logs

echo "Compiling sl_solve_ell_parallel.f90 ..."
$FC $FFLAGS SL_metric_grid.f90 -llapack -lblas -o sl_solve_ell_parallel

echo "Running ell=${ELL_MIN}...${ELL_MAX} with NPROC=${NPROC}"
seq "$ELL_MIN" "$ELL_MAX" | xargs -I{} -P "$NPROC" sh -c '
    ell="$1"
    prefix="$2"
    echo "[ell=${ell}] start"
    ./sl_solve_ell_parallel "$ell" "$ell" "$prefix" > "logs/${prefix}_ell_${ell}.log" 2>&1
    echo "[ell=${ell}] done"
' _ {} "$PREFIX"

echo "Merging files into ${OUTPUT} ..."
python3 merge_sl_ell_files.py --prefix "$PREFIX" --ell-min "$ELL_MIN" --ell-max "$ELL_MAX" --output "$OUTPUT" --overwrite

echo "Done: ${OUTPUT}"
