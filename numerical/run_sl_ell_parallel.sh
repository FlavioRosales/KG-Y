#!/usr/bin/env bash
set -euo pipefail

ELL_MIN=${1:-0}
ELL_MAX=${2:-80}
JOBS=${3:-4}

# Para paralelizar por ell, evita sobre-suscripción de BLAS/OpenMP.
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}
export VECLIB_MAXIMUM_THREADS=${VECLIB_MAXIMUM_THREADS:-1}
export NUMEXPR_NUM_THREADS=${NUMEXPR_NUM_THREADS:-1}

seq "${ELL_MIN}" "${ELL_MAX}" | xargs -n 1 -P "${JOBS}" ./sl_solver_one_ell