#!/usr/bin/env bash
set -euo pipefail

export LC_NUMERIC=C

TEMPLATE="k_05.nml"

N=100
KMIN=0.0
KMAX=10.0

export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_btl_tcp_if_include="lo"

# Paso uniforme en doble precisión (via awk)
dk=$(awk -v kmin="$KMIN" -v kmax="$KMAX" -v n="$N" 'BEGIN{
  if (n<2) { print 0; exit }
  printf "%.17g", (kmax-kmin)/(n-1)
}')

for i in $(seq 0 $((N))); do
  # k0 = KMIN + i*dk (uniforme)
  k0=$(awk -v kmin="$KMIN" -v dk="$dk" -v i="$i" 'BEGIN{
    printf "%.12f", kmin + i*dk
  }')

  # NOMBRE SIN COLISIONES: usa el índice i (monótono y único)
  tag=$(printf "%03d" "$i")
  param_file="k_m22_${tag}.nml"

  echo ">>> [$i/$((N-1))] Generando ${param_file} con k0 = ${k0}"

  # Sustituye k0 en el namelist
  sed -E "s/^([[:space:]]*k0[[:space:]]*=[[:space:]]*).*/\1${k0}/" \
    "$TEMPLATE" > "$param_file"

  echo ">>> Ejecutando simulación con ${param_file}"
  make run NPROC=8 OMP_NUM_THREADS=1 PARAMS="$param_file"
  echo ">>> done"
done

echo "Todas las simulaciones han terminado."
