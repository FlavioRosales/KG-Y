#!/usr/bin/env bash
set -euo pipefail

# Fuerza punto decimal (evita coma por locale)
export LC_NUMERIC=C

# Archivo de plantilla
TEMPLATE="k_05.nml"

# Resoluciones a probar
NR_LIST=(8000 12000 16000)

# Clave en el namelist
NR_KEY="nr"

# MPI: fuerza shared memory y restringe TCP a loopback
export OMPI_MCA_btl="self,vader,tcp"
export OMPI_MCA_btl_tcp_if_include="lo"

for NR in "${NR_LIST[@]}"; do
    # Tag para el nombre del archivo
    tag=$(printf "%05d" "$NR")
    param_file="k_05_nr_${tag}.nml"

    echo ">>> Generando ${param_file} con ${NR_KEY} = ${NR}"

    # Sustituye solo la línea de nr en la plantilla (respeta indentación)
    sed -E "s/^([[:space:]]*${NR_KEY}[[:space:]]*=[[:space:]]*).*/\1${NR}/" \
        "$TEMPLATE" > "$param_file"

    echo ">>> Ejecutando simulación con ${param_file}"
    make run NPROC=8 PARAMS="$param_file"
    echo
done

echo "Todas las simulaciones han terminado."
