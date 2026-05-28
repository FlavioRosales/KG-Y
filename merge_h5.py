#!/usr/bin/env python3

import glob
import os
import sys
import h5py
from tqdm import tqdm

REQUIRED_DATASETS = ("t", "phi", "pi", "N", "F")

def list_mode_groups(fin):
    return [name for name, obj in fin.items()
            if isinstance(obj, h5py.Group) and name.startswith("mode_")]

def verify_mode(fin, fout, gname):
    gin = fin[gname]
    gout = fout[gname]

    for dset in REQUIRED_DATASETS:
        if gin[dset].shape != gout[dset].shape:
            raise RuntimeError(f"Shape distinta en {gname}/{dset}")
        if gin[dset].dtype != gout[dset].dtype:
            raise RuntimeError(f"Dtype distinto en {gname}/{dset}")

def main():
    in_paths = sorted(glob.glob("modes_rank*.h5"))
    if not in_paths:
        print("[ERROR] No encontré archivos modes_rank*.h5", file=sys.stderr)
        sys.exit(1)

    out_path = "modes.h5"
    if os.path.exists(out_path):
        print(f"[ERROR] {out_path} ya existe. Bórralo primero.", file=sys.stderr)
        sys.exit(1)

    total_files = len(in_paths)

    with h5py.File(out_path, "w") as fout:

        for idx, p in enumerate(tqdm(in_paths, desc="Merging ranks", unit="file")):

            with h5py.File(p, "r") as fin:

                # Copiar r solo una vez
                if "r" in fin and "r" not in fout:
                    fin.copy("r", fout)

                mode_groups = list_mode_groups(fin)

                # Copiar modos
                for gname in mode_groups:
                    if gname in fout:
                        raise RuntimeError(f"Duplicado: {gname}")
                    fin.copy(gname, fout)

                fout.flush()

                # Verificación ligera
                for gname in mode_groups:
                    verify_mode(fin, fout, gname)

            # Eliminar archivo tras verificación
            os.remove(p)

    print("\n[OK] Merge completo -> modes.h5")

if __name__ == "__main__":
    main()
