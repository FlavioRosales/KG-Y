#!/usr/bin/env python3
"""
Merge all per-ell Sturm-Liouville HDF5 files in the current folder.

Run inside the folder that contains the partial files:

    python3 merge_sl_folder.py

Expected partial-file structure:
    /r
    /params
    /ell_####/k
    /ell_####/modes
    /ell_####/n_save

The script scans *.h5 files in the current directory, keeps only files that
contain exactly one ell_#### group, sorts them by ell, and writes one merged file:

    sl_spectrum.h5

No command-line parameters are required.
"""

from __future__ import annotations

import re
from pathlib import Path

import h5py
import numpy as np

OUTPUT_NAME = "sl_spectrum.h5"
GRID_ATOL = 1.0e-13
ELL_GROUP_RE = re.compile(r"^ell_(\d{4,})$")


def ell_from_group(group_name: str) -> int:
    match = ELL_GROUP_RE.match(group_name)
    if match is None:
        raise ValueError(f"Invalid ell group name: {group_name}")
    return int(match.group(1))


def get_single_ell_group(path: Path) -> tuple[int, str] | None:
    """
    Return (ell, group_name) if the file is a per-ell file.
    Return None if it is not a per-ell file, for example a merged file.
    """
    try:
        with h5py.File(path, "r") as h5:
            groups = [key for key in h5.keys() if ELL_GROUP_RE.match(key)]
            if len(groups) != 1:
                return None
            if "r" not in h5 or "params" not in h5:
                return None
            group_name = groups[0]
            ell = ell_from_group(group_name)
            return ell, group_name
    except OSError:
        return None


def collect_partial_files(folder: Path) -> list[tuple[int, Path, str]]:
    candidates: list[tuple[int, Path, str]] = []

    for path in sorted(folder.glob("*.h5")):
        if path.name == OUTPUT_NAME:
            continue

        found = get_single_ell_group(path)
        if found is None:
            continue

        ell, group_name = found
        candidates.append((ell, path, group_name))

    if not candidates:
        raise FileNotFoundError(
            "No per-ell HDF5 files found in the current directory. "
            "Expected files with exactly one ell_#### group."
        )

    candidates.sort(key=lambda item: item[0])

    # Detect duplicate ell files before writing anything.
    seen: dict[int, Path] = {}
    for ell, path, _ in candidates:
        if ell in seen:
            raise RuntimeError(
                f"Duplicate files for ell={ell}: {seen[ell].name} and {path.name}. "
                "Remove one of them before merging."
            )
        seen[ell] = path

    return candidates


def read_reference_metadata(path: Path):
    with h5py.File(path, "r") as h5:
        r_ref = h5["r"][:]

        params_real = h5["params/real"][:]
        params_real_extra = (
            h5["params/real_extra"][:]
            if "real_extra" in h5["params"]
            else None
        )
        params_int = h5["params/integer"][:]

        Nr = int(params_int[0])
        N_intervals = int(params_int[1])

        # In the per-ell solver this entry is n_modes_cap.
        n_modes_cap = int(params_int[3]) if len(params_int) >= 4 else -1

    return r_ref, params_real, params_real_extra, Nr, N_intervals, n_modes_cap


def verify_grid(path: Path, r_ref: np.ndarray) -> None:
    with h5py.File(path, "r") as h5:
        r_here = h5["r"][:]

    if r_here.shape != r_ref.shape:
        raise ValueError(
            f"Grid shape mismatch in {path.name}: {r_here.shape} != {r_ref.shape}"
        )

    if not np.allclose(r_here, r_ref, rtol=0.0, atol=GRID_ATOL):
        max_diff = float(np.max(np.abs(r_here - r_ref)))
        raise ValueError(
            f"Grid values mismatch in {path.name}: max |dr| = {max_diff:.3e}"
        )


def main() -> None:
    folder = Path.cwd()
    output = folder / OUTPUT_NAME

    partials = collect_partial_files(folder)
    ell_values = [ell for ell, _, _ in partials]
    ell_min = min(ell_values)
    ell_max = max(ell_values)

    print("Found per-ell files:")
    for ell, path, group_name in partials:
        print(f"  ell={ell:4d}  group={group_name}  file={path.name}")

    r_ref, params_real, params_real_extra, Nr, N_intervals, n_modes_cap = read_reference_metadata(partials[0][1])

    # Allocate up to ell_max so indexing by ell is direct.
    n_save_by_ell = np.full(ell_max + 1, -1, dtype=np.int64)
    k_min_by_ell = np.full(ell_max + 1, np.nan, dtype=np.float64)
    k_max_by_ell = np.full(ell_max + 1, np.nan, dtype=np.float64)

    if output.exists():
        print(f"Removing existing output: {output.name}")
        output.unlink()

    with h5py.File(output, "w") as fout:
        fout.create_dataset("r", data=r_ref)

        gpar = fout.create_group("params")
        gpar.create_dataset("real", data=params_real)
        if params_real_extra is not None:
            gpar.create_dataset("real_extra", data=params_real_extra)

        gpar.create_dataset(
            "integer",
            data=np.array([Nr, N_intervals, ell_max, n_modes_cap, -1], dtype=np.int32),
        )
        gpar.create_dataset("ell_values", data=np.array(ell_values, dtype=np.int32))

        for ell, path, group_name in partials:
            verify_grid(path, r_ref)

            expected = f"ell_{ell:04d}"
            if group_name != expected:
                raise ValueError(
                    f"Group mismatch in {path.name}: found {group_name}, expected {expected}"
                )

            with h5py.File(path, "r") as fin:
                fin.copy(group_name, fout, name=expected)

                k = fin[f"{group_name}/k"][:]
                n_save_by_ell[ell] = len(k)
                if len(k) > 0:
                    k_min_by_ell[ell] = k[0]
                    k_max_by_ell[ell] = k[-1]

        gpar.create_dataset("n_save_by_ell", data=n_save_by_ell)
        gpar.create_dataset("k_min_by_ell", data=k_min_by_ell)
        gpar.create_dataset("k_max_by_ell", data=k_max_by_ell)

    print()
    print(f"Done. Wrote: {output.name}")
    print(f"Included ell range: {ell_min} ... {ell_max}")
    print(f"Number of ell files: {len(partials)}")


if __name__ == "__main__":
    main()
