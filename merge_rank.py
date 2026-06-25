#!/usr/bin/env python3
"""Merge the per-rank HDF5 output files produced by KG-Y.

The code writes one HDF5 file per MPI rank.  In the Fortran layout used by
KG-Y, h5py sees the relevant datasets as

  /diagnostics/N, /diagnostics/F : (time, local_mode)
  /fields/<field>                : (time, local_mode, radial)
  /planes/<field>                : (time, angular, radial)

This script constructs a single, analysis-friendly file:

  * /modes is assembled and sorted by the global nmode index;
  * /diagnostics and /fields are concatenated along their mode axis;
  * /planes is summed element-by-element across ranks, because every rank
    stores a partial modal reconstruction of the same physical plane;
  * /grid and any other static top-level items are copied from rank 0.

The source rank files are deleted only after the merged file has been
written successfully and moved into its final location.
"""

from __future__ import annotations

import argparse
import re
import sys
from contextlib import ExitStack
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

import h5py
import numpy as np


MERGED_TOP_LEVEL = {"modes", "diagnostics", "fields", "planes"}
STATIC_PLANE_NAMES = {"t", "r", "theta", "varphi", "x", "y", "z"}


@dataclass
class RankInfo:
    path: Path
    nmode: np.ndarray
    start: int = -1

    @property
    def nlocal(self) -> int:
        return int(self.nmode.size)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge KG-Y per-rank HDF5 files in RUN_DIR. The merged file is "
            "written as RUN_DIR/RUN_DIR_NAME.h5 and rank files are deleted "
            "after a successful merge."
        )
    )
    parser.add_argument(
        "run_dir",
        help="Directory containing output_rank_*.h5 files.",
    )
    parser.add_argument(
        "--pattern",
        default="output_rank_*.h5",
        help="Glob used to locate rank files (default: %(default)s).",
    )
    parser.add_argument(
        "--mode-block",
        type=int,
        default=64,
        help="Number of modes copied per hyperslab block (default: %(default)s).",
    )
    parser.add_argument(
        "--time-block",
        type=int,
        default=1,
        help="Number of time samples copied per hyperslab block (default: %(default)s).",
    )
    parser.add_argument(
        "--plane-time-block",
        type=int,
        default=1,
        help="Number of plane snapshots accumulated at once (default: %(default)s).",
    )
    parser.add_argument(
        "--gzip",
        type=int,
        choices=range(0, 10),
        default=None,
        metavar="LEVEL",
        help="Optionally gzip-compress the merged datasets, with LEVEL from 0 to 9.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace RUN_DIR/RUN_DIR_NAME.h5 if it already exists.",
    )
    parser.add_argument(
        "--keep-ranks",
        action="store_true",
        help="Keep output_rank_*.h5 files after a successful merge.",
    )
    return parser.parse_args()


def rank_key(path: Path) -> tuple[int, str]:
    match = re.search(r"output_rank_(\d+)\.h5$", path.name)
    if match:
        return int(match.group(1)), path.name
    return sys.maxsize, path.name


def collect_rank_files(run_dir: Path, pattern: str, output: Path) -> list[Path]:
    if not run_dir.is_dir():
        raise NotADirectoryError(f"No existe el directorio de corrida: {run_dir}")

    output_resolved = output.resolve()
    files = []
    for path in run_dir.glob(pattern):
        if path.resolve() != output_resolved and path.is_file():
            files.append(path)

    ordered = sorted(files, key=rank_key)
    if not ordered:
        raise FileNotFoundError(
            f"No se encontraron archivos {pattern!r} en {run_dir}."
        )
    return ordered


def copy_attrs(src: h5py.Group | h5py.Dataset | h5py.File,
               dst: h5py.Group | h5py.Dataset | h5py.File) -> None:
    for key, value in src.attrs.items():
        dst.attrs[key] = value


def require_same_array(reference: np.ndarray, candidate: np.ndarray, label: str) -> None:
    if reference.shape != candidate.shape:
        raise ValueError(
            f"{label}: shape inconsistente {candidate.shape}; se esperaba {reference.shape}."
        )
    if np.issubdtype(reference.dtype, np.floating):
        equal = np.allclose(reference, candidate, rtol=0.0, atol=1.0e-13)
    else:
        equal = np.array_equal(reference, candidate)
    if not equal:
        raise ValueError(f"{label}: los valores no coinciden entre los ranks.")


def check_static_dataset(files: list[h5py.File], path: str) -> None:
    reference = files[0][path][...]
    for fin in files[1:]:
        require_same_array(reference, fin[path][...], f"{path} en {fin.filename}")


def create_dataset_like(
    parent: h5py.Group | h5py.File,
    name: str,
    source: h5py.Dataset,
    shape: tuple[int, ...] | None = None,
    *,
    chunks: tuple[int, ...] | None = None,
    gzip: int | None = None,
) -> h5py.Dataset:
    kwargs: dict[str, object] = {}
    if chunks is not None:
        kwargs["chunks"] = chunks
    if gzip is not None:
        kwargs["compression"] = "gzip"
        kwargs["compression_opts"] = gzip
        kwargs["shuffle"] = True

    if shape is None:
        data = source[...]
        dst = parent.create_dataset(name, data=data, dtype=source.dtype, **kwargs)
    else:
        dst = parent.create_dataset(name, shape=shape, dtype=source.dtype, **kwargs)

    copy_attrs(source, dst)
    return dst


def copy_static_dataset(
    files: list[h5py.File],
    source: h5py.Dataset,
    parent: h5py.Group | h5py.File,
    name: str,
) -> None:
    path = source.name
    for fin in files[1:]:
        require_same_array(source[...], fin[path][...], f"{path} en {fin.filename}")
    create_dataset_like(parent, name, source)


def group_with_attrs(parent: h5py.Group | h5py.File,
                     name: str,
                     source_group: h5py.Group) -> h5py.Group:
    group = parent.create_group(name)
    copy_attrs(source_group, group)
    return group


def modal_chunks(shape: tuple[int, ...], mode_block: int, time_block: int) -> tuple[int, ...]:
    """Chunk layout compatible with h5py's view of Fortran-written datasets."""
    if len(shape) == 1:
        return (min(mode_block, shape[0]),)
    if len(shape) >= 2:
        return (min(time_block, shape[0]), min(mode_block, shape[1]), *shape[2:])
    return shape


def copy_modal_dataset(
    files: list[h5py.File],
    infos: list[RankInfo],
    group_name: str,
    dataset_name: str,
    out_group: h5py.Group,
    total_modes: int,
    mode_block: int,
    time_block: int,
    gzip: int | None,
) -> None:
    source = files[0][group_name][dataset_name]
    if source.ndim < 2:
        raise ValueError(
            f"/{group_name}/{dataset_name} no tiene un eje de modos en la posición 1."
        )

    expected_shape = list(source.shape)
    expected_shape[1] = total_modes
    out_shape = tuple(expected_shape)

    for fin, info in zip(files, infos):
        current = fin[group_name][dataset_name]
        if current.ndim != source.ndim:
            raise ValueError(f"/{group_name}/{dataset_name}: número de dimensiones inconsistente.")
        if current.shape[0] != source.shape[0] or current.shape[2:] != source.shape[2:]:
            raise ValueError(
                f"/{group_name}/{dataset_name}: shape inconsistente en {fin.filename}: "
                f"{current.shape}."
            )
        if current.shape[1] != info.nlocal:
            raise ValueError(
                f"/{group_name}/{dataset_name}: el eje modal no coincide con /modes/nmode "
                f"en {fin.filename}."
            )

    chunks = modal_chunks(out_shape, mode_block, time_block) if gzip is not None else None
    destination = create_dataset_like(
        out_group,
        dataset_name,
        source,
        shape=out_shape,
        chunks=chunks,
        gzip=gzip,
    )

    nt = source.shape[0]
    for fin, info in zip(files, infos):
        current = fin[group_name][dataset_name]
        for t0 in range(0, nt, time_block):
            t1 = min(t0 + time_block, nt)
            for local0 in range(0, info.nlocal, mode_block):
                local1 = min(local0 + mode_block, info.nlocal)
                src_sel = (slice(t0, t1), slice(local0, local1)) + (slice(None),) * (source.ndim - 2)
                dst_sel = (
                    slice(t0, t1),
                    slice(info.start + local0, info.start + local1),
                ) + (slice(None),) * (source.ndim - 2)
                destination[dst_sel] = current[src_sel]


def merge_modes(
    files: list[h5py.File],
    infos: list[RankInfo],
    out_file: h5py.File,
) -> int:
    ref_group = files[0]["modes"]
    if "nmode" not in ref_group:
        raise KeyError("Falta /modes/nmode en el primer archivo.")

    all_nmode = np.concatenate([info.nmode for info in infos])
    if np.unique(all_nmode).size != all_nmode.size:
        raise ValueError("Hay índices nmode repetidos entre los archivos de rank.")

    global_nmode = np.sort(all_nmode)
    total_modes = int(global_nmode.size)

    # The current main.f90 assigns a contiguous, ascending nmode block to each rank.
    # Validate that invariant, then data can be copied through fast hyperslab slices.
    for info in infos:
        if info.nlocal == 0:
            raise ValueError(f"{info.path}: un rank sin modos no está soportado por este formato.")
        expected_local = np.arange(info.nmode[0], info.nmode[0] + info.nlocal, dtype=info.nmode.dtype)
        if not np.array_equal(info.nmode, expected_local):
            raise ValueError(
                f"{info.path}: /modes/nmode no es un bloque contiguo y ascendente. "
                "El script espera la distribución por bloques que usa el main actual."
            )
        info.start = int(np.searchsorted(global_nmode, info.nmode[0]))
        if not np.array_equal(global_nmode[info.start:info.start + info.nlocal], info.nmode):
            raise ValueError(f"{info.path}: nmode no se puede colocar de manera consistente.")

    out_group = group_with_attrs(out_file, "modes", ref_group)

    dataset_names = [
        name for name, obj in ref_group.items()
        if isinstance(obj, h5py.Dataset)
    ]
    for fin in files[1:]:
        names = {name for name, obj in fin["modes"].items() if isinstance(obj, h5py.Dataset)}
        if names != set(dataset_names):
            raise ValueError(f"/modes no contiene los mismos datasets en {fin.filename}.")

    for name in dataset_names:
        source = ref_group[name]
        if source.ndim == 1 and source.shape[0] == infos[0].nlocal:
            merged = np.empty(total_modes, dtype=source.dtype)
            for fin, info in zip(files, infos):
                current = fin["modes"][name]
                if current.shape != (info.nlocal,):
                    raise ValueError(f"/modes/{name}: shape inconsistente en {fin.filename}.")
                merged[info.start:info.start + info.nlocal] = current[...]
            destination = out_group.create_dataset(name, data=merged, dtype=source.dtype)
            copy_attrs(source, destination)
        else:
            copy_static_dataset(files, source, out_group, name)

    require_same_array(global_nmode, out_group["nmode"][...], "/modes/nmode reconstruido")
    return total_modes


def merge_modal_group(
    files: list[h5py.File],
    infos: list[RankInfo],
    out_file: h5py.File,
    group_name: str,
    time_name: str,
    total_modes: int,
    mode_block: int,
    time_block: int,
    gzip: int | None,
) -> None:
    if group_name not in files[0]:
        return
    if any(group_name not in fin for fin in files[1:]):
        raise ValueError(f"/{group_name} no está presente en todos los archivos de rank.")

    ref_group = files[0][group_name]
    out_group = group_with_attrs(out_file, group_name, ref_group)

    ref_names = {name for name, obj in ref_group.items() if isinstance(obj, h5py.Dataset)}
    for fin in files[1:]:
        current_names = {name for name, obj in fin[group_name].items() if isinstance(obj, h5py.Dataset)}
        if current_names != ref_names:
            raise ValueError(f"/{group_name} no contiene los mismos datasets en {fin.filename}.")

    for name in sorted(ref_names):
        source = ref_group[name]
        if name == time_name:
            copy_static_dataset(files, source, out_group, name)
        elif source.ndim >= 2 and source.shape[1] == infos[0].nlocal:
            copy_modal_dataset(
                files, infos, group_name, name, out_group, total_modes,
                mode_block, time_block, gzip,
            )
        else:
            copy_static_dataset(files, source, out_group, name)


def dataset_paths(group: h5py.Group, prefix: str = "") -> list[str]:
    paths: list[str] = []
    for name, obj in group.items():
        relative = f"{prefix}/{name}" if prefix else name
        if isinstance(obj, h5py.Dataset):
            paths.append(relative)
        elif isinstance(obj, h5py.Group):
            paths.extend(dataset_paths(obj, relative))
    return paths


def clone_group_tree(source: h5py.Group, destination: h5py.Group) -> None:
    copy_attrs(source, destination)
    for name, obj in source.items():
        if isinstance(obj, h5py.Group):
            child = destination.create_group(name)
            clone_group_tree(obj, child)


def destination_parent(root: h5py.Group, relative_path: str) -> tuple[h5py.Group, str]:
    pieces = relative_path.split("/")
    parent = root
    for piece in pieces[:-1]:
        parent = parent[piece]
    return parent, pieces[-1]


def merge_planes(
    files: list[h5py.File],
    out_file: h5py.File,
    plane_time_block: int,
    gzip: int | None,
) -> None:
    present = ["planes" in fin for fin in files]
    if not any(present):
        return
    if not all(present):
        raise ValueError("/planes está presente solo en una parte de los archivos de rank.")

    ref_group = files[0]["planes"]
    ref_paths = set(dataset_paths(ref_group))
    for fin in files[1:]:
        current_paths = set(dataset_paths(fin["planes"]))
        if current_paths != ref_paths:
            raise ValueError(f"/planes no contiene los mismos datasets en {fin.filename}.")

    out_group = out_file.create_group("planes")
    clone_group_tree(ref_group, out_group)

    for relative_path in sorted(ref_paths):
        source = ref_group[relative_path]
        parent, name = destination_parent(out_group, relative_path)
        base_name = Path(relative_path).name

        if base_name in STATIC_PLANE_NAMES:
            copy_static_dataset(files, source, parent, name)
            continue

        for fin in files[1:]:
            current = fin["planes"][relative_path]
            if current.shape != source.shape:
                raise ValueError(
                    f"/planes/{relative_path}: shape inconsistente en {fin.filename}: "
                    f"{current.shape}; se esperaba {source.shape}."
                )
            if current.dtype != source.dtype:
                raise ValueError(f"/planes/{relative_path}: dtype inconsistente en {fin.filename}.")

        shape = source.shape
        chunks = None
        if gzip is not None:
            chunks = (min(plane_time_block, shape[0]), *shape[1:]) if shape else None
        destination = create_dataset_like(
            parent,
            name,
            source,
            shape=shape,
            chunks=chunks,
            gzip=gzip,
        )

        if source.ndim == 0:
            value = np.zeros((), dtype=source.dtype)
            for fin in files:
                value[...] += fin["planes"][relative_path][()]
            destination[()] = value
            continue

        for t0 in range(0, shape[0], plane_time_block):
            t1 = min(t0 + plane_time_block, shape[0])
            selection = (slice(t0, t1),) + (slice(None),) * (source.ndim - 1)
            partial_sum = np.zeros((t1 - t0, *shape[1:]), dtype=source.dtype)
            for fin in files:
                partial_sum += fin["planes"][relative_path][selection]
            destination[selection] = partial_sum


def copy_static_top_level(files: list[h5py.File], out_file: h5py.File) -> None:
    reference = files[0]
    for name, obj in reference.items():
        if name in MERGED_TOP_LEVEL:
            continue
        reference.copy(name, out_file)


def merge(files_paths: list[Path], output: Path, args: argparse.Namespace) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    if output.exists() and not args.overwrite:
        raise FileExistsError(f"El archivo de salida ya existe: {output}. Usa --overwrite para reemplazarlo.")

    with ExitStack() as stack:
        files = [stack.enter_context(h5py.File(path, "r")) for path in files_paths]
        for fin in files:
            for required in ("grid", "modes", "diagnostics", "fields"):
                if required not in fin:
                    raise KeyError(f"{fin.filename}: falta el grupo /{required}.")

        check_static_dataset(files, "/grid/r")

        pairs = [
            (fin, RankInfo(path=path, nmode=fin["modes/nmode"][...]))
            for path, fin in zip(files_paths, files)
        ]
        pairs.sort(key=lambda item: int(item[1].nmode[0]))
        files = [item[0] for item in pairs]
        infos = [item[1] for item in pairs]

        with h5py.File(output, "w") as out_file:
            copy_attrs(files[0], out_file)
            out_file.attrs["merged_rank_count"] = len(files)
            out_file.attrs["source_rank_files"] = np.asarray(
                [info.path.name for info in infos], dtype=h5py.string_dtype(encoding="utf-8")
            )

            copy_static_top_level(files, out_file)
            total_modes = merge_modes(files, infos, out_file)
            merge_modal_group(
                files, infos, out_file, "diagnostics", "t", total_modes,
                args.mode_block, args.time_block, args.gzip,
            )
            merge_modal_group(
                files, infos, out_file, "fields", "t", total_modes,
                args.mode_block, args.time_block, args.gzip,
            )
            merge_planes(files, out_file, args.plane_time_block, args.gzip)

            out_file.attrs["merged_mode_count"] = total_modes


def main() -> None:
    args = parse_args()
    if args.mode_block < 1 or args.time_block < 1 or args.plane_time_block < 1:
        raise ValueError("Los tamaños de bloque deben ser enteros positivos.")

    run_dir = Path(args.run_dir).expanduser().resolve()
    output = run_dir / f"{run_dir.name}.h5"
    temporary = run_dir / f".{run_dir.name}.merge.tmp.h5"

    if output.exists() and not args.overwrite:
        raise FileExistsError(
            f"El archivo final ya existe: {output}. Usa --overwrite para reemplazarlo."
        )
    if temporary.exists():
        temporary.unlink()

    files = collect_rank_files(run_dir, args.pattern, output)

    print("Archivos de entrada:")
    for path in files:
        print(f"  {path.name}")
    print(f"Salida: {output}")

    try:
        merge(files, temporary, args)
        temporary.replace(output)
    except Exception:
        if temporary.exists():
            temporary.unlink()
        raise

    if args.keep_ranks:
        print("Merge terminado correctamente. Los archivos por rank se conservaron.")
        return

    for path in files:
        path.unlink()
    print("Merge terminado correctamente. Se eliminaron los archivos por rank.")


if __name__ == "__main__":
    main()
