#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from __future__ import annotations

import os
import argparse
from dataclasses import dataclass
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import h5py


# ============================================================
# Geometría EF-Schwarzschild
# ============================================================

def A_of_r(r: np.ndarray, Rs: float) -> np.ndarray:
    return 1.0 + Rs / r


def p_of_r(r: np.ndarray, Rs: float) -> np.ndarray:
    return r**2 / np.sqrt(A_of_r(r, Rs))


def w_of_r(r: np.ndarray, Rs: float) -> np.ndarray:
    return r**2 * np.sqrt(A_of_r(r, Rs))


def q_of_r(r: np.ndarray, ell: int) -> np.ndarray:
    return ell * (ell + 1.0) / (r**2)


# ============================================================
# Configuración
# ============================================================

@dataclass
class SLConfig:
    rmin: float = 0.8
    rmax: float = 800.8
    Nr: int = 8000
    Rs: float = 0.5
    ell_min: int = 0
    ell_max: int = 80
    n_modes: int | None = None
    output_file: str = "sl_basis.h5"
    eig_sigma: float = 0.0
    max_workers: int | None = None
    compression: str | None = None
    compression_level: int = 0

    def resolved_n_modes(self) -> int:
        return self.Nr // 2 if self.n_modes is None else self.n_modes


# ============================================================
# Construcción matrices
# ============================================================

def build_matrices_uniform(
    rmin: float,
    rmax: float,
    Nr: int,
    Rs: float,
    ell: int,
    kappa: float | None = None,
):
    if Nr < 3:
        raise ValueError("Nr debe ser al menos 3.")
    if rmin <= 0.0:
        raise ValueError("rmin debe ser positivo.")
    if rmax <= rmin:
        raise ValueError("Debe cumplirse rmax > rmin.")
    if ell < 0:
        raise ValueError("ell debe ser no negativo.")

    r = np.linspace(rmin, rmax, Nr, dtype=np.float64)
    dr = (rmax - rmin) / (Nr - 1)

    if kappa is None:
        kappa = (ell + 1.0) / rmax

    p = p_of_r(r, Rs)
    w = w_of_r(r, Rs)
    q = q_of_r(r, ell)

    p_face = 2.0 * p[:-1] * p[1:] / (p[:-1] + p[1:] + 1e-300)
    coeff = p_face / dr

    main = np.zeros(Nr, dtype=np.float64)
    off = np.zeros(Nr - 1, dtype=np.float64)

    main[:-1] += coeff
    main[1:] += coeff
    off[:] -= coeff

    K = sp.diags([off, main, off], offsets=[-1, 0, 1], format="csr")
    M = sp.diags(w * dr, 0, format="csr")
    V = sp.diags(w * q * dr, 0, format="csr")

    keep = np.arange(1, Nr)
    K = K[keep][:, keep]
    M = M[keep][:, keep]
    V = V[keep][:, keep]

    p_rmax = float(p[-1])
    ndof = K.shape[0]
    robin = sp.csr_matrix(([p_rmax * kappa], ([ndof - 1], [ndof - 1])), shape=K.shape)
    K = K + robin

    return r, dr, (K + V).tocsr(), M.tocsr()


# ============================================================
# Problema espectral
# ============================================================

def solve_basis_uniform(
    rmin: float,
    rmax: float,
    Nr: int,
    Rs: float,
    ell: int,
    n_modes: int,
    sigma: float = 0.0,
):
    if n_modes <= 0:
        raise ValueError("n_modes debe ser positivo.")
    if n_modes >= Nr - 1:
        raise ValueError("n_modes debe satisfacer n_modes < Nr-1.")

    r, dr, Aop, M = build_matrices_uniform(rmin, rmax, Nr, Rs, ell)

    evals, evecs = spla.eigsh(
        Aop,
        k=n_modes,
        M=M,
        sigma=sigma,
        which="LM",
    )

    idx = np.argsort(evals)
    evals = np.asarray(evals[idx], dtype=np.float64)
    evecs = np.asarray(evecs[:, idx], dtype=np.float64)

    for j in range(n_modes):
        norm2 = float(evecs[:, j].T @ (M @ evecs[:, j]))
        evecs[:, j] /= np.sqrt(norm2 + 1e-300)

    modes = np.zeros((Nr, n_modes), dtype=np.float64)
    modes[1:, :] = evecs

    return r, evals, modes


# ============================================================
# Worker
# ============================================================

def _worker_compute_one_ell(args):
    (
        ell,
        rmin,
        rmax,
        Nr,
        Rs,
        n_modes,
        sigma,
    ) = args

    # Evita oversubscription dentro de cada proceso
    os.environ.setdefault("OMP_NUM_THREADS", "1")
    os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
    os.environ.setdefault("MKL_NUM_THREADS", "1")
    os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

    r, lam, R = solve_basis_uniform(
        rmin=rmin,
        rmax=rmax,
        Nr=Nr,
        Rs=Rs,
        ell=ell,
        n_modes=n_modes,
        sigma=sigma,
    )

    # Guardamos R_T para lectura cómoda en Fortran
    return ell, r, lam, R.T.copy()


# ============================================================
# Guardado HDF5
# ============================================================

def save_basis_group(
    h5f: h5py.File,
    ell: int,
    r: np.ndarray,
    lam: np.ndarray,
    RT: np.ndarray,
    compression: str | None = None,
    compression_level: int = 0,
) -> None:
    gname = f"sl_basis_ell_{ell:03d}"

    if gname in h5f:
        del h5f[gname]

    grp = h5f.create_group(gname)

    Nm = lam.shape[0]
    Nr = r.shape[0]

    grp.attrs["ell"] = int(ell)
    grp.attrs["Nr"] = int(Nr)
    grp.attrs["Nm"] = int(Nm)

    kwargs = {}
    if compression is not None:
        kwargs["compression"] = compression
        kwargs["compression_opts"] = compression_level
        kwargs["shuffle"] = True

    grp.create_dataset("r", data=r, dtype="f8", **kwargs)
    grp.create_dataset("lam", data=lam, dtype="f8", **kwargs)
    grp.create_dataset("R_T", data=RT, dtype="f8", **kwargs)


# ============================================================
# Pipeline principal
# ============================================================

def run_precompute_parallel(config: SLConfig) -> None:
    n_modes = config.resolved_n_modes()

    if config.max_workers is None:
        max_workers = 8
    else:
        max_workers = config.max_workers

    print(f"{n_modes} modos por ℓ")
    print(f"ℓ desde {config.ell_min} hasta {config.ell_max}")
    print(f"max_workers = {max_workers}")
    print(f"Archivo de salida: {config.output_file}")
    print()

    out_dir = os.path.dirname(os.path.abspath(config.output_file))
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    ells = list(range(config.ell_min, config.ell_max + 1))

    with h5py.File(config.output_file, "w") as h5f:
        h5f.attrs["rmin"] = config.rmin
        h5f.attrs["rmax"] = config.rmax
        h5f.attrs["Nr"] = config.Nr
        h5f.attrs["Rs"] = config.Rs
        h5f.attrs["ell_min"] = config.ell_min
        h5f.attrs["ell_max"] = config.ell_max
        h5f.attrs["n_modes"] = n_modes
        h5f.attrs["eig_sigma"] = config.eig_sigma
        h5f.attrs["format"] = "EF-Schwarzschild Sturm-Liouville basis for Fortran"
        h5f.attrs["storage"] = "datasets: r, lam, R_T"

        jobs = [
            (
                ell,
                config.rmin,
                config.rmax,
                config.Nr,
                config.Rs,
                n_modes,
                config.eig_sigma,
            )
            for ell in ells
        ]

        with ProcessPoolExecutor(max_workers=max_workers) as ex:
            futures = [ex.submit(_worker_compute_one_ell, job) for job in jobs]

            for fut in as_completed(futures):
                ell, r, lam, RT = fut.result()
                save_basis_group(
                    h5f,
                    ell=ell,
                    r=r,
                    lam=lam,
                    RT=RT,
                    compression=config.compression,
                    compression_level=config.compression_level,
                )
                print(f"[OK] ℓ={ell:03d} guardado")

    print("\nTodo terminado correctamente.")


# ============================================================
# CLI
# ============================================================

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Construcción paralela de bases SL en HDF5."
    )

    parser.add_argument("--rmin", type=float, default=0.8)
    parser.add_argument("--rmax", type=float, default=600.8)
    parser.add_argument("--Nr", type=int, default=8000)
    parser.add_argument("--Rs", type=float, default=0.5)
    parser.add_argument("--ell-min", type=int, default=0)
    parser.add_argument("--ell-max", type=int, default=80)
    parser.add_argument("--n-modes", type=int, default=None)
    parser.add_argument("--output-file", type=str, default="sl_basis.h5")
    parser.add_argument("--eig-sigma", type=float, default=0.0)
    parser.add_argument("--max-workers", type=int, default=None)
    parser.add_argument("--compression", type=str, default="none")
    parser.add_argument("--compression-level", type=int, default=0)

    return parser


def config_from_args(args: argparse.Namespace) -> SLConfig:
    comp = args.compression
    if comp is not None and comp.lower() == "none":
        comp = None

    return SLConfig(
        rmin=args.rmin,
        rmax=args.rmax,
        Nr=args.Nr,
        Rs=args.Rs,
        ell_min=args.ell_min,
        ell_max=args.ell_max,
        n_modes=args.n_modes,
        output_file=args.output_file,
        eig_sigma=args.eig_sigma,
        max_workers=args.max_workers,
        compression=comp,
        compression_level=args.compression_level,
    )


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    config = config_from_args(args)
    run_precompute_parallel(config)


if __name__ == "__main__":
    main()
