#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import gc
import sys
import shutil
import subprocess
from pathlib import Path

import numpy as np
import h5py
import matplotlib.pyplot as plt


# =========================================================
# CONFIGURACIÓN GENERAL
# =========================================================

SIM_DIR = "sigma_all_4"

# Control fino del pipeline
RUN_PREPROCESS = False # Pon True solo si quieres rehacer merge_h5.py y sum_lm.py
RUN_SPECTRA = True
RUN_XY = True
RUN_XZ = True
RUN_MOVIES = True
FORCE_PREPROCESS = False # Si RUN_PREPROCESS=True, fuerza rehacer aunque ya existan outputs

CONFIG = {
    "work_dir": SIM_DIR,

    "merge_script_src": "merge_h5.py",
    "sum_script_src": "sum_lm_fixed.py",
    "merge_script_dst": "merge_h5.py",
    "sum_script_dst": "sum_lm.py",

    "python_exec": sys.executable,

    "modes_h5": f"{SIM_DIR}/modes.h5",
    "basis_h5": "/home/flavio/Codes/KG-Y/sl_basis.h5",
    "spectra_out_dir": f"{SIM_DIR}/spectra_all_times",
    "use_spectra_cache": True,
    "spectra_cache_file": None,
    "rebuild_spectra_cache": False,

    "planes_h5": f"{SIM_DIR}/planes_sum.h5",
    "xy_out_dir": f"{SIM_DIR}/frames_phi_xy",
    "xz_out_dir": f"{SIM_DIR}/frames_phi_xz",

    "movie_pk": f"{SIM_DIR}/Pk_movie.mp4",
    "movie_heat": f"{SIM_DIR}/heatmap_movie.mp4",
    "movie_xy": f"{SIM_DIR}/phi_xy_movie.mp4",
    "movie_xz": f"{SIM_DIR}/phi_xz_movie.mp4",

    "field_name_modes": "phi",
    "field_name_xy": "phi_xy",
    "field_name_xz": "phi_xz",

    "RMIN": 0.8,
    "RMAX": 600.8,
    "NR": 8000,
    "RS": 0.5,

    "NBINS": 35,
    "T_FACTOR": 5.0,
    "ELL_MAX": 80,
    "TBLOCK": 8,
    "EPS": 1e-300,

    "r_cut": 500.0,
    "nx": 700,
    "ny": 700,
    "interp_method": "linear",
    "fill_with_nearest": True,
    "clip_frac": 0.9,
    "dpi": 150,
    "cmap": "jet",

    "fps": 10,
    "ffmpeg_crf": 18,
    "ffmpeg_preset": "medium",
}

MODE_RE = re.compile(r"^mode_(\d+)$")


# =========================================================
# UTILIDADES GENERALES
# =========================================================


def log(msg):
    print(msg, flush=True)


def ensure_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def run_cmd(cmd, cwd=None):
    log(f"[RUN] {' '.join(cmd)}")
    subprocess.run(cmd, cwd=cwd, check=True)


def copy_file(src, dst, work_dir="."):
    src_path = Path(src).resolve()
    dst_path = Path(work_dir) / dst

    if not src_path.exists():
        raise FileNotFoundError(f"No existe el archivo fuente: {src_path}")

    ensure_parent(dst_path)
    shutil.copy2(src_path, dst_path)
    log(f"[COPY] {src_path} -> {dst_path}")


# =========================================================
# PASO 1: COPIAR Y EJECUTAR merge_h5.py y sum_lm.py
# =========================================================


def copy_and_run_preprocessing(cfg, force=False):
    sim_dir = Path(cfg["work_dir"]).resolve()
    sim_dir.mkdir(parents=True, exist_ok=True)

    modes_h5 = sim_dir / Path(cfg["modes_h5"]).name
    planes_h5 = sim_dir / Path(cfg["planes_h5"]).name

    if (not force) and modes_h5.exists() and planes_h5.exists():
        log(f"[SKIP] Preprocessing omitido: ya existen {modes_h5.name} y {planes_h5.name}")
        return

    parent_dir = sim_dir.parent

    merge_src = parent_dir / cfg["merge_script_src"]
    sum_src = parent_dir / cfg["sum_script_src"]

    merge_dst = sim_dir / cfg["merge_script_dst"]
    sum_dst = sim_dir / cfg["sum_script_dst"]

    if not merge_src.exists():
        raise FileNotFoundError(f"No existe el archivo fuente: {merge_src}")
    if not sum_src.exists():
        raise FileNotFoundError(f"No existe el archivo fuente: {sum_src}")

    shutil.copy2(merge_src, merge_dst)
    shutil.copy2(sum_src, sum_dst)

    log(f"[COPY] {merge_src} -> {merge_dst}")
    log(f"[COPY] {sum_src} -> {sum_dst}")

    run_cmd([cfg["python_exec"], merge_dst.name], cwd=str(sim_dir))
    run_cmd([cfg["python_exec"], sum_dst.name], cwd=str(sim_dir))


# =========================================================
# PARTE A: ESPECTRO P(l,n) y P(k)
# =========================================================


def A_of_r(r, Rs):
    return 1.0 + Rs / r


def w_of_r(r, Rs):
    return r**2 * np.sqrt(A_of_r(r, Rs))


def n_to_ellm(n: int):
    ell = int(np.sqrt(n))
    rr = n - ell * ell
    m = int(rr) - ell
    return ell, m


def load_sl_basis_h5(basis_h5: str, ell: int):
    if not os.path.exists(basis_h5):
        raise FileNotFoundError(f"No encontré el archivo de bases SL: {basis_h5}")

    gname = f"sl_basis_ell_{ell:03d}"

    with h5py.File(basis_h5, "r") as f:
        if gname not in f:
            raise FileNotFoundError(f"No encontré la base SL para ell={ell} en {basis_h5}")

        grp = f[gname]

        if "lam" not in grp:
            raise KeyError(f"Falta dataset 'lam' en {gname}")
        if "R_T" not in grp:
            raise KeyError(f"Falta dataset 'R_T' en {gname}")

        lam = grp["lam"][...].astype(np.float64, copy=False)
        RT = grp["R_T"][...].astype(np.float64, copy=False)

    if RT.ndim != 2:
        raise RuntimeError(f"R_T en ell={ell} no es 2D, shape={RT.shape}")
    if lam.ndim != 1:
        raise RuntimeError(f"lam en ell={ell} no es 1D, shape={lam.shape}")

    if RT.shape[0] != lam.shape[0]:
        raise RuntimeError(
            f"Inconsistencia en ell={ell}: R_T.shape[0]={RT.shape[0]} != len(lam)={lam.shape[0]}"
        )

    R = np.ascontiguousarray(RT.T, dtype=np.float64)
    lam = np.ascontiguousarray(lam, dtype=np.float64)

    return lam, R


def build_P_ln_block(h5, s0, s1, cfg, ells, ell_to_modes, basis_cache, Nm_ref):
    tblock = s1 - s0
    P_ln_block = np.zeros((len(ells), Nm_ref, tblock), dtype=np.float64)
    field_name = cfg["field_name_modes"]

    for iell, ell in enumerate(ells):
        WRt = basis_cache[ell]["WRt"]
        modes = ell_to_modes[ell]
        n_m = len(modes)

        Phi = np.empty((cfg["NR"], n_m, tblock), dtype=np.complex128)
        for j, (_, gname) in enumerate(modes):
            arr = h5[gname][field_name][:, :, s0:s1]
            Phi[:, j, :] = arr[0, :, :] + 1j * arr[1, :, :]

        Phi2D = Phi.reshape(cfg["NR"], n_m * tblock)
        A2D = WRt @ np.conjugate(Phi2D)
        A = A2D.reshape(Nm_ref, n_m, tblock)

        P_ln_block[iell, :, :] = np.sum(np.abs(A) ** 2, axis=1)

        del Phi, Phi2D, A2D, A
        gc.collect()

    return P_ln_block


def _default_spectra_cache_file(modes_h5, field_name, ell_max, nbins):
    src = Path(modes_h5)
    ell_tag = "all" if ell_max is None else str(int(ell_max))
    name = f"{src.stem}__spectra_cache_{field_name}_ell{ell_tag}_nb{int(nbins)}.h5"
    return str(src.with_name(name))


def _spectra_cache_is_valid(cache_file, cfg):
    if not os.path.exists(cache_file):
        return False

    try:
        modes_stat = os.stat(cfg["modes_h5"])
        basis_stat = os.stat(cfg["basis_h5"])

        with h5py.File(cache_file, "r") as h:
            attrs = h.attrs
            required_attrs = [
                "modes_h5", "modes_size", "modes_mtime_ns",
                "basis_h5", "basis_size", "basis_mtime_ns",
                "field_name_modes", "ELL_MAX", "NBINS",
                "RMIN", "RMAX", "NR", "RS",
                "HMIN", "HMAX", "PK_MIN", "PK_MAX",
            ]
            for key in required_attrs:
                if key not in attrs:
                    return False

            required_dsets = ["P_ln", "P_bin", "P_mask", "k_centers", "k_edges", "ells", "snap_index"]
            if not all(name in h for name in required_dsets):
                return False

            def _same_path(a, b):
                return os.path.abspath(str(a)) == os.path.abspath(str(b))

            if not _same_path(attrs["modes_h5"], cfg["modes_h5"]):
                return False
            if int(attrs["modes_size"]) != int(modes_stat.st_size):
                return False
            if int(attrs["modes_mtime_ns"]) != int(modes_stat.st_mtime_ns):
                return False

            if not _same_path(attrs["basis_h5"], cfg["basis_h5"]):
                return False
            if int(attrs["basis_size"]) != int(basis_stat.st_size):
                return False
            if int(attrs["basis_mtime_ns"]) != int(basis_stat.st_mtime_ns):
                return False

            cache_ell_max = attrs["ELL_MAX"]
            cache_ell_max = None if str(cache_ell_max) == "None" else int(cache_ell_max)
            if cache_ell_max != cfg["ELL_MAX"]:
                return False

            if str(attrs["field_name_modes"]) != str(cfg["field_name_modes"]):
                return False
            if int(attrs["NBINS"]) != int(cfg["NBINS"]):
                return False
            if not np.isclose(float(attrs["RMIN"]), float(cfg["RMIN"])):
                return False
            if not np.isclose(float(attrs["RMAX"]), float(cfg["RMAX"])):
                return False
            if int(attrs["NR"]) != int(cfg["NR"]):
                return False
            if not np.isclose(float(attrs["RS"]), float(cfg["RS"])):
                return False

        return True
    except Exception:
        return False


def _prepare_spectra_cache(cfg):
    modes_h5 = cfg["modes_h5"]
    basis_h5 = cfg["basis_h5"]
    cache_file = cfg.get("spectra_cache_file") or _default_spectra_cache_file(
        modes_h5, cfg["field_name_modes"], cfg["ELL_MAX"], cfg["NBINS"]
    )

    use_cache = bool(cfg.get("use_spectra_cache", True))
    rebuild = bool(cfg.get("rebuild_spectra_cache", False))

    if use_cache and (not rebuild) and _spectra_cache_is_valid(cache_file, cfg):
        log(f"[CACHE] Reutilizando cache de espectros: {cache_file}")
        return cache_file

    if not os.path.exists(modes_h5):
        raise FileNotFoundError(f"No existe modes.h5: {modes_h5}")
    if not os.path.exists(basis_h5):
        raise FileNotFoundError(f"No existe basis_h5: {basis_h5}")

    ensure_parent(cache_file)
    tmp_cache = cache_file + ".tmp"
    if os.path.exists(tmp_cache):
        os.remove(tmp_cache)

    r = np.linspace(cfg["RMIN"], cfg["RMAX"], cfg["NR"])
    dr = (cfg["RMAX"] - cfg["RMIN"]) / (cfg["NR"] - 1)
    weights = (w_of_r(r, cfg["RS"]) * dr).astype(np.float64, copy=False)

    with h5py.File(modes_h5, "r") as h5:
        mode_names = [k for k in h5.keys() if MODE_RE.match(k)]
        mode_names.sort()

        if not mode_names:
            raise RuntimeError("No encontré grupos mode_XXXXX en el HDF5.")

        dset0 = h5[mode_names[0]][cfg["field_name_modes"]]
        if dset0.ndim != 3 or dset0.shape[0] != 2:
            raise RuntimeError(
                f"Dataset {cfg['field_name_modes']} no tiene forma (2,Nr,Nt). Forma: {dset0.shape}"
            )

        if dset0.shape[1] != cfg["NR"]:
            raise RuntimeError(
                f"NR inconsistente: HDF5 Nr={dset0.shape[1]} vs NR={cfg['NR']}"
            )
        Nt = int(dset0.shape[2])

    ell_to_modes = {}
    for gname in mode_names:
        n_lin = int(MODE_RE.match(gname).group(1))
        ell, m = n_to_ellm(n_lin)
        if (cfg["ELL_MAX"] is None) or (ell <= cfg["ELL_MAX"]):
            ell_to_modes.setdefault(ell, []).append((m, gname))

    ells = sorted(ell_to_modes.keys())
    for ell in ells:
        ell_to_modes[ell].sort(key=lambda x: x[0])

    log(f"[INFO] Nt={Nt} | n_ell={len(ells)} | TBLOCK={cfg['TBLOCK']}")

    basis_cache = {}
    Nm_ref = None

    for ell in ells:
        lam, Rm = load_sl_basis_h5(basis_h5, ell)

        if Rm.shape[0] != cfg["NR"]:
            raise RuntimeError(
                f"NR inconsistente en base ell={ell}: {Rm.shape[0]} vs {cfg['NR']}"
            )

        Nm = int(Rm.shape[1])
        if Nm_ref is None:
            Nm_ref = Nm
        elif Nm != Nm_ref:
            raise RuntimeError(f"Nm inconsistente: ell={ell} Nm={Nm} vs ref={Nm_ref}")

        k_ell = np.sqrt(np.maximum(lam, 0.0))
        WRt = np.ascontiguousarray(Rm.T * weights[None, :], dtype=np.float64)

        basis_cache[ell] = {"k": k_ell, "WRt": WRt}

        del lam, Rm
        gc.collect()

    k_all_static = np.concatenate([basis_cache[ell]["k"] for ell in ells])
    k_pos = k_all_static[k_all_static > 0]
    if k_pos.size == 0:
        raise RuntimeError("No hay valores positivos de k para construir los bins logarítmicos.")

    kmin = np.percentile(k_pos, 1)
    kmax = np.percentile(k_pos, 99)
    k_edges = np.logspace(np.log10(kmin), np.log10(kmax), cfg["NBINS"] + 1)
    k_centers = np.sqrt(k_edges[:-1] * k_edges[1:])

    k_all_flat = np.concatenate([basis_cache[ell]["k"] for ell in ells]).astype(np.float64, copy=False)
    idx_bins = np.searchsorted(k_edges, k_all_flat, side="right") - 1
    idx_valid = (idx_bins >= 0) & (idx_bins < cfg["NBINS"])

    global_hmin = +np.inf
    global_hmax = -np.inf
    global_pk_min = +np.inf
    global_pk_max = -np.inf

    modes_stat = os.stat(modes_h5)
    basis_stat = os.stat(basis_h5)

    with h5py.File(tmp_cache, "w") as hc:
        hc.attrs["modes_h5"] = os.path.abspath(modes_h5)
        hc.attrs["modes_size"] = int(modes_stat.st_size)
        hc.attrs["modes_mtime_ns"] = int(modes_stat.st_mtime_ns)
        hc.attrs["basis_h5"] = os.path.abspath(basis_h5)
        hc.attrs["basis_size"] = int(basis_stat.st_size)
        hc.attrs["basis_mtime_ns"] = int(basis_stat.st_mtime_ns)
        hc.attrs["field_name_modes"] = str(cfg["field_name_modes"])
        hc.attrs["ELL_MAX"] = "None" if cfg["ELL_MAX"] is None else int(cfg["ELL_MAX"])
        hc.attrs["NBINS"] = int(cfg["NBINS"])
        hc.attrs["RMIN"] = float(cfg["RMIN"])
        hc.attrs["RMAX"] = float(cfg["RMAX"])
        hc.attrs["NR"] = int(cfg["NR"])
        hc.attrs["RS"] = float(cfg["RS"])
        hc.attrs["Nt"] = int(Nt)
        hc.attrs["Nm_ref"] = int(Nm_ref)
        hc.attrs["n_ell"] = int(len(ells))

        hc.create_dataset("ells", data=np.asarray(ells, dtype=np.int32))
        hc.create_dataset("k_edges", data=k_edges.astype(np.float64, copy=False))
        hc.create_dataset("k_centers", data=k_centers.astype(np.float64, copy=False))
        hc.create_dataset("snap_index", data=np.arange(Nt, dtype=np.int32))

        P_ln_ds = hc.create_dataset(
            "P_ln",
            shape=(Nt, len(ells), Nm_ref),
            dtype=np.float32,
            chunks=(1, len(ells), Nm_ref),
        )
        P_bin_ds = hc.create_dataset(
            "P_bin",
            shape=(Nt, cfg["NBINS"]),
            dtype=np.float64,
            chunks=(max(1, min(cfg["TBLOCK"], Nt)), cfg["NBINS"]),
        )
        P_mask_ds = hc.create_dataset(
            "P_mask",
            shape=(Nt, cfg["NBINS"]),
            dtype=np.bool_,
            chunks=(max(1, min(cfg["TBLOCK"], Nt)), cfg["NBINS"]),
        )

        with h5py.File(modes_h5, "r") as h5:
            for s0 in range(0, Nt, cfg["TBLOCK"]):
                s1 = min(s0 + cfg["TBLOCK"], Nt)
                P_ln_block = build_P_ln_block(
                    h5, s0, s1, cfg, ells, ell_to_modes, basis_cache, Nm_ref
                )

                P_ln_ds[s0:s1, :, :] = np.moveaxis(P_ln_block, 2, 0).astype(np.float32, copy=False)

                for dt, s in enumerate(range(s0, s1)):
                    P_ln = P_ln_block[:, :, dt]
                    L = np.log10(P_ln + cfg["EPS"])
                    global_hmin = min(global_hmin, float(np.nanmin(L)))
                    global_hmax = max(global_hmax, float(np.nanmax(L)))

                    P_flat = P_ln.reshape(-1).astype(np.float64, copy=False)
                    P_bin = np.zeros(cfg["NBINS"], dtype=np.float64)
                    N_bin = np.zeros(cfg["NBINS"], dtype=np.int64)

                    v = idx_valid & np.isfinite(P_flat)
                    np.add.at(P_bin, idx_bins[v], P_flat[v])
                    np.add.at(N_bin, idx_bins[v], 1)
                    mask = (N_bin > 0) & (P_bin > 0) & np.isfinite(P_bin)

                    P_bin_ds[s, :] = P_bin
                    P_mask_ds[s, :] = mask

                    if np.any(mask):
                        global_pk_min = min(global_pk_min, float(np.min(P_bin[mask])))
                        global_pk_max = max(global_pk_max, float(np.max(P_bin[mask])))

                hc.flush()
                log(
                    f"[CACHE-SCAN] {s0}..{s1-1} | "
                    f"heat=[{global_hmin:.3g},{global_hmax:.3g}] | "
                    f"pk=[{global_pk_min:.3g},{global_pk_max:.3g}]"
                )
                del P_ln_block
                gc.collect()

        hc.attrs["HMIN"] = float(global_hmin)
        hc.attrs["HMAX"] = float(global_hmax)
        hc.attrs["PK_MIN"] = float(global_pk_min)
        hc.attrs["PK_MAX"] = float(global_pk_max)

    os.replace(tmp_cache, cache_file)
    log(f"[CACHE] Cache de espectros generado: {cache_file}")
    return cache_file


def generate_spectra_all_times(cfg):
    log("[STEP] Generando heatmaps P(l,n) y curvas P(k)")

    out_dir = cfg["spectra_out_dir"]
    heat_dir = os.path.join(out_dir, "heatmap")
    pk_dir = os.path.join(out_dir, "pk")
    os.makedirs(heat_dir, exist_ok=True)
    os.makedirs(pk_dir, exist_ok=True)

    cache_file = _prepare_spectra_cache(cfg)

    with h5py.File(cache_file, "r") as hc:
        P_ln_ds = hc["P_ln"]
        P_bin_ds = hc["P_bin"]
        P_mask_ds = hc["P_mask"]
        k_centers = hc["k_centers"][...]
        k_edges = hc["k_edges"][...]
        Nt = int(P_ln_ds.shape[0])

        HMIN = float(hc.attrs["HMIN"])
        HMAX = float(hc.attrs["HMAX"])
        PK_MIN = float(hc.attrs["PK_MIN"])
        PK_MAX = float(hc.attrs["PK_MAX"])

        log(f"[INFO] Heatmap fixed (log10): [{HMIN:.6g},{HMAX:.6g}]")
        log(f"[INFO] P(k) fixed: [{PK_MIN:.6g},{PK_MAX:.6g}]")

        for s in range(Nt):
            P_ln = P_ln_ds[s, :, :].astype(np.float64, copy=False)
            P_bin = P_bin_ds[s, :].astype(np.float64, copy=False)
            mask = P_mask_ds[s, :].astype(bool, copy=False)

            fig = plt.figure(figsize=(7.2, 4.0))
            plt.imshow(
                np.log10(P_ln + cfg["EPS"]),
                origin="lower",
                aspect="auto",
                vmin=HMIN,
                vmax=HMAX,
            )
            plt.colorbar(label=r"$\log_{10} P(\ell,n)$")
            plt.xlabel(r"$n$")
            plt.ylabel(r"$\ell$")
            plt.tight_layout()
            fig.savefig(os.path.join(heat_dir, f"P_ln_t{s:04d}.png"), dpi=150)
            plt.close(fig)

            fig = plt.figure(figsize=(6.4, 4.0))
            plt.loglog(k_centers[mask], P_bin[mask], marker="o", linestyle="-")
            plt.xlabel(r"$k$")
            plt.ylabel(r"$P_{\rm bin}(k)$")
            plt.grid(True, which="both")
            plt.xlim(k_edges[0], k_edges[-1])
            plt.ylim(PK_MIN, PK_MAX)
            plt.tight_layout()
            fig.savefig(os.path.join(pk_dir, f"P_k_t{s:04d}.png"), dpi=150)
            plt.close(fig)

            if (s % max(1, cfg["TBLOCK"])) == 0 or s == (Nt - 1):
                log(f"[SAVE] snapshot {s}/{Nt-1}")

        P0 = P_bin_ds[0, :].astype(np.float64, copy=False)
        m0 = P_mask_ds[0, :].astype(bool, copy=False)
        P1 = P_bin_ds[Nt - 1, :].astype(np.float64, copy=False)
        m1 = P_mask_ds[Nt - 1, :].astype(bool, copy=False)

        fig = plt.figure(figsize=(6.8, 4.2))
        plt.loglog(
            k_centers[m0], P0[m0],
            marker="o", linestyle="-",
            label=rf"$t={cfg['T_FACTOR'] * 0:.1f}$",
        )
        plt.loglog(
            k_centers[m1], P1[m1],
            marker="o", linestyle="-",
            label=rf"$t={cfg['T_FACTOR'] * (Nt - 1):.1f}$",
        )
        plt.ylabel(r"$P_{\rm bin}(k)$")
        plt.grid(True, which="both")
        plt.xlim(k_edges[0], k_edges[-1])
        plt.ylim(PK_MIN, PK_MAX)
        plt.legend()
        plt.tight_layout()
        fig.savefig(os.path.join(out_dir, "Pk_first_last.png"), dpi=180)
        plt.close(fig)

    log(f"[OK] Frames en: {heat_dir} y {pk_dir}")
    log(f"[OK] Comparación: {os.path.join(out_dir, 'Pk_first_last.png')}")

# =========================================================
# PARTE B: GRAFICAR XY A TODOS LOS TIEMPOS
# =========================================================


def _build_xy_mesh(r_sel, phi):
    """
    Construye una malla curvilínea (X,Y) a partir de la malla polar (r,phi)
    sin interpolar a una malla cartesiana regular.
    """
    R, PH = np.meshgrid(r_sel, phi, indexing="ij")
    X = R * np.cos(PH)
    Y = R * np.sin(PH)
    return X, Y


def _default_xy_cache_file(h5file, field, r_cut, snap_start, snap_stop, stride):
    src = Path(h5file)
    stop_tag = "all" if snap_stop is None else str(int(snap_stop))
    name = (
        f"{src.stem}__xycache_{field}"
        f"_rcut{r_cut:g}_s{int(snap_start)}_e{stop_tag}_st{int(stride)}.h5"
    )
    return str(src.with_name(name))


def _xy_cache_is_valid(cache_file, source_h5, field, r_cut, snap_start, snap_stop, stride):
    if not os.path.exists(cache_file):
        return False

    try:
        src_stat = os.stat(source_h5)
        with h5py.File(cache_file, "r") as h:
            attrs = h.attrs
            required = [
                "source_h5", "source_size", "source_mtime_ns",
                "field", "r_cut", "snap_start", "snap_stop", "stride",
                "vmax_re", "vmax_im", "vmax_abs",
            ]
            for key in required:
                if key not in attrs:
                    return False

            same_source = os.path.abspath(str(attrs["source_h5"])) == os.path.abspath(source_h5)
            same_size = int(attrs["source_size"]) == int(src_stat.st_size)
            same_mtime = int(attrs["source_mtime_ns"]) == int(src_stat.st_mtime_ns)
            same_field = str(attrs["field"]) == str(field)
            same_rcut = np.isclose(float(attrs["r_cut"]), float(r_cut))
            same_start = int(attrs["snap_start"]) == int(snap_start)
            cache_stop = attrs["snap_stop"]
            cache_stop = None if str(cache_stop) == "None" else int(cache_stop)
            same_stop = cache_stop == (None if snap_stop is None else int(snap_stop))
            same_stride = int(attrs["stride"]) == int(stride)

            if not all([same_source, same_size, same_mtime, same_field,
                        same_rcut, same_start, same_stop, same_stride]):
                return False

            return all(name in h for name in ["frames", "X", "Y"])
    except Exception:
        return False


def _build_xy_cache(
    h5file,
    cache_file,
    field,
    r_cut,
    snap_start=0,
    snap_stop=None,
    stride=1,
    verbose=True,
):
    """
    Hace un solo barrido sobre planes_sum.h5 y guarda en un archivo auxiliar:
      - la malla XY ya transformada
      - los snapshots recortados a r <= r_cut
      - los máximos globales para colorbar fija

    Luego el graficado lee de este cache y evita recalcular / releer todo
    desde el archivo original solo para fijar escalas.
    """
    os.makedirs(os.path.dirname(cache_file) or ".", exist_ok=True)

    with h5py.File(h5file, "r") as f:
        r = np.asarray(f["r"][...], float)
        phi = np.asarray(f["phi_grid"][...], float)
        U = f[field]
        tarr = np.asarray(f["t"][...], float) if "t" in f else None

        Nt = U.shape[0]
        if snap_stop is None or snap_stop > Nt:
            snap_stop = Nt

        mask_r = (r <= r_cut)
        if not np.any(mask_r):
            raise ValueError(f"No hay puntos con r <= {r_cut}. max(r)={r.max()}")

        r_sel = r[mask_r]
        X, Y = _build_xy_mesh(r_sel, phi)

        snap_ids = np.arange(snap_start, snap_stop, stride, dtype=np.int64)
        ns = len(snap_ids)
        if ns == 0:
            raise ValueError("No hay snapshots para cachear con ese rango/stride.")

        vmax_re = 0.0
        vmax_im = 0.0
        vmax_abs = 0.0

        with h5py.File(cache_file, "w") as h:
            h.attrs["source_h5"] = os.path.abspath(h5file)
            st = os.stat(h5file)
            h.attrs["source_size"] = int(st.st_size)
            h.attrs["source_mtime_ns"] = int(st.st_mtime_ns)
            h.attrs["field"] = str(field)
            h.attrs["r_cut"] = float(r_cut)
            h.attrs["snap_start"] = int(snap_start)
            h.attrs["snap_stop"] = "None" if snap_stop is None else int(snap_stop)
            h.attrs["stride"] = int(stride)

            h.create_dataset("r_sel", data=r_sel)
            h.create_dataset("phi", data=phi)
            h.create_dataset("X", data=X)
            h.create_dataset("Y", data=Y)
            h.create_dataset("snap_ids", data=snap_ids)
            if tarr is not None:
                h.create_dataset("t_sel", data=tarr[snap_ids])

            dset = h.create_dataset(
                "frames",
                shape=(ns, len(r_sel), len(phi)),
                dtype=U.dtype,
                chunks=(1, len(r_sel), len(phi)),
            )

            for j, s in enumerate(snap_ids):
                F = np.asarray(U[s, mask_r, :])
                dset[j, :, :] = F

                vmax_re = max(vmax_re, float(np.nanmax(np.abs(np.real(F)))))
                vmax_im = max(vmax_im, float(np.nanmax(np.abs(np.imag(F)))))
                vmax_abs = max(vmax_abs, float(np.nanmax(np.abs(F))))

                if verbose and (j == 0 or j % 20 == 0):
                    log(
                        f"[cache] {j+1}/{ns} | s={s} | "
                        f"vmax_re={vmax_re:.3g} "
                        f"vmax_im={vmax_im:.3g} "
                        f"vmax_abs={vmax_abs:.3g}"
                    )

            h.attrs["vmax_re"] = float(vmax_re)
            h.attrs["vmax_im"] = float(vmax_im)
            h.attrs["vmax_abs"] = float(vmax_abs)

    return cache_file, vmax_re, vmax_im, vmax_abs


def _ensure_xy_cache(
    h5file,
    field,
    r_cut,
    snap_start=0,
    snap_stop=None,
    stride=1,
    cache_file=None,
    rebuild=False,
    verbose=True,
):
    if cache_file is None:
        cache_file = _default_xy_cache_file(h5file, field, r_cut, snap_start, snap_stop, stride)

    if (not rebuild) and _xy_cache_is_valid(cache_file, h5file, field, r_cut, snap_start, snap_stop, stride):
        if verbose:
            log(f"[cache] usando cache existente: {cache_file}")
        with h5py.File(cache_file, "r") as h:
            return (
                cache_file,
                float(h.attrs["vmax_re"]),
                float(h.attrs["vmax_im"]),
                float(h.attrs["vmax_abs"]),
            )

    if verbose:
        log(f"[cache] construyendo cache XY: {cache_file}")
    return _build_xy_cache(
        h5file,
        cache_file,
        field,
        r_cut,
        snap_start=snap_start,
        snap_stop=snap_stop,
        stride=stride,
        verbose=verbose,
    )


def save_all_snapshots_three_panels_fixed_scale(
    h5file="planes_sum.h5",
    out_dir="frames_phi",
    field="phi_xy",
    r_cut=120.0,
    nx=700, ny=700,
    method="linear",
    cmap="jet",
    dpi=150,
    snap_start=0,
    snap_stop=None,
    stride=1,
    fill_with_nearest=True,
    clip_frac=1.0,
    verbose=True,
    use_xy_cache=True,
    xy_cache_file=None,
    rebuild_xy_cache=False,
):
    """
    Versión rápida para XY:
      1) no interpola a una malla cartesiana regular
      2) usa pcolormesh sobre la malla polar transformada
      3) opcionalmente crea un archivo auxiliar con los snapshots recortados
         y los vmax globales, para no recorrer dos veces el HDF5 principal

    Nota:
      - nx, ny, method y fill_with_nearest se conservan solo por compatibilidad.
    """
    del nx, ny, method, fill_with_nearest

    os.makedirs(out_dir, exist_ok=True)

    if use_xy_cache:
        cache_file, vmax_re, vmax_im, vmax_abs = _ensure_xy_cache(
            h5file,
            field,
            r_cut,
            snap_start=snap_start,
            snap_stop=snap_stop,
            stride=stride,
            cache_file=xy_cache_file,
            rebuild=rebuild_xy_cache,
            verbose=verbose,
        )

        with h5py.File(cache_file, "r") as h:
            X = np.asarray(h["X"][...])
            Y = np.asarray(h["Y"][...])
            frames = h["frames"]
            snap_ids = np.asarray(h["snap_ids"][...], dtype=int)
            tsel = np.asarray(h["t_sel"][...], float) if "t_sel" in h else None

            vmax_re *= clip_frac
            vmax_im *= clip_frac
            vmax_abs *= clip_frac

            if verbose:
                log(
                    f"[fixed] vlims: "
                    f"Re=[{-vmax_re:.4g},{vmax_re:.4g}]  "
                    f"Im=[{-vmax_im:.4g},{vmax_im:.4g}]  "
                    f"|.|=[0,{vmax_abs:.4g}]"
                )

            for j, s in enumerate(snap_ids):
                F = np.asarray(frames[j, :, :])
                ReF = np.real(F)
                ImF = np.imag(F)
                AbsF = np.abs(F)

                fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

                im0 = axs[0].pcolormesh(
                    X, Y, ReF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=-vmax_re,
                    vmax=vmax_re,
                    rasterized=True,
                )
                axs[0].set_title(r"$\Re(\Phi)$")
                axs[0].set_xlabel("x")
                axs[0].set_ylabel("y")
                axs[0].set_aspect("equal")
                axs[0].set_xlim(-r_cut, r_cut)
                axs[0].set_ylim(-r_cut, r_cut)

                im1 = axs[1].pcolormesh(
                    X, Y, ImF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=-vmax_im,
                    vmax=vmax_im,
                    rasterized=True,
                )
                axs[1].set_title(r"$\Im(\Phi)$")
                axs[1].set_xlabel("x")
                axs[1].set_ylabel("y")
                axs[1].set_aspect("equal")
                axs[1].set_xlim(-r_cut, r_cut)
                axs[1].set_ylim(-r_cut, r_cut)

                im2 = axs[2].pcolormesh(
                    X, Y, AbsF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=0.0,
                    vmax=vmax_abs,
                    rasterized=True,
                )
                axs[2].set_title(r"$|\Phi|$")
                axs[2].set_xlabel("x")
                axs[2].set_ylabel("y")
                axs[2].set_aspect("equal")
                axs[2].set_xlim(-r_cut, r_cut)
                axs[2].set_ylim(-r_cut, r_cut)

                fig.colorbar(im0, ax=axs[0], shrink=0.85)
                fig.colorbar(im1, ax=axs[1], shrink=0.85)
                fig.colorbar(im2, ax=axs[2], shrink=0.85)

                supt = f"snap={s}" if tsel is None else f"t={tsel[j]:.6g}"
                fig.suptitle(supt, y=1.02)

                out_path = os.path.join(out_dir, f"{field}_xy_{s:05d}.png")
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)

                if verbose and (j == 0 or j % 10 == 0):
                    log(f"[ok] {out_path}")

    else:
        # Ruta compatible, sin cache auxiliar.
        with h5py.File(h5file, "r") as f:
            r = np.asarray(f["r"][...], float)
            phi = np.asarray(f["phi_grid"][...], float)
            U = f[field]
            tarr = np.asarray(f["t"][...], float) if "t" in f else None

            Nt = U.shape[0]
            if snap_stop is None or snap_stop > Nt:
                snap_stop = Nt

            mask_r = (r <= r_cut)
            if not np.any(mask_r):
                raise ValueError(f"No hay puntos con r <= {r_cut}. max(r)={r.max()}")

            r_sel = r[mask_r]
            X, Y = _build_xy_mesh(r_sel, phi)

            vmax_re = 0.0
            vmax_im = 0.0
            vmax_abs = 0.0
            for s in range(snap_start, snap_stop, stride):
                F = U[s, mask_r, :]
                vmax_re = max(vmax_re, float(np.nanmax(np.abs(np.real(F)))))
                vmax_im = max(vmax_im, float(np.nanmax(np.abs(np.imag(F)))))
                vmax_abs = max(vmax_abs, float(np.nanmax(np.abs(F))))

            vmax_re *= clip_frac
            vmax_im *= clip_frac
            vmax_abs *= clip_frac

            for s in range(snap_start, snap_stop, stride):
                F = U[s, mask_r, :]
                ReF = np.real(F)
                ImF = np.imag(F)
                AbsF = np.abs(F)

                fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

                im0 = axs[0].pcolormesh(X, Y, ReF, shading="nearest", cmap=cmap,
                                        vmin=-vmax_re, vmax=vmax_re, rasterized=True)
                im1 = axs[1].pcolormesh(X, Y, ImF, shading="nearest", cmap=cmap,
                                        vmin=-vmax_im, vmax=vmax_im, rasterized=True)
                im2 = axs[2].pcolormesh(X, Y, AbsF, shading="nearest", cmap=cmap,
                                        vmin=0.0, vmax=vmax_abs, rasterized=True)

                axs[0].set_title(r"$\Re(\Phi)$")
                axs[1].set_title(r"$\Im(\Phi)$")
                axs[2].set_title(r"$|\Phi|$")
                for ax in axs:
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    ax.set_aspect("equal")
                    ax.set_xlim(-r_cut, r_cut)
                    ax.set_ylim(-r_cut, r_cut)

                fig.colorbar(im0, ax=axs[0], shrink=0.85)
                fig.colorbar(im1, ax=axs[1], shrink=0.85)
                fig.colorbar(im2, ax=axs[2], shrink=0.85)

                supt = f"snap={s}" if tarr is None else f"t={tarr[s]:.6g}"
                fig.suptitle(supt, y=1.02)
                out_path = os.path.join(out_dir, f"{field}_xy_{s:05d}.png")
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)

                if verbose and (s == snap_start or (s - snap_start) % (10 * stride) == 0):
                    log(f"[ok] {out_path}")

    if verbose:
        log(f"[done] frames en: {out_dir}")


def generate_xy_frames(cfg):
    log("[STEP] Generando frames XY para todos los tiempos")

    if not os.path.exists(cfg["planes_h5"]):
        raise FileNotFoundError(f"No existe planes_sum.h5: {cfg['planes_h5']}")

    save_all_snapshots_three_panels_fixed_scale(
        h5file=cfg["planes_h5"],
        out_dir=cfg["xy_out_dir"],
        field=cfg["field_name_xy"],
        r_cut=cfg["r_cut"],
        nx=cfg["nx"],
        ny=cfg["ny"],
        method=cfg["interp_method"],
        cmap=cfg["cmap"],
        dpi=cfg["dpi"],
        fill_with_nearest=cfg["fill_with_nearest"],
        clip_frac=cfg["clip_frac"],
        stride=1,
        verbose=True,
        use_xy_cache=cfg.get("use_xy_cache", True),
        xy_cache_file=cfg.get("xy_cache_file", None),
        rebuild_xy_cache=cfg.get("rebuild_xy_cache", False),
    )



def _guess_xz_angle_name(f, preferred=None):
    candidates = []
    if preferred:
        candidates.append(preferred)
    candidates.extend(["theta_grid", "theta", "th_grid", "phi_grid"])
    for k in candidates:
        if k in f:
            return k
    raise KeyError(
        "No encontré un dataset angular para XZ. Probé: "
        + ", ".join(candidates)
    )



def _build_xz_mesh(r_sel, ang):
    R, A = np.meshgrid(r_sel, ang, indexing="ij")
    amin = float(np.nanmin(ang))
    amax = float(np.nanmax(ang))
    span = float(amax - amin)

    # Si cubre ~2pi, interpretamos corte polar completo en el plano xz.
    if span > 1.5 * np.pi:
        X = R * np.cos(A)
        Z = R * np.sin(A)
    else:
        # Si cubre ~pi, interpretamos corte meridional usual.
        X = R * np.sin(A)
        Z = R * np.cos(A)
    return X, Z



def _default_xz_cache_file(h5file, field, r_cut, snap_start, snap_stop, stride):
    base = Path(h5file)
    stop_txt = "end" if snap_stop is None else str(int(snap_stop))
    return str(
        base.with_name(
            f"{base.stem}_{field}_xz_cache_r{r_cut:g}_s{snap_start}_{stop_txt}_k{stride}.h5"
        )
    )



def _xz_cache_is_valid(
    cache_file,
    h5file,
    field,
    r_cut,
    snap_start,
    snap_stop,
    stride,
    angle_name=None,
):
    if not os.path.exists(cache_file):
        return False
    if not os.path.exists(h5file):
        return False

    try:
        src_mtime = os.path.getmtime(h5file)
        with h5py.File(cache_file, "r") as h:
            if h.attrs.get("source_file", "") != os.path.abspath(h5file):
                return False
            if h.attrs.get("field", "") != field:
                return False
            if h.attrs.get("plane", "") != "xz":
                return False
            if float(h.attrs.get("r_cut", np.nan)) != float(r_cut):
                return False
            if int(h.attrs.get("snap_start", -1)) != int(snap_start):
                return False
            stop_attr = h.attrs.get("snap_stop", "None")
            stop_attr = None if stop_attr == "None" else int(stop_attr)
            if stop_attr != snap_stop:
                return False
            if int(h.attrs.get("stride", -1)) != int(stride):
                return False
            if abs(float(h.attrs.get("source_mtime", -1.0)) - float(src_mtime)) > 1e-9:
                return False
            if angle_name is not None and h.attrs.get("angle_name", "") != angle_name:
                return False
        return True
    except Exception:
        return False



def _build_xz_cache(
    h5file,
    cache_file,
    field,
    r_cut,
    snap_start=0,
    snap_stop=None,
    stride=1,
    angle_name=None,
    verbose=True,
):
    os.makedirs(os.path.dirname(os.path.abspath(cache_file)) or ".", exist_ok=True)

    with h5py.File(h5file, "r") as f:
        r = np.asarray(f["r"][...], float)
        ang_name = _guess_xz_angle_name(f, preferred=angle_name)
        ang = np.asarray(f[ang_name][...], float)
        U = f[field]
        tarr = np.asarray(f["t"][...], float) if "t" in f else None

        Nt = U.shape[0]
        if snap_stop is None or snap_stop > Nt:
            snap_stop = Nt

        mask_r = (r <= r_cut)
        if not np.any(mask_r):
            raise ValueError(f"No hay puntos con r <= {r_cut}. max(r)={r.max()}")

        r_sel = r[mask_r]
        X, Z = _build_xz_mesh(r_sel, ang)
        snap_ids = np.arange(snap_start, snap_stop, stride, dtype=int)
        ns = len(snap_ids)

        vmax_re = 0.0
        vmax_im = 0.0
        vmax_abs = 0.0

        with h5py.File(cache_file, "w") as h:
            h.attrs["source_file"] = os.path.abspath(h5file)
            h.attrs["source_mtime"] = float(os.path.getmtime(h5file))
            h.attrs["field"] = field
            h.attrs["plane"] = "xz"
            h.attrs["angle_name"] = ang_name
            h.attrs["r_cut"] = float(r_cut)
            h.attrs["snap_start"] = int(snap_start)
            h.attrs["snap_stop"] = "None" if snap_stop is None else int(snap_stop)
            h.attrs["stride"] = int(stride)

            h.create_dataset("r_sel", data=r_sel)
            h.create_dataset("angle", data=ang)
            h.create_dataset("X", data=X)
            h.create_dataset("Z", data=Z)
            h.create_dataset("snap_ids", data=snap_ids)
            if tarr is not None:
                h.create_dataset("t_sel", data=tarr[snap_ids])

            dset = h.create_dataset(
                "frames",
                shape=(ns, len(r_sel), len(ang)),
                dtype=U.dtype,
                chunks=(1, len(r_sel), len(ang)),
            )

            for j, s in enumerate(snap_ids):
                F = np.asarray(U[s, mask_r, :])
                dset[j, :, :] = F

                vmax_re = max(vmax_re, float(np.nanmax(np.abs(np.real(F)))))
                vmax_im = max(vmax_im, float(np.nanmax(np.abs(np.imag(F)))))
                vmax_abs = max(vmax_abs, float(np.nanmax(np.abs(F))))

                if verbose and (j == 0 or j % 20 == 0):
                    log(
                        f"[cache XZ] {j+1}/{ns} | s={s} | "
                        f"vmax_re={vmax_re:.3g} "
                        f"vmax_im={vmax_im:.3g} "
                        f"vmax_abs={vmax_abs:.3g}"
                    )

            h.attrs["vmax_re"] = float(vmax_re)
            h.attrs["vmax_im"] = float(vmax_im)
            h.attrs["vmax_abs"] = float(vmax_abs)

    return cache_file, vmax_re, vmax_im, vmax_abs



def _ensure_xz_cache(
    h5file,
    field,
    r_cut,
    snap_start=0,
    snap_stop=None,
    stride=1,
    cache_file=None,
    rebuild=False,
    angle_name=None,
    verbose=True,
):
    if cache_file is None:
        cache_file = _default_xz_cache_file(h5file, field, r_cut, snap_start, snap_stop, stride)

    if (not rebuild) and _xz_cache_is_valid(
        cache_file, h5file, field, r_cut, snap_start, snap_stop, stride, angle_name=angle_name
    ):
        if verbose:
            log(f"[cache] usando cache XZ existente: {cache_file}")
        with h5py.File(cache_file, "r") as h:
            return (
                cache_file,
                float(h.attrs["vmax_re"]),
                float(h.attrs["vmax_im"]),
                float(h.attrs["vmax_abs"]),
            )

    if verbose:
        log(f"[cache] construyendo cache XZ: {cache_file}")
    return _build_xz_cache(
        h5file,
        cache_file,
        field,
        r_cut,
        snap_start=snap_start,
        snap_stop=snap_stop,
        stride=stride,
        angle_name=angle_name,
        verbose=verbose,
    )



def save_all_snapshots_three_panels_fixed_scale_xz(
    h5file="planes_sum.h5",
    out_dir="frames_phi_xz",
    field="phi_xz",
    r_cut=120.0,
    cmap="jet",
    dpi=150,
    snap_start=0,
    snap_stop=None,
    stride=1,
    clip_frac=1.0,
    verbose=True,
    use_xz_cache=True,
    xz_cache_file=None,
    rebuild_xz_cache=False,
    angle_name=None,
):
    os.makedirs(out_dir, exist_ok=True)

    if use_xz_cache:
        cache_file, vmax_re, vmax_im, vmax_abs = _ensure_xz_cache(
            h5file,
            field,
            r_cut,
            snap_start=snap_start,
            snap_stop=snap_stop,
            stride=stride,
            cache_file=xz_cache_file,
            rebuild=rebuild_xz_cache,
            angle_name=angle_name,
            verbose=verbose,
        )

        with h5py.File(cache_file, "r") as h:
            X = np.asarray(h["X"][...])
            Z = np.asarray(h["Z"][...])
            frames = h["frames"]
            snap_ids = np.asarray(h["snap_ids"][...], dtype=int)
            tsel = np.asarray(h["t_sel"][...], float) if "t_sel" in h else None

            vmax_re *= clip_frac
            vmax_im *= clip_frac
            vmax_abs *= clip_frac

            if verbose:
                log(
                    f"[fixed XZ] vlims: "
                    f"Re=[{-vmax_re:.4g},{vmax_re:.4g}]  "
                    f"Im=[{-vmax_im:.4g},{vmax_im:.4g}]  "
                    f"|.|=[0,{vmax_abs:.4g}]"
                )

            for j, s in enumerate(snap_ids):
                F = np.asarray(frames[j, :, :])
                ReF = np.real(F)
                ImF = np.imag(F)
                AbsF = np.abs(F)

                fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

                im0 = axs[0].pcolormesh(
                    X, Z, ReF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=-vmax_re,
                    vmax=vmax_re,
                    rasterized=True,
                )
                axs[0].set_title(r"$\Re(\Phi)$")
                axs[0].set_xlabel("x")
                axs[0].set_ylabel("z")
                axs[0].set_aspect("equal")
                axs[0].set_xlim(-r_cut, r_cut)
                axs[0].set_ylim(-r_cut, r_cut)

                im1 = axs[1].pcolormesh(
                    X, Z, ImF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=-vmax_im,
                    vmax=vmax_im,
                    rasterized=True,
                )
                axs[1].set_title(r"$\Im(\Phi)$")
                axs[1].set_xlabel("x")
                axs[1].set_ylabel("z")
                axs[1].set_aspect("equal")
                axs[1].set_xlim(-r_cut, r_cut)
                axs[1].set_ylim(-r_cut, r_cut)

                im2 = axs[2].pcolormesh(
                    X, Z, AbsF,
                    shading="nearest",
                    cmap=cmap,
                    vmin=0.0,
                    vmax=vmax_abs,
                    rasterized=True,
                )
                axs[2].set_title(r"$|\Phi|$")
                axs[2].set_xlabel("x")
                axs[2].set_ylabel("z")
                axs[2].set_aspect("equal")
                axs[2].set_xlim(-r_cut, r_cut)
                axs[2].set_ylim(-r_cut, r_cut)

                fig.colorbar(im0, ax=axs[0], shrink=0.85)
                fig.colorbar(im1, ax=axs[1], shrink=0.85)
                fig.colorbar(im2, ax=axs[2], shrink=0.85)

                supt = f"snap={s}" if tsel is None else f"t={tsel[j]:.6g}"
                fig.suptitle(supt, y=1.02)

                out_path = os.path.join(out_dir, f"{field}_xz_{s:05d}.png")
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)

                if verbose and (j == 0 or j % 10 == 0):
                    log(f"[ok] {out_path}")

    else:
        with h5py.File(h5file, "r") as f:
            r = np.asarray(f["r"][...], float)
            ang_name = _guess_xz_angle_name(f, preferred=angle_name)
            ang = np.asarray(f[ang_name][...], float)
            U = f[field]
            tarr = np.asarray(f["t"][...], float) if "t" in f else None

            Nt = U.shape[0]
            if snap_stop is None or snap_stop > Nt:
                snap_stop = Nt

            mask_r = (r <= r_cut)
            if not np.any(mask_r):
                raise ValueError(f"No hay puntos con r <= {r_cut}. max(r)={r.max()}")

            r_sel = r[mask_r]
            X, Z = _build_xz_mesh(r_sel, ang)

            vmax_re = 0.0
            vmax_im = 0.0
            vmax_abs = 0.0
            for s in range(snap_start, snap_stop, stride):
                F = U[s, mask_r, :]
                vmax_re = max(vmax_re, float(np.nanmax(np.abs(np.real(F)))))
                vmax_im = max(vmax_im, float(np.nanmax(np.abs(np.imag(F)))))
                vmax_abs = max(vmax_abs, float(np.nanmax(np.abs(F))))

            vmax_re *= clip_frac
            vmax_im *= clip_frac
            vmax_abs *= clip_frac

            for s in range(snap_start, snap_stop, stride):
                F = U[s, mask_r, :]
                ReF = np.real(F)
                ImF = np.imag(F)
                AbsF = np.abs(F)

                fig, axs = plt.subplots(1, 3, figsize=(16, 5), constrained_layout=True)

                im0 = axs[0].pcolormesh(X, Z, ReF, shading="nearest", cmap=cmap,
                                        vmin=-vmax_re, vmax=vmax_re, rasterized=True)
                im1 = axs[1].pcolormesh(X, Z, ImF, shading="nearest", cmap=cmap,
                                        vmin=-vmax_im, vmax=vmax_im, rasterized=True)
                im2 = axs[2].pcolormesh(X, Z, AbsF, shading="nearest", cmap=cmap,
                                        vmin=0.0, vmax=vmax_abs, rasterized=True)

                axs[0].set_title(r"$\Re(\Phi)$")
                axs[1].set_title(r"$\Im(\Phi)$")
                axs[2].set_title(r"$|\Phi|$")
                for ax in axs:
                    ax.set_xlabel("x")
                    ax.set_ylabel("z")
                    ax.set_aspect("equal")
                    ax.set_xlim(-r_cut, r_cut)
                    ax.set_ylim(-r_cut, r_cut)

                fig.colorbar(im0, ax=axs[0], shrink=0.85)
                fig.colorbar(im1, ax=axs[1], shrink=0.85)
                fig.colorbar(im2, ax=axs[2], shrink=0.85)

                supt = f"snap={s}" if tarr is None else f"t={tarr[s]:.6g}"
                fig.suptitle(supt, y=1.02)
                out_path = os.path.join(out_dir, f"{field}_xz_{s:05d}.png")
                fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
                plt.close(fig)

                if verbose and (s == snap_start or (s - snap_start) % (10 * stride) == 0):
                    log(f"[ok] {out_path}")

    if verbose:
        log(f"[done] frames XZ en: {out_dir}")



def generate_xz_frames(cfg):
    log("[STEP] Generando frames XZ para todos los tiempos")

    if not os.path.exists(cfg["planes_h5"]):
        raise FileNotFoundError(f"No existe planes_sum.h5: {cfg['planes_h5']}")

    save_all_snapshots_three_panels_fixed_scale_xz(
        h5file=cfg["planes_h5"],
        out_dir=cfg["xz_out_dir"],
        field=cfg["field_name_xz"],
        r_cut=cfg["r_cut"],
        cmap=cfg["cmap"],
        dpi=cfg["dpi"],
        clip_frac=cfg["clip_frac"],
        stride=1,
        verbose=True,
        use_xz_cache=cfg.get("use_xz_cache", True),
        xz_cache_file=cfg.get("xz_cache_file", None),
        rebuild_xz_cache=cfg.get("rebuild_xz_cache", False),
        angle_name=cfg.get("xz_angle_name", None),
    )

# =========================================================
# PARTE C: CREAR VIDEOS
# =========================================================


def make_movie_from_frames(frame_pattern, output_file, fps=10, crf=18, preset="medium"):
    ensure_parent(output_file)

    cmd = [
        "ffmpeg",
        "-y",
        "-framerate", str(fps),
        "-i", frame_pattern,
        "-c:v", "libx264",
        "-pix_fmt", "yuv420p",
        "-crf", str(crf),
        "-preset", preset,
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2",
        output_file,
    ]
    run_cmd(cmd)


def _maybe_make_movie(frame_pattern, output_file, fps, crf, preset, label):
    import glob
    first_frame = frame_pattern.replace("%04d", "*").replace("%05d", "*")
    matches = sorted(glob.glob(first_frame))
    if not matches:
        log(f"[SKIP] No hay frames para {label}: {first_frame}")
        return
    make_movie_from_frames(frame_pattern, output_file, fps=fps, crf=crf, preset=preset)


def generate_movies(cfg):
    log("[STEP] Generando películas con ffmpeg")

    pk_pattern = os.path.join(cfg["spectra_out_dir"], "pk", "P_k_t%04d.png")
    heat_pattern = os.path.join(cfg["spectra_out_dir"], "heatmap", "P_ln_t%04d.png")
    xy_pattern = os.path.join(cfg["xy_out_dir"], f"{cfg['field_name_xy']}_xy_%05d.png")
    xz_pattern = os.path.join(cfg["xz_out_dir"], f"{cfg['field_name_xz']}_xz_%05d.png")

    _maybe_make_movie(
        pk_pattern, cfg["movie_pk"],
        fps=cfg["fps"], crf=cfg["ffmpeg_crf"], preset=cfg["ffmpeg_preset"],
        label="P(k)"
    )

    _maybe_make_movie(
        heat_pattern, cfg["movie_heat"],
        fps=cfg["fps"], crf=cfg["ffmpeg_crf"], preset=cfg["ffmpeg_preset"],
        label="heatmap"
    )

    _maybe_make_movie(
        xy_pattern, cfg["movie_xy"],
        fps=cfg["fps"], crf=cfg["ffmpeg_crf"], preset=cfg["ffmpeg_preset"],
        label="XY"
    )

    _maybe_make_movie(
        xz_pattern, cfg["movie_xz"],
        fps=cfg["fps"], crf=cfg["ffmpeg_crf"], preset=cfg["ffmpeg_preset"],
        label="XZ"
    )


# =========================================================
# MAIN
# =========================================================


def main():
    try:
        if RUN_PREPROCESS:
            copy_and_run_preprocessing(CONFIG, force=FORCE_PREPROCESS)
        else:
            log("[SKIP] RUN_PREPROCESS=False")

        if RUN_SPECTRA:
            generate_spectra_all_times(CONFIG)
        else:
            log("[SKIP] RUN_SPECTRA=False")

        if RUN_XY:
            generate_xy_frames(CONFIG)
        else:
            log("[SKIP] RUN_XY=False")

        if RUN_XZ:
            generate_xz_frames(CONFIG)
        else:
            log("[SKIP] RUN_XZ=False")

        if RUN_MOVIES:
            generate_movies(CONFIG)
        else:
            log("[SKIP] RUN_MOVIES=False")

        log("")
        log("[DONE] Pipeline completo terminado.")
        log(f"[OUT] Heatmaps/P(k): {CONFIG['spectra_out_dir']}")
        log(f"[OUT] XY frames     : {CONFIG['xy_out_dir']}")
        log(f"[OUT] XZ frames     : {CONFIG['xz_out_dir']}")
        log(f"[OUT] Video Pk      : {CONFIG['movie_pk']}")
        log(f"[OUT] Video heatmap : {CONFIG['movie_heat']}")
        log(f"[OUT] Video XY      : {CONFIG['movie_xy']}")
        log(f"[OUT] Video XZ      : {CONFIG['movie_xz']}")

    except subprocess.CalledProcessError as e:
        log(f"[ERROR] Falló un comando externo con código {e.returncode}")
        raise
    except Exception as e:
        log(f"[ERROR] {e}")
        raise


if __name__ == "__main__":
    main()
