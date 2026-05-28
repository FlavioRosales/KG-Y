#!/usr/bin/env python3
"""
Reconstruye proyecciones del campo escalar a partir de modos esfericos.

Ejecucion normal:
    python sum_lm.py

Por defecto lee modes.h5 y escribe planes_sum.h5 en la carpeta actual.

Cambio importante de esta version:
    Por defecto reconstruye solo r >= 1.0 R_s. Esto evita que puntos dentro
    del horizonte o cerca de la frontera de excision, donde pueden aparecer
    valores no finitos o desbordamientos numericos, contaminen la proyeccion.

Salida:
    planes_sum.h5
        r, theta, phi_grid, t
        phi_xy, pi_xy      shape (Nt_snap, Nr_out, Nphi)
        phi_xz, pi_xz      shape (Nt_snap, Nr_out, Ntheta)
        N, F               suma en tiempos de snapshot
        N_diagnostic, F_diagnostic si hay series diagnosticas mas finas
"""

import argparse
import re
from pathlib import Path

import h5py
import numpy as np

try:
    from scipy.special import sph_harm_y as _sph_harm_y
except ImportError:  # scipy < 1.15
    _sph_harm_y = None
    from scipy.special import sph_harm as _sph_harm_old


MODE_RE = re.compile(r"^mode_(\d+)$")
COMPLEX_DTYPE = np.complex64
FLOAT32_MAX = np.finfo(np.float32).max
DIAG_CLIP = 1.0e250

TIME_CANDIDATES = (
    "t_N", "tN", "time_N", "timeN", "N_t", "N_time",
    "t_F", "tF", "time_F", "timeF", "F_t", "F_time",
    "t_diag", "tdiag", "diag_t", "time_diag", "diagnostic_t",
    "t_full", "time_full", "time", "times",
)


def n_to_ellm(n: int) -> tuple[int, int]:
    ell = int(np.sqrt(n))
    rem = n - ell * ell
    m = int(rem) - ell
    return ell, m


def ellm_to_n(ell: int, m: int) -> int:
    return ell * ell + (ell + m)


def find_modes(fin: h5py.File) -> list[int]:
    modes = []
    for key, obj in fin.items():
        match = MODE_RE.match(key)
        if match and isinstance(obj, h5py.Group):
            modes.append(int(match.group(1)))
    return sorted(modes)


def nyquist_grids(ell_max: int) -> tuple[np.ndarray, np.ndarray]:
    ntheta = ell_max + 1
    nphi = 2 * ell_max + 1
    theta = (np.arange(ntheta) + 0.5) * np.pi / ntheta
    phi = np.arange(nphi) * 2.0 * np.pi / nphi
    return theta.astype(np.float64), phi.astype(np.float64)


def sph_harm_stable(ell: int, m: int, theta_polar: np.ndarray, phi_azimuth: np.ndarray) -> np.ndarray:
    """Evalua Y_{ell m}(theta, phi) usando armonicos normalizados de SciPy."""
    if _sph_harm_y is not None:
        # scipy.special.sph_harm_y(n, m, theta, phi), theta polar, phi azimutal.
        return _sph_harm_y(ell, m, theta_polar, phi_azimuth)
    # scipy.special.sph_harm(m, n, theta, phi), theta azimutal, phi polar.
    return _sph_harm_old(m, ell, phi_azimuth, theta_polar)


def assert_finite_array(arr: np.ndarray, name: str) -> None:
    bad = ~np.isfinite(arr)
    if np.any(bad):
        idx = np.argwhere(bad)
        first = tuple(int(i) for i in idx[0])
        raise RuntimeError(
            f"{name} contiene {idx.shape[0]} valores no finitos. Primer indice no finito: {first}."
        )


def sanitize_real_pair(real: np.ndarray, imag: np.ndarray, name: str, strict: bool = False) -> tuple[np.ndarray, int]:
    """Limpia NaN/Inf y valores que desbordarian al convertir a float32."""
    bad = (
        ~np.isfinite(real)
        | ~np.isfinite(imag)
        | (np.abs(real) > FLOAT32_MAX)
        | (np.abs(imag) > FLOAT32_MAX)
    )
    nbad = int(np.count_nonzero(bad))

    if nbad:
        first = tuple(int(i) for i in np.argwhere(bad)[0])
        msg = (
            f"{name}: {nbad} valores no finitos o fuera de rango float32. "
            f"Primer indice local: {first}."
        )
        if strict:
            raise RuntimeError(msg)
        print(f"[WARN] {msg} Se reemplazan por 0 para la reconstruccion grafica.")

        real = np.array(real, copy=True)
        imag = np.array(imag, copy=True)
        real[bad] = 0.0
        imag[bad] = 0.0

    z = real.astype(np.float32, copy=False) + 1j * imag.astype(np.float32, copy=False)
    return z.astype(COMPLEX_DTYPE, copy=False), nbad


def Ylm_equator(ell: int, phi_grid: np.ndarray) -> np.ndarray:
    theta = np.full_like(phi_grid, 0.5 * np.pi, dtype=np.float64)
    Y = np.empty((2 * ell + 1, phi_grid.size), dtype=np.complex128)
    for row, m in enumerate(range(-ell, ell + 1)):
        Y[row, :] = sph_harm_stable(ell, m, theta, phi_grid)
    assert_finite_array(Y, f"Ylm_equator(ell={ell})")
    return Y.astype(COMPLEX_DTYPE, copy=False)


def Ylm_meridian(ell: int, theta_grid: np.ndarray, phi0: float = 0.0) -> np.ndarray:
    phi = np.full_like(theta_grid, phi0, dtype=np.float64)
    Y = np.empty((2 * ell + 1, theta_grid.size), dtype=np.complex128)
    for row, m in enumerate(range(-ell, ell + 1)):
        Y[row, :] = sph_harm_stable(ell, m, theta_grid, phi)
    assert_finite_array(Y, f"Ylm_meridian(ell={ell})")
    return Y.astype(COMPLEX_DTYPE, copy=False)


def as_1d_series(ds: h5py.Dataset, name: str = "dataset") -> np.ndarray:
    arr = np.asarray(ds[...])
    arr = np.squeeze(arr)
    if arr.ndim != 1:
        raise RuntimeError(
            f"{name} tiene shape {ds.shape}; despues de squeeze queda {arr.shape}, no es una serie 1D"
        )
    return arr


def find_matching_time_grid(
    fin: h5py.File,
    group: h5py.Group,
    length: int,
    t_snap: np.ndarray,
    label: str,
) -> tuple[np.ndarray, str]:
    for key in TIME_CANDIDATES:
        if key in group and isinstance(group[key], h5py.Dataset):
            arr = as_1d_series(group[key], name=f"{group.name}/{key}")
            if arr.size == length:
                return arr.astype(np.float64), f"{group.name}/{key}"

    for key in TIME_CANDIDATES:
        if key in fin and isinstance(fin[key], h5py.Dataset):
            arr = as_1d_series(fin[key], name=f"/{key}")
            if arr.size == length:
                return arr.astype(np.float64), f"/{key}"

    if length == t_snap.size:
        return t_snap.astype(np.float64), "snapshot t"

    return np.linspace(float(t_snap[0]), float(t_snap[-1]), length, dtype=np.float64), (
        f"linspace(t[0], t[-1], {length}) asumido para {label}"
    )


def series_on_snapshots(series: np.ndarray, series_t: np.ndarray, t_snap: np.ndarray, name: str) -> np.ndarray:
    series = np.asarray(series, dtype=np.float64)
    series_t = np.asarray(series_t, dtype=np.float64)

    if series.size != series_t.size:
        raise RuntimeError(f"{name}: serie y tiempo tienen longitudes distintas: {series.size} vs {series_t.size}")

    if series.size == t_snap.size and np.allclose(series_t, t_snap, rtol=0.0, atol=1.0e-14):
        return series.copy()

    if series_t[0] > series_t[-1]:
        series_t = series_t[::-1]
        series = series[::-1]

    if np.any(np.diff(series_t) < 0):
        raise RuntimeError(f"{name}: la malla temporal diagnostica no es monotona")

    return np.interp(t_snap, series_t, series)


def validate_input(fin: h5py.File, nmodes: list[int], Nr: int, Nt_snap: int, t_snap: np.ndarray) -> None:
    required = ("t", "phi", "pi", "N", "F")
    for nmode in nmodes:
        group_name = f"mode_{nmode:05d}"
        gm = fin[group_name]

        for name in required:
            if name not in gm:
                raise RuntimeError(f"{group_name} no contiene el dataset requerido '{name}'")

        if gm["phi"].shape != (2, Nr, Nt_snap):
            raise RuntimeError(f"{group_name}/phi tiene shape {gm['phi'].shape}; se esperaba {(2, Nr, Nt_snap)}")

        if gm["pi"].shape != (2, Nr, Nt_snap):
            raise RuntimeError(f"{group_name}/pi tiene shape {gm['pi'].shape}; se esperaba {(2, Nr, Nt_snap)}")

        t_mode = as_1d_series(gm["t"], name=f"{group_name}/t")
        if t_mode.size != Nt_snap:
            raise RuntimeError(f"{group_name}/t tiene longitud {t_mode.size}; se esperaba Nt_snap={Nt_snap}")

        if not np.allclose(t_mode, t_snap, rtol=0.0, atol=1.0e-14):
            raise RuntimeError(f"{group_name}/t no coincide con la malla temporal de referencia")

        _ = as_1d_series(gm["N"], name=f"{group_name}/N")
        _ = as_1d_series(gm["F"], name=f"{group_name}/F")


def group_modes_by_ell(nmodes: list[int]) -> tuple[int, dict[int, list[tuple[int, int]]]]:
    ell_to_present: dict[int, list[tuple[int, int]]] = {}
    for nmode in nmodes:
        ell, m = n_to_ellm(nmode)
        expected = ellm_to_n(ell, m)
        if expected != nmode:
            raise RuntimeError(f"Mapeo inconsistente para n={nmode}: ell={ell}, m={m}, n_rec={expected}")
        if abs(m) > ell:
            raise RuntimeError(f"Modo invalido n={nmode}: ell={ell}, m={m}")
        ell_to_present.setdefault(ell, []).append((m, nmode))

    for ell in ell_to_present:
        ell_to_present[ell].sort(key=lambda item: item[0])

    return max(ell_to_present), ell_to_present


def report_missing_modes(ell_max: int, ell_to_present: dict[int, list[tuple[int, int]]], verbose: bool) -> None:
    total_missing = 0
    incomplete_ells = 0

    for ell in range(ell_max + 1):
        present = ell_to_present.get(ell, [])
        present_m = {m for m, _ in present}
        missing_m = [m for m in range(-ell, ell + 1) if m not in present_m]
        if missing_m:
            total_missing += len(missing_m)
            incomplete_ells += 1
            if verbose:
                print(
                    f"[WARN] ell={ell}: faltan {len(missing_m)} valores de m; "
                    f"se sumaran solo los disponibles. Primeros faltantes: {missing_m[:10]}"
                )

    if total_missing > 0:
        print(
            f"[INFO] Hay {incomplete_ells} valores de ell con modos m faltantes "
            f"({total_missing} modos ausentes en total). Se reconstruye solo con los modos disponibles."
        )


def clean_diagnostic_series(series: np.ndarray, name: str) -> np.ndarray:
    series = np.asarray(series, dtype=np.float64)
    bad = ~np.isfinite(series) | (np.abs(series) > DIAG_CLIP)
    nbad = int(np.count_nonzero(bad))
    if nbad:
        first = int(np.argwhere(bad)[0, 0])
        print(
            f"[WARN] {name}: {nbad} valores no finitos o muy grandes en diagnostico. "
            f"Primer indice: {first}. Se reemplazan por 0 para evitar overflow en la suma."
        )
        series = np.array(series, copy=True)
        series[bad] = 0.0
    return series


def sum_diagnostic_series(
    fin: h5py.File,
    nmodes: list[int],
    dataset_name: str,
    t_snap: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, str]:
    first_group = fin[f"mode_{nmodes[0]:05d}"]
    first_series = as_1d_series(first_group[dataset_name], name=f"mode_{nmodes[0]:05d}/{dataset_name}")
    full_len = first_series.size
    full_time, time_source = find_matching_time_grid(fin, first_group, full_len, t_snap, dataset_name)

    total_full = np.zeros(full_len, dtype=np.float64)

    for nmode in nmodes:
        gm = fin[f"mode_{nmode:05d}"]
        series = as_1d_series(gm[dataset_name], name=f"mode_{nmode:05d}/{dataset_name}")
        if series.size != full_len:
            raise RuntimeError(
                f"mode_{nmode:05d}/{dataset_name} tiene longitud {series.size}; "
                f"se esperaba {full_len}, como en mode_{nmodes[0]:05d}/{dataset_name}"
            )

        time_i, source_i = find_matching_time_grid(fin, gm, series.size, t_snap, dataset_name)
        if time_i.size != full_time.size or not np.allclose(time_i, full_time, rtol=0.0, atol=1.0e-12):
            raise RuntimeError(
                f"La malla temporal de {dataset_name} en mode_{nmode:05d} no coincide con la referencia.\n"
                f"Referencia: {time_source}\nActual: {source_i}"
            )

        series = clean_diagnostic_series(series, f"mode_{nmode:05d}/{dataset_name}")
        total_full += series

    total_snap = series_on_snapshots(total_full, full_time, t_snap, name=dataset_name)
    return total_full, full_time, total_snap, time_source


def read_complex_mode_block(
    gm: h5py.Group,
    dataset_name: str,
    r0_index: int,
    s0: int,
    s1: int,
    name: str,
    strict_input: bool,
) -> np.ndarray:
    real = gm[dataset_name][0, r0_index:, s0:s1]
    imag = gm[dataset_name][1, r0_index:, s0:s1]
    z, _ = sanitize_real_pair(real, imag, name=name, strict=strict_input)
    return z


def reconstruct_planes(
    in_file: str,
    out_file: str,
    t_block: int,
    r_min_output: float,
    strict_input: bool = False,
    verbose_missing: bool = False,
) -> None:
    if t_block < 1:
        raise ValueError("t_block debe ser >= 1")

    in_path = Path(in_file)
    if not in_path.exists():
        raise FileNotFoundError(
            f"No encontre {in_file!r} en la carpeta actual: {Path.cwd()}\n"
            "Ejecuta este script dentro de la carpeta que contiene modes.h5 "
            "o usa --in-file /ruta/a/modes.h5"
        )

    with h5py.File(in_path, "r") as fin:
        if "r" not in fin:
            raise RuntimeError("El archivo de entrada no contiene el dataset 'r'")

        r_full = np.asarray(fin["r"][...], dtype=np.float64)
        if r_full.ndim != 1 or r_full.size < 2:
            raise RuntimeError("El dataset r debe ser 1D y tener al menos dos puntos")
        if np.any(np.diff(r_full) <= 0):
            raise RuntimeError("El dataset r debe ser estrictamente creciente")

        Nr_full = r_full.size
        r0_index = int(np.searchsorted(r_full, r_min_output, side="left"))
        if r0_index >= Nr_full:
            raise RuntimeError(f"r_min_output={r_min_output} deja vacia la malla radial")

        r = r_full[r0_index:]
        Nr = r.size

        nmodes = find_modes(fin)
        if not nmodes:
            raise RuntimeError(f"No encontre grupos mode_XXXXX en {in_file}")

        gref = fin[f"mode_{nmodes[0]:05d}"]
        if "t" not in gref:
            raise RuntimeError(f"mode_{nmodes[0]:05d} no contiene el dataset 't'")

        t = as_1d_series(gref["t"], name=f"mode_{nmodes[0]:05d}/t").astype(np.float64)
        Nt = t.size

        validate_input(fin, nmodes, Nr_full, Nt, t)

        ell_max, ell_to_present = group_modes_by_ell(nmodes)
        theta_grid, phi_grid = nyquist_grids(ell_max)
        Nth, Nph = theta_grid.size, phi_grid.size

        print(f"[INFO] Input:  {in_file}")
        print(f"[INFO] Output: {out_file}")
        print(f"[INFO] Nr_full={Nr_full}, Nr_out={Nr}, Nt_snap={Nt}, ell_max={ell_max}, Ntheta={Nth}, Nphi={Nph}")
        print(f"[INFO] r_min_output={r[0]:.12g} R_s, puntos descartados al inicio={r0_index}")
        print(f"[INFO] T_BLOCK={t_block}, dtype={np.dtype(COMPLEX_DTYPE).name}")

        report_missing_modes(ell_max, ell_to_present, verbose=verbose_missing)
        mode_groups = {n: fin[f"mode_{n:05d}"] for n in nmodes}

        print("[INFO] Summing diagnostic series N and F...")
        N_full, t_N_full, N_snap, N_time_source = sum_diagnostic_series(fin, nmodes, "N", t)
        F_full, t_F_full, F_snap, F_time_source = sum_diagnostic_series(fin, nmodes, "F", t)
        print(f"[INFO] N length={N_full.size}; time source: {N_time_source}")
        print(f"[INFO] F length={F_full.size}; time source: {F_time_source}")
        if N_full.size != Nt:
            print("[INFO] N se guardara tambien como N_diagnostic y N se interpolara a los snapshots de phi/pi.")
        if F_full.size != Nt:
            print("[INFO] F se guardara tambien como F_diagnostic y F se interpolara a los snapshots de phi/pi.")

        print("[INFO] Precomputing Ylm caches...")
        Yeq_cache = [None] * (ell_max + 1)
        Yme_cache = [None] * (ell_max + 1)
        for ell in range(ell_max + 1):
            Yeq_cache[ell] = Ylm_equator(ell, phi_grid)
            Yme_cache[ell] = Ylm_meridian(ell, theta_grid, phi0=0.0)

        with h5py.File(out_file, "w") as fout:
            fout.create_dataset("r", data=r)
            fout.create_dataset("r_full_input", data=r_full)
            fout.create_dataset("theta", data=theta_grid)
            fout.create_dataset("phi_grid", data=phi_grid)
            fout.create_dataset("t", data=t)

            fout.attrs["description"] = "Reconstruccion modal en planos xy y xz a partir de modos disponibles."
            fout.attrs["mode_indexing"] = "n = ell**2 + ell + m"
            fout.attrs["xy_plane"] = "theta = pi/2, phi = phi_grid"
            fout.attrs["xz_plane"] = "meridiano phi = 0; para m=0 representa el plano xz axisimetrico"
            fout.attrs["ell_max"] = ell_max
            fout.attrs["complex_dtype"] = np.dtype(COMPLEX_DTYPE).name
            fout.attrs["r_min_output_requested"] = float(r_min_output)
            fout.attrs["r_min_output_actual"] = float(r[0])
            fout.attrs["r0_index_input"] = int(r0_index)
            fout.attrs["N_time_source"] = N_time_source
            fout.attrs["F_time_source"] = F_time_source
            fout.attrs["N_dataset"] = "N is evaluated/interpolated on snapshot times t"
            fout.attrs["F_dataset"] = "F is evaluated/interpolated on snapshot times t"

            chunk_r = min(Nr, 512)
            phi_xy_ds = fout.create_dataset("phi_xy", shape=(Nt, Nr, Nph), dtype=COMPLEX_DTYPE, chunks=(1, chunk_r, Nph))
            pi_xy_ds = fout.create_dataset("pi_xy", shape=(Nt, Nr, Nph), dtype=COMPLEX_DTYPE, chunks=(1, chunk_r, Nph))
            phi_xz_ds = fout.create_dataset("phi_xz", shape=(Nt, Nr, Nth), dtype=COMPLEX_DTYPE, chunks=(1, chunk_r, Nth))
            pi_xz_ds = fout.create_dataset("pi_xz", shape=(Nt, Nr, Nth), dtype=COMPLEX_DTYPE, chunks=(1, chunk_r, Nth))

            for s0 in range(0, Nt, t_block):
                s1 = min(Nt, s0 + t_block)
                bs = s1 - s0
                print(f"[INFO] snapshots {s0}..{s1 - 1}")

                phi_xy_blk = np.zeros((bs, Nr, Nph), dtype=COMPLEX_DTYPE)
                pi_xy_blk = np.zeros_like(phi_xy_blk)
                phi_xz_blk = np.zeros((bs, Nr, Nth), dtype=COMPLEX_DTYPE)
                pi_xz_blk = np.zeros_like(phi_xz_blk)

                for ell in range(ell_max + 1):
                    present = ell_to_present.get(ell, [])
                    if not present:
                        continue

                    row_idx = np.array([m + ell for m, _ in present], dtype=np.int64)
                    Yeq_sel = Yeq_cache[ell][row_idx, :]
                    Yme_sel = Yme_cache[ell][row_idx, :]

                    nm_present = len(present)
                    cphi = np.empty((nm_present, Nr, bs), dtype=COMPLEX_DTYPE)
                    cpi = np.empty((nm_present, Nr, bs), dtype=COMPLEX_DTYPE)

                    for j, (_m, nmode) in enumerate(present):
                        gm = mode_groups[nmode]
                        cphi[j, :, :] = read_complex_mode_block(
                            gm, "phi", r0_index, s0, s1,
                            name=f"mode_{nmode:05d}/phi, snapshots {s0}:{s1}, r>=r_min_output",
                            strict_input=strict_input,
                        )
                        cpi[j, :, :] = read_complex_mode_block(
                            gm, "pi", r0_index, s0, s1,
                            name=f"mode_{nmode:05d}/pi, snapshots {s0}:{s1}, r>=r_min_output",
                            strict_input=strict_input,
                        )

                    assert_finite_array(cphi, f"coeficientes phi ell={ell}, snapshots {s0}:{s1}")
                    assert_finite_array(cpi, f"coeficientes pi ell={ell}, snapshots {s0}:{s1}")

                    tmp_xy_phi = np.tensordot(cphi, Yeq_sel, axes=(0, 0))
                    tmp_xy_pi = np.tensordot(cpi, Yeq_sel, axes=(0, 0))
                    tmp_xz_phi = np.tensordot(cphi, Yme_sel, axes=(0, 0))
                    tmp_xz_pi = np.tensordot(cpi, Yme_sel, axes=(0, 0))

                    phi_xy_blk += np.transpose(tmp_xy_phi, (1, 0, 2)).astype(COMPLEX_DTYPE, copy=False)
                    pi_xy_blk += np.transpose(tmp_xy_pi, (1, 0, 2)).astype(COMPLEX_DTYPE, copy=False)
                    phi_xz_blk += np.transpose(tmp_xz_phi, (1, 0, 2)).astype(COMPLEX_DTYPE, copy=False)
                    pi_xz_blk += np.transpose(tmp_xz_pi, (1, 0, 2)).astype(COMPLEX_DTYPE, copy=False)

                assert_finite_array(phi_xy_blk, f"phi_xy reconstruido snapshots {s0}:{s1}")
                assert_finite_array(pi_xy_blk, f"pi_xy reconstruido snapshots {s0}:{s1}")
                assert_finite_array(phi_xz_blk, f"phi_xz reconstruido snapshots {s0}:{s1}")
                assert_finite_array(pi_xz_blk, f"pi_xz reconstruido snapshots {s0}:{s1}")

                phi_xy_ds[s0:s1, :, :] = phi_xy_blk
                pi_xy_ds[s0:s1, :, :] = pi_xy_blk
                phi_xz_ds[s0:s1, :, :] = phi_xz_blk
                pi_xz_ds[s0:s1, :, :] = pi_xz_blk
                fout.flush()

            fout.create_dataset("N", data=N_snap)
            fout.create_dataset("F", data=F_snap)

            if N_full.size != Nt or not np.allclose(t_N_full, t, rtol=0.0, atol=1.0e-14):
                fout.create_dataset("t_N_diagnostic", data=t_N_full)
                fout.create_dataset("N_diagnostic", data=N_full)
            if F_full.size != Nt or not np.allclose(t_F_full, t, rtol=0.0, atol=1.0e-14):
                fout.create_dataset("t_F_diagnostic", data=t_F_full)
                fout.create_dataset("F_diagnostic", data=F_full)

            fout.flush()

    print(f"[OK] Generado {out_file}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruye proyecciones xy/xz a partir de coeficientes modales. "
            "Ejecucion normal: python sum_lm.py. "
            "Por defecto lee modes.h5 y escribe planes_sum.h5."
        )
    )
    parser.add_argument("--in-file", default="modes.h5", help="Archivo HDF5 de entrada. Default: modes.h5")
    parser.add_argument("--out-file", default="planes_sum.h5", help="Archivo HDF5 de salida. Default: planes_sum.h5")
    parser.add_argument("--t-block", type=int, default=4, help="Snapshots por bloque temporal. Usa 1 o 2 si falta RAM. Default: 4")
    parser.add_argument(
        "--r-min-output",
        type=float,
        default=1.0,
        help="Radio minimo guardado en la reconstruccion. Default: 1.0, el horizonte en unidades Rs.",
    )
    parser.add_argument(
        "--strict-input",
        action="store_true",
        help="Detiene la ejecucion si encuentra NaN/Inf o overflow en phi/pi dentro de r>=r-min-output.",
    )
    parser.add_argument("--verbose-missing", action="store_true", help="Imprime detalle de m faltantes por ell.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    reconstruct_planes(
        in_file=args.in_file,
        out_file=args.out_file,
        t_block=args.t_block,
        r_min_output=args.r_min_output,
        strict_input=args.strict_input,
        verbose_missing=args.verbose_missing,
    )


if __name__ == "__main__":
    main()
