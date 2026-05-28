#!/usr/bin/env python3
import os
import re
import glob
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl


# ============================================================
# Biyectiva n -> (ell,m)
# ============================================================
def n_to_ellm(n):
    ell = int(np.sqrt(n))
    r = n - ell * ell
    m = int(r) - ell
    return ell, m


# ============================================================
# Extraer k0 desde el nombre de la carpeta
# Usa los últimos 3 dígitos y los interpreta como xx.x
# Ejemplo: "..._005" -> 0.5, "..._125" -> 12.5
# ============================================================
def parse_k0_from_folder(folder_name):
    m = re.search(r'(\d{3})$', folder_name)
    if m is None:
        raise ValueError(f"No pude extraer los últimos 3 dígitos de: {folder_name}")
    return int(m.group(1)) / 10.0


# ============================================================
# Leer todos los modos desde modes.h5
# Devuelve:
#   mode_ids : lista de n
#   ells     : lista de ell
#   Ns       : arreglo (nmodes, nt)
#   Fs       : arreglo (nmodes, nt)
#   t        : tiempo lineal entre 0 y tmax
# ============================================================
def read_modes_h5(h5_path, tmax=1800.0):
    mode_ids = []
    ells = []
    Ns = []
    Fs = []

    with h5py.File(h5_path, "r") as f:
        mode_keys = sorted(
            [k for k in f.keys() if k.startswith("mode_")],
            key=lambda s: int(s.split("_")[1])
        )

        if not mode_keys:
            raise RuntimeError(f"No encontré grupos mode_xxxxx en {h5_path}")

        nt_ref = None

        for key in mode_keys:
            n = int(key.split("_")[1])
            ell, mm = n_to_ellm(n)

            g = f[key]
            if "N" not in g or "F" not in g:
                print(f"[WARN] {h5_path}:{key} no tiene N o F, se omite.")
                continue

            N = np.asarray(g["N"][:], dtype=np.float64)
            F = np.asarray(g["F"][:], dtype=np.float64)

            if N.ndim != 1 or F.ndim != 1:
                print(f"[WARN] {h5_path}:{key} tiene N o F con dimensión != 1, se omite.")
                continue

            if len(N) != len(F):
                print(f"[WARN] {h5_path}:{key} tiene longitudes distintas en N y F, se omite.")
                continue

            if nt_ref is None:
                nt_ref = len(N)
            elif len(N) != nt_ref:
                print(f"[WARN] {h5_path}:{key} tiene longitud distinta ({len(N)} vs {nt_ref}), se omite.")
                continue

            mode_ids.append(n)
            ells.append(ell)
            Ns.append(N)
            Fs.append(F)

    if len(Ns) == 0:
        raise RuntimeError(f"No encontré modos válidos en {h5_path}")

    Ns = np.asarray(Ns)
    Fs = np.asarray(Fs)
    t = np.linspace(0.0, tmax, Ns.shape[1])

    return mode_ids, ells, Ns, Fs, t


# ============================================================
# Ordenar por ell
# ============================================================
def sort_by_ell(mode_ids, ells, A):
    order = np.argsort(ells)
    mode_ids_s = np.asarray(mode_ids)[order]
    ells_s = np.asarray(ells)[order]
    A_s = np.asarray(A)[order, :]
    return mode_ids_s, ells_s, A_s


# ============================================================
# Graficar N/N(0)
# ============================================================
def plot_normalized_N(t, ells, Ns, k0, outpath):
    # Normalización robusta
    N0 = Ns[:, 0].copy()
    N0_safe = np.where(np.abs(N0) > 0.0, N0, np.nan)
    Ns_norm = Ns / N0_safe[:, None]

    fig, ax = plt.subplots(figsize=(7, 5))

    norm = mpl.colors.Normalize(vmin=np.min(ells), vmax=np.max(ells))
    cmap = plt.get_cmap("viridis")

    for ell, series in zip(ells, Ns_norm):
        ax.plot(t, series, linewidth=0.7, color=cmap(norm(ell)))

    ax.set_xlabel("t")
    ax.set_ylabel(r"$N_\ell(t)/N_\ell(0)$")
    ax.set_title(rf"Noether charge $k_0={k0}$")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"mode $\ell$")

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Graficar F = dN/dt
# ============================================================
def plot_F(t, ells, Fs, k0, outpath):
    fig, ax = plt.subplots(figsize=(7, 5))

    norm = mpl.colors.Normalize(vmin=np.min(ells), vmax=np.max(ells))
    cmap = plt.get_cmap("viridis")

    for ell, series in zip(ells, Fs):
        ax.plot(t, series, linewidth=0.7, color=cmap(norm(ell)))

    ax.set_xlabel("t")
    ax.set_ylabel(r"$\dot N_\ell(t)$")
    ax.set_title(rf"$\dot{{N}}_\ell(t)$  ($k_0={k0}$)")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label(r"mode $\ell$")

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)


# ============================================================
# Procesar una carpeta
# ============================================================
def process_folder(folder, tmax=1500.0):
    folder = os.path.abspath(folder)
    folder_name = os.path.basename(folder)
    h5_path = os.path.join(folder, "modes.h5")

    if not os.path.exists(h5_path):
        print(f"[WARN] No existe {h5_path}, se omite.")
        return

    try:
        k0 = parse_k0_from_folder(folder_name)
    except Exception as e:
        print(f"[WARN] No pude extraer k0 de {folder_name}: {e}")
        return

    print(f"[INFO] Procesando {folder_name}  ->  k0 = {k0}")

    mode_ids, ells, Ns, Fs, t = read_modes_h5(h5_path, tmax=tmax)

    mode_ids_s, ells_s, Ns_s = sort_by_ell(mode_ids, ells, Ns)
    _,         _,    Fs_s = sort_by_ell(mode_ids, ells, Fs)

    out1 = os.path.join(folder, "N_normalized_by_ell.png")
    out2 = os.path.join(folder, "F_by_ell.png")

    plot_normalized_N(t, ells_s, Ns_s, k0, out1)
    plot_F(t, ells_s, Fs_s, k0, out2)

    print(f"[OK] Guardé:")
    print(f"     {out1}")
    print(f"     {out2}")


# ============================================================
# Main
# ============================================================
def main():
    # Busca carpetas en el directorio actual cuyo nombre termine en 3 dígitos
    folders = sorted([d for d in glob.glob("*") if os.path.isdir(d) and re.search(r"\d{3}$", d)])

    if not folders:
        print("[ERROR] No encontré carpetas que terminen en 3 dígitos.")
        return

    for folder in folders:
        try:
            process_folder(folder, tmax=1500.0)
        except Exception as e:
            print(f"[ERROR] Falló {folder}: {e}")


if __name__ == "__main__":
    main()
