import numpy as np 
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import Normalize
from tools_h5 import read_coords, num_datasets, read_scalar

# Estilo general: usa fuentes LaTeX y ajusta tamaños
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern"],
    "axes.labelsize": 16,
    "axes.titlesize": 18,
    "xtick.labelsize": 14,
    "ytick.labelsize": 14,
    "legend.fontsize": 14,
    "figure.titlesize": 18,
    "text.latex.preamble": r"\usepackage{amsmath}"
})

def sph_cutXZ(r, th, ph, f, title, norm_override=None, savepath=None, show=True,
              vmin=None, vmax=None):
    # 1) Quita posible punto duplicado en 2π
    if np.isclose((ph[-1] - ph[0]) % (2*np.pi), 0.0) and len(ph) > 1:
        ph = ph[:-1]
        f  = f[:, :, :-1]

    # 2) Índices φ≈0 y φ≈π
    k0 = int(np.argmin(np.abs((ph + np.pi) % (2*np.pi) - np.pi)))
    kπ = (k0 + len(ph)//2) % len(ph)

    # 3) Malla física
    R, TH, PH = np.meshgrid(r, th, ph, indexing='ij')
    X = R * np.sin(TH) * np.cos(PH)
    Z = R * np.cos(TH)

    # 4) Normalización
    from matplotlib.colors import Normalize
    if norm_override is not None:
        norm = norm_override
    else:
        if (vmin is not None) or (vmax is not None):
            if vmin is None:
                vmin = np.nanmin(np.stack([f[:, :, k0], f[:, :, kπ]]))
            if vmax is None:
                vmax = np.nanmax(np.stack([f[:, :, k0], f[:, :, kπ]]))
            if not (vmin < vmax):
                raise ValueError(f"vmin({vmin}) debe ser < vmax({vmax}).")
            norm = Normalize(vmin=vmin, vmax=vmax)
        else:
            vmin_ = np.nanmin(np.stack([f[:, :, k0], f[:, :, kπ]]))
            vmax_ = np.nanmax(np.stack([f[:, :, k0], f[:, :, kπ]]))
            norm = Normalize(vmin=vmin_, vmax=vmax_)

    # 5) Figura y ejes con márgenes FIJOS (sin tight_layout, sin bbox_inches='tight')
    fig = plt.figure(figsize=(10, 10), constrained_layout=False)
    ax  = fig.add_subplot(111)

    pcm1 = ax.pcolormesh(X[:, :, k0], Z[:, :, k0], f[:, :, k0],
                         shading='auto', cmap='viridis', norm=norm)
    ax.pcolormesh(X[:, :, kπ], Z[:, :, kπ], f[:, :, kπ],
                  shading='auto', cmap='viridis', norm=norm)

    ax.set_xlabel(r"$x$", labelpad=10)
    ax.set_ylabel(r"$z$", labelpad=10)
    ax.set_title((title), pad=12)
    ax.set_aspect('equal', adjustable='box')
    ax.grid(False)

    # Barra de color con fracción/padding fijos (no altera el tamaño de la figura)
    cbar = fig.colorbar(pcm1, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label((title), fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    # Márgenes fijos
    fig.subplots_adjust(left=0.10, right=0.90, bottom=0.10, top=0.90)

    if savepath is not None:
        # OJO: no usar bbox_inches="tight" para no cambiar dimensiones
        plt.savefig(savepath, dpi=150)
    if show:
        plt.show()
    else:
        plt.close(fig)


import os
from matplotlib.colors import Normalize

def make_xz_video_scalar(h5file, out_mp4="xz_scalar.mp4", fps=20, dpi=150,
                         dt=None, t0=0.0, title_prefix=None, keep_frames=False,
                         vmin=None, vmax=None, symmetric=False,
                         per_snapshot=False, per_snapshot_symmetric=False):
    r, th, ph = read_coords(h5file)

    N = num_datasets(h5file)
    if N <= 0:
        raise RuntimeError("No se encontraron snapshots escalares en el archivo.")

    # ---------- Escala GLOBAL (si per_snapshot=False) ----------
    if not per_snapshot:
        use_scan = True
        if (vmin is not None) and (vmax is not None):
            if not (vmin < vmax):
                raise ValueError(f"vmin({vmin}) debe ser < vmax({vmax}).")
            norm_global = Normalize(vmin=vmin, vmax=vmax)
            use_scan = False

        if use_scan:
            vmin_glob, vmax_glob = np.inf, -np.inf
            for i in range(N):
                f = read_scalar(h5file, i)
                vmin_glob = min(vmin_glob, np.nanmin(f))
                vmax_glob = max(vmax_glob, np.nanmax(f))

            if (vmin is None) and (vmax is None) and symmetric:
                m = max(abs(vmin_glob), abs(vmax_glob))
                vmin_use, vmax_use = -m, +m
            else:
                vmin_use = vmin_glob if vmin is None else vmin
                vmax_use = vmax_glob if vmax is None else vmax
                if not (vmin_use < vmax_use):
                    raise ValueError(f"vmin({vmin_use}) debe ser < vmax({vmax_use}).")
            norm_global = Normalize(vmin=vmin_use, vmax=vmax_use)
    # -----------------------------------------------------------

    frames_dir = os.path.splitext(out_mp4)[0] + "_frames"
    os.makedirs(frames_dir, exist_ok=True)

    if title_prefix is None:
        base = os.path.basename(h5file)
        title_prefix = os.path.splitext(base)[0]

    png_paths = []
    for i in range(N):
        f = read_scalar(h5file, i)
        f = f**2

        # ---------- Escala POR SNAPSHOT ----------
        if per_snapshot:
            if per_snapshot_symmetric:
                m = max(abs(np.nanmin(f)), abs(np.nanmax(f)))
                norm_this = Normalize(vmin=-m, vmax=+m)
            else:
                # deja que sph_cutXZ calcule vmin/vmax del cuadro (no pases norm)
                norm_this = None
        else:
            norm_this = norm_global
        # ----------------------------------------

        if dt is not None:
            t = t0 + i*dt
            core = (title_prefix) if title_prefix else ""
            title = rf"{core} \quad (t={t:.6f})"
        else:
            core = (title_prefix) if title_prefix else ""
            title = rf"{core} \quad (i={i})"

        png_path = os.path.join(frames_dir, f"frame_{i:05d}.png")
        sph_cutXZ(
            r, th, ph, f, title,
            norm_override=norm_this,   # None => escala propia del cuadro
            savepath=png_path,
            show=False
        )
        png_paths.append(png_path)

    # --- MP4 (igual que antes) ---
    made_video = False
    try:
        import imageio.v2 as imageio
        with imageio.get_writer(out_mp4, fps=fps) as writer:
            for p in png_paths:
                writer.append_data(imageio.imread(p))
        made_video = True
    except Exception as e:
        print(f"[Aviso] No se pudo crear MP4 con imageio ({e}). "
              f"Frames en: {frames_dir}\n"
              f"ffmpeg -r {fps} -i {frames_dir}/frame_%05d.png -c:v libx264 -pix_fmt yuv420p {out_mp4}")

    if made_video and not keep_frames:
        for p in png_paths:
            try: os.remove(p)
            except OSError: pass
        try: os.rmdir(frames_dir)
        except OSError: pass

    if made_video:
        print(f"[OK] Video escrito en: {out_mp4}")





def sph_cutYZ(r, th, ph, f):
    k = len(ph)//4  
    R, TH, PH = np.meshgrid(r, th, ph, indexing='ij')
    X = R * np.sin(TH) * np.cos(PH)
    Y = R * np.sin(TH) * np.sin(PH)
    Z = R * np.cos(TH)



    # Figura
    plt.figure(figsize=(10, 10))
    plt.pcolormesh(Y[:, :, k], Z[:, :, k], f[:, :, k], shading='auto', cmap='viridis')
    plt.pcolormesh(Y[:, :, k + len(ph)//2], Z[:, :, k + len(ph)//2], f[:, :, k + len(ph)//2], shading='auto', cmap='viridis')

    # Etiquetas y título
    plt.xlabel(r"$y$", labelpad=10)
    plt.ylabel(r"$z$", labelpad=10)
    plt.title(r"$J^2 = \gamma_{i j} J^i J^j$", pad=12)

    # Barra de color
    cbar = plt.colorbar()
    cbar.set_label(r"$J^2$", fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    # Estética del eje
    plt.axis("equal")
    plt.grid(False)
    plt.tight_layout()
    plt.show()

def sph_cutXY(r, th, ph, f):
    R, TH, PH = np.meshgrid(r, th, ph, indexing='ij')
    X = R * np.sin(TH) * np.cos(PH)
    Y = R * np.sin(TH) * np.sin(PH)
    Z = R * np.cos(TH)

    j= len(th)//2

    # Figura
    plt.figure(figsize=(10, 10))
    plt.pcolormesh(X[:, j, :], Y[:, j, :], f[:, j, :], shading='auto', cmap='viridis')

    # Etiquetas y título
    plt.xlabel(r"$x$", labelpad=10)
    plt.ylabel(r"$y$", labelpad=10)
    plt.title(r"$J^r$", pad=12)

    # Barra de color
    cbar = plt.colorbar()
    cbar.set_label(r"$J^2$", fontsize=16)
    cbar.ax.tick_params(labelsize=14)

    # Estética del eje
    plt.axis("equal")
    plt.grid(False)
    plt.tight_layout()
    plt.show()

import numpy as np
import matplotlib.pyplot as plt

def sph_vecXZ(r, th, ph, vr, vth, vph, arrow_frac=0.1):
    """
    Dibuja el campo vectorial (J^r, J^θ, J^φ) en el plano φ = π
    proyectado a coordenadas cartesianas (x, z), con escala automática.
    """

    k = len(ph) // 2  # corte φ = π    
    dis = len(ph)//2

    # Construcción de malla
    R, TH, PH = np.meshgrid(r, th, ph, indexing='ij')
    X = R * np.sin(TH) * np.cos(PH)
    Z = R * np.cos(TH)

    Jx = vr * np.sin(TH) * np.cos(PH) + vth * np.cos(TH) * np.cos(PH) - vph * np.sin(PH)
    Jz = vr * np.cos(TH) - vth * np.sin(TH)

    # Submuestreo para visualización
    skip_r = slice(0, None, 3)
    skip_th = slice(0, None, 3)
    Xk = X[skip_r, skip_th, k]
    Zk = Z[skip_r, skip_th, k]
    Jxk = Jx[skip_r, skip_th, k]
    Jzk = Jz[skip_r, skip_th, k]
    Xk_m = X[skip_r, skip_th, k + dis]
    Zk_m = Z[skip_r, skip_th, k + dis]
    Jxk_m = Jx[skip_r, skip_th, k + dis]
    Jzk_m = Jz[skip_r, skip_th, k + dis]

    # Magnitud de los vectores
    Jmag = np.sqrt(Jxk**2 + Jzk**2)
    Jmag_m = np.sqrt(Jxk_m**2 + Jzk_m**2)
    Jmax = max(np.max(Jmag),np.max(Jmag_m))

    # Rango espacial del dominio
    dx = Xk.max() - Xk.min()
    dz = Zk.max() - Zk.min()
    domain_diag = np.sqrt(dx**2 + dz**2)

    # Estimación de escala: que la flecha más larga mida arrow_frac del dominio
    desired_length = arrow_frac * domain_diag
    scale = Jmax / desired_length

    # Gráfica
    plt.figure(figsize=(10, 8))
    plt.quiver(Xk, Zk, Jxk, Jzk,
               scale=scale, scale_units='xy', width=0.003,
               color='black')
    plt.quiver(Xk_m, Zk_m, Jxk_m, Jzk_m,
               scale=scale, scale_units='xy', width=0.003,
               color='black')
    plt.xlabel(r'$x$', fontsize=16)
    plt.ylabel(r'$z$', fontsize=16)
    plt.title(r"Campo $\vec{J}$ proyectado en el plano $\phi = \pi$", fontsize=18)
    plt.xticks(fontsize=13)
    plt.yticks(fontsize=13)
    plt.axis('equal')
    plt.tight_layout()
    plt.show()


    

