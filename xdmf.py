#!/usr/bin/env python3
import argparse, os, h5py, numpy as np
from scipy.special import sph_harm  # cámbialo a sph_harm_y cuando actualices SciPy

# --- XDMF: rectilinear con X=r, Y=theta, Z=phi_ext (phi cerrado) ---
def build_xdmf_timeseries_rect_r_as_x(
    xdmf_path,
    h5_path,
    Nr,
    Nth,
    Nph_ext,
    times=None,
    group="/phi2_rx",
    name="phi2",
    frame_indices=None,
):
    """
    Genera un XDMF para phi^2(r,theta,phi_ext).

    frame_indices: índices k (0..Nsnap-1) para los que existen datasets phi2_rx/k.
                   Si es None, se asume 0..Nsnap-1 con Nsnap = len(times).
    """
    h5b = os.path.basename(h5_path)

    if frame_indices is None:
        if times is None:
            raise ValueError("Se requiere 'frame_indices' o 'times' para inferir nsteps.")
        frame_indices = range(len(times))

    frame_indices = list(frame_indices)

    L = [
        '<?xml version="1.0" ?>',
        '<Xdmf Version="2.0">',
        '  <Domain>',
        '    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">',
    ]

    for j, k in enumerate(frame_indices):
        if times is not None and k < len(times):
            tval = float(times[k])
        else:
            tval = float(j)

        ds = f"{h5b}:{group}/{k:05d}"

        L += [
            f'      <Grid Name="step_{k:05d}" GridType="Uniform">',
            f'        <Time Value="{tval}"/>',
            f'        <Topology TopologyType="3DRectMesh" Dimensions="{Nph_ext} {Nth} {Nr}"/>',
            '        <Geometry GeometryType="VXVYVZ">',
            f'          <DataItem Dimensions="{Nr}"      NumberType="Float" Precision="8" Format="HDF">{h5b}:/r</DataItem>',
            f'          <DataItem Dimensions="{Nth}"     NumberType="Float" Precision="8" Format="HDF">{h5b}:/theta</DataItem>',
            f'          <DataItem Dimensions="{Nph_ext}" NumberType="Float" Precision="8" Format="HDF">{h5b}:/phi_ext</DataItem>',
            '        </Geometry>',
            f'        <Attribute Name="{name}" AttributeType="Scalar" Center="Node">',
            f'          <DataItem Dimensions="{Nph_ext} {Nth} {Nr}" NumberType="Float" Precision="8" Format="HDF">{ds}</DataItem>',
            '        </Attribute>',
            '      </Grid>',
        ]

    L += [
        '    </Grid>',
        '  </Domain>',
        '</Xdmf>',
    ]

    with open(xdmf_path, "w") as f:
        f.write("\n".join(L))


def main():
    ap = argparse.ArgumentParser(
        description="Fast |phi|^2(r,theta,phi) → HDF5/XDMF (X=r,Y=theta,Z=phi) usando GEMM."
    )
    ap.add_argument("--re",      default="rePhi.h5")
    ap.add_argument("--im",      default="imPhi.h5")
    ap.add_argument("--lmax",    type=int, required=True)
    ap.add_argument("--nth",     type=int, default=96)
    ap.add_argument("--nph",     type=int, default=192)
    ap.add_argument("--outbase", default="phi2_timeseries_r_as_x_fast")
    ap.add_argument("--gzip",    type=int, default=4)
    ap.add_argument("--chunk-z", type=int, default=32,  help="chunk Nz (phi-ext)")
    ap.add_argument("--chunk-y", type=int, default=32,  help="chunk Ny (theta)")
    ap.add_argument("--chunk-x", type=int, default=128, help="chunk Nx (r)")
    ap.add_argument(
        "--nxdmf",
        type=int,
        default=0,
        help="Número de frames a procesar y listar en el XDMF, espaciados linealmente (0 = usar todos).",
    )
    args = ap.parse_args()

    re_file, im_file = args.re, args.im
    lmax, Nth, Nph = args.lmax, args.nth, args.nph
    out_h5   = f"{args.outbase}.h5"
    out_xdmf = f"{args.outbase}.xdmf"

    # ================== LEER EJES / TIEMPOS ==================
    with h5py.File(re_file, "r") as fre:
        if "r" not in fre:
            raise KeyError("'r' no encontrado en rePhi.h5")
        r = fre["r"][:].astype(np.float64)
        if "time" in fre:
            t = fre["time"][:]
        elif "t" in fre:
            t = fre["t"][:]
        else:
            t = None

    Nr = r.size

    theta  = np.linspace(0.0, np.pi,   Nth, dtype=np.float64)
    phi    = np.linspace(0.0, 2*np.pi, Nph, endpoint=False, dtype=np.float64)
    phi_ext = np.concatenate([phi, [2*np.pi]])
    Nph_ext = Nph + 1
    Nang = Nth * Nph_ext

    pairs = [(ell, m) for ell in range(lmax + 1) for m in range(-ell, ell + 1)]

    # ================== ABRIR Y PREPARAR SALIDA ==================
    with h5py.File(re_file, "r") as fre, \
         h5py.File(im_file, "r") as fim, \
         h5py.File(out_h5, "w") as fout:

        fout.create_dataset("r",       data=r,       dtype="f8")
        fout.create_dataset("theta",   data=theta,   dtype="f8")
        fout.create_dataset("phi_ext", data=phi_ext, dtype="f8")

        # ---- Determinar Nsnap con el primer modo válido ----
        Nsnap = None
        for ell, m in pairs:
            n = ell*(ell+1) + m
            dname = f"mode_{n:05d}"
            if dname not in fre:
                continue
            Re = fre[dname]
            if Re.ndim == 1 and Re.shape[0] == Nr:
                Nsnap = 1
                break
            elif Re.ndim == 2 and Re.shape[0] == Nr:
                Nsnap = Re.shape[1]
                break
            elif Re.ndim == 2 and Re.shape[1] == Nr:
                Nsnap = Re.shape[0]
                break
        if Nsnap is None:
            raise RuntimeError("No se pudo determinar Nsnap: ningún modo válido encontrado.")

        # ---- Elegir frames a procesar ----
        if args.nxdmf > 0 and args.nxdmf < Nsnap:
            frame_indices = np.linspace(0, Nsnap - 1, args.nxdmf, dtype=int)
            frame_indices = np.unique(frame_indices)
            print(f"[INFO] Solo se procesarán {len(frame_indices)} frames (linealmente muestreados).")
        else:
            frame_indices = np.arange(Nsnap, dtype=int)
            print(f"[INFO] Se procesarán todos los {Nsnap} frames.")
        frame_indices = np.asarray(frame_indices, dtype=int)
        nF = frame_indices.size

        # ---- Submuestrear tiempo en salida ----
        if t is not None and t.shape[0] == Nsnap:
            t_sel = np.asarray(t, dtype=np.float64)[frame_indices]
            fout.create_dataset("time", data=t_sel, dtype="f8")
        else:
            t_sel = None
            if t is not None and t.shape[0] != Nsnap:
                print(f"[WARN] Tamaño de 'time' ({t.shape[0]}) != Nsnap({Nsnap}), no se escribe eje tiempo.")

        g = fout.create_group("phi2_rx")
        chunks = (
            min(args.chunk_z, Nph_ext),
            min(args.chunk_y, Nth),
            min(args.chunk_x, Nr),
        )

        # ================== PRECALCULAR Y_lm(θ,φ_ext) ==================
        TH, PH = np.meshgrid(theta, phi_ext, indexing="ij")  # (Nth, Nph_ext)

        valid_pairs = []
        Y_rows = []
        coeffs_list = []

        # ---- Leer TODOS los modos válidos completos (como antes) ----
        for ell, m in pairs:
            n = ell*(ell+1) + m
            dname = f"mode_{n:05d}"

            if dname not in fre or dname not in fim:
                print(f"[WARN] {dname} no encontrado, se omite.")
                continue

            Re = fre[dname][:]
            Im = fim[dname][:]

            if Re.shape != Im.shape:
                print(f"[WARN] {dname} Re/Im shapes {Re.shape} vs {Im.shape}, se omite.")
                continue

            # Normalizar a (Nsnap, Nr)
            if Re.ndim == 1 and Re.shape[0] == Nr:
                if Nsnap != 1:
                    print(f"[WARN] {dname} tiene 1 snapshot pero Nsnap={Nsnap}, se omite.")
                    continue
                Re = Re[None, :]      # (1,Nr)
                Im = Im[None, :]
            elif Re.ndim == 2 and Re.shape[0] == Nr:
                if Re.shape[1] != Nsnap:
                    print(f"[WARN] {dname} Nsnap inconsistente {Re.shape[1]} != {Nsnap}, se omite.")
                    continue
                Re = Re.T             # (Nsnap,Nr)
                Im = Im.T
            elif Re.ndim == 2 and Re.shape[1] == Nr:
                if Re.shape[0] != Nsnap:
                    print(f"[WARN] {dname} Nsnap inconsistente {Re.shape[0]} != {Nsnap}, se omite.")
                    continue
                # ya está (Nsnap,Nr)
            else:
                print(f"[WARN] {dname} shape inesperado {Re.shape}, se omite.")
                continue

            # Registrar modo
            valid_pairs.append((ell, m))
            coeffs_list.append(Re + 1j*Im)

            # Y_lm para este modo
            Ylm = sph_harm(m, ell, PH, TH)  # (Nth,Nph_ext)
            Y_rows.append(Ylm.reshape(-1))

        if len(valid_pairs) == 0:
            raise RuntimeError("No se encontró ningún modo válido para construir phi^2.")

        coeffs = np.stack(coeffs_list, axis=0)   # (M_valid, Nsnap, Nr)
        Y_flat = np.stack(Y_rows,    axis=0)     # (M_valid, Nang)
        M_valid = coeffs.shape[0]

        print(f"[INFO] Modos válidos usados: {M_valid}")
        print(f"[INFO] Nsnap total: {Nsnap}")
        print(f"[INFO] Frames efectivamente computados: {nF}")

        # ================== LOOP SOLO SOBRE frame_indices ==================
        for j, k in enumerate(frame_indices, start=1):
            # C_k: (Nr, M_valid)
            C = coeffs[:, k, :].T   # (M_valid,Nr) -> (Nr,M_valid)

            # GEMM: (Nr,M) @ (M,Nang) -> (Nr,Nang)
            Phi_flat = C @ Y_flat

            # |phi|^2
            phi2 = np.abs(Phi_flat)**2           # (Nr,Nang)
            phi2 = phi2.reshape(Nr, Nth, Nph_ext)
            phi2_out = np.ascontiguousarray(phi2.transpose(2,1,0))  # (Nph_ext,Nth,Nr)

            g.create_dataset(
                f"{k:05d}",
                data=phi2_out,
                dtype="f8",
                compression=("gzip" if args.gzip > 0 else None),
                compression_opts=(int(args.gzip) if args.gzip > 0 else None),
                chunks=chunks,
            )

            print(f"[OK] frame {j}/{nF} (k={k}) procesado")

    # ================== XDMF ==================
    build_xdmf_timeseries_rect_r_as_x(
        out_xdmf,
        out_h5,
        Nr,
        Nth,
        Nph_ext,
        times=t if t_sel is None else t,  # usa tiempos originales; XDMF se indexa con frame_indices
        group="/phi2_rx",
        name="phi2",
        frame_indices=frame_indices,
    )

    print(f"[OK] HDF5: {out_h5}")
    print(f"[OK] XDMF: {out_xdmf}")


if __name__ == "__main__":
    main()
