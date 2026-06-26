program kg_y_main

  use kg_mpi
  use run_control
  use mesh
  use geometry
  use state
  use initial_data
  use sl_spectrum
  use cfl
  use time_integrators
  use hdf5_lib
  use plane_projection
  use sl_projection, only: sl_projector_t, sl_projector_load, &
                           sl_projector_power_phi, sl_projector_free
  use hdf5, only: HID_T

  implicit none

  integer, parameter :: RE = 1, IM = 2
  integer, parameter :: IO_BLOCK_MODES       = 64
  integer, parameter :: PROJ_BLOCK_MODES     = 16
  integer, parameter :: SPECTRUM_BLOCK_MODES = 16
  integer, parameter :: FIELD_PHI = 1, FIELD_PI = 2

  type :: sl_cache_slot_t
    integer :: ell = -999999
    real(kind=8), allocatable :: lambda(:)
    real(kind=8), allocatable :: eigenR(:,:)
  end type sl_cache_slot_t

  type :: local_mode_t
    integer :: nmode = -1
    integer :: ell   = -1
    integer :: mm    = 0
    type(state_t) :: S
  end type local_mode_t

  type(mesh_t)            :: Mesh
  type(geometry_t)        :: G
  type(plane_projector_t) :: Projector
  type(sl_projector_t)    :: SpectrumProjector

  type(local_mode_t), allocatable :: Modes(:)

  integer :: ell_min, ell_max
  integer :: n_total, nmode
  integer :: Nr, Nt, istep
  integer :: nlocal, q
  integer :: ierr_mpi, ierr_h5
  integer :: ell, mm
  integer :: k, slot, Nm_out
  integer :: cache_ptr
  integer :: n_ini, n_fin, base, rem

  integer :: nout_0D, nout_1D, nout_planes
  integer :: iout_0D, iout_1D, iout_planes
  integer :: nout_spectrum
  integer :: nphi_plane, ntheta_plane

  integer(HID_T) :: fid_out
  integer(HID_T) :: gid_grid, gid_modes, gid_diag, gid_fields
  integer(HID_T) :: gid_planes, gid_xy, gid_xz
  integer(HID_T) :: gid_spectrum, did_t_spectrum
  integer(HID_T), allocatable :: gid_spectrum_ell(:)
  integer(HID_T), allocatable :: did_spectrum_power(:)
  integer(HID_T) :: did_t0D, did_N, did_F, did_Finf
  integer(HID_T) :: did_t1D, did_phi_re, did_phi_im, did_pi_re, did_pi_im
  integer(HID_T) :: did_t_planes
  integer(HID_T) :: did_phi_xy_re, did_phi_xy_im, did_pi_xy_re, did_pi_xy_im
  integer(HID_T) :: did_phi_xz_plus_re, did_phi_xz_plus_im
  integer(HID_T) :: did_pi_xz_plus_re,  did_pi_xz_plus_im
  integer(HID_T) :: did_phi_xz_minus_re, did_phi_xz_minus_im
  integer(HID_T) :: did_pi_xz_minus_re,  did_pi_xz_minus_im

  integer, allocatable :: mode_nmode(:), mode_ell(:), mode_mm(:)

  real(kind=8) :: dt
  real(kind=8) :: t_sim0, t_sim1
  real(kind=8) :: N_loc, N_star, N_tot, Amp

  real(kind=8), allocatable :: N_local(:), F_local(:), Finf_local(:)
  real(kind=8), allocatable :: io_buffer(:,:), proj_im_buffer(:,:), plane(:,:)

  real(kind=8), allocatable :: spectrum_phi_re(:,:), spectrum_phi_im(:,:)
  real(kind=8), allocatable :: spectrum_power(:), spectrum_power_block(:)

  logical :: use_sl_ic
  logical, allocatable :: has_spectrum_ell(:)

  character(len=256) :: sl_filename
  character(len=256) :: output_filename

  real(kind=8), allocatable :: lambda(:), eigenR(:,:)

  integer, parameter :: CACHE_K = 1
  type(sl_cache_slot_t) :: sl_cache(CACHE_K)

  ! ============================================================
  ! INICIALIZACIÓN
  ! ============================================================
  call MPI_CREATE()
  call h5_init(ierr_h5)

  call sim_setup()
  call set_sigma_k(p0_p)
  call build_mesh(Mesh)

  if (rank == 0) call print_mesh_info(Mesh)

  ell_min = ell_min_p
  ell_max = ell_max_p

  call G%init(M_bh=M_bh_p, Mesh=Mesh)

  Nr = Mesh%Nr

  dt = dt_cfl(Mesh%dr_min, G%alpha, G%beta, G%guu, cfl_p)
  Nt = int(floor((t_end_p + 1.0d-12) / dt))

  if (Nt < 1) Nt = 1

  dt = t_end_p / dble(Nt)

  if (dt_out_0D_p /= -1.0d0) then
    every_0D_p = max(1, nint(dt_out_0D_p / dt))
  end if

  if (dt_out_1D_p /= -1.0d0) then
    every_1D_p = max(1, nint(dt_out_1D_p / dt))
  end if

  if (dt_out_planes_p /= -1.0d0) then
    every_planes_p = max(1, nint(dt_out_planes_p / dt))
  end if

  if (rank == 0) then
    write(*,'("=== Run parameters ===")')
    write(*,'(A,1X,ES12.5)') 'dt       =', dt
    write(*,'(A,1X,I0)')     'Nt       =', Nt
    write(*,'(A,1X,ES12.5)') 't_end    =', t_end_p
    write(*,'(A,1X,ES12.5)') 'mu       =', mu_p
    write(*,'(A,1X,A)')      'IC       =', trim(initial_conditions_p)
    write(*,'(A,1X,I0)')     'nmax     =', nmax_p
    write(*,'(A,1X,I0)')     'every_0D =', every_0D_p
    write(*,'(A,1X,I0)')     'every_1D =', every_1D_p

    if (save_planes_p) then
      write(*,'(A,1X,I0)') 'every_planes =', every_planes_p
    end if

    write(*,'("=== End run parameters ===")')
  end if

  ! ============================================================
  ! MODOS GLOBALES Y DISTRIBUCIÓN LOCAL
  ! ============================================================
  if (only_m0_p) then
    n_total = ell_max - ell_min + 1
  else
    n_total = modes_total(ell_max)
  end if

  if (rank == 0) then
    if (only_m0_p) then
      write(*,'("=== Solo m=0, ell=",I0,"...",I0," | n_total=",I0," ===")') &
        ell_min, ell_max, n_total
    else
      write(*,'("=== Todos los modos hasta ell_max=",I0," | n_total=",I0," ===")') &
        ell_max, n_total
    end if
  end if

  base = n_total / nproc
  rem  = mod(n_total, nproc)

  if (rank < rem) then
    n_ini = rank * (base + 1)
    n_fin = n_ini + base
  else
    n_ini = rem * (base + 1) + (rank - rem) * base
    n_fin = n_ini + base - 1
  end if

  nlocal = n_fin - n_ini + 1

  allocate(Modes(nlocal))
  allocate(mode_nmode(nlocal), mode_ell(nlocal), mode_mm(nlocal))

  q = 0

  do nmode = n_ini, n_fin
    q = q + 1

    if (only_m0_p) then
      ell = ell_min + nmode
      mm  = 0
    else
      call lm_of_n(nmode, ell, mm)
    end if

    Modes(q)%nmode = nmode
    Modes(q)%ell   = ell
    Modes(q)%mm    = mm

    mode_nmode(q) = nmode
    mode_ell(q)   = ell
    mode_mm(q)    = mm
  end do

  write(*,'("Rank ",I0," | nmode range = [",I0,",",I0,"] | nlocal = ",I0)') &
    rank, n_ini, n_fin, nlocal

  flush(6)

  ! ============================================================
  ! CACHE DE LA BASE STURM-LIOUVILLE
  ! ============================================================
  cache_ptr = 0

  do k = 1, CACHE_K
    sl_cache(k)%ell = -999999
  end do

  sl_filename = '/home/flavio/Codes/KG-Y/numerical/sl_spectrum.h5'

  use_sl_ic = (trim(initial_conditions_p) == 'mode_from_SL'     .or. &
               trim(initial_conditions_p) == 'mode_from_sl'     .or. &
               trim(initial_conditions_p) == 'mode_from_SL_all' .or. &
               trim(initial_conditions_p) == 'mode_from_sl_all')

  ! ============================================================
  ! CONDICIONES INICIALES
  ! ============================================================
  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

  if (rank == 0) then
    print*, "======================================================"
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

  N_loc  = 0.0d0
  N_star = 100.0d0

  do q = 1, nlocal

    call Modes(q)%S%init(Nr=Nr, Nt=0, mu_field=mu_p)

    slot = -1

    if (use_sl_ic) then
      call get_sl_cache_slot(Modes(q)%ell, slot, Nm_out)
    end if

    select case (trim(initial_conditions_p))

    case ('gaussian_positive_charge')
      call set_IC_incoming_gaussian_positive_charge(k0_p, Modes(q)%S, G, Mesh)

    case ('mode_from_SL', 'mode_from_sl')
      call set_IC_mode_from_SL(sl_cache(slot)%lambda, &
                               sl_cache(slot)%eigenR, &
                               Modes(q)%S, G, Mesh)

    case ('mode_from_SL_all', 'mode_from_sl_all')
      call mode_from_SL_all(mu_p, sl_cache(slot)%lambda, &
                            sl_cache(slot)%eigenR, &
                            Modes(q)%S, Mesh)

    end select

    N_loc = N_loc + Modes(q)%S%NC(G, Mesh)

  end do

  call MPI_Allreduce(N_loc, N_tot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, &
                     MPI_COMM_WORLD, ierr_mpi)

  Amp = dsqrt(N_star / N_tot)

  do concurrent (q = 1:nlocal)
    Modes(q)%S%phi = Amp * Modes(q)%S%phi
    Modes(q)%S%psi = Amp * Modes(q)%S%psi
    Modes(q)%S%pi  = Amp * Modes(q)%S%pi
  end do

  print*, "Initial data was setting for rank", rank

  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

  if (rank == 0) then
    print*, "======================================================"
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

  ! ============================================================
  ! ARCHIVO LOCAL DEL RANK Y DATASETS
  ! ============================================================
  nout_0D = count_output_samples(Nt, every_0D_p)
  nout_1D = count_output_samples(Nt, every_1D_p)

  allocate(N_local(nlocal), F_local(nlocal), Finf_local(nlocal))
  allocate(io_buffer(Nr, min(IO_BLOCK_MODES, nlocal)))

  write(output_filename,'("output_rank_",I5.5,".h5")') rank

  call h5_open_new(trim(output_filename), fid_out, ierr_h5)

  call h5_create_group(fid_out, 'grid',        gid_grid,   ierr_h5)
  call h5_create_group(fid_out, 'modes',       gid_modes,  ierr_h5)
  call h5_create_group(fid_out, 'diagnostics', gid_diag,   ierr_h5)
  call h5_create_group(fid_out, 'fields',      gid_fields, ierr_h5)

  ! ============================================================
  ! ESPECTRO SL
  !
  ! Se guarda en los mismos tiempos que los campos 1D:
  !
  !   /spectrum/t
  !   /spectrum/ell_xxxx/k
  !   /spectrum/ell_xxxx/lambda
  !   /spectrum/ell_xxxx/P_phi_local(n, it)
  !
  ! Cada rank escribe únicamente los ell que contiene.
  ! ============================================================
  if (use_sl_ic) then

    nout_spectrum = nout_1D

    allocate(gid_spectrum_ell(0:ell_max))
    allocate(did_spectrum_power(0:ell_max))
    allocate(has_spectrum_ell(0:ell_max))

    gid_spectrum_ell   = -1_HID_T
    did_spectrum_power = -1_HID_T
    has_spectrum_ell   = .false.

    allocate(spectrum_phi_re(Nr, min(SPECTRUM_BLOCK_MODES, nlocal)))
    allocate(spectrum_phi_im(Nr, min(SPECTRUM_BLOCK_MODES, nlocal)))

    call h5_create_group(fid_out, 'spectrum', gid_spectrum, ierr_h5)

    call h5_create_real_time_series(gid_spectrum, 't', &
                                    nout_spectrum, did_t_spectrum, ierr_h5)

    call spectrum_create_layout()

  end if

  call h5_write_real_1d(gid_grid,  'r',     Mesh%r,     ierr_h5)
  call h5_write_int_1d (gid_modes, 'nmode', mode_nmode, ierr_h5)
  call h5_write_int_1d (gid_modes, 'ell',   mode_ell,   ierr_h5)
  call h5_write_int_1d (gid_modes, 'mm',    mode_mm,    ierr_h5)

  call h5_create_real_time_series(gid_diag, 't', nout_0D, did_t0D, ierr_h5)
  call h5_create_real_0d_series(gid_diag, 'N',    nlocal, nout_0D, did_N,    ierr_h5)
  call h5_create_real_0d_series(gid_diag, 'F',    nlocal, nout_0D, did_F,    ierr_h5)
  call h5_create_real_0d_series(gid_diag, 'Finf', nlocal, nout_0D, did_Finf, ierr_h5)

  call h5_create_real_time_series(gid_fields, 't', nout_1D, did_t1D, ierr_h5)

  call h5_create_real_1d_series(gid_fields, 'phi_re', Nr, nlocal, nout_1D, did_phi_re, ierr_h5)
  call h5_create_real_1d_series(gid_fields, 'phi_im', Nr, nlocal, nout_1D, did_phi_im, ierr_h5)
  call h5_create_real_1d_series(gid_fields, 'pi_re',  Nr, nlocal, nout_1D, did_pi_re,  ierr_h5)
  call h5_create_real_1d_series(gid_fields, 'pi_im',  Nr, nlocal, nout_1D, did_pi_im,  ierr_h5)

  if (save_planes_p) then

    nphi_plane   = nphi_plane_p
    ntheta_plane = ntheta_plane_p

    if (nphi_plane < 1) then
      nphi_plane = 2 * ell_max + 1
    end if

    if (ntheta_plane < 1) then
      ntheta_plane = ell_max + 1
    end if

    call plane_projector_init(Projector, Mesh, mode_ell, mode_mm, r_plane_max_p, &
                              nphi_plane, ntheta_plane)

    nout_planes = count_output_samples(Nt, every_planes_p)

    allocate(proj_im_buffer(Projector%nr, min(PROJ_BLOCK_MODES, nlocal)))
    allocate(plane(Projector%nr, max(Projector%nphi, Projector%ntheta)))

    call h5_create_group(fid_out,    'planes', gid_planes, ierr_h5)
    call h5_create_group(gid_planes, 'xy',     gid_xy,     ierr_h5)
    call h5_create_group(gid_planes, 'xz',     gid_xz,     ierr_h5)

    call h5_write_real_1d(gid_xy, 'r',      Projector%r,      ierr_h5)
    call h5_write_real_1d(gid_xy, 'varphi', Projector%varphi, ierr_h5)
    call h5_write_real_1d(gid_xz, 'r',      Projector%r,      ierr_h5)
    call h5_write_real_1d(gid_xz, 'theta',  Projector%theta,  ierr_h5)

    call h5_create_real_time_series(gid_planes, 't', nout_planes, did_t_planes, ierr_h5)

    call h5_create_real_2d_series(gid_xy, 'phi_re', Projector%nr, Projector%nphi, &
                                  nout_planes, did_phi_xy_re, ierr_h5)

    call h5_create_real_2d_series(gid_xy, 'phi_im', Projector%nr, Projector%nphi, &
                                  nout_planes, did_phi_xy_im, ierr_h5)

    call h5_create_real_2d_series(gid_xy, 'pi_re', Projector%nr, Projector%nphi, &
                                  nout_planes, did_pi_xy_re, ierr_h5)

    call h5_create_real_2d_series(gid_xy, 'pi_im', Projector%nr, Projector%nphi, &
                                  nout_planes, did_pi_xy_im, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'phi_re_plus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_phi_xz_plus_re, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'phi_im_plus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_phi_xz_plus_im, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'pi_re_plus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_pi_xz_plus_re, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'pi_im_plus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_pi_xz_plus_im, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'phi_re_minus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_phi_xz_minus_re, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'phi_im_minus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_phi_xz_minus_im, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'pi_re_minus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_pi_xz_minus_re, ierr_h5)

    call h5_create_real_2d_series(gid_xz, 'pi_im_minus', Projector%nr, Projector%ntheta, &
                                  nout_planes, did_pi_xz_minus_im, ierr_h5)

  end if

  ! ============================================================
  ! SNAPSHOT INICIAL
  ! ============================================================
  iout_0D = 1

  do q = 1, nlocal
    N_local(q)    = Modes(q)%S%NC(G, Mesh)
    F_local(q)    = Modes(q)%S%Noether_rate_hor(G, Mesh)
    Finf_local(q) = Modes(q)%S%Noether_rate_inf(G, Mesh)
  end do

  call h5_write_real_time_sample(did_t0D, iout_0D, 0.0d0, ierr_h5)
  call h5_write_real_0d_sample(did_N,    iout_0D, N_local,    ierr_h5)
  call h5_write_real_0d_sample(did_F,    iout_0D, F_local,    ierr_h5)
  call h5_write_real_0d_sample(did_Finf, iout_0D, Finf_local, ierr_h5)

  iout_1D = 1

  call h5_write_real_time_sample(did_t1D, iout_1D, 0.0d0, ierr_h5)
  call write_fields_snapshot(iout_1D)

  if (save_planes_p) then
    iout_planes = 1
    call write_planes_snapshot(iout_planes, 0.0d0)
  end if

  if (use_sl_ic) then
    call write_spectrum_snapshot(iout_1D, 0.0d0)
  end if

  ! ============================================================
  ! LOOP TEMPORAL GLOBAL
  ! ============================================================
  t_sim0 = MPI_Wtime()

  do istep = 1, Nt

    do q = 1, nlocal
      call ssprk3_step_state(dt, Modes(q)%S, G, Mesh, Modes(q)%ell, mu_p)
    end do

    if (mod(istep, every_0D_p) == 0 .or. istep == Nt) then

      iout_0D = iout_0D + 1

      do q = 1, nlocal
        N_local(q)    = Modes(q)%S%NC(G, Mesh)
        F_local(q)    = Modes(q)%S%Noether_rate_hor(G, Mesh)
        Finf_local(q) = Modes(q)%S%Noether_rate_inf(G, Mesh)
      end do

      call h5_write_real_time_sample(did_t0D, iout_0D, dt*dble(istep), ierr_h5)
      call h5_write_real_0d_sample(did_N,    iout_0D, N_local,    ierr_h5)
      call h5_write_real_0d_sample(did_F,    iout_0D, F_local,    ierr_h5)
      call h5_write_real_0d_sample(did_Finf, iout_0D, Finf_local, ierr_h5)

    end if

    if (mod(istep, every_1D_p) == 0 .or. istep == Nt) then

      iout_1D = iout_1D + 1

      call h5_write_real_time_sample(did_t1D, iout_1D, dt*dble(istep), ierr_h5)
      call write_fields_snapshot(iout_1D)

      if (use_sl_ic) then
        call write_spectrum_snapshot(iout_1D, dt*dble(istep))
      end if

    end if

    if (save_planes_p) then
      if (mod(istep, every_planes_p) == 0 .or. istep == Nt) then

        iout_planes = iout_planes + 1

        call write_planes_snapshot(iout_planes, dt*dble(istep))

        if (rank == 0) then
          write(*,'("t = ",ES12.5)') dt*dble(istep)
        end if

      end if
    end if

  end do

  call MPI_Barrier(MPI_COMM_WORLD, ierr_mpi)

  t_sim1 = MPI_Wtime()

  if (rank == 0) then
    write(*,'("=== Total wall time: ",F12.4," s ===")') t_sim1 - t_sim0
  end if

  ! ============================================================
  ! CIERRE DE HDF5
  ! ============================================================
  call h5_close_dataset(did_t0D,    ierr_h5)
  call h5_close_dataset(did_N,      ierr_h5)
  call h5_close_dataset(did_F,      ierr_h5)
  call h5_close_dataset(did_Finf,   ierr_h5)
  call h5_close_dataset(did_t1D,    ierr_h5)
  call h5_close_dataset(did_phi_re, ierr_h5)
  call h5_close_dataset(did_phi_im, ierr_h5)
  call h5_close_dataset(did_pi_re,  ierr_h5)
  call h5_close_dataset(did_pi_im,  ierr_h5)

  if (save_planes_p) then

    call h5_close_dataset(did_t_planes,         ierr_h5)
    call h5_close_dataset(did_phi_xy_re,        ierr_h5)
    call h5_close_dataset(did_phi_xy_im,        ierr_h5)
    call h5_close_dataset(did_pi_xy_re,         ierr_h5)
    call h5_close_dataset(did_pi_xy_im,         ierr_h5)
    call h5_close_dataset(did_phi_xz_plus_re,   ierr_h5)
    call h5_close_dataset(did_phi_xz_plus_im,   ierr_h5)
    call h5_close_dataset(did_pi_xz_plus_re,    ierr_h5)
    call h5_close_dataset(did_pi_xz_plus_im,    ierr_h5)
    call h5_close_dataset(did_phi_xz_minus_re,  ierr_h5)
    call h5_close_dataset(did_phi_xz_minus_im,  ierr_h5)
    call h5_close_dataset(did_pi_xz_minus_re,   ierr_h5)
    call h5_close_dataset(did_pi_xz_minus_im,   ierr_h5)

    call h5_close_group(gid_xy,     ierr_h5)
    call h5_close_group(gid_xz,     ierr_h5)
    call h5_close_group(gid_planes, ierr_h5)

  end if

  if (use_sl_ic) then

    do ell = 0, ell_max
      if (.not. has_spectrum_ell(ell)) cycle

      call h5_close_dataset(did_spectrum_power(ell), ierr_h5)
      call h5_close_group(gid_spectrum_ell(ell), ierr_h5)
    end do

    call h5_close_dataset(did_t_spectrum, ierr_h5)
    call h5_close_group(gid_spectrum, ierr_h5)

  end if

  call h5_close_group(gid_grid,   ierr_h5)
  call h5_close_group(gid_modes,  ierr_h5)
  call h5_close_group(gid_diag,   ierr_h5)
  call h5_close_group(gid_fields, ierr_h5)

  call h5_close_file(fid_out, ierr_h5)

  ! ============================================================
  ! LIMPIEZA
  ! ============================================================
  do q = 1, nlocal
    call Modes(q)%S%free()
  end do

  if (allocated(io_buffer))            deallocate(io_buffer)
  if (allocated(proj_im_buffer))       deallocate(proj_im_buffer)
  if (allocated(spectrum_phi_re))      deallocate(spectrum_phi_re)
  if (allocated(spectrum_phi_im))      deallocate(spectrum_phi_im)
  if (allocated(spectrum_power))       deallocate(spectrum_power)
  if (allocated(spectrum_power_block)) deallocate(spectrum_power_block)
  if (allocated(gid_spectrum_ell))     deallocate(gid_spectrum_ell)
  if (allocated(did_spectrum_power))   deallocate(did_spectrum_power)
  if (allocated(has_spectrum_ell))     deallocate(has_spectrum_ell)
  if (allocated(plane))                deallocate(plane)
  if (allocated(N_local))              deallocate(N_local)
  if (allocated(F_local))              deallocate(F_local)
  if (allocated(Finf_local))           deallocate(Finf_local)
  if (allocated(mode_nmode))           deallocate(mode_nmode)
  if (allocated(mode_ell))             deallocate(mode_ell)
  if (allocated(mode_mm))              deallocate(mode_mm)
  if (allocated(Modes))                deallocate(Modes)

  call plane_projector_free(Projector)
  call sl_projector_free(SpectrumProjector)

  do k = 1, CACHE_K
    if (allocated(sl_cache(k)%lambda)) deallocate(sl_cache(k)%lambda)
    if (allocated(sl_cache(k)%eigenR)) deallocate(sl_cache(k)%eigenR)
  end do

  call h5_finalize(ierr_h5)
  call MPI_Finalize(ierr_mpi)

contains

  integer function count_output_samples(nsteps, every) result(nout)
    integer, intent(in) :: nsteps, every

    nout = 1 + nsteps / every

    if (mod(nsteps, every) /= 0) then
      nout = nout + 1
    end if
  end function count_output_samples


  subroutine write_fields_snapshot(isnap)
    integer, intent(in) :: isnap

    integer :: ibeg, iend, nb, j, qloc

    do ibeg = 1, nlocal, IO_BLOCK_MODES

      iend = min(ibeg + IO_BLOCK_MODES - 1, nlocal)
      nb   = iend - ibeg + 1

      do j = 1, nb
        qloc = ibeg + j - 1
        io_buffer(:,j) = Modes(qloc)%S%phi(:,RE)
      end do

      call h5_write_real_1d_block(did_phi_re, isnap, ibeg, &
                                  io_buffer(:,1:nb), ierr_h5)

      do j = 1, nb
        qloc = ibeg + j - 1
        io_buffer(:,j) = Modes(qloc)%S%phi(:,IM)
      end do

      call h5_write_real_1d_block(did_phi_im, isnap, ibeg, &
                                  io_buffer(:,1:nb), ierr_h5)

      do j = 1, nb
        qloc = ibeg + j - 1
        io_buffer(:,j) = Modes(qloc)%S%pi(:,RE)
      end do

      call h5_write_real_1d_block(did_pi_re, isnap, ibeg, &
                                  io_buffer(:,1:nb), ierr_h5)

      do j = 1, nb
        qloc = ibeg + j - 1
        io_buffer(:,j) = Modes(qloc)%S%pi(:,IM)
      end do

      call h5_write_real_1d_block(did_pi_im, isnap, ibeg, &
                                  io_buffer(:,1:nb), ierr_h5)

    end do
  end subroutine write_fields_snapshot


  subroutine write_planes_snapshot(isnap, time)
    integer,      intent(in) :: isnap
    real(kind=8), intent(in) :: time

    call h5_write_real_time_sample(did_t_planes, isnap, time, ierr_h5)

    call write_plane_component(did_phi_xy_re,       isnap, PLANE_XY,  1, FIELD_PHI, PART_RE)
    call write_plane_component(did_phi_xy_im,       isnap, PLANE_XY,  1, FIELD_PHI, PART_IM)
    call write_plane_component(did_pi_xy_re,        isnap, PLANE_XY,  1, FIELD_PI,  PART_RE)
    call write_plane_component(did_pi_xy_im,        isnap, PLANE_XY,  1, FIELD_PI,  PART_IM)

    call write_plane_component(did_phi_xz_plus_re,  isnap, PLANE_XZ,  1, FIELD_PHI, PART_RE)
    call write_plane_component(did_phi_xz_plus_im,  isnap, PLANE_XZ,  1, FIELD_PHI, PART_IM)
    call write_plane_component(did_pi_xz_plus_re,   isnap, PLANE_XZ,  1, FIELD_PI,  PART_RE)
    call write_plane_component(did_pi_xz_plus_im,   isnap, PLANE_XZ,  1, FIELD_PI,  PART_IM)

    call write_plane_component(did_phi_xz_minus_re, isnap, PLANE_XZ, -1, FIELD_PHI, PART_RE)
    call write_plane_component(did_phi_xz_minus_im, isnap, PLANE_XZ, -1, FIELD_PHI, PART_IM)
    call write_plane_component(did_pi_xz_minus_re,  isnap, PLANE_XZ, -1, FIELD_PI,  PART_RE)
    call write_plane_component(did_pi_xz_minus_im,  isnap, PLANE_XZ, -1, FIELD_PI,  PART_IM)

  end subroutine write_planes_snapshot


  subroutine write_plane_component(did, isnap, which_plane, side, which_field, part)
    integer(HID_T), intent(in) :: did
    integer,        intent(in) :: isnap, which_plane, side, which_field, part

    integer :: ibeg, iend, nb, j, qloc
    real(kind=8) :: azimuth_sign

    plane = 0.0d0

    do ibeg = 1, nlocal, PROJ_BLOCK_MODES

      iend = min(ibeg + PROJ_BLOCK_MODES - 1, nlocal)
      nb   = iend - ibeg + 1

      do j = 1, nb

        qloc = ibeg + j - 1

        azimuth_sign = 1.0d0

        if (which_plane == PLANE_XZ .and. side < 0) then
          if (mod(abs(Modes(qloc)%mm), 2) == 1) azimuth_sign = -1.0d0
        end if

        if (which_field == FIELD_PHI) then
          io_buffer(1:Projector%nr,j) = azimuth_sign * &
                                         Modes(qloc)%S%phi(1:Projector%nr,RE)

          proj_im_buffer(:,j) = azimuth_sign * &
                                 Modes(qloc)%S%phi(1:Projector%nr,IM)
        else
          io_buffer(1:Projector%nr,j) = azimuth_sign * &
                                         Modes(qloc)%S%pi(1:Projector%nr,RE)

          proj_im_buffer(:,j) = azimuth_sign * &
                                 Modes(qloc)%S%pi(1:Projector%nr,IM)
        end if

      end do

      if (which_plane == PLANE_XY) then
        call plane_projector_add_xy(Projector, ibeg, &
                                    io_buffer(1:Projector%nr,1:nb), &
                                    proj_im_buffer(:,1:nb), part, &
                                    plane(:,1:Projector%nphi))
      else
        call plane_projector_add_xz(Projector, ibeg, &
                                    io_buffer(1:Projector%nr,1:nb), &
                                    proj_im_buffer(:,1:nb), part, &
                                    plane(:,1:Projector%ntheta))
      end if

    end do

    if (which_plane == PLANE_XY) then
      call h5_write_real_2d_sample(did, isnap, &
                                   plane(:,1:Projector%nphi), ierr_h5)
    else
      call h5_write_real_2d_sample(did, isnap, &
                                   plane(:,1:Projector%ntheta), ierr_h5)
    end if

  end subroutine write_plane_component


  subroutine write_spectrum_snapshot(isnap, time)
    !!
    !! Escribe el espectro en una salida 1D:
    !!
    !!   P_{ell n}^{local}(t)
    !!   =
    !!   sum_{m en rank} |a_{ell m n}|^2 .
    !!
    integer,      intent(in) :: isnap
    real(kind=8), intent(in) :: time

    integer :: q_first, q_last
    integer :: ell_here
    integer :: ibeg, iend, nb
    integer :: j, qloc

    call h5_write_real_time_sample(did_t_spectrum, isnap, time, ierr_h5)

    q_first = 1

    do while (q_first <= nlocal)

      ell_here = Modes(q_first)%ell
      q_last   = q_first

      do while (q_last < nlocal)
        if (Modes(q_last+1)%ell /= ell_here) exit
        q_last = q_last + 1
      end do

      call sl_projector_load(SpectrumProjector, trim(sl_filename), &
                             ell_here, Mesh)

      allocate(spectrum_power(SpectrumProjector%nspec))
      allocate(spectrum_power_block(SpectrumProjector%nspec))

      spectrum_power = 0.0d0

      do ibeg = q_first, q_last, SPECTRUM_BLOCK_MODES

        iend = min(ibeg + SPECTRUM_BLOCK_MODES - 1, q_last)
        nb   = iend - ibeg + 1

        do j = 1, nb
          qloc = ibeg + j - 1

          spectrum_phi_re(:,j) = Modes(qloc)%S%phi(:,RE)
          spectrum_phi_im(:,j) = Modes(qloc)%S%phi(:,IM)
        end do

        call sl_projector_power_phi(SpectrumProjector, &
                                    spectrum_phi_re(:,1:nb), &
                                    spectrum_phi_im(:,1:nb), &
                                    Mesh, G, spectrum_power_block)

        spectrum_power = spectrum_power + spectrum_power_block

      end do

      call h5_write_real_0d_sample( &
        did_spectrum_power(ell_here), isnap, spectrum_power, ierr_h5)

      deallocate(spectrum_power)
      deallocate(spectrum_power_block)

      call sl_projector_free(SpectrumProjector)

      q_first = q_last + 1

    end do

  end subroutine write_spectrum_snapshot


    subroutine spectrum_create_layout()

    integer :: q_first, q_last
    integer :: ell_here
    character(len=32) :: ell_name

    q_first = 1

    do while (q_first <= nlocal)

      ell_here = Modes(q_first)%ell
      q_last   = q_first

      do while (q_last < nlocal)
        if (Modes(q_last + 1)%ell /= ell_here) exit
        q_last = q_last + 1
      end do

      write(ell_name, '("ell_", I4.4)') ell_here

      call h5_create_group(gid_spectrum, trim(ell_name), &
                           gid_spectrum_ell(ell_here), ierr_h5)

      call sl_projector_load(SpectrumProjector, trim(sl_filename), &
                             ell_here, Mesh)

      call h5_write_real_1d(gid_spectrum_ell(ell_here), 'k', &
                            SpectrumProjector%k, ierr_h5)

      call h5_write_real_1d(gid_spectrum_ell(ell_here), 'lambda', &
                            SpectrumProjector%lambda, ierr_h5)

      call h5_create_real_0d_series( &
        gid_spectrum_ell(ell_here), 'P_phi_local', &
        SpectrumProjector%nspec, nout_spectrum, &
        did_spectrum_power(ell_here), ierr_h5)

      has_spectrum_ell(ell_here) = .true.

      call sl_projector_free(SpectrumProjector)

      q_first = q_last + 1

    end do

  end subroutine spectrum_create_layout


  subroutine get_sl_cache_slot(ell_in, slot_out, nmodes_out)
    integer, intent(in)  :: ell_in
    integer, intent(out) :: slot_out
    integer, intent(out) :: nmodes_out

    integer :: kk

    do kk = 1, CACHE_K
      if (sl_cache(kk)%ell == ell_in) then
        slot_out   = kk
        nmodes_out = size(sl_cache(kk)%lambda)
        return
      end if
    end do

    if (allocated(lambda)) deallocate(lambda)
    if (allocated(eigenR)) deallocate(eigenR)

    call read_sl_modes_h5(trim(sl_filename), ell_in, lambda, eigenR, Mesh)

    cache_ptr = cache_ptr + 1

    if (cache_ptr > CACHE_K) then
      cache_ptr = 1
    end if

    slot_out = cache_ptr

    if (allocated(sl_cache(slot_out)%lambda)) then
      deallocate(sl_cache(slot_out)%lambda)
    end if

    if (allocated(sl_cache(slot_out)%eigenR)) then
      deallocate(sl_cache(slot_out)%eigenR)
    end if

    sl_cache(slot_out)%ell = ell_in

    call move_alloc(lambda, sl_cache(slot_out)%lambda)
    call move_alloc(eigenR, sl_cache(slot_out)%eigenR)

    nmodes_out = size(sl_cache(slot_out)%lambda)

  end subroutine get_sl_cache_slot

end program kg_y_main