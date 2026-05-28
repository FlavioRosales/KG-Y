program SL_conservative_nonuniform
    use hdf5
    implicit none

    real(kind=8) :: r_max
    real(kind=8) :: dr_min
    real(kind=8) :: b_i, q_i
    real(kind=8) :: norm
    real(kind=8) :: p_minus, p_plus
    real(kind=8) :: g
    integer :: Nr
    integer :: ell, ell_max
    integer :: i, j, info, lwork
    integer :: n_modes, n_save, n_points
    integer :: NH, Nout, N_intervals
    character(len=256) :: h5_file
    character(len=32)  :: group_name

    real(kind=8), allocatable :: rb(:)
    real(kind=8), allocatable :: r(:), omega(:), mdiag(:)
    real(kind=8), allocatable :: dr_edge(:), drcell(:)
    real(kind=8), allocatable :: D(:), E(:), Z(:,:)
    real(kind=8), allocatable :: work(:)

    real(kind=8), allocatable :: lambda_out(:), k_out(:)
    real(kind=8), allocatable :: modes_out(:,:)
    real(kind=8), allocatable :: params_real(:)
    integer, allocatable :: params_int(:)

    integer(hid_t) :: file_id
    integer(hid_t) :: params_group_id
    integer(hid_t) :: ell_group_id
    integer :: hdferr

    external dstev

    !------------------------------------------------------------
    ! Parámetros
    !------------------------------------------------------------
    r_max  = 10000.0d0
    dr_min = 0.1d0

    ell_max = 5
    n_modes = 800

    NH   = nint(1.0d0/dr_min)
    Nout = 10000

    N_intervals = NH + Nout
    Nr = N_intervals - 1

    n_save   = min(n_modes, Nr)
    n_points = Nr + 2

    h5_file = "sl_spectrum_rmax_10000_nout_10000.h5"

    !------------------------------------------------------------
    ! Factor geométrico exterior
    !------------------------------------------------------------
    g = find_g(r_max, dr_min, Nout)

    write(*,*) 'Grid parameters:'
    write(*,*) '  r_max       = ', r_max
    write(*,*) '  dr_min      = ', dr_min
    write(*,*) '  NH          = ', NH
    write(*,*) '  Nout        = ', Nout
    write(*,*) '  Nr          = ', Nr
    write(*,*) '  n_save      = ', n_save
    write(*,*) '  g           = ', g
    write(*,*) '  HDF5 file   = ', trim(h5_file)

    !------------------------------------------------------------
    ! Arreglos
    !------------------------------------------------------------
    allocate(rb(0:Nr+1))
    allocate(r(Nr), omega(Nr), mdiag(Nr))
    allocate(dr_edge(0:Nr), drcell(Nr))
    allocate(D(Nr), E(Nr-1), Z(Nr,Nr))

    allocate(lambda_out(n_save), k_out(n_save))
    allocate(modes_out(n_points, n_save))

    lwork = max(1, 2*Nr - 2)
    allocate(work(lwork))

    call build_grid(Nr, NH, Nout, r_max, dr_min, g, rb, dr_edge)

    do i = 1, Nr
        r(i) = rb(i)
        drcell(i) = 0.5d0*(dr_edge(i-1) + dr_edge(i))
    end do

    omega = r**2 * dsqrt(1.0d0 + 1.0d0/r)

    do i = 1, Nr
        mdiag(i) = omega(i)*drcell(i)
    end do

    !------------------------------------------------------------
    ! Chequeos de la malla
    !------------------------------------------------------------
    if (abs(rb(0)) > 1.0d-14) then
        write(*,*) 'ERROR: rb(0) no es cero.'
        stop
    end if

    if (abs(rb(NH) - 1.0d0) > 1.0d-12) then
        write(*,*) 'ERROR: la malla no contiene r=1 exactamente.'
        write(*,*) 'rb(NH) = ', rb(NH)
        stop
    end if

    if (abs(rb(Nr+1) - r_max) > 1.0d-8) then
        write(*,*) 'ERROR: la malla no termina en r_max.'
        write(*,*) 'rb(Nr+1) = ', rb(Nr+1)
        stop
    end if

    write(*,*) '  first interior point = ', r(1)
    write(*,*) '  last interior point  = ', r(Nr)
    write(*,*) '  first dr             = ', dr_edge(0)
    write(*,*) '  last dr              = ', dr_edge(Nr)

    !------------------------------------------------------------
    ! Crear archivo HDF5
    !------------------------------------------------------------
    call h5open_f(hdferr)
    call h5fcreate_f(trim(h5_file), H5F_ACC_TRUNC_F, file_id, hdferr)

    ! Malla y pesos
    call write_real_1d(file_id, "r",       rb)
    call write_real_1d(file_id, "r_int",   r)
    call write_real_1d(file_id, "dr_edge", dr_edge)
    call write_real_1d(file_id, "drcell",  drcell)
    call write_real_1d(file_id, "omega",   omega)
    call write_real_1d(file_id, "mdiag",   mdiag)

    ! Parámetros
    allocate(params_real(3))
    params_real(1) = r_max
    params_real(2) = dr_min
    params_real(3) = g

    allocate(params_int(6))
    params_int(1) = Nr
    params_int(2) = NH
    params_int(3) = Nout
    params_int(4) = ell_max
    params_int(5) = n_modes
    params_int(6) = n_save

    call h5gcreate_f(file_id, "params", params_group_id, hdferr)
    call write_real_1d(params_group_id, "real", params_real)
    call write_int_1d(params_group_id, "integer", params_int)
    call h5gclose_f(params_group_id, hdferr)

    !------------------------------------------------------------
    ! Loop en ell
    !------------------------------------------------------------
    do ell = 0, ell_max

        !------------------------------------------------------------
        ! Construcción de diagonal D de B = M^{-1/2} A M^{-1/2}
        !------------------------------------------------------------
        do i = 1, Nr

            q_i = dble(ell)*(dble(ell)+1.0d0) * &
                  dsqrt(1.0d0 + 1.0d0/r(i))

            p_minus = p(0.5d0*(rb(i-1) + rb(i)))
            p_plus  = p(0.5d0*(rb(i)   + rb(i+1)))

            b_i = p_plus/dr_edge(i) + &
                  p_minus/dr_edge(i-1) + &
                  q_i*drcell(i)

            ! Frontera izquierda
            if (ell == 0 .and. i == 1) then
                b_i = b_i - p_minus/dr_edge(i-1)
            end if

            ! Frontera derecha
            if (i == Nr) then
                b_i = b_i - p_plus/dr_edge(i)
            end if

            D(i) = b_i/mdiag(i)

        end do

        !------------------------------------------------------------
        ! Sub/superdiagonal E de B
        !------------------------------------------------------------
        do i = 1, Nr-1

            p_plus = p(0.5d0*(rb(i) + rb(i+1)))

            E(i) = -p_plus/dr_edge(i)
            E(i) = E(i)/dsqrt(mdiag(i)*mdiag(i+1))

        end do

        !------------------------------------------------------------
        ! Resolver B y = lambda y
        !------------------------------------------------------------
        call dstev('V', Nr, D, E, Z, Nr, work, info)

        if (info /= 0) then
            write(*,*) 'DSTEV failed for ell = ', ell, ' info = ', info
            stop
        end if

        !------------------------------------------------------------
        ! Recuperar modos físicos R = M^{-1/2} y
        !------------------------------------------------------------
        do j = 1, n_save
            do i = 1, Nr
                Z(i,j) = Z(i,j)/dsqrt(mdiag(i))
            end do
        end do

        !------------------------------------------------------------
        ! Normalización discreta:
        !
        !   sum_i omega_i R_i^2 drcell_i = 1
        !------------------------------------------------------------
        do j = 1, n_save

            norm = 0.0d0

            do i = 1, Nr
                norm = norm + omega(i)*Z(i,j)*Z(i,j)*drcell(i)
            end do

            norm = dsqrt(norm)

            if (norm > 0.0d0) then
                do i = 1, Nr
                    Z(i,j) = Z(i,j)/norm
                end do
            end if

        end do

        !------------------------------------------------------------
        ! Preparar datasets para HDF5
        !------------------------------------------------------------
        do j = 1, n_save
            lambda_out(j) = D(j)
            k_out(j) = dsqrt(abs(D(j)))
        end do

        ! Frontera izquierda
        do j = 1, n_save
            if (ell == 0) then
                modes_out(1,j) = Z(1,j)
            else
                modes_out(1,j) = 0.0d0
            end if
        end do

        ! Puntos interiores
        do j = 1, n_save
            do i = 1, Nr
                modes_out(i+1,j) = Z(i,j)
            end do
        end do

        ! Frontera derecha: Neumann R(r_max)=R_Nr
        do j = 1, n_save
            modes_out(n_points,j) = Z(Nr,j)
        end do

        !------------------------------------------------------------
        ! Guardar grupo /ell_XXXX
        !------------------------------------------------------------
        write(group_name,'("ell_",i4.4)') ell

        call h5gcreate_f(file_id, trim(group_name), ell_group_id, hdferr)

        call write_real_1d(ell_group_id, "lambda", lambda_out)
        call write_real_1d(ell_group_id, "k",      k_out)
        call write_real_2d(ell_group_id, "modes",  modes_out)

        call h5gclose_f(ell_group_id, hdferr)

        write(*,*) 'Finished ell = ', ell

    end do

    call h5fclose_f(file_id, hdferr)
    call h5close_f(hdferr)

    write(*,*) 'Done. Wrote ', trim(h5_file)

contains

    real(kind=8) function p(r_in)
        implicit none
        real(kind=8), intent(in) :: r_in

        p = r_in**2 / dsqrt(1.0d0 + 1.0d0/r_in)

    end function p

    subroutine build_grid(Nr, NH, Nout, r_max, dr_min, g, rb, dr_edge)
        implicit none

        integer, intent(in) :: Nr, NH, Nout
        real(kind=8), intent(in) :: r_max, dr_min, g
        real(kind=8), intent(out) :: rb(0:Nr+1)
        real(kind=8), intent(out) :: dr_edge(0:Nr)

        integer :: i, j

        rb(0) = 0.0d0

        do i = 1, NH
            rb(i) = dble(i)*dr_min
        end do

        rb(NH) = 1.0d0

        do j = 0, Nout-1
            rb(NH+j+1) = rb(NH+j) + dr_min*g**j
        end do

        rb(Nr+1) = r_max

        do i = 0, Nr
            dr_edge(i) = rb(i+1) - rb(i)

            if (dr_edge(i) <= 0.0d0) then
                write(*,*) 'ERROR: dr_edge <= 0 en i = ', i
                write(*,*) 'rb(i)   = ', rb(i)
                write(*,*) 'rb(i+1) = ', rb(i+1)
                stop
            end if
        end do

    end subroutine build_grid

    real(kind=8) function find_g(r_max, dr_min, Nout)
        implicit none

        real(kind=8), intent(in) :: r_max, dr_min
        integer, intent(in) :: Nout

        real(kind=8) :: target
        real(kind=8) :: glo, ghi, gmid
        real(kind=8) :: smid
        integer :: it, itmax

        target = (r_max - 1.0d0)/dr_min

        if (abs(target - dble(Nout)) < 1.0d-12*dble(Nout)) then
            find_g = 1.0d0
            return
        end if

        if (target < dble(Nout)) then
            write(*,*) 'ERROR: con estos parámetros no existe g >= 1.'
            write(*,*) 'target = ', target
            write(*,*) 'Nout   = ', Nout
            stop
        end if

        glo = 1.0d0
        ghi = 1.01d0

        do while (geom_sum_limited(ghi, Nout, target) < target)
            ghi = 1.0d0 + 2.0d0*(ghi - 1.0d0)

            if (ghi > 2.0d0) then
                write(*,*) 'ERROR: no se pudo encerrar g.'
                write(*,*) 'ghi    = ', ghi
                write(*,*) 'target = ', target
                stop
            end if
        end do

        itmax = 200

        do it = 1, itmax

            gmid = 0.5d0*(glo + ghi)
            smid = geom_sum_limited(gmid, Nout, target)

            if (smid < target) then
                glo = gmid
            else
                ghi = gmid
            end if

        end do

        find_g = 0.5d0*(glo + ghi)

    end function find_g

    real(kind=8) function geom_sum_limited(g, n, limit)
        implicit none

        real(kind=8), intent(in) :: g
        integer, intent(in) :: n
        real(kind=8), intent(in) :: limit

        integer :: j
        real(kind=8) :: term

        geom_sum_limited = 0.0d0
        term = 1.0d0

        do j = 0, n-1

            geom_sum_limited = geom_sum_limited + term

            if (geom_sum_limited >= limit) then
                return
            end if

            term = term*g

        end do

    end function geom_sum_limited

    subroutine write_real_1d(parent_id, name, x)
        use hdf5
        implicit none

        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: x(:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t), dimension(1) :: dims
        integer :: hdferr

        dims(1) = int(size(x), kind=hsize_t)

        call h5screate_simple_f(1, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, x, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine write_real_1d

    subroutine write_int_1d(parent_id, name, x)
        use hdf5
        implicit none

        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: x(:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t), dimension(1) :: dims
        integer :: hdferr

        dims(1) = int(size(x), kind=hsize_t)

        call h5screate_simple_f(1, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine write_int_1d

    subroutine write_real_2d(parent_id, name, a)
        use hdf5
        implicit none

        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: a(:,:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t), dimension(2) :: dims
        integer :: hdferr

        dims(1) = int(size(a,1), kind=hsize_t)
        dims(2) = int(size(a,2), kind=hsize_t)

        call h5screate_simple_f(2, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, a, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)

    end subroutine write_real_2d

end program SL_conservative_nonuniform