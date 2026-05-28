program SL_solve_ell_parallel
    !====================================================================
    ! Per-ell Sturm-Liouville solver for the KG radial basis.
    !
    ! Purpose:
    !   Run independent ell jobs in parallel, each writing its own HDF5 file:
    !       <output_prefix>_ell_####.h5
    !   Then merge the per-ell files with merge_sl_ell_files.py.
    !
    ! Usage:
    !   ./sl_solve_ell_parallel ell
    !   ./sl_solve_ell_parallel ell_start ell_end
    !   ./sl_solve_ell_parallel ell_start ell_end output_prefix
    !
    ! Example:
    !   ./sl_solve_ell_parallel 0 0 sl_part
    !   ./sl_solve_ell_parallel 12 12 sl_part
    !
    ! This version preserves the FV conservative tridiagonal operator,
    ! uses the discrete Nyquist-safe cutoff, and saves a dynamic number
    ! of modes per ell.
    !====================================================================
    use hdf5
    implicit none

    ! ---------- parameters ---------------------------------------------
    real(kind=8), parameter :: r_target = 10000.0d0
    real(kind=8), parameter :: dr_min   =     0.1d0
    real(kind=8), parameter :: dr_max   =     1.0d0
    real(kind=8), parameter :: r_flat   =  1000.0d0

    integer,      parameter :: n_modes_cap = 5200

    ! Safety fraction in the continuous phase Nyquist q_N = pi/dr_max.
    ! The saved cutoff is converted to the discrete second-order k_h.
    real(kind=8), parameter :: nyquist_safety = 0.80d0
    real(kind=8), parameter :: pi = acos(-1.0d0)
    real(kind=8), parameter :: q_safe = nyquist_safety*pi/dr_max
    real(kind=8), parameter :: k_safe = (2.0d0/dr_max) * dsin(0.5d0*q_safe*dr_max)

    real(kind=8), parameter :: lambda_tol = 1.0d-12
    ! -------------------------------------------------------------------

    integer :: Nr, N_intervals
    integer :: ell_start, ell_end, ell
    integer :: i, j, lwork, info, n_save_ell
    character(len=256) :: output_prefix

    real(kind=8) :: r_outer, drc_i, norm, kj

    real(kind=8), allocatable :: r(:), dr_edge(:), mdiag(:)
    real(kind=8), allocatable :: D(:), E(:), Z(:,:), work(:), k_save(:)

    integer :: hdferr

    external dstev

    call parse_args(ell_start, ell_end, output_prefix)

    call count_grid(r_target, N_intervals, r_outer)
    Nr = N_intervals - 1

    call print_grid_info(Nr, N_intervals, r_target, r_outer, ell_start, ell_end, output_prefix)

    allocate(r(0:Nr+1), dr_edge(0:Nr), mdiag(Nr))
    call build_grid(Nr, r, dr_edge)
    r_outer = r(Nr+1)

    do i = 1, Nr
        drc_i = 0.5d0*(dr_edge(i-1) + dr_edge(i))
        mdiag(i) = w(r(i)) * drc_i
        if (mdiag(i) <= 0.0d0) then
            write(*,*) 'ERROR: nonpositive mass matrix entry at i=', i
            write(*,*) 'r(i)=', r(i), 'mdiag=', mdiag(i)
            stop
        end if
    end do

    call print_mesh_quality(Nr, r, dr_edge)

    allocate(D(Nr), E(Nr-1), Z(Nr,Nr))
    lwork = max(1, 2*Nr - 2)
    allocate(work(lwork), k_save(Nr))

    call h5open_f(hdferr)

    do ell = ell_start, ell_end

        call build_B_matrix(Nr, ell, r, dr_edge, mdiag, D, E)

        call dstev('V', Nr, D, E, Z, Nr, work, info)
        if (info /= 0) then
            write(*,'(A,I6,A,I8)') 'ERROR: dstev failed for ell=', ell, ', info=', info
            stop
        end if

        n_save_ell = 0
        do j = 1, Nr
            if (D(j) < -lambda_tol) then
                write(*,'(A,I6,A,I8,A,ES24.16)') &
                    'ERROR: negative eigenvalue for ell=', ell, ', j=', j, ', lambda=', D(j)
                stop
            else if (D(j) < 0.0d0) then
                kj = 0.0d0
            else
                kj = dsqrt(D(j))
            end if

            if (kj <= k_safe .and. n_save_ell < n_modes_cap) then
                n_save_ell = n_save_ell + 1
                k_save(n_save_ell) = kj
            end if
        end do

        if (n_save_ell <= 0) then
            write(*,'(A,I6,A,ES12.4)') &
                'ERROR: no modes below k_safe for ell=', ell, ', k_safe=', k_safe
            stop
        end if

        do j = 1, n_save_ell
            do i = 1, Nr
                Z(i,j) = Z(i,j) / dsqrt(mdiag(i))
            end do
        end do

        do j = 1, n_save_ell
            norm = 0.0d0
            do i = 1, Nr
                norm = norm + Z(i,j)*Z(i,j)*mdiag(i)
            end do
            norm = dsqrt(norm)
            if (norm <= 0.0d0) then
                write(*,*) 'ERROR: zero norm for ell,j=', ell, j
                stop
            end if
            do i = 1, Nr
                Z(i,j) = Z(i,j)/norm
            end do
        end do

        call write_single_ell_file(trim(output_prefix), ell, Nr, N_intervals, r_outer, &
                                   r, n_save_ell, k_save(1:n_save_ell), Z(:,1:n_save_ell))

        write(*,'(A,I6,A,I8,A,ES12.4,A,ES12.4,A,ES12.4)') &
            '  ell=', ell, '  n_save=', n_save_ell, &
            '  k_min=', k_save(1), '  k_max=', k_save(n_save_ell), &
            '  k_safe=', k_safe

    end do

    call h5close_f(hdferr)

    write(*,'(A)') 'Done.'

contains

    !=================================================================
    ! Command-line parsing
    !=================================================================
    subroutine parse_args(ell_start, ell_end, output_prefix)
        integer, intent(out) :: ell_start, ell_end
        character(len=*), intent(out) :: output_prefix

        integer :: nargs, ios
        character(len=256) :: arg

        nargs = command_argument_count()

        if (nargs < 1) then
            write(*,'(A)') 'Usage:'
            write(*,'(A)') '  ./sl_solve_ell_parallel ell'
            write(*,'(A)') '  ./sl_solve_ell_parallel ell_start ell_end'
            write(*,'(A)') '  ./sl_solve_ell_parallel ell_start ell_end output_prefix'
            stop
        end if

        call get_command_argument(1, arg)
        read(arg,*,iostat=ios) ell_start
        if (ios /= 0) then
            write(*,*) 'ERROR: could not parse ell_start from argument 1: ', trim(arg)
            stop
        end if

        if (nargs >= 2) then
            call get_command_argument(2, arg)
            read(arg,*,iostat=ios) ell_end
            if (ios /= 0) then
                write(*,*) 'ERROR: could not parse ell_end from argument 2: ', trim(arg)
                stop
            end if
        else
            ell_end = ell_start
        end if

        if (nargs >= 3) then
            call get_command_argument(3, output_prefix)
        else
            output_prefix = 'sl_part'
        end if

        if (ell_start < 0 .or. ell_end < 0) then
            write(*,*) 'ERROR: ell values must be nonnegative.'
            stop
        end if
        if (ell_end < ell_start) then
            write(*,*) 'ERROR: ell_end must be >= ell_start.'
            stop
        end if
    end subroutine parse_args

    !=================================================================
    ! Metric-derived coefficients
    !=================================================================
    real(kind=8) function gamma_rr(x)
        real(kind=8), intent(in) :: x
        gamma_rr = 1.0d0 + 1.0d0/x
    end function gamma_rr

    real(kind=8) function p(x)
        real(kind=8), intent(in) :: x
        p = x*x / dsqrt(gamma_rr(x))
    end function p

    real(kind=8) function w(x)
        real(kind=8), intent(in) :: x
        w = x*x * dsqrt(gamma_rr(x))
    end function w

    real(kind=8) function Q_func(ell, x)
        integer, intent(in) :: ell
        real(kind=8), intent(in) :: x
        Q_func = dble(ell) * dble(ell + 1) * dsqrt(gamma_rr(x))
    end function Q_func

    !=================================================================
    ! Build B = M^{-1/2} A M^{-1/2}
    !=================================================================
    subroutine build_B_matrix(Nr, ell, r, dr_edge, mdiag, D, E)
        integer, intent(in) :: Nr, ell
        real(kind=8), intent(in) :: r(0:Nr+1), dr_edge(0:Nr), mdiag(Nr)
        real(kind=8), intent(out) :: D(Nr), E(Nr-1)

        integer :: i
        real(kind=8) :: p_l, p_r, drc_i, A_ii

        do i = 1, Nr
            drc_i = 0.5d0*(dr_edge(i-1) + dr_edge(i))
            p_l = p(0.5d0*(r(i-1) + r(i)))
            p_r = p(0.5d0*(r(i) + r(i+1)))

            A_ii = p_r/dr_edge(i) + p_l/dr_edge(i-1) + Q_func(ell,r(i))*drc_i

            if (ell == 0 .and. i == 1) A_ii = A_ii - p_l/dr_edge(0)
            if (i == Nr) A_ii = A_ii - p_r/dr_edge(Nr)

            D(i) = A_ii / mdiag(i)
        end do

        do i = 1, Nr-1
            E(i) = -p(0.5d0*(r(i) + r(i+1))) &
                   / (dr_edge(i) * dsqrt(mdiag(i)*mdiag(i+1)))
        end do
    end subroutine build_B_matrix

    !=================================================================
    ! Grid
    !=================================================================
    real(kind=8) function local_dr(x)
        real(kind=8), intent(in) :: x
        real(kind=8) :: xi

        if (x <= 1.0d0) then
            local_dr = dr_min
        else if (x >= r_flat) then
            local_dr = dr_max
        else
            xi = dlog(0.5d0*x*(x + 1.0d0)) / dlog(0.5d0*r_flat*(r_flat + 1.0d0))
            local_dr = dr_min + (dr_max - dr_min)*xi
        end if
    end function local_dr

    real(kind=8) function next_dr(x)
        real(kind=8), intent(in) :: x
        real(kind=8), parameter :: eps = 1.0d-12

        next_dr = local_dr(x)

        if (x < 1.0d0 - eps .and. x + next_dr > 1.0d0) then
            next_dr = 1.0d0 - x
        end if
    end function next_dr

    subroutine count_grid(r_target, N_intervals, r_outer)
        real(kind=8), intent(in) :: r_target
        integer, intent(out) :: N_intervals
        real(kind=8), intent(out) :: r_outer

        real(kind=8) :: x, dx
        real(kind=8), parameter :: eps = 1.0d-12

        x = 0.0d0
        N_intervals = 0
        do while (x < r_target - eps)
            dx = next_dr(x)
            if (dx <= 0.0d0) then
                write(*,*) 'ERROR in count_grid: dx <= 0 at x=', x
                stop
            end if
            x = x + dx
            N_intervals = N_intervals + 1
        end do
        r_outer = x
    end subroutine count_grid

    subroutine build_grid(Nr, r, dr_edge)
        integer, intent(in) :: Nr
        real(kind=8), intent(out) :: r(0:Nr+1), dr_edge(0:Nr)

        integer :: i
        real(kind=8) :: dx

        r(0) = 0.0d0
        do i = 0, Nr
            dx = next_dr(r(i))
            if (dx <= 0.0d0) then
                write(*,*) 'ERROR in build_grid: dx <= 0 at i=', i, ' r=', r(i)
                stop
            end if
            r(i+1) = r(i) + dx
        end do

        do i = 0, Nr
            dr_edge(i) = r(i+1) - r(i)
            if (dr_edge(i) <= 0.0d0) then
                write(*,'(A,I8,2ES24.16)') 'ERROR: dr_edge <= 0 at i=', i, r(i), r(i+1)
                stop
            end if
        end do
    end subroutine build_grid

    !=================================================================
    ! HDF5 output for one ell
    !=================================================================
    subroutine write_single_ell_file(prefix, ell, Nr, N_intervals, r_outer, r, n_save_ell, k_save, modes_save)
        character(len=*), intent(in) :: prefix
        integer, intent(in) :: ell, Nr, N_intervals, n_save_ell
        real(kind=8), intent(in) :: r_outer
        real(kind=8), intent(in) :: r(0:Nr+1)
        real(kind=8), intent(in) :: k_save(n_save_ell)
        real(kind=8), intent(in) :: modes_save(Nr,n_save_ell)

        character(len=512) :: h5_file
        character(len=32)  :: group_name
        integer(hid_t) :: file_id, grp_id
        integer :: hdferr

        write(h5_file,'(A,"_ell_",I4.4,".h5")') trim(prefix), ell

        call h5fcreate_f(trim(h5_file), H5F_ACC_TRUNC_F, file_id, hdferr)
        call write_real_1d(file_id, 'r', r)

        call h5gcreate_f(file_id, 'params', grp_id, hdferr)
        call write_real_1d(grp_id, 'real',       [r_target, dr_min, dr_max, r_flat])
        call write_real_1d(grp_id, 'real_extra', [r_outer, nyquist_safety, k_safe, q_safe])
        call write_int_1d (grp_id, 'integer',    [Nr, N_intervals, ell, ell, n_modes_cap, n_save_ell])
        call h5gclose_f(grp_id, hdferr)

        write(group_name,'("ell_",I4.4)') ell
        call h5gcreate_f(file_id, trim(group_name), grp_id, hdferr)
        call write_real_1d(grp_id, 'k', k_save)
        call write_real_2d(grp_id, 'modes', modes_save)
        call write_int_scalar(grp_id, 'n_save', n_save_ell)
        call write_int_scalar(grp_id, 'ell', ell)
        call write_real_1d(grp_id, 'spectral_cut', [nyquist_safety, q_safe, k_safe])
        call h5gclose_f(grp_id, hdferr)

        call h5fclose_f(file_id, hdferr)
    end subroutine write_single_ell_file

    !=================================================================
    ! Diagnostics
    !=================================================================
    subroutine print_grid_info(Nr, N_intervals, r_target, r_outer, ell_start, ell_end, prefix)
        integer, intent(in) :: Nr, N_intervals, ell_start, ell_end
        real(kind=8), intent(in) :: r_target, r_outer
        character(len=*), intent(in) :: prefix

        write(*,'(A)')          'SL per-ell grid parameters:'
        write(*,'(A,F12.4)')    '  r_target        = ', r_target
        write(*,'(A,F12.4)')    '  r_outer         = ', r_outer
        write(*,'(A,ES12.4)')   '  overshoot       = ', r_outer - r_target
        write(*,'(A,F12.4)')    '  dr_min          = ', dr_min
        write(*,'(A,F12.4)')    '  dr_max          = ', dr_max
        write(*,'(A,F12.4)')    '  r_flat          = ', r_flat
        write(*,'(A,I10)')      '  N_intervals     = ', N_intervals
        write(*,'(A,I10)')      '  Nr              = ', Nr
        write(*,'(A,I10)')      '  n_modes_cap     = ', n_modes_cap
        write(*,'(A,F12.4)')    '  nyquist_safety  = ', nyquist_safety
        write(*,'(A,ES12.4)')   '  q_safe          = ', q_safe
        write(*,'(A,ES12.4)')   '  k_safe          = ', k_safe
        write(*,'(A,I6,A,I6)')  '  ell range       = ', ell_start, ' ... ', ell_end
        write(*,'(A,A)')        '  output_prefix   = ', trim(prefix)
    end subroutine print_grid_info

    subroutine print_mesh_quality(Nr, r, dr_edge)
        integer, intent(in) :: Nr
        real(kind=8), intent(in) :: r(0:Nr+1), dr_edge(0:Nr)

        integer :: i, imax
        real(kind=8) :: ratio, ratio_max, a, b

        ratio_max = 0.0d0
        imax = 0
        do i = 1, Nr
            a = dr_edge(i-1)
            b = dr_edge(i)
            ratio = max(a/b, b/a)
            if (ratio > ratio_max) then
                ratio_max = ratio
                imax = i
            end if
        end do

        write(*,'(A)')          'Mesh quality:'
        write(*,'(A,ES12.4)')   '  min(dr)      = ', minval(dr_edge)
        write(*,'(A,ES12.4)')   '  max(dr)      = ', maxval(dr_edge)
        write(*,'(A,ES12.4)')   '  last dr      = ', dr_edge(Nr)
        write(*,'(A,ES12.4)')   '  max ratio    = ', ratio_max
        write(*,'(A,I10)')      '  ratio index  = ', imax
        write(*,'(A,ES12.4)')   '  r(index)     = ', r(imax)
    end subroutine print_mesh_quality

    !=================================================================
    ! HDF5 helpers
    !=================================================================
    subroutine write_real_1d(parent_id, name, x)
        use hdf5
        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: x(:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(1)
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
        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: x(:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(1)
        integer :: hdferr

        dims(1) = int(size(x), kind=hsize_t)
        call h5screate_simple_f(1, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, x, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_1d

    subroutine write_int_scalar(parent_id, name, x)
        use hdf5
        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: x

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(1)
        integer :: hdferr
        integer :: buffer(1)

        buffer(1) = x
        dims(1) = 1
        call h5screate_simple_f(1, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, buffer, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_int_scalar

    subroutine write_real_2d(parent_id, name, a)
        use hdf5
        integer(hid_t), intent(in) :: parent_id
        character(len=*), intent(in) :: name
        real(kind=8), intent(in) :: a(:,:)

        integer(hid_t) :: dset_id, space_id
        integer(hsize_t) :: dims(2)
        integer :: hdferr

        dims(1) = int(size(a,1), kind=hsize_t)
        dims(2) = int(size(a,2), kind=hsize_t)
        call h5screate_simple_f(2, dims, space_id, hdferr)
        call h5dcreate_f(parent_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dset_id, hdferr)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, a, dims, hdferr)
        call h5dclose_f(dset_id, hdferr)
        call h5sclose_f(space_id, hdferr)
    end subroutine write_real_2d

end program SL_solve_ell_parallel