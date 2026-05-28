module SL_config
    implicit none

    integer, parameter :: dp = kind(1.0d0)

    real(dp), parameter :: R_s = 1.0d0

    integer, parameter :: ell_min = 0
    integer, parameter :: ell_max = 80

    integer, parameter :: N_nodes = 18000
    integer, parameter :: N_elem  = N_nodes - 1

    real(dp), parameter :: r_min = 0.5d0
    real(dp), parameter :: r_max = 1200.5d0

    real(dp), parameter :: nyq_factor = 0.5d0
    real(dp), parameter :: lambda_min = -1.0d-12

contains

    real(dp) function p_coeff(r)
        implicit none
        real(dp), intent(in) :: r

        p_coeff = r**2 / dsqrt(1.0d0 + R_s/r)
    end function p_coeff


    real(dp) function w_coeff(r)
        implicit none
        real(dp), intent(in) :: r

        w_coeff = r**2 * dsqrt(1.0d0 + R_s/r)
    end function w_coeff


    real(dp) function q_coeff(r, ell)
        implicit none
        real(dp), intent(in) :: r
        integer, intent(in) :: ell

        q_coeff = dsqrt(1.0d0 + R_s/r) * dble(ell*(ell + 1))
    end function q_coeff

end module SL_config


module SL_hdf5_io
    use hdf5
    use SL_config, only: dp
    implicit none

contains

    subroutine write_real_1d(loc_id, name, data, n)
        implicit none

        integer(HID_T), intent(in) :: loc_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        real(dp), intent(in) :: data(n)

        integer :: error
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims = (/ int(n, HSIZE_T) /)

        call h5screate_simple_f(1, dims, space_id, error)
        call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

    end subroutine write_real_1d


    subroutine write_int_1d(loc_id, name, data, n)
        implicit none

        integer(HID_T), intent(in) :: loc_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: n
        integer, intent(in) :: data(n)

        integer :: error
        integer(HID_T) :: dset_id, space_id
        integer(HSIZE_T), dimension(1) :: dims

        dims = (/ int(n, HSIZE_T) /)

        call h5screate_simple_f(1, dims, space_id, error)
        call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_INTEGER, space_id, dset_id, error)
        call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, dims, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

    end subroutine write_int_1d


    subroutine write_real_2d(loc_id, name, data, n1, n2)
        implicit none

        integer(HID_T), intent(in) :: loc_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: n1, n2
        real(dp), intent(in) :: data(n1,n2)

        integer :: error
        integer(HID_T) :: dset_id, space_id
        integer(HID_T) :: dcpl_id
        integer(HSIZE_T), dimension(2) :: dims
        integer(HSIZE_T), dimension(2) :: chunk_dims

        dims = (/ int(n1, HSIZE_T), int(n2, HSIZE_T) /)

        chunk_dims(1) = int(n1, HSIZE_T)
        chunk_dims(2) = int(min(n2, 16), HSIZE_T)

        call h5screate_simple_f(2, dims, space_id, error)

        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_id, error)
        call h5pset_chunk_f(dcpl_id, 2, chunk_dims, error)

        call h5dcreate_f(loc_id, trim(name), H5T_NATIVE_DOUBLE, space_id, dset_id, error, dcpl_id=dcpl_id)
        call h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, data, dims, error)

        call h5pclose_f(dcpl_id, error)
        call h5dclose_f(dset_id, error)
        call h5sclose_f(space_id, error)

    end subroutine write_real_2d

end module SL_hdf5_io


program SL_solver
    use SL_config
    use SL_hdf5_io
    use hdf5
    implicit none

    integer :: i, j, ell
    integer :: n
    real(dp) :: dr

    real(dp), allocatable :: r(:)

    real(dp), allocatable :: K_diag(:), K_off(:)
    real(dp), allocatable :: M_diag(:), M_off(:)

    real(dp) :: r_mid, h
    real(dp) :: p_mid, q_mid, w_mid
    real(dp) :: kg, kp, mm

    real(dp) :: drho_i, max_drho
    real(dp) :: k_nyq, k_cut, lambda_cut

    integer :: ka, kb, ldab, ldbb, ldq, ldz, info
    integer :: il, iu, m_found

    character(len=1) :: jobz, range, uplo

    real(dp) :: VL, VU, abstol

    real(dp), allocatable :: AB(:,:), BB(:,:)
    real(dp), allocatable :: Q(:,:)
    real(dp), allocatable :: Q_dummy(:,:)
    real(dp), allocatable :: Z(:,:)
    real(dp), allocatable :: Z_dummy(:,:)

    real(dp), allocatable :: eig(:)
    real(dp), allocatable :: eig_save(:)
    real(dp), allocatable :: k_save(:)

    real(dp), allocatable :: work(:)
    integer, allocatable :: iwork(:)
    integer, allocatable :: ifail(:)

    integer :: error
    integer(HID_T) :: file_id
    integer(HID_T) :: group_id

    character(len=32) :: group_name
    character(len=*), parameter :: filename = "SL_modes_cut_ell.h5"

    real(dp), dimension(8) :: params_real
    integer, dimension(5) :: params_int

    real(dp), dimension(3) :: spectral_cut
    integer, dimension(2) :: ell_info

    external dsbgvx

    n  = N_nodes
    dr = (r_max - r_min) / dble(N_nodes - 1)

    allocate(r(n))

    do i = 1, n
        r(i) = r_min + dble(i-1)*dr
    end do

    max_drho = 0.0d0

    do i = 1, N_elem
        h = r(i+1) - r(i)
        r_mid = 0.5d0*(r(i) + r(i+1))

        drho_i = dsqrt(1.0d0 + R_s/r_mid) * h
        max_drho = max(max_drho, drho_i)
    end do

    k_nyq     = acos(-1.0d0) / max_drho
    k_cut     = nyq_factor * k_nyq
    lambda_cut = k_cut**2

    params_real(1) = R_s
    params_real(2) = r_min
    params_real(3) = r_max
    params_real(4) = dr
    params_real(5) = k_nyq
    params_real(6) = k_cut
    params_real(7) = lambda_cut
    params_real(8) = nyq_factor

    params_int(1) = N_nodes
    params_int(2) = N_elem
    params_int(3) = ell_min
    params_int(4) = ell_max
    params_int(5) = n

    call h5open_f(error)
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, error)

    call write_real_1d(file_id, "r", r, n)
    call write_real_1d(file_id, "params_real", params_real, 8)
    call write_int_1d(file_id, "params_int", params_int, 5)

    allocate(K_diag(n), K_off(n-1))
    allocate(M_diag(n), M_off(n-1))

    ka = 1
    kb = 1

    ldab = ka + 1
    ldbb = kb + 1

    allocate(AB(ldab,n), BB(ldbb,n))
    allocate(eig(n))
    allocate(work(7*n))
    allocate(iwork(5*n))
    allocate(ifail(n))

    allocate(Q_dummy(1,1))
    allocate(Z_dummy(1,1))

    jobz  = 'N'
    range = 'V'
    uplo  = 'U'

    VL = lambda_min
    VU = lambda_cut

    il = 1
    iu = n

    abstol = 0.0d0

    do ell = ell_min, ell_max

        print *, "Counting modes for ell =", ell

        K_diag = 0.0d0
        K_off  = 0.0d0

        M_diag = 0.0d0
        M_off  = 0.0d0

        do i = 1, N_elem

            h     = r(i+1) - r(i)
            r_mid = 0.5d0*(r(i) + r(i+1))

            p_mid = p_coeff(r_mid)
            q_mid = q_coeff(r_mid, ell)
            w_mid = w_coeff(r_mid)

            kg = p_mid / h
            kp = q_mid * h / 6.0d0
            mm = w_mid * h / 6.0d0

            K_diag(i)   = K_diag(i)   + kg + 2.0d0*kp
            K_diag(i+1) = K_diag(i+1) + kg + 2.0d0*kp

            K_off(i)    = K_off(i)    - kg + kp

            M_diag(i)   = M_diag(i)   + 2.0d0*mm
            M_diag(i+1) = M_diag(i+1) + 2.0d0*mm

            M_off(i)    = M_off(i)    + mm

        end do

        AB = 0.0d0
        BB = 0.0d0

        do j = 1, n
            AB(2,j) = K_diag(j)
            BB(2,j) = M_diag(j)
        end do

        AB(1,1) = 0.0d0
        BB(1,1) = 0.0d0

        do j = 2, n
            AB(1,j) = K_off(j-1)
            BB(1,j) = M_off(j-1)
        end do

        call dsbgvx(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                    Q_dummy, 1, VL, VU, il, iu, abstol, m_found, eig, &
                    Z_dummy, 1, work, iwork, ifail, info)

        if (info /= 0) then
            print *, "DSBGVX count failed for ell =", ell
            print *, "INFO =", info
            stop
        end if

        print *, "Solving ell =", ell, "modes found =", m_found

        allocate(Q(n,n))
        allocate(Z(n,m_found))

        AB = 0.0d0
        BB = 0.0d0

        do j = 1, n
            AB(2,j) = K_diag(j)
            BB(2,j) = M_diag(j)
        end do

        AB(1,1) = 0.0d0
        BB(1,1) = 0.0d0

        do j = 2, n
            AB(1,j) = K_off(j-1)
            BB(1,j) = M_off(j-1)
        end do

        jobz = 'V'
        ldq  = n
        ldz  = n

        call dsbgvx(jobz, range, uplo, n, ka, kb, AB, ldab, BB, ldbb, &
                    Q, ldq, VL, VU, il, iu, abstol, m_found, eig, &
                    Z, ldz, work, iwork, ifail, info)

        if (info /= 0) then
            print *, "DSBGVX solve failed for ell =", ell
            print *, "INFO =", info
            stop
        end if

        allocate(eig_save(m_found))
        allocate(k_save(m_found))

        do j = 1, m_found

            eig_save(j) = eig(j)

            if (eig(j) > 0.0d0) then
                k_save(j) = dsqrt(eig(j))
            else
                k_save(j) = 0.0d0
            end if

        end do

        write(group_name,'("ell_",I4.4)') ell

        call h5gcreate_f(file_id, trim(group_name), group_id, error)

        call write_real_1d(group_id, "eigenvalues", eig_save, m_found)
        call write_real_1d(group_id, "k", k_save, m_found)
        call write_real_2d(group_id, "modes", Z, n, m_found)

        spectral_cut(1) = k_nyq
        spectral_cut(2) = k_cut
        spectral_cut(3) = lambda_cut

        ell_info(1) = ell
        ell_info(2) = m_found

        call write_real_1d(group_id, "spectral_cut", spectral_cut, 3)
        call write_int_1d(group_id, "ell_info", ell_info, 2)

        call h5gclose_f(group_id, error)

        deallocate(Q)
        deallocate(Z)
        deallocate(eig_save)
        deallocate(k_save)

        jobz = 'N'

    end do

    call h5fclose_f(file_id, error)
    call h5close_f(error)

    deallocate(r)

    deallocate(K_diag, K_off)
    deallocate(M_diag, M_off)

    deallocate(AB, BB)
    deallocate(eig)
    deallocate(work)
    deallocate(iwork)
    deallocate(ifail)

    deallocate(Q_dummy)
    deallocate(Z_dummy)

    print *, "Done."
    print *, "Saved file: ", trim(filename)

end program SL_solver
