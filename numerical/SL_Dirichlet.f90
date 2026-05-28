program SL_Dirichlet 
    implicit none 

    real(kind=8) :: dr, r_max 
    integer :: Nr
    real(kind=8), allocatable :: r(:), omega(:)
    real(kind=8), allocatable :: Pw(:), Pw_prime(:)
    real(kind=8), allocatable :: C_above(:), C_below(:), C_diag(:)
    real(kind=8), allocatable :: C(:,:)
    real(kind=8), allocatable :: wr(:), wi(:), vr(:,:), vl(:,:)
    real(kind=8), allocatable :: work(:)
    real(kind=8) :: norm
    integer :: ell, ell_max
    integer :: i, j, info, lwork
    integer :: n_modes
    integer :: unit_lam, unit_modes
    character(len=128) :: file_lam, file_modes

    external dgeev

    ! Set parameters
    r_max = 800.0d0
    Nr = 8000
    dr = r_max / dble(Nr+1)
    ell_max = 5
    n_modes = 80

    allocate(r(Nr), omega(Nr))
    allocate(Pw(Nr), Pw_prime(Nr))
    allocate(C_above(Nr), C_below(Nr), C_diag(Nr))
    allocate(C(Nr,Nr))
    allocate(wr(Nr), wi(Nr), vr(Nr,Nr), vl(1,1))

    lwork = 4*Nr
    allocate(work(lwork))

    do i=1, Nr 
        r(i) = dble(i) * dr
    end do

    omega = r**2 * sqrt(1.0d0 + 1.0d0/r)

    Pw = 1.0d0 / (1.0d0 + 1.0d0/r)
    Pw_prime = (2.0d0 + 0.5d0/(r + 1.0d0)) / (r + 1.0d0)

    C_below = -Pw/dr**2 + Pw_prime/(2.0d0*dr)
    C_above = -Pw/dr**2 - Pw_prime/(2.0d0*dr)

    do ell=0, ell_max

        C_diag = 2.0d0*Pw/dr**2 + dble(ell)*(dble(ell)+1.0d0)/r**2

        C = 0.0d0

        do i=1, Nr
            C(i,i) = C_diag(i)

            if (i < Nr) then
                C(i,i+1) = C_above(i)
            end if

            if (i > 1) then
                C(i,i-1) = C_below(i)
            end if
        end do

        ! Condición de frontera izquierda:
        ! ell = 0  -> Neumann R'(0)=0, equivalente a R_0 = R_1
        ! ell > 0  -> Dirichlet R(0)=0, equivalente a R_0 = 0
        if (ell == 0) then
            C(1,1) = C_diag(1) + C_below(1)
        end if

        ! Condición de frontera derecha:
        ! Todos los ell -> Neumann R'(r_max)=0,
        ! equivalente a R_{Nr+1} = R_Nr
        C(Nr,Nr) = C(Nr,Nr) + C_above(Nr)


        call dgeev('N','V', Nr, C, Nr, wr, wi, vl, 1, vr, Nr, work, lwork, info)

        if (info /= 0) then
            write(*,*) 'DGEEV failed for ell = ', ell, ' info = ', info
            stop
        end if

        call sort_eigenpairs(Nr, wr, wi, vr)

        do j=1, min(n_modes,Nr)

            norm = 0.0d0

            do i=1, Nr
                norm = norm + omega(i)*vr(i,j)*vr(i,j)*dr
            end do

            norm = sqrt(norm)

            if (norm > 0.0d0) then
                do i=1, Nr
                    vr(i,j) = vr(i,j)/norm
                end do
            end if

        end do

        write(file_lam,'("lambda_ell_rmax_800",i4.4,".dat")') ell
        write(file_modes,'("modes_ell_rmax_800",i4.4,".dat")') ell

        open(newunit=unit_lam, file=file_lam, status='replace', action='write')

        do j=1, min(n_modes,Nr)
            write(unit_lam,'(i8,3es24.16)') j, wr(j), wi(j), sqrt(abs(wr(j)))
        end do

        close(unit_lam)

        open(newunit=unit_modes, file=file_modes, status='replace', action='write')

        do i=1, Nr
            write(unit_modes,'(es24.16)', advance='no') r(i)

            do j=1, min(n_modes,Nr)
                write(unit_modes,'(1x,es24.16)', advance='no') vr(i,j)
            end do

            write(unit_modes,*)
        end do

        close(unit_modes)

        write(*,*) 'Finished ell = ', ell

    end do

contains

    subroutine sort_eigenpairs(N, wr, wi, vr)
        implicit none

        integer, intent(in) :: N
        real(kind=8), intent(inout) :: wr(N), wi(N), vr(N,N)

        integer :: i, j, k
        real(kind=8) :: tmp
        real(kind=8), allocatable :: tmpvec(:)

        allocate(tmpvec(N))

        do i=1, N-1

            k = i

            do j=i+1, N
                if (wr(j) < wr(k)) then
                    k = j
                end if
            end do

            if (k /= i) then

                tmp = wr(i)
                wr(i) = wr(k)
                wr(k) = tmp

                tmp = wi(i)
                wi(i) = wi(k)
                wi(k) = tmp

                tmpvec = vr(:,i)
                vr(:,i) = vr(:,k)
                vr(:,k) = tmpvec

            end if

        end do

        deallocate(tmpvec)

    end subroutine sort_eigenpairs

end program SL_Dirichlet