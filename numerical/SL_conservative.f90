program SL_conservative 
    implicit none 

    real(kind=8) :: dr, r_max
    real(kind=8) :: b_i, q_i
    real(kind=8) :: norm
    integer :: Nr
    integer :: ell, ell_max
    integer :: i, j, info, lwork
    integer :: n_modes
    integer :: unit_lam, unit_modes
    character(len=128) :: file_lam, file_modes

    real(kind=8), allocatable :: r(:), omega(:)
    real(kind=8), allocatable :: D(:), E(:), Z(:,:)
    real(kind=8), allocatable :: work(:)

    external dstev

    ! Set parameters
    r_max = 800.0d0
    Nr = 8000
    dr = r_max / dble(Nr+1)
    ell_max = 5
    n_modes = 80

    allocate(r(Nr), omega(Nr))
    allocate(D(Nr), E(Nr-1), Z(Nr,Nr))

    lwork = max(1, 2*Nr - 2)
    allocate(work(lwork))

    do i=1, Nr 
        r(i) = dble(i) * dr
    end do

    omega = r**2 * dsqrt(1.0d0 + 1.0d0/r)

    do ell=0, ell_max

        !------------------------------------------------------------
        ! Construcción de la diagonal D de B = W^{-1/2} A W^{-1/2}
        !------------------------------------------------------------
        do i=1, Nr

            q_i = dble(ell)*(dble(ell)+1.0d0) * dsqrt(1.0d0 + 1.0d0/r(i))

            b_i = (p(r(i)+0.5d0*dr) + p(r(i)-0.5d0*dr))/dr**2 + q_i

            ! Frontera izquierda:
            ! ell = 0 -> Neumann R'(0)=0, equivalente a R_0 = R_1.
            ! Esto suma A_below(1) a la diagonal de A.
            if (ell == 0 .and. i == 1) then
                b_i = b_i - p(r(i)-0.5d0*dr)/dr**2
            end if

            ! Frontera derecha:
            ! Todos los ell -> Neumann R'(r_max)=0,
            ! equivalente a R_{Nr+1} = R_Nr.
            ! Esto suma A_above(Nr) a la diagonal de A.
            if (i == Nr) then
                b_i = b_i - p(r(i)+0.5d0*dr)/dr**2
            end if

            D(i) = b_i / omega(i)

        end do

        !------------------------------------------------------------
        ! Construcción de la sub/superdiagonal E de B
        !------------------------------------------------------------
        do i=1, Nr-1
            E(i) = -p(r(i)+0.5d0*dr) / (dr**2 * dsqrt(omega(i)*omega(i+1)))
        end do

        !------------------------------------------------------------
        ! Resolver B y = lambda y
        ! D contiene la diagonal de entrada y sale con eigenvalores.
        ! E contiene la off-diagonal de entrada y es destruida.
        ! Z sale con los eigenvectores y.
        !------------------------------------------------------------
        call dstev('V', Nr, D, E, Z, Nr, work, info)

        if (info /= 0) then
            write(*,*) 'DSTEV failed for ell = ', ell, ' info = ', info
            stop
        end if

        !------------------------------------------------------------
        ! Recuperar modos físicos R = W^{-1/2} y
        !------------------------------------------------------------
        do j=1, min(n_modes,Nr)
            do i=1, Nr
                Z(i,j) = Z(i,j) / dsqrt(omega(i))
            end do
        end do

        !------------------------------------------------------------
        ! Normalización continua:
        ! sum_i omega_i R_i^2 dr = 1
        !------------------------------------------------------------
        do j=1, min(n_modes,Nr)

            norm = 0.0d0

            do i=1, Nr
                norm = norm + omega(i)*Z(i,j)*Z(i,j)*dr
            end do

            norm = dsqrt(norm)

            if (norm > 0.0d0) then
                do i=1, Nr
                    Z(i,j) = Z(i,j)/norm
                end do
            end if

        end do

        write(file_lam,'("lambda_cons_ell_rmax_800",i4.4,".dat")') ell
        write(file_modes,'("modes_cons_ell_rmax_800",i4.4,".dat")') ell

        open(newunit=unit_lam, file=file_lam, status='replace', action='write')

        do j=1, min(n_modes,Nr)
            write(unit_lam,'(i8,3es24.16)') j, D(j), 0.0d0, dsqrt(abs(D(j)))
        end do

        close(unit_lam)

        open(newunit=unit_modes, file=file_modes, status='replace', action='write')

        !------------------------------------------------------------
        ! Escribir frontera izquierda normalizada
        !
        ! ell = 0  -> R(0) = R_1
        ! ell > 0  -> R(0) = 0
        !------------------------------------------------------------
        write(unit_modes,'(es24.16)', advance='no') 0.0d0

        do j=1, min(n_modes,Nr)
            if (ell == 0) then
                write(unit_modes,'(1x,es24.16)', advance='no') Z(1,j)
            else
                write(unit_modes,'(1x,es24.16)', advance='no') 0.0d0
            end if
        end do

        write(unit_modes,*)

        ! Puntos interiores
        do i=1, Nr
            write(unit_modes,'(es24.16)', advance='no') r(i)

            do j=1, min(n_modes,Nr)
                write(unit_modes,'(1x,es24.16)', advance='no') Z(i,j)
            end do

            write(unit_modes,*)
        end do

        !------------------------------------------------------------
        ! Escribir frontera derecha normalizada
        !
        ! Neumann en r_max -> R(r_max) = R_Nr
        !------------------------------------------------------------
        write(unit_modes,'(es24.16)', advance='no') r_max

        do j=1, min(n_modes,Nr)
            write(unit_modes,'(1x,es24.16)', advance='no') Z(Nr,j)
        end do

        write(unit_modes,*)

        close(unit_modes)

        write(*,*) 'Finished ell = ', ell

    end do

contains

    real(kind=8) function p(r_in)
        implicit none
        real(kind=8), intent(in) :: r_in

        p = r_in**2 / dsqrt(1.0d0 + 1.0d0/r_in)

    end function p

end program SL_conservative