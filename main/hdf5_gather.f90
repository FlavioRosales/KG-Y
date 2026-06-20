module hdf5_gather
  use iso_fortran_env, only: real64, int64
  use hdf5_lib,  only: create_hdf5_file, h5_write_rgrid, &
                       h5_read_3d, h5_create_3d_dset, h5_write_3d_block, &
                       h5_read_1d
  use kg_mpi,    only: modes_total, n_range_for_rank
  implicit none

  private
  public :: assemble_global_files, assemble_ncharge_global

  ! Nombres de archivos globales (resultado final)
  character(len=*), parameter :: file_phire_global = "rePhi.h5"
  character(len=*), parameter :: file_phiim_global = "imPhi.h5"
  character(len=*), parameter :: file_psire_global = "rePsi.h5"
  character(len=*), parameter :: file_psiim_global = "imPsi.h5"
  character(len=*), parameter :: file_pire_global  = "rePi.h5"
  character(len=*), parameter :: file_piim_global  = "imPi.h5"

contains

  !---------------------------------------------------------------------------
  ! Interfaz que llama el main: ensambla TODOS los campos (φ, ψ, π).
  !   - ell_max   -> para saber n_total = (ell_max+1)^2
  !   - Nsnap     -> snapshots sin contar t=0 (Nt3 = Nsnap+1)
  !   - Nr        -> puntos radiales
  !   - nproc     -> nº de ranks
  !   - rgrid(:)  -> malla radial global
  !---------------------------------------------------------------------------
  subroutine assemble_global_files(ell_max, Nsnap, Nr, nproc, rgrid)
    integer,      intent(in) :: ell_max, Nsnap, Nr, nproc
    real(real64), intent(in) :: rgrid(:)

    integer(int64) :: n_total

    n_total = modes_total(ell_max)

    call assemble_one_field_3d("rePhi", file_phire_global, n_total, nproc, Nsnap, Nr, rgrid)
    call assemble_one_field_3d("imPhi", file_phiim_global, n_total, nproc, Nsnap, Nr, rgrid)
    call assemble_one_field_3d("rePsi", file_psire_global, n_total, nproc, Nsnap, Nr, rgrid)
    call assemble_one_field_3d("imPsi", file_psiim_global, n_total, nproc, Nsnap, Nr, rgrid)
    call assemble_one_field_3d("rePi",  file_pire_global,  n_total, nproc, Nsnap, Nr, rgrid)
    call assemble_one_field_3d("imPi",  file_piim_global,  n_total, nproc, Nsnap, Nr, rgrid)
  end subroutine assemble_global_files

  !---------------------------------------------------------------------------
  ! Ensambla un solo campo 3D:
  !   local (por rank):  data(Nt3, Nr, Nmodes_rank)
  !   global:            data(Nt3, Nr, n_total)
  !---------------------------------------------------------------------------
  subroutine assemble_one_field_3d(base_name, global_filename, n_total, nproc, Nsnap, Nr, rgrid)
    character(len=*), intent(in) :: base_name, global_filename
    integer(int64),   intent(in) :: n_total
    integer,          intent(in) :: nproc, Nsnap, Nr
    real(real64),     intent(in) :: rgrid(:)

    integer           :: r, ierr_h5, ierr_cmd
    integer(int64)    :: n0r, n1r
    integer           :: Nt3_global, Nmodes_rank
    character(len=256) :: srcfile, cmd
    real(real64), allocatable :: buf(:,:,:)

    Nt3_global = Nsnap + 1   ! t=0..Nsnap -> Nt3_global capas temporales

    ! Crear/truncar archivo global + escribir r-grid y dataset 3D vacío
    call create_hdf5_file(trim(global_filename))
    call h5_write_rgrid(trim(global_filename), rgrid, ierr_h5)
    call h5_create_3d_dset(trim(global_filename), "data", Nt3_global, Nr, &
                           int(n_total,kind=4), ierr_h5)

    do r = 0, nproc-1
      ! Archivo local de este rank para este campo
      write(srcfile, '(A,"_rank",I4.4,".h5")') trim(base_name), r

      ! Rango de modos globales que le tocaban a este rank
      call n_range_for_rank(r, nproc, n_total, n0r, n1r)
      Nmodes_rank = int(n1r - n0r + 1, kind=4)

      ! Leer dataset "data" de este rank: buf(Nt3_global, Nr, Nmodes_rank)
      call h5_read_3d(trim(srcfile), "data", buf, ierr_h5)
      if (ierr_h5 /= 0) then
        write(*,'("[gather] WARNING: no pude leer data de ",A," (ierr=",I0,")")') &
              trim(srcfile), ierr_h5
      else
        if ( size(buf,1) /= Nt3_global .or. size(buf,2) /= Nr .or. size(buf,3) /= Nmodes_rank ) then
          write(*,*) "[gather] WARNING: dims inesperadas en ",trim(srcfile), &
                     " leí (", size(buf,1),",",size(buf,2),",",size(buf,3), &
                     ") pero esperaba (", Nt3_global,",",Nr,",",Nmodes_rank,")"
        end if

        ! Escribir bloque de modos [n0r..n1r] en el dataset global
        ! El eje 3 es el índice de modo (1-based): j = n0r+1 .. n1r+1
        call h5_write_3d_block(trim(global_filename), "data", buf, &
                               int(n0r,kind=4)+1, ierr_h5)
      end if

      if (allocated(buf)) deallocate(buf)

      ! Eliminar archivo temporal de este rank
      write(cmd, '(A,A)') 'rm -f ', trim(srcfile)
      call execute_command_line(trim(cmd), exitstat=ierr_cmd)
    end do

  end subroutine assemble_one_field_3d

    !---------------------------------------------------------------------------
  ! Ensambla los diagnósticos 1D (t, Noether, Flujo) en un archivo global
  !   ncharge.h5, copiando los datasets
  !   "mode_XXXXX_t", "mode_XXXXX_N", "mode_XXXXX_F"
  !   desde ncharge_rankXXXX.h5 (uno por rank).
  !
  ! Se asume que:
  !   - Cada rank r escribe en:  ncharge_rankrrrr.h5
  !   - Los modos que maneja cada rank se obtienen con n_range_for_rank.
  !---------------------------------------------------------------------------
  subroutine assemble_ncharge_global(ell_max, nproc)
    use iso_fortran_env, only: real64, int64
    implicit none
    integer,          intent(in) :: ell_max, nproc

    integer(int64) :: n_total, n0r, n1r, nmode
    integer        :: r, ierr_h5, ierr_cmd
    integer        :: nlen
    character(len=256) :: srcfile, cmd, base, dname_t, dname_N, dname_F
    real(real64), allocatable :: v(:)

    ! Número total de modos globales
    n_total = modes_total(ell_max)

    ! Crear/truncar archivo global de diagnósticos
    call create_hdf5_file("ncharge.h5")

    do r = 0, nproc-1
      ! Archivo fuente de este rank
      write(srcfile, '(A,"_rank",I4.4,".h5")') "ncharge", r

      ! Rango global de modos que le tocaron a este rank
      call n_range_for_rank(r, nproc, n_total, n0r, n1r)

      do nmode = n0r, n1r
        write(base, '("mode_",I5.5)') nmode
        dname_t = trim(base)//"_t"
        dname_N = trim(base)//"_N"
        dname_F = trim(base)//"_F"

        ! --- Copiar t ---
        call h5_read_1d(trim(srcfile), trim(dname_t), v, nlen, ierr_h5)
        if (ierr_h5 == 0) then
          call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_t))
        else
          write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
                trim(dname_t), trim(srcfile), ierr_h5
        end if
        if (allocated(v)) deallocate(v)

        ! --- Copiar Noether ---
        call h5_read_1d(trim(srcfile), trim(dname_N), v, nlen, ierr_h5)
        if (ierr_h5 == 0) then
          call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_N))
        else
          write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
                trim(dname_N), trim(srcfile), ierr_h5
        end if
        if (allocated(v)) deallocate(v)

        ! --- Copiar Flujo ---
        call h5_read_1d(trim(srcfile), trim(dname_F), v, nlen, ierr_h5)
        if (ierr_h5 == 0) then
          call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_F))
        else
          write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
                trim(dname_F), trim(srcfile), ierr_h5
        end if
        if (allocated(v)) deallocate(v)

      end do  ! nmode

      ! Eliminar archivo temporal de este rank
      write(cmd, '(A,A)') 'rm -f ', trim(srcfile)
      call execute_command_line(trim(cmd), exitstat=ierr_cmd)

    end do  ! r = 0..nproc-1

  end subroutine assemble_ncharge_global


! end module hdf5_gather
! module hdf5_gather
!   use iso_fortran_env, only: real64, int64
!   use hdf5_lib,  only: create_hdf5_file, h5_write_rgrid, &
!                        h5_read_3d, h5_create_3d_dset, h5_write_3d_block, &
!                        h5_read_1d
!   use kg_mpi,    only: n_range_for_rank   ! <<< ya no necesitas modes_total
!   implicit none

!   private
!   public :: assemble_global_files, assemble_ncharge_global

!   character(len=*), parameter :: file_phire_global = "rePhi.h5"
!   character(len=*), parameter :: file_phiim_global = "imPhi.h5"
!   character(len=*), parameter :: file_psire_global = "rePsi.h5"
!   character(len=*), parameter :: file_psiim_global = "imPsi.h5"
!   character(len=*), parameter :: file_pire_global  = "rePi.h5"
!   character(len=*), parameter :: file_piim_global  = "imPi.h5"

! contains

!   !---------------------------------------------------------------------------
!   ! Ensambla TODOS los campos (φ, ψ, π) asumiendo SOLO ell:
!   !   - total global = Nell = ell_max + 1   (ell = 0..ell_max)
!   !---------------------------------------------------------------------------
!   subroutine assemble_global_files(ell_max, Nsnap, Nr, nproc, rgrid)
!     integer,      intent(in) :: ell_max, Nsnap, Nr, nproc
!     real(real64), intent(in) :: rgrid(:)

!     integer(int64) :: n_total

!     ! Total global de "modos" ahora = número de ells
!     n_total = int(ell_max + 1, kind=int64)

!     call assemble_one_field_3d("rePhi", file_phire_global, n_total, nproc, Nsnap, Nr, rgrid)
!     call assemble_one_field_3d("imPhi", file_phiim_global, n_total, nproc, Nsnap, Nr, rgrid)
!     call assemble_one_field_3d("rePsi", file_psire_global, n_total, nproc, Nsnap, Nr, rgrid)
!     call assemble_one_field_3d("imPsi", file_psiim_global, n_total, nproc, Nsnap, Nr, rgrid)
!     call assemble_one_field_3d("rePi",  file_pire_global,  n_total, nproc, Nsnap, Nr, rgrid)
!     call assemble_one_field_3d("imPi",  file_piim_global,  n_total, nproc, Nsnap, Nr, rgrid)
!   end subroutine assemble_global_files

!   !---------------------------------------------------------------------------
!   ! Ensambla un campo 3D:
!   !   local:  data(Nt3, Nr, Nell_rank)
!   !   global: data(Nt3, Nr, Nell_total)
!   !---------------------------------------------------------------------------
!   subroutine assemble_one_field_3d(base_name, global_filename, n_total, nproc, Nsnap, Nr, rgrid)
!     character(len=*), intent(in) :: base_name, global_filename
!     integer(int64),   intent(in) :: n_total
!     integer,          intent(in) :: nproc, Nsnap, Nr
!     real(real64),     intent(in) :: rgrid(:)

!     integer            :: r, ierr_h5, ierr_cmd
!     integer(int64)     :: n0r, n1r
!     integer            :: Nt3_global, Nmodes_rank
!     character(len=256) :: srcfile, cmd
!     real(real64), allocatable :: buf(:,:,:)

!     Nt3_global = Nsnap + 1

!     call create_hdf5_file(trim(global_filename))
!     call h5_write_rgrid(trim(global_filename), rgrid, ierr_h5)
!     call h5_create_3d_dset(trim(global_filename), "data", Nt3_global, Nr, &
!                            int(n_total,kind=4), ierr_h5)

!     do r = 0, nproc-1

!       write(srcfile, '(A,"_rank",I4.4,".h5")') trim(base_name), r

!       ! Rango global [n0r..n1r] ahora es rango de ell (contiguo)
!       call n_range_for_rank(r, nproc, n_total, n0r, n1r)
!       Nmodes_rank = int(n1r - n0r + 1_int64, kind=4)

!       call h5_read_3d(trim(srcfile), "data", buf, ierr_h5)
!       if (ierr_h5 /= 0) then
!         write(*,'("[gather] WARNING: no pude leer data de ",A," (ierr=",I0,")")') &
!               trim(srcfile), ierr_h5
!       else
!         if ( size(buf,1) /= Nt3_global .or. size(buf,2) /= Nr .or. size(buf,3) /= Nmodes_rank ) then
!           write(*,*) "[gather] WARNING: dims inesperadas en ",trim(srcfile), &
!                      " leí (", size(buf,1),",",size(buf,2),",",size(buf,3), &
!                      ") pero esperaba (", Nt3_global,",",Nr,",",Nmodes_rank,")"
!         end if

!         ! Escribe el bloque a partir del índice global (1-based)
!         call h5_write_3d_block(trim(global_filename), "data", buf, &
!                                int(n0r,kind=4)+1, ierr_h5)
!       end if

!       if (allocated(buf)) deallocate(buf)

!       write(cmd, '(A,A)') 'rm -f ', trim(srcfile)
!       call execute_command_line(trim(cmd), exitstat=ierr_cmd)

!     end do

!   end subroutine assemble_one_field_3d

!   !---------------------------------------------------------------------------
!   ! Ensambla diagnósticos 1D en ncharge.h5.
!   ! IMPORTANTE:
!   !   Esto asume que el main escribió datasets como:
!   !     "mode_XXXXX_t", "mode_XXXXX_N", "mode_XXXXX_F"
!   !   donde XXXXX ahora corresponde a ell (porque estás evolucionando solo m=0).
!   !---------------------------------------------------------------------------
!   subroutine assemble_ncharge_global(ell_max, nproc)
!     implicit none
!     integer, intent(in) :: ell_max, nproc

!     integer(int64) :: n_total, n0r, n1r, nmode
!     integer        :: r, ierr_h5, ierr_cmd
!     integer        :: nlen
!     character(len=256) :: srcfile, cmd, base, dname_t, dname_N, dname_F
!     real(real64), allocatable :: v(:)

!     n_total = int(ell_max + 1, kind=int64)

!     call create_hdf5_file("ncharge.h5")

!     do r = 0, nproc-1

!       write(srcfile, '(A,"_rank",I4.4,".h5")') "ncharge", r

!       call n_range_for_rank(r, nproc, n_total, n0r, n1r)

!       do nmode = n0r, n1r
!         ! nmode aquí es ell global
!         write(base, '("mode_",I5.5)') nmode
!         dname_t = trim(base)//"_t"
!         dname_N = trim(base)//"_N"
!         dname_F = trim(base)//"_F"

!         call h5_read_1d(trim(srcfile), trim(dname_t), v, nlen, ierr_h5)
!         if (ierr_h5 == 0) then
!           call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_t))
!         else
!           write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
!                 trim(dname_t), trim(srcfile), ierr_h5
!         end if
!         if (allocated(v)) deallocate(v)

!         call h5_read_1d(trim(srcfile), trim(dname_N), v, nlen, ierr_h5)
!         if (ierr_h5 == 0) then
!           call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_N))
!         else
!           write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
!                 trim(dname_N), trim(srcfile), ierr_h5
!         end if
!         if (allocated(v)) deallocate(v)

!         call h5_read_1d(trim(srcfile), trim(dname_F), v, nlen, ierr_h5)
!         if (ierr_h5 == 0) then
!           call h5_write_rgrid("ncharge.h5", v, ierr_h5, trim(dname_F))
!         else
!           write(*,'("[gather-ncharge] WARNING: no pude leer ",A," de ",A," (ierr=",I0,")")') &
!                 trim(dname_F), trim(srcfile), ierr_h5
!         end if
!         if (allocated(v)) deallocate(v)

!       end do

!       write(cmd, '(A,A)') 'rm -f ', trim(srcfile)
!       call execute_command_line(trim(cmd), exitstat=ierr_cmd)

! !     end do

! !   end subroutine assemble_ncharge_global

end module hdf5_gather
