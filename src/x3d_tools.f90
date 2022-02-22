!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!##################################################################
function rl(complexnumber)

  !use param
  use decomp_2d, only : mytype

  implicit none

  real(mytype) :: rl
  complex(mytype) :: complexnumber

  rl = real(complexnumber, kind=mytype)

end function rl
!##################################################################

!##################################################################
function iy(complexnumber)

  !use param
  use decomp_2d, only : mytype

  implicit none

  real(mytype) :: iy
  complex(mytype) :: complexnumber

  iy = aimag(complexnumber)

end function iy
!##################################################################

!##################################################################
function cx(realpart,imaginarypart)

  !use param
  use decomp_2d, only : mytype

  implicit none

  complex(mytype) :: cx
  real(mytype) :: realpart, imaginarypart

  cx = cmplx(realpart, imaginarypart, kind=mytype)

end function cx
!################################################################################
! SUBROUTINE DEFINITION
!##################################################################
subroutine boot_xcompact3d()
  
  use MPI
  use decomp_2d, only : nrank, nproc, decomp_2d_abort
  
  implicit none

  integer :: code

  !! Initialise MPI
  call MPI_INIT(code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_INIT")
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_RANK")
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,code)
  if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_COMM_SIZE")


endsubroutine boot_xcompact3d
!########################################################################

!########################################################################
subroutine init_xcompact3d(ndt_max)

  use decomp_2d, only : decomp_2d_init, decomp_info_init
  use decomp_2d, only : init_coarser_mesh_statS, &
                        init_coarser_mesh_statV, &
                        init_coarser_mesh_statP
  use decomp_2d, only : ph1, ph2, ph3, phG
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use x3d_operator_x_data, only : x3d_operator_x_data_init
  use x3d_operator_y_data, only : x3d_operator_y_data_init
  use x3d_operator_z_data, only : x3d_operator_z_data_init
  use x3d_operator_1d, only : x3d_operator_1d_init
  use x3d_derive, only : x3d_derive_init
  use parameters
  use case

  use var

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : test_mode
  use variables, only : nstat, nprobe, nvisu

  implicit none

  !real, intent(inout) :: trun
  integer, intent(inout) :: ndt_max

  integer :: nargin, arg, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()

  !! Don't want to read input files - just basic numbers necessary for compute
  ! 1) nx = 16
  ! 2) ny = 16
  ! 3) nz = 16
  ! 4) p_row = 0
  ! 5) p_col = 0
  nx = 32; ny = 32; nz = 32
  p_row = 0; p_col = 0
  !trun = 5.0
  ndt_max = 5
  test_mode = .false. 
  do arg = 1, nargin
     call get_command_argument(arg, InputFN, FNLength, status)
     read(InputFN, *, iostat=status) DecInd
     if (arg.eq.1) then
        nx = DecInd
     elseif (arg.eq.2) then
        ny = DecInd
     elseif (arg.eq.3) then
        nz = DecInd
     elseif (arg.eq.4) then
        p_row = DecInd
     elseif (arg.eq.5) then
        p_col = DecInd
     elseif (arg.eq.6) then
        !trun = real(DecInd)
        ndt_max = DecInd
     elseif (arg.eq.7) then
        if (DecInd.eq.0) then
           test_mode = .false.
        else
           test_mode = .true.
        end if
        write(*,*) 'Test mode ', test_mode
     else
        print *, "Error: Too many arguments!"
        print *, "  x3div accepts"
        print *, "  1) nx (default=16)"
        print *, "  2) ny (default=16)"
        print *, "  3) nz (default=16)"
        print *, "  4) p_row (default=0)"
        print *, "  5) p_col (default=0)"
        print *, "  6) ndt_max (default=10)"
        print *, "  7) test_mode logical 0/1 (default=0)"
     endif
  enddo

  call parameter()
  call case_boot()
  call listing()
  call case_listing()

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call var_init()
  call x3d_operator_x_data_init(nx, nzm)
  call x3d_operator_y_data_init(ny, nym)
  call x3d_operator_z_data_init(nz, nzm)
  call x3d_operator_1d_init()
  call x3d_derive_init()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine finalise_xcompact3d(flag)

  use MPI
  use decomp_2d, only : decomp_2d_finalize, decomp_info_finalize, &
                        ph1, ph2, ph3, phG
  use decomp_2d_poisson, only : decomp_2d_poisson_finalize
  use x3d_operator_x_data, only : x3d_operator_x_data_finalize
  use x3d_operator_y_data, only : x3d_operator_y_data_finalize
  use x3d_operator_z_data, only : x3d_operator_z_data_finalize
  use x3d_operator_1d, only : x3d_operator_1d_finalize
  use x3d_derive, only : x3d_derive_finalize
  use var, only : var_finalize
  use variables, only : nz

  implicit none

  logical, intent(in) :: flag
  integer :: ierr

  call decomp_info_finalize(ph1)
  call decomp_info_finalize(ph2)
  call decomp_info_finalize(ph3)
  call decomp_info_finalize(phG)
  call decomp_2d_poisson_finalize()

  call x3d_derive_finalize()
  call x3d_operator_1d_finalize()
  call x3d_operator_x_data_finalize()
  call x3d_operator_y_data_finalize()
  call x3d_operator_z_data_finalize(nz)
  call var_finalize()

  call decomp_2d_finalize()
  if (flag) then
    CALL MPI_FINALIZE(ierr)
  endif

endsubroutine finalise_xcompact3d
