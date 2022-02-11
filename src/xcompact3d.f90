!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use MPI

  use var

  use transeq, only : calculate_transeq_rhs
  use navier, only : solve_poisson, cor_vel

  implicit none

  double precision :: tstart, tend, telapsed, tmin
  real :: trun
  integer :: ndt
  integer :: ierr
  integer :: mpi_real_type

  call init_xcompact3d(trun)

  telapsed = 0
  tmin = telapsed

  ndt = 0

  if (kind(telapsed) == kind(0.0d0)) then
     mpi_real_type = MPI_DOUBLE
  else
     mpi_real_type = MPI_FLOAT
  end if
  
  do while(tmin < trun)
     call init_flowfield()

     tstart = MPI_Wtime()

     call calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)

     ux1(:,:,:) = ux1(:,:,:) + dt * dux1(:,:,:,1)
     uy1(:,:,:) = uy1(:,:,:) + dt * duy1(:,:,:,1)
     uz1(:,:,:) = uz1(:,:,:) + dt * duz1(:,:,:,1)
     
     divu3(:,:,:) = zero
     call solve_poisson(pp3,px1,py1,pz1,ux1,uy1,uz1)
     call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

     tend = MPI_Wtime()
     telapsed = telapsed + (tend - tstart)

     call MPI_Allreduce(telapsed, tmin, 1, mpi_real_type, MPI_MIN, MPI_COMM_WORLD, ierr)
     if (nrank == 0) then
        print *, "Tmin = ", tmin, " of ", trun
     end if

     ndt = ndt + 1
  end do

  if (nrank == 0) then
     print *, "Elapsed time [s]: ", telapsed
     print *, "Timesteps completed: ", ndt
     print *, "Compute rate [dt/s]: ", ndt / telapsed
  end if

  call finalise_xcompact3d()

end program xcompact3d
!########################################################################
!########################################################################
subroutine init_xcompact3d(trun)

  use MPI
  use decomp_2d
  USE decomp_2d_poisson, ONLY : decomp_2d_poisson_init
  use case

  use var

  use variables, only : nx, ny, nz, nxm, nym, nzm
  use variables, only : p_row, p_col
  use variables, only : test_mode

  implicit none

  real, intent(inout) :: trun

  integer :: ierr

  integer :: nargin, arg, FNLength, status, DecInd
  logical :: back
  character(len=80) :: InputFN, FNBase

  !! Initialise MPI
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,nrank,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  ! Handle input file like a boss -- GD
  nargin=command_argument_count()

  !! Don't want to read input files - just basic numbers necessary for compute
  ! 1) nx = 16
  ! 2) ny = 16
  ! 3) nz = 16
  ! 4) p_row = 0
  ! 5) p_col = 0
  nx = 16; ny = 16; nz = 16
  p_row = 0; p_col = 0
  trun = 5.0
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
        trun = real(DecInd)
     elseif (arg.eq.7) then
        if (DecInd.eq.0) then
           test_mode = .false.
        else
           test_mode = .true.
        end if
     else
        print *, "Error: Too many arguments!"
        print *, "  x3div accepts"
        print *, "  1) nx (default=16)"
        print *, "  2) ny (default=16)"
        print *, "  3) nz (default=16)"
        print *, "  4) p_row (default=0)"
        print *, "  5) p_col (default=0)"
        print *, "  6) trun (default=5)"
        print *, "  7) test_mode logical 0/1 (default=0)"
     endif
  enddo

  call parameter()

  call decomp_2d_init(nx,ny,nz,p_row,p_col)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()

  call init_flowfield()

endsubroutine init_xcompact3d
!########################################################################
!########################################################################
subroutine init_flowfield()
  
  use case

  use var

  call init(rho1,ux1,uy1,uz1,ep1,phi1,drho1,dux1,duy1,duz1,dphi1,pp3,px1,py1,pz1)
  itime = 0

  divu3(:,:,:) = zero

end subroutine
!########################################################################
!########################################################################
subroutine finalise_xcompact3d()

  use MPI
  use decomp_2d

  implicit none

  integer :: ierr
  
  call decomp_2d_finalize
  CALL MPI_FINALIZE(ierr)

endsubroutine finalise_xcompact3d
