!################################################################################
!This file is part of Xcompact3d.
!
!Xcompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Xcompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Xcompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Xcompact3d/Incompact3d in your
!    publications and presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for
!    incompressible flows: a simple and efficient method with the quasi-spectral
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

program xcompact3d

  use MPI

  use var

  use transeq, only : calculate_transeq_rhs
  use navier, only : solve_poisson, cor_vel

  implicit none

  real :: tstart, tend, telapsed, trun, tmin
  integer :: ierr

  call init_xcompact3d(trun)

  telapsed = 0
  tmin = telapsed

  do while(tmin < trun)
     call init_flowfield()

     call cpu_time(tstart)

     call calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)

     ux1(:,:,:) = ux1(:,:,:) + dt * dux1(:,:,:,1)
     uy1(:,:,:) = uy1(:,:,:) + dt * duy1(:,:,:,1)
     uz1(:,:,:) = uz1(:,:,:) + dt * duz1(:,:,:,1)
     
     divu3(:,:,:) = zero
     call solve_poisson(pp3,px1,py1,pz1,ux1,uy1,uz1)
     call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

     call cpu_time(tend)
     telapsed = telapsed + (tend - tstart)

     call MPI_Allreduce(telapsed, tmin, 1, MPI_INT, MPI_MIN, MPI_COMM_WORLD, ierr)
     if (nrank == 0) then
        print *, "Tmin = ", tmin, " of ", trun
     end if
  end do

  if (nrank == 0) then
     print *, "Elapsed time: ", telapsed
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
     else
        print *, "Error: Too many arguments!"
        print *, "  x3div accepts"
        print *, "  1) nx"
        print *, "  2) ny"
        print *, "  3) nz"
        print *, "  4) p_row"
        print *, "  5) p_col"
     endif
  enddo

  call parameter(InputFN)

  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)    !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)    !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)
  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)
  call decomp_info_init(nxm, nym, nz, ph3)

  call init_variables()

  call schemes()

  call decomp_2d_poisson_init()
  call decomp_info_init(nxm,nym,nzm,phG)

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
