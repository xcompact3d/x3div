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
  use decomp_2d, only : nrank, xsize, real_type, decomp_2d_warning 
  use param,   only : dt, zero, itr
  use transeq, only : calculate_transeq_rhs
  use navier,  only : solve_poisson, cor_vel
  use mom,     only : test_du, test_dv, test_dw
  use nvtx

  implicit none

  double precision :: tstart, tend, telapsed, tmin, tmax
  !real :: trun
  integer :: i, j, k
  integer :: ndt, ndt_max
  integer :: code

  call boot_xcompact3d()

  call init_xcompact3d(ndt_max)

  telapsed = 0
  tmin = telapsed

  ndt = 1

  do while(ndt < ndt_max)
     itr = 1 ! no inner iterations
     !call init_flowfield()

     tstart = MPI_Wtime()

     call nvtxStartRange("transeq")
     call calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)
     call nvtxEndRange
     do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
       ux1(i,j,k) = ux1(i,j,k) + dt * dux1(i,j,k,1)
       uy1(i,j,k) = uy1(i,j,k) + dt * duy1(i,j,k,1)
       uz1(i,j,k) = uz1(i,j,k) + dt * duz1(i,j,k,1)
     enddo
     
     !do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
     !  divu3(:,:,:) = zero
     !enddo
     call nvtxStartRange("solve_poisson")
     call solve_poisson(pp3,px1,py1,pz1,ux1,uy1,uz1)
     call nvtxEndRange
     call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

     tend = MPI_Wtime()
     telapsed = telapsed + (tend - tstart)
     tmin = telapsed
     tmax = telapsed

     !call MPI_Allreduce(telapsed, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, code)
     !if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_Allreduce")
     !call MPI_Allreduce(telapsed, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, code)
     !if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_Allreduce")
     if (nrank == 0) then
        print *, "Elapse time min ", tmin, " max ", tmax
     end if

     ndt = ndt + 1
  end do

  if (nrank == 0) then
     print *, "Elapsed time (min-max) [s]: ", tmin, tmax
     print *, "Timesteps completed: ", ndt
     print *, "Compute rate (min-max)[dt/s]: ", ndt / tmin, ndt /  tmax
  end if

  call finalise_xcompact3d(.true.)

end program xcompact3d
!########################################################################
!########################################################################

