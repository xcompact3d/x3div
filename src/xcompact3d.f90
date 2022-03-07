!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use MPI

  use var
  use decomp_2d, only : nrank, xsize, real_type, decomp_2d_warning
  use param,   only : dt, zero, itr, itime, ioutput
  use transeq, only : calculate_transeq_rhs
  use navier,  only : solve_poisson, cor_vel
  use case
  use time_integrators, only : int_time
  use visu, only : write_snapshot
  
  implicit none

  double precision :: tstart, tend, telapsed, tmin, tmax
  !real :: trun
  integer :: i, j, k
  integer :: ndt, ndt_max
  integer :: code
  character(len=32) :: num

  call boot_xcompact3d()

  call init_xcompact3d(ndt_max)

  call case_init(ux1, uy1, uz1)

  telapsed = 0
  tmin = telapsed
  ndt = 1

  ioutput = 5
  print *, "=== WARNING ==="
  print *, " Writing with ioutput=",ioutput,", edit xcompact3d.f90 and recompile to change"
  print *, "=== WARNING ==="
  
  itime = 0
  do while(ndt < ndt_max)

     itr = 1 ! no inner iterations
     !call init_flowfield()

     tstart = MPI_Wtime()

     call case_bc(ux1, uy1, uz1)

     call calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)
     call int_time(ux1,uy1,uz1,dux1,duy1,duz1)

     !do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
     !  divu3(:,:,:) = zero
     !enddo
     call solve_poisson(pp3,px1,py1,pz1,ux1,uy1,uz1)
     call cor_vel(ux1,uy1,uz1,px1,py1,pz1)

     tend = MPI_Wtime()
     telapsed = telapsed + (tend - tstart)
     tmin = telapsed
     tmax = telapsed

     call MPI_Allreduce(telapsed, tmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, code)
     if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_Allreduce")
     call MPI_Allreduce(telapsed, tmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, code)
     if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_Allreduce")
     if (nrank == 0) then
        print *, "Elapse time min ", tmin, " max ", tmax
     end if

     !! End of timestep
     itime = itime + 1
     ndt = ndt + 1
     if (mod(itime,ioutput) == 0) then
        call write_snapshot(ux1, uy1, uz1, itime, num)
     endif
     call case_postprocess(ux1, uy1, uz1, ndt)

  end do

  if (nrank == 0) then
     print *, "Elapsed time (min-max) [s]: ", tmin, tmax
     print *, "Timesteps completed: ", ndt
     print *, "Compute rate (min-max)[dt/s]: ", ndt / tmin, ndt /  tmax
  end if

  call case_finalize()
  call finalise_xcompact3d(.true.)

end program xcompact3d
!########################################################################
!########################################################################

