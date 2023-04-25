!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

program xcompact3d

  use MPI

  use var
  use decomp_2d_constants, only : real_type 
  use decomp_2d_mpi, only : nrank, decomp_2d_warning 
  use decomp_2d, only : xsize 
  use param,   only : dt, zero, itr
  use param,   only : gdt
  use transeq, only : calculate_transeq_rhs
  use navier,  only : solve_poisson, cor_vel
  use decomp_2d_poisson, only : ax, bx, ay, by, az, bz, kxyz
  use mom,     only : test_du, test_dv, test_dw
  use time_integrators, only : int_time
  !use nvtx
  use x3d_operator_1d
  use variables

  implicit none

  double precision :: tstart, tend, telapsed, tmin, tmax
  !real :: trun
  integer :: i, j, k
  integer :: ndt, ndt_max
  integer :: code

  call boot_xcompact3d()

  !call nvtxStartRange("init_xcompact3d")
  call init_xcompact3d(ndt_max)
  !call nvtxEndRange

  telapsed = 0
  tmin = telapsed

  ndt = 1

  !$acc data copyin(gdt) async 
  !$acc data copyin(x3d_op_derx%f,x3d_op_derx%s,x3d_op_derx%w,x3d_op_derx%periodic) async 
  !$acc data copyin(x3d_op_dery%f,x3d_op_dery%s,x3d_op_dery%w,x3d_op_dery%periodic) async
  !$acc data copyin(x3d_op_derz%f,x3d_op_derz%s,x3d_op_derz%w,x3d_op_derz%periodic) async
  !
  !$acc data copyin(x3d_op_derxp%f,x3d_op_derxp%s,x3d_op_derxp%w,x3d_op_derxp%periodic) async
  !$acc data copyin(x3d_op_deryp%f,x3d_op_deryp%s,x3d_op_deryp%w,x3d_op_deryp%periodic) async
  !$acc data copyin(x3d_op_derzp%f,x3d_op_derzp%s,x3d_op_derzp%w,x3d_op_derzp%periodic) async
  !
  !$acc data copyin(x3d_op_derxvp%f,x3d_op_derxvp%s,x3d_op_derxvp%w,x3d_op_derxvp%periodic) async 
  !$acc data copyin(x3d_op_deryvp%f,x3d_op_deryvp%s,x3d_op_deryvp%w,x3d_op_deryvp%periodic) async 
  !$acc data copyin(x3d_op_derzvp%f,x3d_op_derzvp%s,x3d_op_derzvp%w,x3d_op_derzvp%periodic) async 
  !!
  !$acc data copyin(x3d_op_derxpv%f,x3d_op_derxpv%s,x3d_op_derxpv%w,x3d_op_derxpv%periodic) async 
  !$acc data copyin(x3d_op_derypv%f,x3d_op_derypv%s,x3d_op_derypv%w,x3d_op_derypv%periodic) async 
  !$acc data copyin(x3d_op_derzpv%f,x3d_op_derzpv%s,x3d_op_derzpv%w,x3d_op_derzpv%periodic) async 
  !
  !$acc data copyin(x3d_op_intxvp%f,x3d_op_intxvp%s,x3d_op_intxvp%w,x3d_op_intxvp%periodic) async 
  !$acc data copyin(x3d_op_intyvp%f,x3d_op_intyvp%s,x3d_op_intyvp%w,x3d_op_intyvp%periodic) async 
  !$acc data copyin(x3d_op_intzvp%f,x3d_op_intzvp%s,x3d_op_intzvp%w,x3d_op_intzvp%periodic) async 
  !
  !$acc data copyin(x3d_op_intxpv%f,x3d_op_intxpv%s,x3d_op_intxpv%w,x3d_op_intxpv%periodic) async 
  !$acc data copyin(x3d_op_intypv%f,x3d_op_intypv%s,x3d_op_intypv%w,x3d_op_intypv%periodic) async 
  !$acc data copyin(x3d_op_intzpv%f,x3d_op_intzpv%s,x3d_op_intzpv%w,x3d_op_intzpv%periodic) async 
  !$acc data copyin(kxyz,az,bz,ay,by,ax,bx) async 
  !
  !$acc data copy(ux1,uy1,uz1) async 
  !$acc data copy(dux1,duy1,duz1) async
  !$acc data copy(pp3,px1,py1,pz1) async
  !$acc wait
  call init_flowfield()
  do while(ndt <= ndt_max)
     itr = 1 ! no inner iterations
     !call init_flowfield()

     tstart = MPI_Wtime()
     !call nvtxStartRange("calculate_transeq_rhs")
     call calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)
     !call nvtxEndRange

     !call nvtxStartRange("int_time")
     call int_time(ux1,uy1,uz1,dux1,duy1,duz1)
     !call nvtxEndRange
     !
     !!do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
     !!  divu3(:,:,:) = zero
     !!enddo
     !call nvtxStartRange("solve_poisson")
     call solve_poisson(pp3,px1,py1,pz1,ux1,uy1,uz1)
     !call nvtxEndRange
     !
     !call nvtxStartRange("cor_vel")
     call cor_vel(ux1,uy1,uz1,px1,py1,pz1)
     !call nvtxEndRange

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

     ndt = ndt + 1
  end do
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data
  !$acc end data

  if (nrank == 0) then
     print *, "Elapsed time (min-max) [s]: ", tmin, tmax
     print *, "Timesteps completed: ", ndt
     print *, "Compute rate (min-max)[dt/s]: ", ndt / tmin, ndt /  tmax
  end if

  call finalise_xcompact3d(.true.)

end program xcompact3d
!########################################################################
!########################################################################

