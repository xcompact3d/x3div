!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module transeq

  use decomp_2d_constants, only : mytype
  
  private
  public :: calculate_transeq_rhs

contains
  !############################################################################
  !!  SUBROUTINE: calculate_transeq_rhs
  !! DESCRIPTION: Calculates the right hand sides of all transport
  !!              equations - momentum, scalar transport, etc.
  !############################################################################
  subroutine calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)

    use decomp_2d, only : xsize, zsize
    use param, only : ntime

    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1

    !! Outputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    !! Momentum equations
    call momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1)

  end subroutine calculate_transeq_rhs
  !############################################################################
  !############################################################################
  !!
  !!  subroutine: momentum_rhs_eq
  !! DESCRIPTION: Calculation of convective and diffusion terms of momentum
  !!              equation
  !!
  !############################################################################
  !############################################################################
  subroutine momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1)

    use param
    use variables
    use x3d_operator_1d
    use x3d_transpose
    use x3d_derive
    use decomp_2d , only : xsize, ysize, zsize
    use var, only : ta1,tb1,tc1,td1,te1,tf1
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3
    use mom, only : test_du, test_dv, test_dw
    !use nvtx 

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1


    integer :: i,j,k,is
    integer :: xsz1, xsz2, xsz3
    integer :: ysz1, ysz2, ysz3
    integer :: zsz1, zsz2, zsz3
    xsz1=xsize(1)
    xsz2=xsize(2)
    xsz3=xsize(3)
    ysz1=ysize(1)
    ysz2=ysize(2)
    ysz3=ysize(3)
    zsz1=zsize(1)
    zsz2=zsize(2)
    zsz3=zsize(3)

    !$acc data create(ux2,uy2,uz2,ux3,uy3,uz3) present(ux1,uy1,uz1,dux1,duy1,duz1)
    !$acc enter data create(ta1,td1,ta2,td2,tg2,ta3,td3,tg3)
    ! Transpose the velocity in all 3 system of refence
    !call nvtxStartRange("Tran XY")
    call x3d_transpose_x_to_y(ux1,ux2)
    call x3d_transpose_x_to_y(uy1,uy2)
    call x3d_transpose_x_to_y(uz1,uz2)
    !call nvtxEndRange
    
    !call nvtxStartRange("Tran YZ")
    call x3d_transpose_y_to_z(ux2,ux3)
    call x3d_transpose_y_to_z(uy2,uy3)
    call x3d_transpose_y_to_z(uz2,uz3)
    !call nvtxEndRange

    ! Compute dux
    !SKEW SYMMETRIC FORM
    !WORK X-PENCIL
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      ta1(i,j,k) = ux1(i,j,k) * ux1(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der X")
    call derx (td1,ta1,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (ta1,ux1,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    !call nvtxEndRange
    call test_du(ta1)

    ! Convective terms of x-pencil are stored directly in dux-y-z1
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      dux1(i,j,k,1) = td1(i,j,k) + ux1(i,j,k) * ta1(i,j,k) !duudx+u*dudx
    enddo
    !$acc end kernels

    !WORK Y-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      td2(i,j,k) = ux2(i,j,k) * uy2(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Y")
    call dery (tg2,td2,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (td2,ux2,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    !call nvtxEndRange

    ! Convective terms of y-pencil in tg2,th2,ti2
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      tg2(i,j,k) = tg2(i,j,k) + uy2(i,j,k) * td2(i,j,k)
    enddo
    !$acc end kernels
    
    !WORK Z-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      td3(i,j,k) = ux3(i,j,k) * uz3(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Z")
    call derz (tg3,td3,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (td3,ux3,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    !call nvtxEndRange

    ! Convective terms of z-pencil in ta3,tb3,tc3
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      ta3(i,j,k) = tg3(i,j,k) + uz3(i,j,k) * td3(i,j,k)
    enddo 
    !$acc end kernels

    ! Send z-rhs (ta3,tb3,tc3) to y-pencil (td2,te2,tf2)
    !call nvtxStartRange("Tran ZY")
    call x3d_transpose_z_to_y(ta3,td2)
    !call nvtxEndRange

    ! Combine y-rhs with z-rhs
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      tg2(i,j,k) = tg2(i,j,k) + td2(i,j,k)
    enddo
    !$acc end kernels

    !WORK X-PENCILS
    !call nvtxStartRange("Tran YX")
    call x3d_transpose_y_to_x(tg2,ta1)
    !call nvtxEndRange

    !FINAL SUM: DIFF TERMS + CONV TERMS
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      dux1(i,j,k,1) = dux1(i,j,k,1) + ta1(i,j,k)
    enddo
    !$acc end kernels
    !$acc exit data delete(ta1,td1,ta2,td2,tg2,ta3,td3,tg3)

    ! Compute duy
    !$acc enter data create(tb1,te1,tb2,te2,th2,tb3,te3,th3)
    !SKEW SYMMETRIC FORM
    !WORK X-PENCIL
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      tb1(i,j,k) = ux1(i,j,k) * uy1(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der X")
    call derx (te1,tb1,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tb1,uy1,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    !call nvtxEndRange

    ! Convective terms of x-pencil are stored directly in dux-y-z1
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      duy1(i,j,k,1) = te1(i,j,k) + ux1(i,j,k) * tb1(i,j,k)
    enddo
    !$acc end kernels

    !WORK Y-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      te2(i,j,k) = uy2(i,j,k) * uy2(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Y")
    call dery (th2,te2,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call dery (te2,uy2,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    !call nvtxEndRange
    call test_dv(te2)

    ! Convective terms of y-pencil in tg2,th2,ti2
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      th2(i,j,k) = th2(i,j,k) + uy2(i,j,k) * te2(i,j,k)
    enddo
    !$acc end kernels
    
    !WORK Z-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      te3(i,j,k) = uy3(i,j,k) * uz3(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Z")
    call derz (th3,te3,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (te3,uy3,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    !call nvtxEndRange

    ! Convective terms of z-pencil in ta3,tb3,tc3
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      tb3(i,j,k) = th3(i,j,k) + uz3(i,j,k) * te3(i,j,k)
    enddo 
    !$acc end kernels

    ! Send z-rhs (ta3,tb3,tc3) to y-pencil (td2,te2,tf2)
    !call nvtxStartRange("Tran ZY")
    call x3d_transpose_z_to_y(tb3,te2)
    !call nvtxEndRange

    ! Combine y-rhs with z-rhs
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      th2(i,j,k) = th2(i,j,k) + te2(i,j,k)
    enddo
    !$acc end kernels

    !WORK X-PENCILS
    !call nvtxStartRange("Tran YX")
    call x3d_transpose_y_to_x(th2,tb1)
    !call nvtxEndRange

    !FINAL SUM: DIFF TERMS + CONV TERMS
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      duy1(i,j,k,1) = duy1(i,j,k,1) + tb1(i,j,k)
    enddo
    !$acc end kernels
    !$acc exit data delete(tb1,te1,tb2,te2,th2,tb3,te3,th3)

    ! Compute duz
    !$acc enter data create(tc1,tf1,tc2,tf2,ti2,tc3,tf3,ti3)
    !SKEW SYMMETRIC FORM
    !WORK X-PENCIL
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      tc1(i,j,k) = ux1(i,j,k) * uz1(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der X")
    call derx (tf1,tc1,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tc1,uz1,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    !call nvtxEndRange

    ! Convective terms of x-pencil are stored directly in dux-y-z1
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      duz1(i,j,k,1) = tf1(i,j,k) + ux1(i,j,k) * tc1(i,j,k)
    enddo
    !$acc end kernels

    !WORK Y-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      tf2(i,j,k) = uz2(i,j,k) * uy2(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Y")
    call dery (ti2,tf2,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (tf2,uz2,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    !call nvtxEndRange

    ! Convective terms of y-pencil in tg2,th2,ti2
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      ti2(i,j,k) = ti2(i,j,k) + uy2(i,j,k) * tf2(i,j,k)
    enddo
    !$acc end kernels
    
    !WORK Z-PENCILS
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      tf3(i,j,k) = uz3(i,j,k) * uz3(i,j,k)
    enddo
    !$acc end kernels

    !call nvtxStartRange("Group Der Z")
    call derz (ti3,tf3,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (tf3,uz3,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    !call nvtxEndRange
    call test_dw(tf3)

    ! Convective terms of z-pencil in ta3,tb3,tc3
    !$acc kernels default(present)
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      tc3(i,j,k) = ti3(i,j,k) + uz3(i,j,k) * tf3(i,j,k)
    enddo 
    !$acc end kernels

    ! Send z-rhs (ta3,tb3,tc3) to y-pencil (td2,te2,tf2)
    !call nvtxStartRange("Tran ZY")
    call x3d_transpose_z_to_y(tc3,tf2)
    !call nvtxEndRange

    ! Combine y-rhs with z-rhs
    !$acc kernels default(present)
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      ti2(i,j,k) = ti2(i,j,k) + tf2(i,j,k)
    enddo
    !$acc end kernels

    !WORK X-PENCILS
    !call nvtxStartRange("Tran YX")
    call x3d_transpose_y_to_x(ti2,tc1) !diff+conv. terms
    !call nvtxEndRange

    !FINAL SUM: DIFF TERMS + CONV TERMS
    !$acc kernels default(present)
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      duz1(i,j,k,1) = duz1(i,j,k,1) + tc1(i,j,k)
    enddo
    !$acc end kernels
    !$acc exit data delete(tc1,tf1,tc2,tf2,ti2,tc3,tf3,ti3)

    !$acc end data

  end subroutine momentum_rhs_eq
  !############################################################################
  !############################################################################
end module transeq
