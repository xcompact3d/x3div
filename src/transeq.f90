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
module transeq

  private
  public :: calculate_transeq_rhs

contains
  !############################################################################
  !!  SUBROUTINE: calculate_transeq_rhs
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Calculates the right hand sides of all transport
  !!              equations - momentum, scalar transport, etc.
  !############################################################################
  subroutine calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)

    use decomp_2d, only : mytype
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
  !!      AUTHOR: ?
  !!    MODIFIED: Kay Schäfer
  !! DESCRIPTION: Calculation of convective and diffusion terms of momentum
  !!              equation
  !!
  !############################################################################
  !############################################################################
  subroutine momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1)

    use param
    use variables
    use x3d_operator_1d
    use decomp_2d, only : mytype
    use x3d_transpose
    use x3d_derive
    use decomp_2d , only : xsize, ysize, zsize
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3
    use nvtx
    use mom, only : test_du, test_dv, test_dw
    
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

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    
    call nvtxStartRange("Trans collapse 256")
    !!!!$acc kernels vector_length(256)
    !do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
    !$acc parallel loop gang vector collapse(3)
    do k=1,xsz3
    do j=1,xsz2
    do i=1,xsz3
      ta1(i,j,k) = ux1(i,j,k) * ux1(i,j,k)
      tb1(i,j,k) = ux1(i,j,k) * uy1(i,j,k)
      tc1(i,j,k) = ux1(i,j,k) * uz1(i,j,k)
    enddo
    enddo
    enddo
    !!!!$acc end kernel
    call nvtxEndRange
    
    call nvtxStartRange("Trans Do concurr standard")
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      ta1(i,j,k) = ux1(i,j,k) * ux1(i,j,k)
      tb1(i,j,k) = ux1(i,j,k) * uy1(i,j,k)
      tc1(i,j,k) = ux1(i,j,k) * uz1(i,j,k)
    enddo
    call nvtxEndRange

    call derx (td1,ta1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (te1,tb1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tf1,tc1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (ta1,ux1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tb1,uy1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (tc1,uz1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))

    ! Convective terms of x-pencil are stored in tg1,th1,ti1
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      tg1(i,j,k) = td1(i,j,k) + ux1(i,j,k) * ta1(i,j,k)
      th1(i,j,k) = te1(i,j,k) + ux1(i,j,k) * tb1(i,j,k)
      ti1(i,j,k) = tf1(i,j,k) + ux1(i,j,k) * tc1(i,j,k)
    enddo
    ! TODO: save the x-convective terms already in dux1, duy1, duz1

    call test_du(ta1)
    
    call x3d_transpose_x_to_y(ux1,ux2)
    call x3d_transpose_x_to_y(uy1,uy2)
    call x3d_transpose_x_to_y(uz1,uz2)

    !WORK Y-PENCILS
    
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      td2(i,j,k) = ux2(i,j,k) * uy2(i,j,k)
      te2(i,j,k) = uy2(i,j,k) * uy2(i,j,k)
      tf2(i,j,k) = uz2(i,j,k) * uy2(i,j,k)
    enddo

    call dery (tg2,td2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (th2,te2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call dery (ti2,tf2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (td2,ux2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call dery (te2,uy2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (tf2,uz2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))

    ! Convective terms of y-pencil in tg2,th2,ti2
    do concurrent (k=1:ysz3, j=1:ysz2, i=1:ysz1)
      tg2(i,j,k) = tg2(i,j,k) + uy2(i,j,k) * td2(i,j,k)
      th2(i,j,k) = th2(i,j,k) + uy2(i,j,k) * te2(i,j,k)
      ti2(i,j,k) = ti2(i,j,k) + uy2(i,j,k) * tf2(i,j,k)
    enddo
    
    call test_dv(te2)
    
    call x3d_transpose_y_to_z(ux2,ux3)
    call x3d_transpose_y_to_z(uy2,uy3)
    call x3d_transpose_y_to_z(uz2,uz3)

    !WORK Z-PENCILS
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      td3(i,j,k) = ux3(i,j,k) * uz3(i,j,k)
      te3(i,j,k) = uy3(i,j,k) * uz3(i,j,k)
      tf3(i,j,k) = uz3(i,j,k) * uz3(i,j,k)
    enddo

    call derz (tg3,td3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (th3,te3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (ti3,tf3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (td3,ux3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (te3,uy3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (tf3,uz3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))

    ! Convective terms of z-pencil in ta3,tb3,tc3
    do concurrent (k=1:zsz3, j=1:zsz2, i=1:zsz1)
      ta3(i,j,k) = tg3(i,j,k) + uz3(i,j,k) * td3(i,j,k)
      tb3(i,j,k) = th3(i,j,k) + uz3(i,j,k) * te3(i,j,k)
      tc3(i,j,k) = ti3(i,j,k) + uz3(i,j,k) * tf3(i,j,k)
    enddo 

    call test_dw(tf3)
    
    !WORK Y-PENCILS
    call x3d_transpose_z_to_y(ta3,td2)
    call x3d_transpose_z_to_y(tb3,te2)
    call x3d_transpose_z_to_y(tc3,tf2)

    !WORK X-PENCILS
    call x3d_transpose_y_to_x(td2,ta1)
    call x3d_transpose_y_to_x(te2,tb1)
    call x3d_transpose_y_to_x(tf2,tc1) !diff+conv. terms

    !FINAL SUM: DIFF TERMS + CONV TERMS
    do concurrent (k=1:xsz3, j=1:xsz2, i=1:xsz1)
      dux1(i,j,k,1) = ta1(i,j,k)
      duy1(i,j,k,1) = tb1(i,j,k)
      duz1(i,j,k,1) = tc1(i,j,k)
    enddo

  end subroutine momentum_rhs_eq
  !############################################################################
  !############################################################################
end module transeq
