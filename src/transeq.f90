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

    use x3dprecision, only : mytype
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
    use x3dprecision, only : mytype
    use decomp_2d , only : xsize, ysize, zsize
    use decomp_2d, only : transpose_x_to_y, &
                          transpose_y_to_z, &
                          transpose_z_to_y, &
                          transpose_y_to_x
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3

    use mom, only : test_du, test_dv, test_dw
    
    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1


    integer :: i,j,k,is

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    ta1(:,:,:) = ux1(:,:,:) * ux1(:,:,:)
    tb1(:,:,:) = ux1(:,:,:) * uy1(:,:,:)
    tc1(:,:,:) = ux1(:,:,:) * uz1(:,:,:)

    call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

    ! Convective terms of x-pencil are stored in tg1,th1,ti1
    tg1(:,:,:) = td1(:,:,:) + ux1(:,:,:) * ta1(:,:,:)
    th1(:,:,:) = te1(:,:,:) + ux1(:,:,:) * tb1(:,:,:)
    ti1(:,:,:) = tf1(:,:,:) + ux1(:,:,:) * tc1(:,:,:)
    ! TODO: save the x-convective terms already in dux1, duy1, duz1

    call test_du(ta1)
    
    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)

    !WORK Y-PENCILS
    td2(:,:,:) = ux2(:,:,:) * uy2(:,:,:)
    te2(:,:,:) = uy2(:,:,:) * uy2(:,:,:)
    tf2(:,:,:) = uz2(:,:,:) * uy2(:,:,:)

    call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

    ! Convective terms of y-pencil in tg2,th2,ti2
    tg2(:,:,:) = tg2(:,:,:) + uy2(:,:,:) * td2(:,:,:)
    th2(:,:,:) = th2(:,:,:) + uy2(:,:,:) * te2(:,:,:)
    ti2(:,:,:) = ti2(:,:,:) + uy2(:,:,:) * tf2(:,:,:)

    call test_dv(te2)
    
    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)

    !WORK Z-PENCILS
    td3(:,:,:) = ux3(:,:,:) * uz3(:,:,:)
    te3(:,:,:) = uy3(:,:,:) * uz3(:,:,:)
    tf3(:,:,:) = uz3(:,:,:) * uz3(:,:,:)

    call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

    ! Convective terms of z-pencil in ta3,tb3,tc3
    ta3(:,:,:) = tg3(:,:,:) + uz3(:,:,:) * td3(:,:,:)
    tb3(:,:,:) = th3(:,:,:) + uz3(:,:,:) * te3(:,:,:)
    tc3(:,:,:) = ti3(:,:,:) + uz3(:,:,:) * tf3(:,:,:)

    call test_dw(tf3)
    
    !WORK Y-PENCILS
    call transpose_z_to_y(ta3,td2)
    call transpose_z_to_y(tb3,te2)
    call transpose_z_to_y(tc3,tf2)

    !WORK X-PENCILS
    call transpose_y_to_x(td2,ta1)
    call transpose_y_to_x(te2,tb1)
    call transpose_y_to_x(tf2,tc1) !diff+conv. terms

    !FINAL SUM: DIFF TERMS + CONV TERMS
    dux1(:,:,:,1) = ta1(:,:,:)
    duy1(:,:,:,1) = tb1(:,:,:)
    duz1(:,:,:,1) = tc1(:,:,:)

  end subroutine momentum_rhs_eq
  !############################################################################
  !############################################################################
end module transeq
