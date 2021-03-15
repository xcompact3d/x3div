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
  subroutine calculate_transeq_rhs(drho1,dux1,duy1,duz1,dphi1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    use decomp_2d, only : mytype, xsize, zsize
    use variables, only : numscalar
    use param, only : ntime, ilmn, nrhotime, ilmn_solve_temp

    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), nrhotime), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), numscalar), intent(in) :: phi1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ep1
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: divu3

    !! Outputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: drho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime, numscalar) :: dphi1

    !! Momentum equations
    call momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

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
  subroutine momentum_rhs_eq(dux1,duy1,duz1,rho1,ux1,uy1,uz1,ep1,phi1,divu3)

    use param
    use variables
    use decomp_2d
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,mu1,mu2,mu3
    use var, only : rho2,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
    use var, only : rho3,ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
    use var, only : sgsx1,sgsy1,sgsz1

    use case, only : momentum_forcing

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),intent(in),dimension(zsize(1),zsize(2),zsize(3)) :: divu3

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1


    integer :: i,j,k,is

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    if (ilmn) then
      ta1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * ux1(:,:,:)
      tb1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uy1(:,:,:)
      tc1(:,:,:) = rho1(:,:,:,1) * ux1(:,:,:) * uz1(:,:,:)
    else
      ta1(:,:,:) = ux1(:,:,:) * ux1(:,:,:)
      tb1(:,:,:) = ux1(:,:,:) * uy1(:,:,:)
      tc1(:,:,:) = ux1(:,:,:) * uz1(:,:,:)
    endif

    call derx (td1,ta1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (te1,tb1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tf1,tc1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
    call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
    call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)

    ! Convective terms of x-pencil are stored in tg1,th1,ti1
    if (ilmn) then
      tg1(:,:,:) = td1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * ta1(:,:,:)
      th1(:,:,:) = te1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tb1(:,:,:)
      ti1(:,:,:) = tf1(:,:,:) + rho1(:,:,:,1) * ux1(:,:,:) * tc1(:,:,:)
    else
      tg1(:,:,:) = td1(:,:,:) + ux1(:,:,:) * ta1(:,:,:)
      th1(:,:,:) = te1(:,:,:) + ux1(:,:,:) * tb1(:,:,:)
      ti1(:,:,:) = tf1(:,:,:) + ux1(:,:,:) * tc1(:,:,:)
    endif
    ! TODO: save the x-convective terms already in dux1, duy1, duz1

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derx (td1,rho1(:,:,:,1),di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
       tg1(:,:,:) = tg1(:,:,:) + ux1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       th1(:,:,:) = th1(:,:,:) + uy1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
       ti1(:,:,:) = ti1(:,:,:) + uz1(:,:,:) * ux1(:,:,:) * td1(:,:,:)
    endif

    call transpose_x_to_y(ux1,ux2)
    call transpose_x_to_y(uy1,uy2)
    call transpose_x_to_y(uz1,uz2)

    if (ilmn) then
       call transpose_x_to_y(rho1(:,:,:,1),rho2)
       call transpose_x_to_y(mu1,mu2)
    else
       rho2(:,:,:) = one
    endif

    !WORK Y-PENCILS
    if (ilmn) then
      td2(:,:,:) = rho2(:,:,:) * ux2(:,:,:) * uy2(:,:,:)
      te2(:,:,:) = rho2(:,:,:) * uy2(:,:,:) * uy2(:,:,:)
      tf2(:,:,:) = rho2(:,:,:) * uz2(:,:,:) * uy2(:,:,:)
    else
      td2(:,:,:) = ux2(:,:,:) * uy2(:,:,:)
      te2(:,:,:) = uy2(:,:,:) * uy2(:,:,:)
      tf2(:,:,:) = uz2(:,:,:) * uy2(:,:,:)
    endif

    call dery (tg2,td2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (th2,te2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (ti2,tf2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (td2,ux2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
    call dery (te2,uy2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
    call dery (tf2,uz2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)

    ! Convective terms of y-pencil in tg2,th2,ti2
    if (ilmn) then
      tg2(:,:,:) = tg2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * td2(:,:,:)
      th2(:,:,:) = th2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
      ti2(:,:,:) = ti2(:,:,:) + rho2(:,:,:) * uy2(:,:,:) * tf2(:,:,:)
    else
      tg2(:,:,:) = tg2(:,:,:) + uy2(:,:,:) * td2(:,:,:)
      th2(:,:,:) = th2(:,:,:) + uy2(:,:,:) * te2(:,:,:)
      ti2(:,:,:) = ti2(:,:,:) + uy2(:,:,:) * tf2(:,:,:)
    endif

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call dery (te2,rho2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
       tg2(:,:,:) = tg2(:,:,:) + ux2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       th2(:,:,:) = th2(:,:,:) + uy2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
       ti2(:,:,:) = ti2(:,:,:) + uz2(:,:,:) * uy2(:,:,:) * te2(:,:,:)
    endif

    call transpose_y_to_z(ux2,ux3)
    call transpose_y_to_z(uy2,uy3)
    call transpose_y_to_z(uz2,uz3)

    !WORK Z-PENCILS
    if (ilmn) then
       call transpose_y_to_z(rho2,rho3)
       call transpose_y_to_z(mu2,mu3)

       td3(:,:,:) = rho3(:,:,:) * ux3(:,:,:) * uz3(:,:,:)
       te3(:,:,:) = rho3(:,:,:) * uy3(:,:,:) * uz3(:,:,:)
       tf3(:,:,:) = rho3(:,:,:) * uz3(:,:,:) * uz3(:,:,:)
    else
       td3(:,:,:) = ux3(:,:,:) * uz3(:,:,:)
       te3(:,:,:) = uy3(:,:,:) * uz3(:,:,:)
       tf3(:,:,:) = uz3(:,:,:) * uz3(:,:,:)
    endif

    call derz (tg3,td3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (th3,te3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
    call derz (ti3,tf3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (td3,ux3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (te3,uy3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
    call derz (tf3,uz3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)

    ! Convective terms of z-pencil in ta3,tb3,tc3
    if (ilmn) then
      ta3(:,:,:) = tg3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * td3(:,:,:)
      tb3(:,:,:) = th3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * te3(:,:,:)
      tc3(:,:,:) = ti3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
    else
      ta3(:,:,:) = tg3(:,:,:) + uz3(:,:,:) * td3(:,:,:)
      tb3(:,:,:) = th3(:,:,:) + uz3(:,:,:) * te3(:,:,:)
      tc3(:,:,:) = ti3(:,:,:) + uz3(:,:,:) * tf3(:,:,:)
    endif

    if (ilmn) then
       !! Quasi-skew symmetric terms
       call derz (tf3,rho3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
       ta3(:,:,:) = ta3(:,:,:) + ux3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + uy3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + uz3(:,:,:) * uz3(:,:,:) * tf3(:,:,:)

       !! Add the additional divu terms
       ta3(:,:,:) = ta3(:,:,:) + rho3(:,:,:) * ux3(:,:,:) * divu3(:,:,:)
       tb3(:,:,:) = tb3(:,:,:) + rho3(:,:,:) * uy3(:,:,:) * divu3(:,:,:)
       tc3(:,:,:) = tc3(:,:,:) + rho3(:,:,:) * uz3(:,:,:) * divu3(:,:,:)
    endif

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
