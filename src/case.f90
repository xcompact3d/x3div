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

module case

  use mom
  use param
  use decomp_2d
  use variables

  use var, only : nzmsize

  implicit none

  private ! All functions/subroutines private by default
  public :: init, boundary_conditions, &
            momentum_forcing, scalar_forcing, set_fluid_properties

contains
  !##################################################################
  subroutine init (rho1, ux1, uy1, uz1, ep1, phi1, drho1, dux1, duy1, duz1, dphi1, &
       pp3, px1, py1, pz1)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1,drho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime,numscalar) :: dphi1
    real(mytype), dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    INTEGER :: it, is

    !! Zero out the pressure field
    pp3(:,:,:,1) = zero
    px1(:,:,:) = zero
    py1(:,:,:) = zero
    pz1(:,:,:) = zero

    !! Default density and pressure0 to one
    pressure0 = one
    rho1(:,:,:,:) = one

    call vel(ux1, uy1, uz1)
    
    !! Setup old arrays
    do it = 1, ntime
       drho1(:,:,:,it) = rho1(:,:,:,1)
       dux1(:,:,:,it)=ux1(:,:,:)
       duy1(:,:,:,it)=uy1(:,:,:)
       duz1(:,:,:,it)=uz1(:,:,:)
    enddo

    do it = 2, nrhotime
       rho1(:,:,:,it) = rho1(:,:,:,1)
    enddo

    do is = 1, numscalar
       do it = 1, ntime
          dphi1(:,:,:,it,is) = phi1(:,:,:,is)
       enddo
    enddo

  end subroutine init
  !##################################################################
  !##################################################################
  subroutine boundary_conditions (rho,ux,uy,uz,phi,ep)

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz,ep
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),numscalar) :: phi
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),nrhotime) :: rho

  end subroutine boundary_conditions
  !##################################################################
  !##################################################################
  !!
  !!  SUBROUTINE: momentum_forcing
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Calls case-specific forcing functions for the
  !!              momentum equations.
  !!
  !##################################################################
  subroutine momentum_forcing(dux1, duy1, duz1, rho1, ux1, uy1, uz1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

  end subroutine momentum_forcing
  !##################################################################
  !##################################################################
  !!
  !!  SUBROUTINE: scalar_forcing
  !!      AUTHOR: Kay Schäfer
  !! DESCRIPTION: Calls case-specific forcing functions for the
  !!              scalar transport equations.
  !!
  !##################################################################
  subroutine scalar_forcing(dphi1, rho1, ux1, uy1, uz1, phi1)

    implicit none

    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1, phi1
    real(mytype), intent(in), dimension(xsize(1), xsize(2), xsize(3), nrhotime) :: rho1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dphi1


  end subroutine scalar_forcing
  !##################################################################
  !##################################################################
  subroutine set_fluid_properties(rho1, mu1)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: rho1
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: mu1

  endsubroutine set_fluid_properties
  !##################################################################
  !##################################################################
end module case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
