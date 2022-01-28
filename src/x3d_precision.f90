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
!##################################################################
! This is the main 2D pencil decomposition module

module x3dprecision

  use MPI
  use, intrinsic :: iso_fortran_env, only : real32, real64
  use decomp_2d, only : mytype
  
  implicit none

#ifdef DOUBLE_PREC
  real(mytype),parameter, public :: pi=dacos(-1._real64)
  real(mytype),parameter, public :: twopi=2._real64*dacos(-1._real64)
#else
  real(mytype),parameter, public :: pi=acos(-1._real32)
  real(mytype),parameter, public :: twopi=2._real32*acos(-1._real32)
#endif

  ! Make everything private unless declared public
  private

  public ::  sin_prec,  cos_prec,  tan_prec, &
            asin_prec, acos_prec, atan_prec, &
            sinh_prec, cosh_prec, tanh_prec, &
             exp_prec,  log_prec,log10_prec, &
            sqrt_prec,  abs_prec
contains

  !##################################################################
  !********************************************************************
  ! Math functions for Single/double precision
  !-------------------------------------------
  elemental function sin_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dsin(x)
#else
    y = sin(x)
#endif
  end function sin_prec
  !-------------------------------------------
  elemental function cos_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dcos(x)
#else
    y = cos(x)
#endif
  end function cos_prec
  !-------------------------------------------
  elemental function tan_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dtan(x)
#else
    y = tan(x)
#endif
  end function tan_prec
  !-------------------------------------------
  elemental function asin_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dasin(x)
#else
    y = asin(x)
#endif
  end function asin_prec
  !-------------------------------------------
  elemental function acos_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dacos(x)
#else
    y = acos(x)
#endif
  end function acos_prec
  !-------------------------------------------
  elemental function atan_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = datan(x)
#else
    y = atan(x)
#endif
  end function atan_prec
  !-------------------------------------------
  elemental function sinh_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dsinh(x)
#else
    y = sinh(x)
#endif
  end function sinh_prec
  !-------------------------------------------
  elemental function cosh_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dcosh(x)
#else
    y = cosh(x)
#endif
  end function cosh_prec
  !-------------------------------------------
  elemental function tanh_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dtanh(x)
#else
    y = tanh(x)
#endif
  end function tanh_prec
  !-------------------------------------------
  elemental function exp_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dexp(x)
#else
    y = exp(x)
#endif
  end function exp_prec
  !-------------------------------------------
  elemental function log_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dlog(x)
#else
    y = alog(x)
#endif
  end function log_prec
  !-------------------------------------------
  elemental function log10_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dlog10(x)
#else
    y = alog10(x)
#endif
  end function log10_prec
  !-------------------------------------------
  elemental function sqrt_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dsqrt(x)
#else
    y = sqrt(x)
#endif
  end function sqrt_prec
  !-------------------------------------------
  elemental function abs_prec(x) result(y)
    USE decomp_2d, only : mytype
    real(mytype), intent(in) :: x
    real(mytype) :: y
#ifdef DOUBLE_PREC
    y = dabs(x)
#else
    y = abs(x)
#endif
  end function abs_prec

  !********************************************************************


end module x3dprecision
