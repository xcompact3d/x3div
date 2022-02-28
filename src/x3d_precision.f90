!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_precision

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

end module x3d_precision
