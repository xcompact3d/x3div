!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_filters

   use decomp_2d, only: mytype
   use x3d_operator_1d, only: x3doperator1d
   use param
   use thomas

   implicit none

   ! Make everything public unless declared private
   public

   ABSTRACT INTERFACE
      SUBROUTINE FILTER_X(t, u, ff, fs, fw, nx, ny, nz, npaire)
         use decomp_2d, only: mytype
         integer, intent(in) :: nx, ny, nz, npaire
         real(mytype), intent(out), dimension(nx, ny, nz) :: t
         real(mytype), intent(in), dimension(nx, ny, nz) :: u
         real(mytype), intent(in), dimension(nx):: ff, fs, fw
      END SUBROUTINE FILTER_X
      SUBROUTINE FILTER_Y(t, u, ff, fs, fw, nx, ny, nz, npaire)
         use decomp_2d, only: mytype
         integer, intent(in) :: nx, ny, nz, npaire
         real(mytype), intent(out), dimension(nx, ny, nz) :: t
         real(mytype), intent(in), dimension(nx, ny, nz) :: u
         real(mytype), intent(in), dimension(ny):: ff, fs, fw
      END SUBROUTINE FILTER_Y
      SUBROUTINE FILTER_Z(t, u, ff, fs, fw, nx, ny, nz, npaire)
         use decomp_2d, only: mytype
         integer, intent(in) :: nx, ny, nz, npaire
         real(mytype), intent(out), dimension(nx, ny, nz) :: t
         real(mytype), intent(in), dimension(nx, ny, nz) :: u
         real(mytype), intent(in), dimension(nz):: ff, fs, fw
      END SUBROUTINE FILTER_Z
   END INTERFACE

   PROCEDURE(FILTER_X) :: filx_00, filx_11, filx_12, filx_21, filx_22
   PROCEDURE(FILTER_Y) :: fily_00, fily_11, fily_12, fily_21, fily_22
   PROCEDURE(FILTER_Z) :: filz_00, filz_11, filz_12, filz_21, filz_22
   PROCEDURE(FILTER_X), POINTER :: filx => null(), filxS => null()
   PROCEDURE(FILTER_Y), POINTER :: fily => null(), filyS => null()
   PROCEDURE(FILTER_Z), POINTER :: filz => null(), filzS => null()

contains

end module x3d_filters
