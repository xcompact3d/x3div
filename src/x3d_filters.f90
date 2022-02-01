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

module x3d_filters

  use decomp_2d, only : mytype
  use x3d_operator_1d, only : x3doperator1d
  use param
  use thomas
  
  implicit none

  ! Make everything public unless declared private
  public

  ABSTRACT INTERFACE
     SUBROUTINE FILTER_X(t,u,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer, intent(in) :: nx,ny,nz,npaire
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(ny,nz):: s
       real(mytype), intent(in), dimension(nx):: ff,fs,fw
     END SUBROUTINE FILTER_X
     SUBROUTINE FILTER_Y(t,u,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer, intent(in) :: nx,ny,nz,npaire
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(nx,nz):: s
       real(mytype), intent(in), dimension(ny):: ff,fs,fw
     END SUBROUTINE FILTER_Y
     SUBROUTINE FILTER_Z(t,u,s,ff,fs,fw,nx,ny,nz,npaire)
       use decomp_2d, only : mytype
       integer, intent(in) :: nx,ny,nz,npaire
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(nx,ny):: s
       real(mytype), intent(in), dimension(nz):: ff,fs,fw
     END SUBROUTINE FILTER_Z
  END INTERFACE

  PROCEDURE (FILTER_X) :: filx_00, filx_11, filx_12, filx_21, filx_22
  PROCEDURE (FILTER_Y) :: fily_00, fily_11, fily_12, fily_21, fily_22
  PROCEDURE (FILTER_Z) :: filz_00, filz_11, filz_12, filz_21, filz_22
  PROCEDURE (FILTER_X), POINTER :: filx=>null(), filxS=>null()
  PROCEDURE (FILTER_Y), POINTER :: fily=>null(), filyS=>null()
  PROCEDURE (FILTER_Z), POINTER :: filz=>null(), filzS=>null()

contains



end module x3d_filters
