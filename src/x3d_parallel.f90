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

module x3dparallel

  use decomp_2d, only : mytype, nproc

  implicit none

  private        ! Make everything private unless declared public


  public ::  x3d_transpose_x_to_y, &
             x3d_transpose_y_to_z, &
             x3d_transpose_z_to_y, &
             x3d_transpose_y_to_x
contains

  !############################################################################
  !!  SUBROUTINE: x3d_transpose_x_to_y
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_x_to_y(data_in, data_out)
    
    use decomp_2d, only : xsize, ysize
    use decomp_2d, only : transpose_x_to_y

    implicit none

    !! Input/Output
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: data_in
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    !! Local
    integer :: i, j, k

    if (nproc == 1) then
      do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_x_to_y(data_in,data_out)
    endif

  end subroutine x3d_transpose_x_to_y
  !############################################################################
  !!  SUBROUTINE: x3d_transpose_y_to_z
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_y_to_z(data_in, data_out)
    
    use decomp_2d, only : ysize, zsize
    use decomp_2d, only : transpose_y_to_z

    implicit none

    !! Input/Output
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(out) :: data_out

    !! Local
    integer :: i, j, k

    if (nproc == 1) then
      do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_z(data_in,data_out)
    endif

  end subroutine x3d_transpose_y_to_z
  !############################################################################
  !!  SUBROUTINE: x3d_transpose_z_to_y
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_z_to_y(data_in, data_out)
    
    use decomp_2d, only : ysize, zsize
    use decomp_2d, only : transpose_z_to_y

    implicit none

    !! Input/Output
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: data_in
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    !! Local
    integer :: i, j, k

    if (nproc == 1) then
      do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_z_to_y(data_in,data_out)
    endif

  end subroutine x3d_transpose_z_to_y
  !############################################################################
  !!  SUBROUTINE: x3d_transpose_y_to_z
  !!      AUTHOR: Stefano Rolfo
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_y_to_x(data_in, data_out)
    
    use decomp_2d, only : ysize, xsize
    use decomp_2d, only : transpose_y_to_x

    implicit none

    !! Input/Output
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: data_out

    !! Local
    integer :: i, j, k

    if (nproc == 1) then
      do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_x(data_in,data_out)
    endif

  end subroutine x3d_transpose_y_to_x
  !############################################################################

end module x3dparallel
