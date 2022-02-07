!################################################################################
!This file is part of Xcompact3d
!
!Xcompact3d
!Copyright (c) 2012-2022, Xcompact3d
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!BSD 3-Clause License
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.!
!
!3. Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

module x3d_transpose

  use decomp_2d, only : mytype, decomp_info, decomp_main
  use variables, only : p_row, p_col

  implicit none

  ! Make everything private unless declared public
  private

  public ::  x3d_transpose_x_to_y, &
             x3d_transpose_y_to_z, &
             x3d_transpose_z_to_y, &
             x3d_transpose_y_to_x

  interface x3d_transpose_x_to_y
    module procedure x3d_transpose_x_to_y_real_short
    module procedure x3d_transpose_x_to_y_real
    module procedure x3d_transpose_x_to_y_cplx_short
    module procedure x3d_transpose_x_to_y_cplx
  end interface x3d_transpose_x_to_y

  interface x3d_transpose_y_to_z
    module procedure x3d_transpose_y_to_z_real_short
    module procedure x3d_transpose_y_to_z_real
    module procedure x3d_transpose_y_to_z_cplx_short
    module procedure x3d_transpose_y_to_z_cplx
  end interface x3d_transpose_y_to_z

  interface x3d_transpose_z_to_y
    module procedure x3d_transpose_z_to_y_real_short
    module procedure x3d_transpose_z_to_y_real
    module procedure x3d_transpose_z_to_y_cplx_short
    module procedure x3d_transpose_z_to_y_cplx
  end interface x3d_transpose_z_to_y

  interface x3d_transpose_y_to_x
    module procedure x3d_transpose_y_to_x_real_short
    module procedure x3d_transpose_y_to_x_real
    module procedure x3d_transpose_y_to_x_cplx_short
    module procedure x3d_transpose_y_to_x_cplx
  end interface x3d_transpose_y_to_x


contains


  !############################################################################
  !!  SUBROUTINE: x3d_transpose_x_to_y
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_x_to_y_real(data_in, data_out, decomp)

    use decomp_2d, only : transpose_x_to_y

    implicit none

    !! Input/Output
    real(mytype), dimension(:,:,:), intent(in) :: data_in
    real(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_row == 1) then
      do concurrent (k=1:decomp%xsz(3), j=1:decomp%xsz(2), i=1:decomp%xsz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_x_to_y(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_x_to_y_real
  !############################################################################
  subroutine x3d_transpose_x_to_y_cplx(data_in, data_out, decomp)

    use decomp_2d, only : transpose_x_to_y

    implicit none

    !! Input/Output
    complex(mytype), dimension(:,:,:), intent(in) :: data_in
    complex(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_row == 1) then
      do concurrent (k=1:decomp%xsz(3), j=1:decomp%xsz(2), i=1:decomp%xsz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_x_to_y(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_x_to_y_cplx


  !############################################################################
  !!  SUBROUTINE: x3d_transpose_y_to_z
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_y_to_z_real(data_in, data_out, decomp)

    use decomp_2d, only : transpose_y_to_z

    implicit none

    !! Input/Output
    real(mytype), dimension(:,:,:), intent(in) :: data_in
    real(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_col == 1) then
      do concurrent (k=1:decomp%ysz(3), j=1:decomp%ysz(2), i=1:decomp%ysz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_z(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_y_to_z_real
  !############################################################################
  subroutine x3d_transpose_y_to_z_cplx(data_in, data_out, decomp)

    use decomp_2d, only : transpose_y_to_z

    implicit none

    !! Input/Output
    complex(mytype), dimension(:,:,:), intent(in) :: data_in
    complex(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_col == 1) then
      do concurrent (k=1:decomp%ysz(3), j=1:decomp%ysz(2), i=1:decomp%ysz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_z(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_y_to_z_cplx


  !############################################################################
  !!  SUBROUTINE: x3d_transpose_z_to_y
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_z_to_y_real(data_in, data_out, decomp)

    use decomp_2d, only : transpose_z_to_y

    implicit none

    !! Input/Output
    real(mytype), dimension(:,:,:), intent(in) :: data_in
    real(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_col == 1) then
      do concurrent (k=1:decomp%zsz(3), j=1:decomp%zsz(2), i=1:decomp%zsz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_z_to_y(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_z_to_y_real
  !############################################################################
  subroutine x3d_transpose_z_to_y_cplx(data_in, data_out, decomp)

    use decomp_2d, only : transpose_z_to_y

    implicit none

    !! Input/Output
    complex(mytype), dimension(:,:,:), intent(in) :: data_in
    complex(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_col == 1) then
      do concurrent (k=1:decomp%zsz(3), j=1:decomp%zsz(2), i=1:decomp%zsz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_z_to_y(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_z_to_y_cplx


  !############################################################################
  !!  SUBROUTINE: x3d_transpose_y_to_z
  !! DESCRIPTION: Wrapper around decomp2d_transpose to avoid MPI in case of 
  !!              single core calculation
  !############################################################################
  subroutine x3d_transpose_y_to_x_real(data_in, data_out, decomp)
    
    use decomp_2d, only : transpose_y_to_x

    implicit none

    !! Input/Output
    real(mytype), dimension(:,:,:), intent(in) :: data_in
    real(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_row == 1) then
      do concurrent (k=1:decomp%ysz(3), j=1:decomp%ysz(2), i=1:decomp%ysz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_x(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_y_to_x_real
  !############################################################################
  subroutine x3d_transpose_y_to_x_cplx(data_in, data_out, decomp)
    
    use decomp_2d, only : transpose_y_to_x

    implicit none

    !! Input/Output
    complex(mytype), dimension(:,:,:), intent(in) :: data_in
    complex(mytype), dimension(:,:,:), intent(out) :: data_out
    type(decomp_info), intent(in) :: decomp

    !! Local
    integer :: i, j, k

    if (p_row == 1) then
      do concurrent (k=1:decomp%ysz(3), j=1:decomp%ysz(2), i=1:decomp%ysz(1))
        data_out(i,j,k) = data_in(i,j,k)
      enddo
    else 
      call transpose_y_to_x(data_in,data_out,decomp)
    endif

  end subroutine x3d_transpose_y_to_x_cplx


  !############################################################################
  !!  SUBROUTINE: x3d_transpose_*_to_*_short
  !! DESCRIPTION: Call the x3d_transpose_*_to_* with the decomp_info object
  !############################################################################
  subroutine x3d_transpose_x_to_y_real_short(data_in, data_out)

    use decomp_2d, only : xsize, ysize

    implicit none

    !! Input/Output
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: data_in
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    call x3d_transpose_x_to_y(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_x_to_y_real_short
  !############################################################################
  subroutine x3d_transpose_x_to_y_cplx_short(data_in, data_out)

    use decomp_2d, only : xsize, ysize

    implicit none

    !! Input/Output
    complex(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: data_in
    complex(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    call x3d_transpose_x_to_y(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_x_to_y_cplx_short
  !############################################################################
  subroutine x3d_transpose_y_to_z_real_short(data_in, data_out)

    use decomp_2d, only : ysize, zsize

    implicit none

    !! Input/Output
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(out) :: data_out

    call x3d_transpose_y_to_z(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_y_to_z_real_short
  !############################################################################
  subroutine x3d_transpose_y_to_z_cplx_short(data_in, data_out)

    use decomp_2d, only : ysize, zsize

    implicit none

    !! Input/Output
    complex(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    complex(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(out) :: data_out

    call x3d_transpose_y_to_z(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_y_to_z_cplx_short
  !############################################################################
  subroutine x3d_transpose_z_to_y_real_short(data_in, data_out)

    use decomp_2d, only : zsize, ysize

    implicit none

    !! Input/Output
    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: data_in
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    call x3d_transpose_z_to_y(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_z_to_y_real_short
  !############################################################################
  subroutine x3d_transpose_z_to_y_cplx_short(data_in, data_out)

    use decomp_2d, only : zsize, ysize

    implicit none

    !! Input/Output
    complex(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: data_in
    complex(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(out) :: data_out

    call x3d_transpose_z_to_y(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_z_to_y_cplx_short
  !############################################################################
  subroutine x3d_transpose_y_to_x_real_short(data_in, data_out)

    use decomp_2d, only : ysize, xsize

    implicit none

    !! Input/Output
    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: data_out

    call x3d_transpose_y_to_x(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_y_to_x_real_short
  !############################################################################
  subroutine x3d_transpose_y_to_x_cplx_short(data_in, data_out)

    use decomp_2d, only : ysize, xsize

    implicit none

    !! Input/Output
    complex(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: data_in
    complex(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: data_out

    call x3d_transpose_y_to_x(data_in, data_out, decomp_main)

  end subroutine x3d_transpose_y_to_x_cplx_short

end module x3d_transpose
