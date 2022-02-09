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

module mom

  use MPI
 
  use decomp_2d, only : mytype, real_type, decomp_2d_warning
  use decomp_2d, only : nrank, nproc
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use x3dprecision, only : twopi
  use param    , only : dx, dy, dz, xlx, yly, zlz
  use variables, only : test_mode
  use variables, only : nx, ny, nz
  use param    , only : zero, one, two

  implicit none
  
  private
  public :: vel, test_tgv2d

contains

  subroutine vel(u, v, w)

    implicit none 

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: u, v, w

    real(mytype) :: x, y, z
    integer :: i, j, k

    do k = 1, xsize(3)
       z = (k + xstart(3) - 2) * dz
       do j = 1, xsize(2)
          y = (j + xstart(2) - 2) * dy
          do i = 1, xsize(1)
             x = (i + xstart(1) - 2) * dx
             u(i, j, k) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly))
             v(i, j, k) = -cos(twopi * (x / xlx)) * sin(twopi * (y / yly))
             w(i, j, k) = zero
          enddo
       enddo
    enddo
             
  endsubroutine vel

  !
  ! Compare the velocity field with the analytical solution
  !
  subroutine test_tgv2d(u, v, w, ndt)

    use param, only : dt, xnu

    implicit none

    ! Arguments
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: u, v, w
    integer, intent(in) :: ndt

    ! Local variables
    real(mytype) :: x, y, z, factor, kxtgv, kytgv, kx2tgv, ky2tgv
    real(mytype), dimension(3) :: sol, err_l1, err_l2, err_linf, uvw_max, uvw_min
    integer :: i, j, k, it, code

    ! Temporal variation (explicit Euler time scheme)
    kxtgv = twopi / xlx
    kytgv = twopi / yly
    call compute_kx2(kxtgv, kx2tgv)
    call compute_ky2(kytgv, ky2tgv)
    factor = one
    do it = 2, ndt
       factor = factor * (one - dt*xnu*(kx2tgv+ky2tgv))
    enddo

    ! Init
    err_l1(:) = zero
    err_l2(:) = zero
    err_linf(:) = zero

    do k = 1, xsize(3)
       z = (k + xstart(3) - 2)*dz
       do j = 1, xsize(2)
          y = (j + xstart(2) - 2)*dy
          do i = 1, xsize(1)
             x = (i + xstart(1) - 2)*dx

             ! Initial condition
             sol(1) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly))
             sol(2) = -cos(twopi * (x / xlx)) * sin(twopi * (y / yly))
             sol(3) = zero

             ! Analytical discrete solution
             sol(:) = sol(:) * factor

             ! Update the errors
             err_l1(1) = err_l1(1) + abs(u(i,j,k) - sol(1))
             err_l1(2) = err_l1(2) + abs(v(i,j,k) - sol(2))
             err_l1(3) = err_l1(3) + abs(w(i,j,k) - sol(3))
             err_l2(1) = err_l2(1) + (u(i,j,k) - sol(1))**2
             err_l2(2) = err_l2(2) + (v(i,j,k) - sol(2))**2
             err_l2(3) = err_l2(3) + (w(i,j,k) - sol(3))**2
             err_linf(1) = max(err_linf(1), abs(u(i,j,k)-sol(1)))
             err_linf(2) = max(err_linf(2), abs(v(i,j,k)-sol(2)))
             err_linf(3) = max(err_linf(3), abs(w(i,j,k)-sol(3)))

          enddo
       enddo
    enddo

    uvw_max(1) = maxval(u)
    uvw_max(2) = maxval(v)
    uvw_max(3) = maxval(w)
    uvw_min(1) = minval(u)
    uvw_min(2) = minval(v)
    uvw_min(3) = minval(w)

    ! Compute global errors if needed
    if (nproc > 1) then
       call MPI_ALLREDUCE(MPI_IN_PLACE, err_l1, 3, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       call MPI_ALLREDUCE(MPI_IN_PLACE, err_l2, 3, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       call MPI_ALLREDUCE(MPI_IN_PLACE, err_linf, 3, real_type, MPI_MAX, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       call MPI_ALLREDUCE(MPI_IN_PLACE, uvw_max, 3, real_type, MPI_MAX, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       call MPI_ALLREDUCE(MPI_IN_PLACE, uvw_min, 3, real_type, MPI_MIN, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
    endif

    ! Rescale L1 and L2 errors
    err_l1(:) = err_l1(:) / dble(nx*ny*nz)
    err_l2(:) = sqrt(err_l2(:) / dble(nx*ny*nz))

    ! Print the error for each velocity component
    if (nrank == 0) then
       write(*,*) "Max: ", uvw_max
       write(*,*) "Min: ", uvw_min
       write(*,*) "Amplitude: ", real(uvw_max-uvw_min)
       write(*,*) "TGV 2D error, u, ", err_l1(1), err_l2(1), err_linf(1)
       write(*,*) "TGV 2D error, v, ", err_l1(2), err_l2(2), err_linf(2)
       write(*,*) "TGV 2D error, w, ", err_l1(3), err_l2(3), err_linf(3)
    endif

  end subroutine test_tgv2d
 
  ! Compute the modified wavenumber for the second derivative in x
  subroutine compute_kx2(kin,k2out)

    use param, only : dx2, three, four, half, nine, eight
    use x3d_operator_x_data
    use x3dprecision, only : pi

    implicit none

    real(mytype), intent(in) :: kin
    real(mytype), intent(out) :: k2out

    if (kin.lt.zero .or. kin.gt.pi/min(dx,dy)) then
      if (nrank==0) write(*,*) "TGV2D: Warning, incorrect wavenumber provided."
    endif

    k2out = asix * two * (one - cos(kin*dx)) &
          + four * bsix * half * (one - cos(two*kin*dx)) &
          + nine * csix * (two / nine) * (one - cos(three*kin*dx)) &
          + 16._mytype * dsix * (one / eight) * (one - cos(four*kin*dx))
    k2out = k2out / (one + two * alsaix * cos(kin*dx))

  end subroutine compute_kx2

  ! Compute the modified wavenumber for the second derivative in y
  subroutine compute_ky2(kin,k2out)

    use param, only : dx2, three, four, half, nine, eight
    use x3d_operator_y_data
    use x3dprecision, only : pi

    implicit none

    real(mytype), intent(in) :: kin
    real(mytype), intent(out) :: k2out

    if (kin.lt.zero .or. kin.gt.pi/min(dx,dy)) then
      if (nrank==0) write(*,*) "TGV2D: Warning, incorrect wavenumber provided."
    endif

    k2out = asjy * two * (one - cos(kin*dx)) &
          + four * bsjy * half * (one - cos(two*kin*dx)) &
          + nine * csjy * (two / nine) * (one - cos(three*kin*dx)) &
          + 16._mytype * dsjy * (one / eight) * (one - cos(four*kin*dx))
    k2out = k2out / (one + two * alsajy * cos(kin*dx))

  end subroutine compute_ky2

end module mom
