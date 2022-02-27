!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module bc_tgv

  use MPI

  use decomp_2d, only : mytype, real_type, decomp_2d_warning
  use decomp_2d, only : nrank, nproc, xsize, xstart
  use x3d_precision, only : pi, twopi
  use variables, only : nx, ny, nz
  use param, only : dx, dy, dz, xlx, yly, zlz, dt, xnu
  use param, only : zero, one, two

  implicit none

  integer, save :: tgv_iounit
  integer, parameter :: tgv_verbose = 0
  character(len=*), parameter :: tgv_file = "tgv.dat"

  ! Make everything private unless declared public
  private
  public :: tgv_boot, &
            tgv_listing, &
            tgv_init, &
            tgv_postprocess, &
            tgv_finalize

contains

  !
  ! Initialize case-specific parameters and IO
  !
  subroutine tgv_boot()

    use variables, only : nxm, nym, nzm
    use param, only : dt, dx, dy, dz, dx2, dy2, dz2, &
                      xlx, yly, zlz, re, xnu, one

    implicit none

    ! Local variable
    logical :: fexists

    ! Default domain size is 2 pi x 2 pi x 1
    ! This should be inside input.i3d
    xlx = twopi
    yly = twopi
    zlz = twopi
    dx = xlx / real(nxm, mytype)
    dy = yly / real(nym, mytype)
    dz = zlz / real(nzm, mytype)
    dx2 = dx * dx
    dy2 = dy * dy
    dz2 = dz * dz

    ! Default time step : CFL = 0.2 and U = 1
    ! This should be inside input.i3d ?
    dt = 0.2_mytype * dx

    ! Default Re is 1600 when Re = 0
    ! This should be inside input.i3d ?
    if (abs(re) <= epsilon(re)) then
      re = 1600._mytype
      xnu = one / re
    endif

    ! Check if the file is present and get IO unit
    if (nrank == 0) then

      inquire(file=tgv_file, exist=fexists)
      if (fexists) then
        if (tgv_verbose > 0) call decomp_2d_warning(1, "TGV: file "//tgv_file//" replaced.")
        open(newunit=tgv_iounit, file=tgv_file, action='write', status='replace')
      else
        open(newunit=tgv_iounit, file=tgv_file, action='write', status='new')
      endif    

    endif

  end subroutine tgv_boot

  !
  ! Case-specific parameters in the listing
  !
  subroutine tgv_listing()

    implicit none

    if (nrank == 0) then
      write(*,*)' 3D TGV test case'
      write(*,*)'==========================================================='
    endif

  end subroutine tgv_listing

  !
  ! Initialize 3D fields
  !
  subroutine tgv_init(ux1, uy1, uz1)

    implicit none 

    ! Arguments
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    ! Local variables
    real(mytype) :: x, y, z
    integer :: i, j, k, ip, it

    ! Initial velocity
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
       z = (k + xstart(3) - 2) * dz
       y = (j + xstart(2) - 2) * dy
       x = (i + xstart(1) - 2) * dx
       ux1(i, j, k) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
       uy1(i, j, k) = -cos(twopi * (x / xlx)) * sin(twopi * (y / yly)) * cos(twopi * (z / zlz))
       uz1(i, j, k) = zero
    enddo

    ! Check initial error
    call tgv_postprocess(ux1, uy1, uz1, 1)

  endsubroutine tgv_init

  !
  ! Compare the velocity field with the analytical solution
  !
  subroutine tgv_postprocess(ux1, uy1, uz1, ndt)

    implicit none

    ! Arguments
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    integer, intent(in) :: ndt

    ! Local variables
    real(mytype) :: x, y, z
    integer :: i, j, k, code

  end subroutine tgv_postprocess

  !
  ! Finalize case-specific IO
  !
  subroutine tgv_finalize()

    implicit none

    if (nrank == 0) close(tgv_iounit)

  end subroutine tgv_finalize

end module bc_tgv
