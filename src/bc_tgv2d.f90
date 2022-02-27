!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module bc_tgv2d

  use MPI

  use decomp_2d, only : mytype, real_type, decomp_2d_warning
  use decomp_2d, only : nrank, nproc, xsize, xstart
  use x3d_precision, only : pi, twopi
  use variables, only : nx, ny, nz
  use param, only : dx, dy, dz, xlx, yly, zlz, dt, xnu
  use param, only : zero, one, two

  implicit none

  integer, save :: tgv2d_iounit
  integer, parameter :: tgv2d_verbose = 0
  character(len=*), parameter :: tgv2d_file = "tgv2d.dat"

  ! Make everything private unless declared public
  private
  public :: tgv2d_boot, &
            tgv2d_listing, &
            tgv2d_init, &
            tgv2d_postprocess, &
            tgv2d_finalize

contains

  !
  ! Initialize case-specific parameters and IO
  !
  subroutine tgv2d_boot()

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
    zlz = one
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

      inquire(file=tgv2d_file, exist=fexists)
      if (fexists) then
        if (tgv2d_verbose > 0) call decomp_2d_warning(1, "TGV2D: file "//tgv2d_file//" replaced.")
        open(newunit=tgv2d_iounit, file=tgv2d_file, action='write', status='replace')
      else
        open(newunit=tgv2d_iounit, file=tgv2d_file, action='write', status='new')
      endif    

    endif

  end subroutine tgv2d_boot

  !
  ! Case-specific parameters in the listing
  !
  subroutine tgv2d_listing()

    implicit none

    if (nrank == 0) then
      write(*,*)' 2D TGV test case'
      write(*,*)' Error estimator valid for explicit Euler only'
      write(*,*)'==========================================================='
    endif

  end subroutine tgv2d_listing

  !
  ! Initialize 3D fields
  !
  subroutine tgv2d_init(ux1, uy1, uz1)

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
       ux1(i, j, k) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly))
       uy1(i, j, k) = -cos(twopi * (x / xlx)) * sin(twopi * (y / yly))
       uz1(i, j, k) = zero
    enddo

    ! Check initial error
    call tgv2d_postprocess(ux1, uy1, uz1, 1)

  endsubroutine tgv2d_init

  !
  ! Compare the velocity field with the analytical solution
  !
  subroutine tgv2d_postprocess(ux1, uy1, uz1, ndt)

    implicit none

    ! Arguments
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
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
             err_l1(1) = err_l1(1) + abs(ux1(i,j,k) - sol(1))
             err_l1(2) = err_l1(2) + abs(uy1(i,j,k) - sol(2))
             err_l1(3) = err_l1(3) + abs(uz1(i,j,k) - sol(3))
             err_l2(1) = err_l2(1) + (ux1(i,j,k) - sol(1))**2
             err_l2(2) = err_l2(2) + (uy1(i,j,k) - sol(2))**2
             err_l2(3) = err_l2(3) + (uz1(i,j,k) - sol(3))**2
             err_linf(1) = max(err_linf(1), abs(ux1(i,j,k)-sol(1)))
             err_linf(2) = max(err_linf(2), abs(uy1(i,j,k)-sol(2)))
             err_linf(3) = max(err_linf(3), abs(uz1(i,j,k)-sol(3)))

          enddo
       enddo
    enddo

    uvw_max(1) = maxval(ux1)
    uvw_max(2) = maxval(uy1)
    uvw_max(3) = maxval(uz1)
    uvw_min(1) = minval(ux1)
    uvw_min(2) = minval(uy1)
    uvw_min(3) = minval(uz1)

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
       ! Listing
       write(*,*) "Max: ", uvw_max
       write(*,*) "Min: ", uvw_min
       write(*,*) "Amplitude: ", real(uvw_max-uvw_min)
       write(*,*) "TGV 2D error, u, ", err_l1(1), err_l2(1), err_linf(1)
       write(*,*) "TGV 2D error, v, ", err_l1(2), err_l2(2), err_linf(2)
       write(*,*) "TGV 2D error, w, ", err_l1(3), err_l2(3), err_linf(3)
       ! tgv2d.dat
       write(tgv2d_iounit,*) "Max: ", ndt, uvw_max
       write(tgv2d_iounit,*) "Min: ", ndt, uvw_min
       write(tgv2d_iounit,*) "Amplitude: ", ndt, real(uvw_max-uvw_min)
       write(tgv2d_iounit,*) "TGV 2D error, u, ", ndt, err_l1(1), err_l2(1), err_linf(1)
       write(tgv2d_iounit,*) "TGV 2D error, v, ", ndt, err_l1(2), err_l2(2), err_linf(2)
       write(tgv2d_iounit,*) "TGV 2D error, w, ", ndt, err_l1(3), err_l2(3), err_linf(3)
    endif

  end subroutine tgv2d_postprocess

  !
  ! Finalize case-specific IO
  !
  subroutine tgv2d_finalize()

    implicit none

    if (nrank == 0) close(tgv2d_iounit)

  end subroutine tgv2d_finalize

  !
  ! Compute the modified wavenumber for the second derivative in x
  !
  subroutine compute_kx2(kin,k2out)

    use param, only : dx2, three, four, half, nine, eight
    use x3d_operator_x_data

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

  !
  ! Compute the modified wavenumber for the second derivative in y
  !
  subroutine compute_ky2(kin,k2out)

    use param, only : dx2, three, four, half, nine, eight
    use x3d_operator_y_data

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

end module bc_tgv2d
