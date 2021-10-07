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

module mom

  use MPI
 
  use decomp_2d, only : nrank 
  use decomp_2d, only : mytype, real_type
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use param, only : dx, dy, dz, xlx, yly, zlz
  use var, only : nx, ny, nz
  use var, only : zero, two, pi

  implicit none
  
  private
  public :: vel, test_du, test_dv, test_dw

contains

  subroutine vel(u, v, w)

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: u, v, w

    real(mytype) :: x, y, z
    integer :: i, j, k

    do k = 1, xsize(3)
       z = (k + xstart(3) - 2) * dz
       do j = 1, xsize(2)
          y = (j + xstart(2) - 2) * dy
          do i = 1, xsize(1)
             x = (i + xstart(1) - 2) * dx

             u(i, j, k) = sin(two * pi * (x / xlx)) * cos(two * pi * (y / yly)) * cos(two * pi * (z / zlz))
             v(i, j, k) = cos(two * pi * (x / xlx)) * sin(two * pi * (y / yly)) * cos(two * pi * (z / zlz))
             w(i, j, k) = -two * cos(two * pi * (x / xlx)) * cos(two * pi * (y / yly)) * sin(two * pi * (z / zlz))
          enddo
       enddo
    enddo
             
  endsubroutine vel

  subroutine test_du(du)

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: du

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ierr

    real(mytype) :: err, errloc
    real(mytype) :: du_ana

    if (nrank .eq. 0 ) then
       open(101, file="du.dat", action="write", status="unknown")
    endif
    
    errloc = zero
    do k = 1, xsize(3)
       z = (k + xstart(3) - 2) * dz
       do j = 1, xsize(2)
          y = (j + xstart(2) - 2) * dy
          do i = 1, xsize(1)
             x = (i + xstart(1) - 2) * dx

             du_ana = cos(two * pi * (x / xlx)) * cos(two * pi * (y / yly)) * cos(two * pi * (z / zlz))
             du_ana = two * pi * du_ana / xlx
             errloc = errloc + (du(i, j, k) - du_ana)**2

             if (nrank .eq. 0) then
                if ((j.eq.1) .and. (k.eq.1)) then
                   write(101, *) x, du(i, j, k), du_ana, errloc
                endif
             endif
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
    err = sqrt(err / nx / ny / nz)

    if (nrank .eq. 0) then
       print *, "RMS error in dudx: ", err
    end if

    if (nrank .eq. 0) then
       close(101)
    endif
    
  endsubroutine test_du

  subroutine test_dv(dv)

    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: dv

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ierr

    real(mytype) :: err, errloc
    real(mytype) :: dv_ana
    
    if (nrank .eq. 0) then
       open(102, file="dv.dat", action="write", status="unknown")
    endif

    errloc = zero
    do k = 1, ysize(3)
       z = (k + ystart(3) - 2) * dz
       do j = 1, ysize(2)
          y = (j + ystart(2) - 2) * dy
          do i = 1, ysize(1)
             x = (i + ystart(1) - 2) * dx

             dv_ana = cos(two * pi * (x / xlx)) * cos(two * pi * (y / yly)) * cos(two * pi * (z / zlz))
             dv_ana = two * pi * dv_ana / yly
             errloc = errloc + (dv(i, j, k) - dv_ana)**2
             if (nrank .eq. 0) then
                if ((i.eq.1) .and. (k.eq.1)) then
                   write(102, *) y, dv(i, j, k), dv_ana, errloc
                endif
             endif
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
    err = sqrt(err / nx / ny / nz)

    if (nrank .eq. 0) then
       print *, "RMS error in dvdy: ", err
    end if

    if (nrank .eq. 0) then
       close(102)
    endif

  endsubroutine test_dv

  subroutine test_dw(dw)

    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: dw

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ierr

    real(mytype) :: err, errloc
    real(mytype) :: dw_ana
    
    if (nrank .eq. 0) then
       open(103, file="dw.dat", action="write", status="unknown")
    endif
    
    errloc = zero
    do k = 1, zsize(3)
       z = (k + zstart(3) - 2) * dz
       do j = 1, zsize(2)
          y = (j + zstart(2) - 2) * dy
          do i = 1, zsize(1)
             x = (i + zstart(1) - 2) * dx

             dw_ana = -two * cos(two * pi * (x / xlx)) * cos(two * pi * (y / yly)) * cos(two * pi * (z / zlz))
             dw_ana = two * pi * dw_ana / zlz
             errloc = errloc + (dw(i, j, k) - dw_ana)**2

             if (nrank .eq. 0) then
                if ((i.eq.1) .and. (j.eq.1)) then
                   write(103, *) z, dw(i, j, k), dw_ana, errloc
                endif
             endif
          enddo
       enddo
    enddo
    call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, ierr)
    err = sqrt(err / nx / ny / nz)

    if (nrank .eq. 0) then
       print *, "RMS error in dwdz: ", err
    end if

    if (nrank .eq. 0) then
       close(103)
    endif

  endsubroutine test_dw
  
end module mom
