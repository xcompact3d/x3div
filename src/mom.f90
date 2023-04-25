!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module mom

  use MPI
 
  use decomp_2d_constants, only : mytype, real_type
  use decomp_2d_mpi, only : nrank, nproc, decomp_2d_warning, decomp_2d_abort
  use decomp_2d, only : xsize, ysize, zsize
  use decomp_2d, only : xstart, ystart, zstart
  use x3d_precision, only : twopi
  use param    , only : dx, dy, dz, xlx, yly, zlz
  use variables, only : test_mode
  use variables, only : nx, ny, nz
  use param    , only : zero, two

  implicit none
  
  private
  public :: vel, test_du, test_dv, test_dw

contains

  subroutine vel(u, v, w)

    implicit none 

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(out) :: u, v, w

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: xe1, xe2, xe3
    integer :: xs1, xs2, xs3

    xe1 = xsize(1)
    xe2 = xsize(2)
    xe3 = xsize(3)
    xs1 = xstart(1)
    xs2 = xstart(2)
    xs3 = xstart(3)

    ! Initial velocity
    !$acc kernels default(present)
    do concurrent (k=1:xe3, j=1:xe2, i=1:xe1)
       z = (k + xs3 - 2) * dz
       y = (j + xs2 - 2) * dy
       x = (i + xs1 - 2) * dx
       u(i, j, k) = sin(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
       v(i, j, k) = cos(twopi * (x / xlx)) * sin(twopi * (y / yly)) * cos(twopi * (z / zlz))
       w(i, j, k) = -two * cos(twopi * (x / xlx)) * cos(twopi * (y / yly)) * sin(twopi * (z / zlz))
    enddo
    !$acc end kernels

             
  endsubroutine vel

  subroutine test_du(du)

    implicit none

    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: du

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ne1, ne2, ne3
    integer :: ns1, ns2, ns3
    integer :: code

    real(mytype) :: err, errloc
    real(mytype) :: du_ana

    if (test_mode) then
       !if (nrank .eq. 0 ) then
       !   open(101, file="du.dat", action="write", status="unknown")
       !endif
      
       ne1 = xsize(1)
       ne2 = xsize(2)
       ne3 = xsize(3)
       ns1 = xstart(1)
       ns2 = xstart(2)
       ns3 = xstart(3)
    
       errloc = zero
       !$acc parallel loop default(present) reduction(+:errloc)
       do k = 1, ne3
          z = (k + ns3 - 2) * dz
          do j = 1, ne2
             y = (j + ns2 - 2) * dy
             do i = 1, ne1
                x = (i + ns1 - 2) * dx
                
                du_ana = cos(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
                du_ana = twopi * du_ana / xlx
                errloc = errloc + (du(i, j, k) - du_ana)**2
                
                !if (nrank .eq. 0) then
                !   if ((j.eq.1) .and. (k.eq.1)) then
                !      write(101, *) x, du(i, j, k), du_ana, errloc
                !   endif
                !endif
             enddo
          enddo
       enddo
       !$acc end loop
       call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       err = sqrt(err / nx / ny / nz)

       if (nrank .eq. 0) then
          print *, "RMS error in dudx: ", err
       end if

       !if (nrank .eq. 0) then
       !   close(101)
       !endif
    end if
    
  endsubroutine test_du

  subroutine test_dv(dv)

    implicit none

    real(mytype), dimension(ysize(1), ysize(2), ysize(3)), intent(in) :: dv

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ne1, ne2, ne3
    integer :: ns1, ns2, ns3
    integer :: code

    real(mytype) :: err, errloc
    real(mytype) :: dv_ana

    if (test_mode) then
       !if (nrank .eq. 0) then
       !   open(102, file="dv.dat", action="write", status="unknown")
       !endif
       ne1 = ysize(1)
       ne2 = ysize(2)
       ne3 = ysize(3)
       ns1 = ystart(1)
       ns2 = ystart(2)
       ns3 = ystart(3)

       errloc = zero
       !$acc parallel loop default(present) reduction(+:errloc)
       do k = 1, ne3
          z = (k + ns3 - 2) * dz
          do j = 1, ne2
             y = (j + ns2 - 2) * dy
             do i = 1, ne1
                x = (i + ns1 - 2) * dx
                
                dv_ana = cos(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
                dv_ana = twopi * dv_ana / yly
                errloc = errloc + (dv(i, j, k) - dv_ana)**2
                !if (nrank .eq. 0) then
                !   if ((i.eq.1) .and. (k.eq.1)) then
                !      write(102, *) y, dv(i, j, k), dv_ana, errloc
                !   endif
                !endif
             enddo
          enddo
       enddo
       !$acc end loop
       call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       err = sqrt(err / nx / ny / nz)

       if (nrank .eq. 0) then
          print *, "RMS error in dvdy: ", err
       end if
       
       !if (nrank .eq. 0) then
       !   close(102)
       !endif
    end if
    
  endsubroutine test_dv

  subroutine test_dw(dw)

    real(mytype), dimension(zsize(1), zsize(2), zsize(3)), intent(in) :: dw

    real(mytype) :: x, y, z
    integer :: i, j, k
    integer :: ne1, ne2, ne3
    integer :: ns1, ns2, ns3
    integer :: code

    real(mytype) :: err, errloc
    real(mytype) :: dw_ana

    if (test_mode) then
       !if (nrank .eq. 0) then
       !   open(103, file="dw.dat", action="write", status="unknown")
       !endif
       ne1 = zsize(1)
       ne2 = zsize(2)
       ne3 = zsize(3)
       ns1 = zstart(1)
       ns2 = zstart(2)
       ns3 = zstart(3)
       errloc = zero
       !$acc parallel loop default(present) reduction(+:errloc)
       do k = 1, ne3
          z = (k + ns3 - 2) * dz
          do j = 1, ne2
             y = (j + ns2 - 2) * dy
             do i = 1, ne1
                x = (i + ns1 - 2) * dx

                dw_ana = -two * cos(twopi * (x / xlx)) * cos(twopi * (y / yly)) * cos(twopi * (z / zlz))
                dw_ana = twopi * dw_ana / zlz
                errloc = errloc + (dw(i, j, k) - dw_ana)**2

                !if (nrank .eq. 0) then
                !   if ((i.eq.1) .and. (j.eq.1)) then
                !      write(103, *) z, dw(i, j, k), dw_ana, errloc
                !   endif
                !endif
             enddo
          enddo
       enddo
       !$acc end loop
       call MPI_ALLREDUCE(errloc, err, 1, real_type, MPI_SUM, MPI_COMM_WORLD, code)
       if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
       err = sqrt(err / nx / ny / nz)

       if (nrank .eq. 0) then
          print *, "RMS error in dwdz: ", err
       end if

       !if (nrank .eq. 0) then
       !   close(103)
       !endif
    end if

  endsubroutine test_dw
  
end module mom
