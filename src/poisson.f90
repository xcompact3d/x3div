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

module decomp_2d_poisson

  use decomp_2d, only : mytype
  use decomp_2d, only : DECOMP_INFO
  use decomp_2d, only : decomp_info_init, &
                        decomp_info_finalize
  use decomp_2d_fft, only : decomp_2d_fft_init, &
                            decomp_2d_fft_3d,   &
                            decomp_2d_fft_finalize 
  use x3d_transpose
  use param
  use variables

  implicit none

  private        ! Make everything private unless declared public

  !  real(mytype), private, parameter :: PI = 3.14159265358979323846_mytype

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16_mytype
#else
  real(mytype), parameter :: epsilon = 1.e-8_mytype
#endif

  ! boundary conditions
  integer, save :: bcx, bcy, bcz

  ! decomposition object for physical space
  type(DECOMP_INFO), save :: ph

  ! decomposition object for spectral space
  type(DECOMP_INFO), save :: sp

  ! store sine/cosine factors
  real(mytype), save, allocatable, dimension(:) :: az,bz
  real(mytype), save, allocatable, dimension(:) :: ay,by
  real(mytype), save, allocatable, dimension(:) :: ax,bx

  ! wave numbers
  complex(mytype), save, allocatable, dimension(:,:,:) :: kxyz
  !wave numbers for stretching in a pentadiagonal matrice
  complex(mytype), save, allocatable, dimension(:,:,:,:) :: a,a2,a3
  ! work arrays, 
  ! naming convention: cw (complex); rw (real); 
  !                    1 = X-pencil; 2 = Y-pencil; 3 = Z-pencil
  real(mytype), allocatable, dimension(:,:,:) :: rw1,rw1b,rw2,rw2b,rw3
  complex(mytype), allocatable, dimension(:,:,:) :: cw1,cw1b,cw2,cw22,cw2b,cw2c

  ! underlying FFT library only needs to be initialised once
  logical, save :: fft_initialised = .false.

  abstract interface
     subroutine poisson_xxx(rhs)
       use decomp_2d, only : mytype
       real(mytype), dimension(:,:,:), intent(inout) :: rhs
     end subroutine poisson_xxx
  end interface
  procedure (poisson_xxx), pointer :: poisson

  public :: decomp_2d_poisson_init,decomp_2d_poisson_finalize,poisson
contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise Poisson solver for given boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_init()

    use decomp_2d, only : nrank, nx_global, ny_global, nz_global

    implicit none

    integer :: nx, ny, nz, i
    
    real(mytype) :: rl, iy
    external  rl, iy

    if (nclx) then
       bcx=0
    else
       bcx=1
    endif
    if (ncly) then
       bcy=0
    else
       bcy=1
    endif
    if (nclz) then
       bcz=0
    else
       bcz=1
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Top level wrapper
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       poisson => poisson_000
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       poisson => poisson_100
    else
       stop 'boundary condition not supported'
    end if

    nx = nx_global
    ny = ny_global
    nz = nz_global

    ! pressure-grid having 1 fewer point for non-periodic directions
    if (bcx==1) nx=nx-1
    if (bcy==1) ny=ny-1
    if (bcz==1) nz=nz-1

#ifdef DEBUG 
    if (nrank .eq. 0) write(*,*)'# decomp_2d_poisson_init start'
#endif

    allocate(ax(nx),bx(nx))
    allocate(ay(ny),by(ny))
    allocate(az(nz),bz(nz))
    call abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)

#ifdef DEBUG 
    if (nrank .eq. 0) write(*,*)'# decomp_2d_poisson_init decomp_info_init'
#endif

    call decomp_info_init(nx, ny, nz, ph)
    call decomp_info_init(nx, ny, nz/2+1, sp)

#ifdef DEBUG 
    if (nrank .eq. 0) write(*,*)'# decomp_2d_poisson_init decomp_info_init ok'
#endif

    ! allocate work space
    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(kxyz(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
    else if (bcx==1 .and. bcy==1) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(cw1b(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw22(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2b(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(cw2c(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw1b(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       allocate(rw2b(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       if (bcz==1) then  
          allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       end if
       allocate(kxyz(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))    
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a2(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))      
    end if

#ifdef DEBUG 
    if (nrank .eq. 0) write(*,*)'# decomp_2d_poisson_init before waves'
#endif

    call waves()

#ifdef DEBUG 
    if (nrank .eq. 0) write(*,*)'# decomp_2d_poisson_init end'
#endif

    return
  end subroutine decomp_2d_poisson_init


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Release memory used by Poisson solver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_finalize

    implicit none

    deallocate(ax,bx,ay,by,az,bz)

    call decomp_info_finalize(ph)
    call decomp_info_finalize(sp)

    call decomp_2d_fft_finalize
    fft_initialised = .false.

    deallocate(kxyz)

    if (bcx==0 .and. bcy==0 .and. bcz==0) then
       deallocate(cw1)
       deallocate(a,a2,a3)
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       deallocate(cw1,cw1b,rw1,rw1b,rw2)
       deallocate(a,a2,a3)
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       deallocate(cw1,cw2,cw2b,rw2,rw2b)
       deallocate(a,a2,a3)
    else if (bcx==1 .and. bcy==1) then
       deallocate(cw1,cw1b,cw2,cw2b,rw1,rw1b,rw2,rw2b)
       deallocate(a,a2,a3)
       if (bcz==1) then
          deallocate(rw3)
       end if
    end if

    return
  end subroutine decomp_2d_poisson_finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Solving 3D Poisson equation with periodic B.C in all 3 dimensions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine poisson_000(rhs)

    use x3d_operator_x_data
    use x3d_operator_y_data
    use x3d_operator_z_data
    use decomp_2d, only : nx_global, ny_global, nz_global
    use decomp_2d_fft, only : PHYSICAL_IN_Z

    ! right-hand-side of Poisson as input
    ! solution of Poisson as output
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    integer, dimension(3) :: fft_start, fft_end, fft_size

    complex(mytype) :: xyzk

    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


    real(mytype) :: tmp1, tmp2,x ,y, z

    integer :: nx,ny,nz, i,j,k

    complex(mytype) :: cx
    real(mytype) :: rl, iy
    external cx, rl, iy

    nx = nx_global
    ny = ny_global
    nz = nz_global

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z)
       fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)

    do concurrent(k=sp%xst(3):sp%xen(3), j=sp%xst(2):sp%xen(2),i=sp%xst(1):sp%xen(1))
             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * bz(k) + tmp2 * az(k), &
                             tmp2 * bz(k) - tmp1 * az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * by(j) + tmp2 * ay(j), &
                             tmp2 * by(j) - tmp1 * ay(j), kind=mytype)
             if (j > (ny/2+1)) cw1(i,j,k) = -cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * bx(i) + tmp2 * ax(i), &
                             tmp2 * bx(i) - tmp1 * ax(i), kind=mytype)
             if (i > (nx/2+1)) cw1(i,j,k) = -cw1(i,j,k)

             ! Solve Poisson
             tmp1 = real(kxyz(i,j,k), kind=mytype)
             tmp2 = aimag(kxyz(i,j,k))
             ! CANNOT DO A DIVISION BY ZERO
             if ((tmp1 < epsilon).or.(tmp2 < epsilon)) then
                cw1(i,j,k) = zero
             else
                cw1(i,j,k) = cmplx(real(cw1(i,j,k),kind=mytype) / (-tmp1), &
                                  aimag(cw1(i,j,k)) / (-tmp2), kind=mytype)
             end if

             !Print result in spectal space after Poisson
             !     if (abs(out(i,j,k)) > 1.0e-4) then
             !        write(*,*) 'AFTER',i,j,k,out(i,j,k),xyzk
             !     end if

             ! post-processing backward

             ! POST PROCESSING IN Z
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * bz(k) - tmp2 * az(k), &
                            -tmp2 * bz(k) - tmp1 * az(k), kind=mytype)

             ! POST PROCESSING IN Y
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * by(j) + tmp2 * ay(j), &
                             tmp2 * by(j) - tmp1 * ay(j), kind=mytype)
             if (j > (ny/2 + 1)) cw1(i,j,k) = -cw1(i,j,k)

             ! POST PROCESSING IN X
             tmp1 = real(cw1(i,j,k), kind=mytype)
             tmp2 = aimag(cw1(i,j,k))
             cw1(i,j,k) = cmplx(tmp1 * bx(i) + tmp2 * ax(i), &
                            -tmp2 * bx(i) + tmp1 * ax(i), kind=mytype)
             if (i > (nx/2+1)) cw1(i,j,k) = -cw1(i,j,k)
       ! post-processing in spectral space


    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    !   call decomp_2d_fft_finalize

    return
  end subroutine poisson_000


  subroutine poisson_100(rhs)

    use decomp_2d, only : nx_global, ny_global, nz_global
    use decomp_2d_fft, only : PHYSICAL_IN_Z
    
    implicit none

    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    complex(mytype) :: xyzk
    real(mytype) :: tmp1, tmp2, tmp3, tmp4
    real(mytype) :: xx1,xx2,xx3,xx4,xx5,xx6,xx7,xx8

    integer :: nx,ny,nz, i,j,k, itmp

    complex(mytype) :: cx
    real(mytype) :: rl, iy
    external cx, rl, iy

100 format(1x,a8,3I4,2F12.6)

    nx = nx_global - 1
    ny = ny_global
    nz = nz_global

    ! rhs is in Z-pencil but requires global operations in X
    call x3d_transpose_z_to_y(rhs,rw2,ph)
    call x3d_transpose_y_to_x(rw2,rw1,ph)
    do k=ph%xst(3),ph%xen(3)
       do j=ph%xst(2),ph%xen(2)
          do i=1,nx/2
             rw1b(i,j,k)=rw1(2*(i-1)+1,j,k)
          enddo
          do i=nx/2+1,nx
             rw1b(i,j,k)=rw1(2*nx-2*i+2,j,k)
          enddo
       enddo
    end do

    call x3d_transpose_x_to_y(rw1b,rw2,ph)
    call x3d_transpose_y_to_z(rw2,rhs,ph)

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z,nx,ny,nz)
       fft_initialised = .true.
    end if

    ! compute r2c transform 
    call decomp_2d_fft_3d(rhs,cw1)

    ! normalisation
    cw1 = cw1 / real(nx, kind=mytype) /real(ny, kind=mytype) &
         / real(nz, kind=mytype)
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'START',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! post-processing in spectral space

    ! POST PROCESSING IN Z
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          do i = sp%xst(1), sp%xen(1)
             tmp1 = rl(cw1(i,j,k))
             tmp2 = iy(cw1(i,j,k))
             cw1(i,j,k) = cx(tmp1 * bz(k) + tmp2 * az(k), &
                             tmp2 * bz(k) - tmp1 * az(k))
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after z',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Y
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          do i = sp%xst(1), sp%xen(1)
             tmp1 = rl(cw1(i,j,k))
             tmp2 = iy(cw1(i,j,k))
             cw1(i,j,k) = cx(tmp1 * by(j) + tmp2 * ay(j), &
                             tmp2 * by(j) - tmp1 * ay(j))
             if (j > (ny/2+1)) cw1(i,j,k) = -cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'after y',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN X
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          cw1b(1,j,k) = cw1(1,j,k)
          do i = 2, nx
             tmp1 = rl(cw1(i,j,k))
             tmp2 = iy(cw1(i,j,k))
             tmp3 = rl(cw1(nx-i+2,j,k))
             tmp4 = iy(cw1(nx-i+2,j,k))
             xx1=tmp1 * bx(i)
             xx2=tmp1 * ax(i)
             xx3=tmp2 * bx(i)
             xx4=tmp2 * ax(i)
             xx5=tmp3 * bx(i)
             xx6=tmp3 * ax(i)
             xx7=tmp4 * bx(i)
             xx8=tmp4 * ax(i)
             cw1b(i,j,k) = half * cx(xx1 + xx4 + xx5 - xx8, &
                                    -xx2 + xx3 + xx6 + xx7)
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1b(i,j,k)) > 1.0e-4) then
                write(*,100) 'after x',i,j,k,cw1b(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! Solve Poisson
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          do i = sp%xst(1), sp%xen(1)
             tmp1 = rl(kxyz(i,j,k))
             tmp2 = iy(kxyz(i,j,k))
             ! CANNOT DO A DIVISION BY ZERO
             if ((abs(tmp1) < epsilon).and.(abs(tmp2) < epsilon)) then    
                cw1b(i,j,k)=cx(zero, zero)
             end if
             if ((abs(tmp1) < epsilon).and.(abs(tmp2) >= epsilon)) then
                cw1b(i,j,k)=cx(zero, iy(cw1b(i,j,k)) / (-tmp2))
             end if
             if ((abs(tmp1) >= epsilon).and.(abs(tmp2) < epsilon)) then    
                cw1b(i,j,k)=cx(rl(cw1b(i,j,k)) / (-tmp1), zero)
             end if
             if ((abs(tmp1) >= epsilon).and.(abs(tmp2) >= epsilon)) then
                cw1b(i,j,k)=cx(rl(cw1b(i,j,k)) / (-tmp1), iy(cw1b(i,j,k)) / (-tmp2))
             end if
#ifdef DEBUG
             if (abs(cw1b(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER',i,j,k,cw1b(i,j,k)
#endif
          end do
       end do
    end do

    ! post-processing backward

    ! POST PROCESSING IN X
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          cw1(1,j,k) = cw1b(1,j,k)
          do i = 2, nx 
             tmp1 = rl(cw1b(i,j,k))
             tmp2 = iy(cw1b(i,j,k))
             tmp3 = rl(cw1b(nx-i+2,j,k))
             tmp4 = iy(cw1b(nx-i+2,j,k))
             xx1 = tmp1 * bx(i)
             xx2 = tmp1 * ax(i)
             xx3 = tmp2 * bx(i)
             xx4 = tmp2 * ax(i)
             xx5 = tmp3 * bx(i)
             xx6 = tmp3 * ax(i)
             xx7 = tmp4 * bx(i)
             xx8 = tmp4 * ax(i)
             cw1(i,j,k) = cx(xx1-xx4+xx6+xx7, &
                           -(-xx2-xx3+xx5-xx8))
          end do
       end do
    end do
#ifdef DEBUG
    do k = sp%xst(3),sp%xen(3)
       do j = sp%xst(2),sp%xen(2)
          do i = sp%xst(1),sp%xen(1)
             if (abs(cw1(i,j,k)) > 1.0e-4) then
                write(*,100) 'AFTER X',i,j,k,cw1(i,j,k)
             end if
          end do
       end do
    end do
#endif

    ! POST PROCESSING IN Y
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          do i = sp%xst(1), sp%xen(1)
             tmp1 = rl(cw1(i,j,k))
             tmp2 = iy(cw1(i,j,k))
             cw1(i,j,k) = cx(tmp1 * by(j) - tmp2 * ay(j), &
                             tmp2 * by(j) + tmp1 * ay(j))
             if (j > (ny/2+1)) cw1(i,j,k) = -cw1(i,j,k)
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'AFTER Y',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! POST PROCESSING IN Z
    do k = sp%xst(3), sp%xen(3)
       do j = sp%xst(2), sp%xen(2)
          do i = sp%xst(1), sp%xen(1)
             tmp1 = rl(cw1(i,j,k))
             tmp2 = iy(cw1(i,j,k))
             cw1(i,j,k) = cx(tmp1 * bz(k) - tmp2 * az(k), &
                             tmp2 * bz(k) + tmp1 * az(k))
#ifdef DEBUG
             if (abs(cw1(i,j,k)) > 1.0e-4) &
                  write(*,100) 'END',i,j,k,cw1(i,j,k)
#endif
          end do
       end do
    end do

    ! compute c2r transform
    call decomp_2d_fft_3d(cw1,rhs)

    ! rhs is in Z-pencil but requires global operations in X
    call x3d_transpose_z_to_y(rhs,rw2,ph)
    call x3d_transpose_y_to_x(rw2,rw1,ph)
    do k = ph%xst(3), ph%xen(3)
       do j = ph%xst(2), ph%xen(2)
          do i = 1, nx/2
             rw1b(2*i-1,j,k) = rw1(i,j,k)
          enddo
          do i = 1, nx/2
             rw1b(2*i,j,k) = rw1(nx-i+1,j,k)
          enddo
       enddo
    end do
    call x3d_transpose_x_to_y(rw1b,rw2,ph)
    call x3d_transpose_y_to_z(rw2,rhs,ph)

    !  call decomp_2d_fft_finalize

    return
  end subroutine poisson_100

  subroutine abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)

    use param
    use x3dprecision, only : pi, sin_prec, cos_prec

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: bcx,bcy,bcz
    real(mytype), dimension(:), intent(OUT) :: ax,bx
    real(mytype), dimension(:), intent(OUT) :: ay,by
    real(mytype), dimension(:), intent(OUT) :: az,bz

    integer :: i,j,k

    if (bcx==0) then
       do i=1,nx
          ax(i) = sin_prec(real(i-1, kind=mytype)*pi/real(nx, kind=mytype))
          bx(i) = cos_prec(real(i-1, kind=mytype)*pi/real(nx, kind=mytype))
       end do
    else if (bcx==1) then
       do i=1,nx
          ax(i) = sin_prec(real(i-1, kind=mytype)*pi/two/ &
               real(nx, kind=mytype))
          bx(i) = cos_prec(real(i-1, kind=mytype)*pi/two/ &
               real(nx, kind=mytype))
       end do
    end if

    if (bcy==0) then
       do j=1,ny
          ay(j) = sin_prec(real(j-1, kind=mytype)*pi/real(ny, kind=mytype))
          by(j) = cos_prec(real(j-1, kind=mytype)*pi/real(ny, kind=mytype))
       end do
    else if (bcy==1) then
       do j=1,ny
          ay(j) = sin_prec(real(j-1, kind=mytype)*pi/two/ &
               real(ny, kind=mytype))
          by(j) = cos_prec(real(j-1, kind=mytype)*pi/two/ &
               real(ny, kind=mytype))
       end do
    end if

    if (bcz==0) then
       do k=1,nz
          az(k) = sin_prec(real(k-1, kind=mytype)*pi/real(nz, kind=mytype))
          bz(k) = cos_prec(real(k-1, kind=mytype)*pi/real(nz, kind=mytype))
       end do
    else if (bcz==1) then
       do k=1,nz
          az(k) = sin_prec(real(k-1, kind=mytype)*pi/two/ &
               real(nz, kind=mytype))
          bz(k) = cos_prec(real(k-1, kind=mytype)*pi/two/ &
               real(nz, kind=mytype))
       end do
    end if

    return
  end subroutine abxyz

  ! ***********************************************************
  !
  subroutine waves ()
    !
    !***********************************************************

    use x3d_operator_x_data 
    use x3d_operator_y_data
    use x3d_operator_z_data
    use param
    use variables
    use decomp_2d_fft
    use x3dprecision, only: sin_prec, cos_prec, pi, twopi

    implicit none

    integer :: i,j,k
    real(mytype) :: w,wp,w1,w1p 
    complex(mytype) :: xyzk
    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2,tmp1,tmp2,tmp3
    complex(mytype) :: tmp4,tmp5,tmp6

    real(mytype) :: rlexs
    real(mytype) :: rleys
    real(mytype) :: rlezs, iyezs

    real(mytype) :: ytt_rl,xtt_rl,ztt_rl,yt1_rl,xt1_rl,zt1_rl
    real(mytype) :: xtt1_rl,ytt1_rl,ztt1_rl

    complex(mytype) :: cx
    real(mytype) :: rl, iy
    external cx, rl, iy

    xkx = zero
    xk2 = zero
    yky = zero
    yk2 = zero
    zkz = zero
    zk2 = zero

    !WAVE NUMBER IN X
    if (bcx == 0) then
       do i = 1, nx/2 + 1
          w = twopi * (i-1) / nx
          wp = acix6 * two * dx * sin_prec(w * half) + bcix6 * two * dx * sin_prec(three * half * w)
          wp = wp / (one + two * alcaix6 * cos_prec(w))
!
          xkx(i) = cx_one_one * (nx * wp / xlx)
          exs(i) = cx_one_one * (nx * w / xlx)
          xk2(i) = cx_one_one * (nx * wp / xlx)**2
!
       enddo
       do i = nx/2 + 2, nx
          xkx(i) = xkx(nx-i+2)
          exs(i) = exs(nx-i+2)
          xk2(i) = xk2(nx-i+2)
       enddo
    else
       do i = 1, nx
          w = twopi * half * (i-1) / nxm
          wp = acix6 * two * dx * sin_prec(w * half) +(bcix6 * two * dx) * sin_prec(three * half * w)
          wp = wp / (one + two * alcaix6 * cos_prec(w))
!
          xkx(i) = cx_one_one * nxm * wp / xlx
          exs(i) = cx_one_one * nxm * w / xlx
          xk2(i) = cx_one_one * (nxm * wp / xlx)**2
!      
       enddo
       xkx(1) = zero
       exs(1) = zero
       xk2(1) = zero
    endif
!
    !WAVE NUMBER IN Y
    if (bcy == 0) then
       do j = 1, ny/2 + 1
          w = twopi * (j-1) / ny
          wp = aciy6 * two * dy * sin_prec(w * half) + bciy6 * two * dy * sin_prec(three * half * w)
          wp = wp / (one + two * alcaiy6 * cos_prec(w))
!
          if (istret == 0) yky(j) = cx_one_one * (ny * wp / yly)
          if (istret /= 0) yky(j) = cx_one_one * (ny * wp)
          eys(j) = cx_one_one * (ny * w / yly)
          yk2(j) = cx_one_one * (ny * wp / yly)**2
!      
       enddo
       do j = ny/2 + 2, ny
          yky(j) = yky(ny-j+2)
          eys(j) = eys(ny-j+2)
          yk2(j) = yk2(ny-j+2)
       enddo
    else
       do j = 1, ny
          w = twopi * half * (j-1) / nym
          wp = aciy6 * two * dy * sin_prec(w * half) +(bciy6 * two *dy) * sin_prec(three * half * w)
          wp = wp / (one + two * alcaiy6 * cos_prec(w))
!
          if (istret == 0) yky(j) = cx_one_one * (nym * wp / yly)
          if (istret /= 0) yky(j) = cx_one_one * (nym * wp)
          eys(j)=cx_one_one * (nym * w / yly)
          yk2(j)=cx_one_one * (nym * wp / yly)**2
!      
       enddo
       yky(1) = zero
       eys(1) = zero
       yk2(1) = zero
    endif

    !WAVE NUMBER IN Z
    if (bcz == 0) then
       do k = 1, nz/2 + 1
          w = twopi * (k-1) / nz
          wp = aciz6 * two * dz * sin_prec(w * half) + (bciz6 * two * dz) * sin_prec(three * half * w)
          wp = wp / (one + two * alcaiz6 * cos_prec(w))
!
          zkz(k) = cx_one_one * (nz * wp / zlz)
          ezs(k) = cx_one_one * (nz * w / zlz)
          zk2(k) = cx_one_one * (nz * wp / zlz)**2
!
       enddo
    else
       do k= 1, nz/2 + 1
          w = pi * (k-1) / nzm
          w1 = pi * (nzm-k+1) / nzm
          wp = aciz6 * two * dz * sin_prec(w * half)+(bciz6 * two * dz) * sin_prec(three * half * w)
          wp = wp / (one + two * alcaiz6 * cos_prec(w))
          w1p = aciz6 * two * dz * sin_prec(w1 * half) + (bciz6 * two * dz) * sin_prec(three * half * w1)
          w1p = w1p / (one + two * alcaiz6 * cos_prec(w1))
!
          zkz(k) = cx(nzm * wp / zlz, -nzm * w1p / zlz)
          ezs(k) = cx(nzm * w / zlz, nzm * w1 / zlz)
          zk2(k) = cx((nzm * wp / zlz)**2, (nzm * w1p / zlz)**2)
!
       enddo
    endif

    if ((bcx == 0).and.(bcz == 0).and.(bcy /= 0)) then
       do k = sp%yst(3), sp%yen(3)
!
          rlezs = rl(ezs(k)) * dz
!
          do j = sp%yst(2), sp%yen(2)
!
             rleys = rl(eys(j)) * dy
!
             do i = sp%yst(1), sp%yen(1)
!
                rlexs = rl(exs(i)) * dx
!
                xtt_rl = two * &
     (bicix6 * cos_prec(rlexs * onepfive) + cicix6 * cos_prec(rlexs * twopfive) + dicix6 * cos_prec(rlexs * threepfive))
!
                ytt_rl = two * &
     (biciy6 * cos_prec(rleys * onepfive) + ciciy6 * cos_prec(rleys * twopfive) + diciy6 * cos_prec(rleys * threepfive))
!
                ztt_rl = two * &
     (biciz6 * cos_prec(rlezs * onepfive) + ciciz6 * cos_prec(rlezs * twopfive) + diciz6 * cos_prec(rlezs * threepfive))
!
                xtt1_rl = two * aicix6 * cos_prec(rlexs * half)
                ytt1_rl = two * aiciy6 * cos_prec(rleys * half)
                ztt1_rl = two * aiciz6 * cos_prec(rlezs * half)
!
                xt1_rl = one + two * ailcaix6 * cos_prec(rlexs)
                yt1_rl = one + two * ailcaiy6 * cos_prec(rleys)
                zt1_rl = one + two * ailcaiz6 * cos_prec(rlezs)
!
                xt2 = xk2(i) * ((((ytt1_rl + ytt_rl) / yt1_rl) * ((ztt1_rl + ztt_rl) / zt1_rl))**2)
                yt2 = yk2(j) * ((((xtt1_rl + xtt_rl) / xt1_rl) * ((ztt1_rl + ztt_rl) / zt1_rl))**2)
                zt2 = zk2(k) * ((((xtt1_rl + xtt_rl) / xt1_rl) * ((ytt1_rl + ytt_rl) / yt1_rl))**2)
!
                xyzk = xt2 + yt2 + zt2
                kxyz(i,j,k) = xyzk
!
             enddo
          enddo
       enddo

    else
       if (bcz==0) then
          do k = sp%xst(3),sp%xen(3)
!
             rlezs = rl(ezs(k)) * dz
!
             do j = sp%xst(2),sp%xen(2)
!
                rleys = rl(eys(j)) * dy
!
                do i = sp%xst(1),sp%xen(1)
!
                   rlexs = rl(exs(i)) * dx
!
                   xtt_rl = two * &  
  (bicix6 * cos_prec(rlexs * onepfive) + cicix6 * cos_prec(rlexs * twopfive) + dicix6 * cos_prec(rlexs * threepfive))
!
                   ytt_rl = two * &
  (biciy6 * cos_prec(rleys * onepfive) + ciciy6 * cos_prec(rleys * twopfive) + diciy6 * cos_prec(rleys * threepfive))
!
                   ztt_rl = two * &
  (biciz6 * cos_prec(rlezs * onepfive) + ciciz6 * cos_prec(rlezs * twopfive) + diciz6 * cos_prec(rlezs * threepfive))
!
                   xtt1_rl = two * aicix6 * cos_prec(rlexs * half)
                   ytt1_rl = two * aiciy6 * cos_prec(rleys * half)
                   ztt1_rl = two * aiciz6 * cos_prec(rlezs * half)
!
                   xt1_rl = one + two * ailcaix6 * cos_prec(rlexs)
                   yt1_rl = one + two * ailcaiy6 * cos_prec(rleys)
                   zt1_rl = one + two * ailcaiz6 * cos_prec(rlezs)
!
                   xt2 = xk2(i) * ((((ytt1_rl + ytt_rl) / yt1_rl) * ((ztt1_rl + ztt_rl) / zt1_rl))**2)
                   yt2 = yk2(j) * ((((xtt1_rl + xtt_rl) / xt1_rl) * ((ztt1_rl + ztt_rl) / zt1_rl))**2)
                   zt2 = zk2(k) * ((((xtt1_rl + xtt_rl) / xt1_rl) * ((ytt1_rl + ytt_rl) / yt1_rl))**2)
!
                   xyzk = xt2 + yt2 + zt2
                   kxyz(i,j,k) = xyzk
!
                enddo
             enddo
          enddo

       else
          do k = sp%xst(3), sp%xen(3)
!
             rlezs = rl(ezs(k)) * dz
             iyezs = iy(ezs(k)) * dz
!
             do j = sp%xst(2), sp%xen(2)
!
                rleys = rl(eys(j)) * dy
!
                do i = sp%xst(1), sp%xen(1)  
!
                   rlexs = rl(exs(i)) * dx
!
                   xtt_rl = two * &
  (bicix6 * cos_prec(rlexs * onepfive) + cicix6 * cos_prec(rlexs * twopfive) + dicix6 * cos_prec(rlexs * threepfive))
!
                   ytt_rl = two * &
  (biciy6 * cos_prec(rleys * onepfive) + ciciy6 * cos_prec(rleys * twopfive) + diciy6 * cos_prec(rleys * threepfive))
!
                   ztt = two * cx( &
  biciz6 * cos_prec(rlezs * onepfive) + ciciz6 * cos_prec(rlezs * twopfive) + diciz6 * cos_prec(rlezs * threepfive),&
  biciz6 * cos_prec(iyezs * onepfive) + ciciz6 * cos_prec(iyezs * twopfive) + diciz6 * cos_prec(iyezs * threepfive))
!
                   xtt1_rl = two * aicix6 * cos_prec(rlexs * half)
                   ytt1_rl = two * aiciy6 * cos_prec(rleys * half)
!
                   ztt1 = two * cx(aiciz6 * cos_prec(rlezs * half),&
                                   aiciz6 * cos_prec(iyezs * half))
!
                   xt1_rl = one + two * ailcaix6 * cos_prec(rlexs)
                   yt1_rl = one + two * ailcaiy6 * cos_prec(rleys)
!
                   zt1 = cx((one + two * ailcaiz6 * cos_prec(rlezs)),&
                            (one + two * ailcaiz6 * cos_prec(iyezs)))
!
                   tmp1 = cx(rl(ztt1 + ztt) / rl(zt1),&
                             iy(ztt1 + ztt) / iy(zt1))
!
                   tmp2 = cx_one_one * (ytt1_rl + ytt_rl) / yt1_rl
!
                   tmp3 = cx_one_one * (xtt1_rl + xtt_rl) / xt1_rl
!
                   tmp4 = rl(tmp2)**2 * cx(rl(tmp1)**2, iy(tmp1)**2)
!
                   tmp5 = rl(tmp3)**2 * cx(rl(tmp1)**2, iy(tmp1)**2)
!
                   tmp6 = (rl(tmp3) * rl(tmp2))**2 * cx_one_one
!
                   tmp1 = cx(rl(tmp4) * rl(xk2(i)), iy(tmp4) * iy(xk2(i)))
!
                   tmp2 = cx(rl(tmp5) * rl(yk2(j)), iy(tmp5) * iy(yk2(j)))
!
                   tmp3 = rl(tmp6) * zk2(k)
!
                   xyzk = tmp1 + tmp2 + tmp3
                   kxyz(i,j,k) = xyzk
!
                enddo
             enddo
          enddo
!
       endif
    endif

  return
  end subroutine waves

end module decomp_2d_poisson

