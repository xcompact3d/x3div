!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module decomp_2d_poisson

  use decomp_2d_constants, only : mytype
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
  real(mytype), save, allocatable, dimension(:),public :: az,bz
  real(mytype), save, allocatable, dimension(:),public :: ay,by
  real(mytype), save, allocatable, dimension(:),public :: ax,bx

  ! wave numbers
  complex(mytype), save, allocatable, dimension(:,:,:), public :: kxyz
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
       use decomp_2d_constants, only : mytype
       real(mytype), dimension(:,:,:), intent(inout) :: rhs
     end subroutine poisson_xxx
  end interface
  procedure (poisson_xxx), pointer :: poisson=>null()

  public :: decomp_2d_poisson_init,decomp_2d_poisson_finalize,poisson
contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialise Poisson solver for given boundary conditions
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine decomp_2d_poisson_init()

    use decomp_2d_mpi, only : nrank
    use decomp_2d, only : nx_global, ny_global, nz_global

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
    ax = zero
    bx = zero
    allocate(ay(ny),by(ny))
    ay = zero
    by = zero
    allocate(az(nz),bz(nz))
    az = zero
    bz = zero
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
       cw1 = zero
       allocate(kxyz, source=cw1)
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       a = zero
       allocate(a2, source=a)
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
       a3 = zero
    else if (bcx==1 .and. bcy==0 .and. bcz==0) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       cw1 = zero
       allocate(cw1b, source=cw1)
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       rw1 = zero
       allocate(rw1b, source=rw1)
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       rw2 = zero
       allocate(kxyz, source=cw1)
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       a = zero
       allocate(a2, source=a)
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
       a3 = zero
    else if (bcx==0 .and. bcy==1 .and. bcz==0) then
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       rw2 = zero
       allocate(rw2b, source=rw2)
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       cw1 = zero
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       cw2 = zero
       allocate(cw22, cw2b, cw2c, kxyz, source=cw2)
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       a = zero
       allocate(a2, source=a)
       allocate(a3(sp%yst(1):sp%yen(1),ny,sp%yst(3):sp%yen(3),5))
       a3 = zero
    else if (bcx==1 .and. bcy==1) then
       allocate(cw1(sp%xst(1):sp%xen(1),sp%xst(2):sp%xen(2), &
            sp%xst(3):sp%xen(3)))
       cw1 = zero
       allocate(cw1b, source=cw1)
       allocate(cw2(sp%yst(1):sp%yen(1),sp%yst(2):sp%yen(2), &
            sp%yst(3):sp%yen(3)))
       cw2 = zero
       allocate(cw22, cw2b, cw2c, source=cw2)
       allocate(rw1(ph%xst(1):ph%xen(1),ph%xst(2):ph%xen(2), &
            ph%xst(3):ph%xen(3)))
       rw1 = zero
       allocate(rw1b, source=rw1)
       allocate(rw2(ph%yst(1):ph%yen(1),ph%yst(2):ph%yen(2), &
            ph%yst(3):ph%yen(3)))
       rw2 = zero
       allocate(rw2b, source=rw2)
       if (bcz==1) then  
          allocate(rw3(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
          rw3 = zero
       end if
       allocate(kxyz, source=cw1)
       allocate(a(sp%yst(1):sp%yen(1),ny/2,sp%yst(3):sp%yen(3),5))
       a = zero
       allocate(a2, source=a)
       allocate(a3(sp%yst(1):sp%yen(1),nym,sp%yst(3):sp%yen(3),5))      
       a3 = zero
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

    nullify(poisson)

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
    use decomp_2d_constants, only : PHYSICAL_IN_Z
    !use nvtx

    ! right-hand-side of Poisson as input
    ! solution of Poisson as output
    real(mytype), dimension(:,:,:), intent(INOUT) :: rhs

    integer, dimension(3) :: fft_start, fft_end, fft_size

    complex(mytype) :: xyzk

    complex(mytype) :: ytt,xtt,ztt,yt1,xt1,yt2,xt2
    complex(mytype) :: xtt1,ytt1,ztt1,zt1,zt2


    real(mytype) :: tmp1, tmp2,x ,y, z

    integer :: nx,ny,nz, i,j,k
    integer :: sp_xst3, sp_xst2, sp_xst1
    integer :: sp_xen3, sp_xen2, sp_xen1

    complex(mytype) :: cw1_tmp
    real(mytype) :: ntot
    integer, dimension(3) :: dim3d

    nx = nx_global
    ny = ny_global
    nz = nz_global
    ntot = real(nx, kind=mytype) &
         * real(ny, kind=mytype) &
         * real(nz, kind=mytype)

    sp_xst1=sp%xst(1);sp_xst2=sp%xst(2);sp_xst3=sp%xst(3);
    sp_xen1=sp%xen(1);sp_xen2=sp%xen(2);sp_xen3=sp%xen(3);

    if (.not. fft_initialised) then
       call decomp_2d_fft_init(PHYSICAL_IN_Z)
       fft_initialised = .true.
    end if
    
    !$acc data create(cw1) present(kxyz,az,bz,ay,by,ax,bx,rhs)
    ! compute r2c transform 
    !call nvtxStartRange("call decomp_fft_r2c")
    call decomp_2d_fft_3d(rhs,cw1)
    !call nvtxEndRange
    !do k = sp%xst(3), sp%xen(3)
    !   do j = sp%xst(2), sp%xen(2)
    !      do i = sp%xst(1), sp%xen(1)
    !         if (abs(cw1(i,j,k)) > 1.0e-4) then
    !            write(*,*) 'After r2c', i, j, k, cw1(i,j,k)
    !         end if
    !      end do
    !   end do
    !end do

    !call nvtxStartRange("call normalisation")
    !$acc kernels default(present)
    do concurrent(k=sp_xst3:sp_xen3, j=sp_xst2:sp_xen2,i=sp_xst1:sp_xen1)
      ! POST PROCESSING IN Z
      cw1_tmp = cw1(i,j,k)
      tmp1 = real(cw1_tmp, kind=mytype)/ntot
      tmp2 = aimag(cw1_tmp)/ntot
      cw1_tmp = cmplx(tmp1 * bz(k) + tmp2 * az(k), &
                      tmp2 * bz(k) - tmp1 * az(k), kind=mytype)
      
      ! POST PROCESSING IN Y
      tmp1 = real(cw1_tmp, kind=mytype)
      tmp2 = aimag(cw1_tmp)
      cw1_tmp = cmplx(tmp1 * by(j) + tmp2 * ay(j), &
                      tmp2 * by(j) - tmp1 * ay(j), kind=mytype)
      if (j > (ny/2+1)) cw1_tmp = -cw1_tmp
      
      ! POST PROCESSING IN X
      tmp1 = real(cw1_tmp, kind=mytype)
      tmp2 = aimag(cw1_tmp)
      cw1_tmp = cmplx(tmp1 * bx(i) + tmp2 * ax(i), &
                      tmp2 * bx(i) - tmp1 * ax(i), kind=mytype)
      if (i > (nx/2+1)) cw1_tmp = -cw1_tmp
      
      ! Solve Poisson
      tmp1 = real(kxyz(i,j,k), kind=mytype)
      tmp2 = aimag(kxyz(i,j,k))
      ! CANNOT DO A DIVISION BY ZERO
      if ((tmp1 < epsilon).or.(tmp2 < epsilon)) then
         cw1_tmp = zero
      else
         cw1_tmp = cmplx(real(cw1_tmp,kind=mytype) / (-tmp1), &
                           aimag(cw1_tmp) / (-tmp2), kind=mytype)
      end if
      
      ! post-processing backward
      
      ! POST PROCESSING IN Z
      tmp1 = real(cw1_tmp, kind=mytype)
      tmp2 = aimag(cw1_tmp)
      cw1_tmp = cmplx(tmp1 * bz(k) - tmp2 * az(k), &
                     -tmp2 * bz(k) - tmp1 * az(k), kind=mytype)
      
      ! POST PROCESSING IN Y
      tmp1 = real(cw1_tmp, kind=mytype)
      tmp2 = aimag(cw1_tmp)
      cw1_tmp = cmplx(tmp1 * by(j) + tmp2 * ay(j), &
                      tmp2 * by(j) - tmp1 * ay(j), kind=mytype)
      if (j > (ny/2 + 1)) cw1_tmp = -cw1_tmp
      
      ! POST PROCESSING IN X
      tmp1 = real(cw1_tmp, kind=mytype)
      tmp2 = aimag(cw1_tmp)
      cw1_tmp = cmplx(tmp1 * bx(i) + tmp2 * ax(i), &
                     -tmp2 * bx(i) + tmp1 * ax(i), kind=mytype)
      if (i > (nx/2+1)) cw1_tmp = -cw1_tmp
      ! post-processing in spectral space
      cw1(i,j,k) = cw1_tmp

    end do
    !$acc end kernels
    !call nvtxEndRange
#ifdef DEBUG
    dim3d = shape(cw1)
    do k = 1, dim3d(3),dim3d(3)/2+1
      do j = 1, dim3d(2),dim3d(2)/2+1
        do i = 1, dim3d(1),dim3d(1)/2+1
          print "(i3,i3,i3,1x,e12.5,1x,e12.5)", i, j, k, real(cw1(i,j,k)),&
                             aimag(cw1(i,j,k))
        enddo
      enddo
    enddo
#endif
    ! compute c2r transform
    !call nvtxStartRange("call decomp_c2r")
    call decomp_2d_fft_3d(cw1,rhs)
    !call nvtxEndRange
    !$acc end data 

    !   call decomp_2d_fft_finalize

    return
  end subroutine poisson_000

  subroutine abxyz(ax,ay,az,bx,by,bz,nx,ny,nz,bcx,bcy,bcz)

    use param
    use x3d_precision, only : pi

    implicit none

    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: bcx,bcy,bcz
    real(mytype), dimension(:), intent(OUT) :: ax,bx
    real(mytype), dimension(:), intent(OUT) :: ay,by
    real(mytype), dimension(:), intent(OUT) :: az,bz

    integer :: i,j,k

    if (bcx==0) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*pi/real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*pi/real(nx, kind=mytype))
       end do
    else if (bcx==1) then
       do i=1,nx
          ax(i) = sin(real(i-1, kind=mytype)*pi/two/ &
               real(nx, kind=mytype))
          bx(i) = cos(real(i-1, kind=mytype)*pi/two/ &
               real(nx, kind=mytype))
       end do
    end if

    if (bcy==0) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*pi/real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*pi/real(ny, kind=mytype))
       end do
    else if (bcy==1) then
       do j=1,ny
          ay(j) = sin(real(j-1, kind=mytype)*pi/two/ &
               real(ny, kind=mytype))
          by(j) = cos(real(j-1, kind=mytype)*pi/two/ &
               real(ny, kind=mytype))
       end do
    end if

    if (bcz==0) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*pi/real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*pi/real(nz, kind=mytype))
       end do
    else if (bcz==1) then
       do k=1,nz
          az(k) = sin(real(k-1, kind=mytype)*pi/two/ &
               real(nz, kind=mytype))
          bz(k) = cos(real(k-1, kind=mytype)*pi/two/ &
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
    use x3d_precision, only: pi, twopi

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
          wp = acix6 * two * dx * sin(w * half) + bcix6 * two * dx * sin(three * half * w)
          wp = wp / (one + two * alcaix6 * cos(w))
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
          wp = acix6 * two * dx * sin(w * half) +(bcix6 * two * dx) * sin(three * half * w)
          wp = wp / (one + two * alcaix6 * cos(w))
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
          wp = aciy6 * two * dy * sin(w * half) + bciy6 * two * dy * sin(three * half * w)
          wp = wp / (one + two * alcaiy6 * cos(w))
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
          wp = aciy6 * two * dy * sin(w * half) +(bciy6 * two *dy) * sin(three * half * w)
          wp = wp / (one + two * alcaiy6 * cos(w))
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
          wp = aciz6 * two * dz * sin(w * half) + (bciz6 * two * dz) * sin(three * half * w)
          wp = wp / (one + two * alcaiz6 * cos(w))
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
          wp = aciz6 * two * dz * sin(w * half)+(bciz6 * two * dz) * sin(three * half * w)
          wp = wp / (one + two * alcaiz6 * cos(w))
          w1p = aciz6 * two * dz * sin(w1 * half) + (bciz6 * two * dz) * sin(three * half * w1)
          w1p = w1p / (one + two * alcaiz6 * cos(w1))
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
     (bicix6 * cos(rlexs * onepfive) + cicix6 * cos(rlexs * twopfive) + dicix6 * cos(rlexs * threepfive))
!
                ytt_rl = two * &
     (biciy6 * cos(rleys * onepfive) + ciciy6 * cos(rleys * twopfive) + diciy6 * cos(rleys * threepfive))
!
                ztt_rl = two * &
     (biciz6 * cos(rlezs * onepfive) + ciciz6 * cos(rlezs * twopfive) + diciz6 * cos(rlezs * threepfive))
!
                xtt1_rl = two * aicix6 * cos(rlexs * half)
                ytt1_rl = two * aiciy6 * cos(rleys * half)
                ztt1_rl = two * aiciz6 * cos(rlezs * half)
!
                xt1_rl = one + two * ailcaix6 * cos(rlexs)
                yt1_rl = one + two * ailcaiy6 * cos(rleys)
                zt1_rl = one + two * ailcaiz6 * cos(rlezs)
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
  (bicix6 * cos(rlexs * onepfive) + cicix6 * cos(rlexs * twopfive) + dicix6 * cos(rlexs * threepfive))
!
                   ytt_rl = two * &
  (biciy6 * cos(rleys * onepfive) + ciciy6 * cos(rleys * twopfive) + diciy6 * cos(rleys * threepfive))
!
                   ztt_rl = two * &
  (biciz6 * cos(rlezs * onepfive) + ciciz6 * cos(rlezs * twopfive) + diciz6 * cos(rlezs * threepfive))
!
                   xtt1_rl = two * aicix6 * cos(rlexs * half)
                   ytt1_rl = two * aiciy6 * cos(rleys * half)
                   ztt1_rl = two * aiciz6 * cos(rlezs * half)
!
                   xt1_rl = one + two * ailcaix6 * cos(rlexs)
                   yt1_rl = one + two * ailcaiy6 * cos(rleys)
                   zt1_rl = one + two * ailcaiz6 * cos(rlezs)
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
  (bicix6 * cos(rlexs * onepfive) + cicix6 * cos(rlexs * twopfive) + dicix6 * cos(rlexs * threepfive))
!
                   ytt_rl = two * &
  (biciy6 * cos(rleys * onepfive) + ciciy6 * cos(rleys * twopfive) + diciy6 * cos(rleys * threepfive))
!
                   ztt = two * cx( &
  biciz6 * cos(rlezs * onepfive) + ciciz6 * cos(rlezs * twopfive) + diciz6 * cos(rlezs * threepfive),&
  biciz6 * cos(iyezs * onepfive) + ciciz6 * cos(iyezs * twopfive) + diciz6 * cos(iyezs * threepfive))
!
                   xtt1_rl = two * aicix6 * cos(rlexs * half)
                   ytt1_rl = two * aiciy6 * cos(rleys * half)
!
                   ztt1 = two * cx(aiciz6 * cos(rlezs * half),&
                                   aiciz6 * cos(iyezs * half))
!
                   xt1_rl = one + two * ailcaix6 * cos(rlexs)
                   yt1_rl = one + two * ailcaiy6 * cos(rleys)
!
                   zt1 = cx((one + two * ailcaiz6 * cos(rlezs)),&
                            (one + two * ailcaiz6 * cos(iyezs)))
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

