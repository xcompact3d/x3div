!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module var

  use decomp_2d, only : mytype

  implicit none

  ! define all major arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3, po3, dv3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: pp3
  real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
  real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
  real(mytype), save, allocatable, dimension(:,:,:) :: divu3
  real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
  real(mytype), save, allocatable, dimension(:,:,:,:) :: dux1,duy1,duz1  ! Output of convdiff


  ! define all work arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
       te1,tf1,tg1,th1,ti1
  real(mytype), save, allocatable, dimension(:,:,:) :: pp1,pgy1,pgz1
  real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
       te2,tf2,tg2,th2,ti2,tj2
  real(mytype), save, allocatable, dimension(:,:,:) :: pp2,ppi2,pgy2,pgz2,pgzi2,duxdxp2,uyp2,uzp2,upi2,duydypi2
  real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
       te3,tf3,tg3,th3,ti3
  real(mytype), save, allocatable, dimension(:,:,:) :: pgz3,ppi3,duxydxyp3,uzp3


  interface var_zero
     module procedure var1_zero
     module procedure cvar1_zero
     module procedure var2_zero
     module procedure var3_zero
     module procedure var4_zero
  end interface var_zero
 
contains


  !
  ! Allocate memory and initialize arrays
  !
  subroutine var_init()

    use variables
    use param
    use decomp_2d, only : DECOMP_INFO
    use decomp_2d , only : alloc_x, alloc_y, alloc_z
    use decomp_2d , only : xsize, ysize, zsize, ph1, ph3
    use decomp_2d , only : nrank

    implicit none

    integer :: i, j , k

#ifdef DEBUG
    if (nrank == 0) write(*,*) '# var_init start'
#endif

    if (nrank == 0) write(*,*) '# Initializing variables...'

    !xsize(i), ysize(i), zsize(i), i=1,2,3 - sizes of the sub-domains held by the current process. The first letter refers to the pencil orientation and the three 1D array elements contain the sub-domain sizes in X, Y and Z directions, respectively. In a 2D pencil decomposition, there is always one dimension which completely resides in local memory. So by definition xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global.

    !xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 - the starting and ending indices for each sub-domain, as in the global coordinate system. Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. It may be convenient for certain applications to use global coordinate (for example when extracting a 2D plane from a 3D domain, it is easier to know which process owns the plane if global index is used).


    !X PENCILS
    call alloc_x(ux1)
    call var_zero(ux1)
    allocate(uy1, uz1, px1, py1, pz1, source=ux1)

    call alloc_x(ta1)
    call var_zero(ta1)
    allocate(tb1, tc1, td1, te1, tf1, tg1, th1, ti1, source=ta1)

    allocate(pp1(nxm,xsize(2),xsize(3)))
    call var_zero(pp1)
    allocate(pgy1, source=pp1)
    allocate(pgz1, source=pp1)

    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)))
    call var_zero(dpdyx1)
    allocate(dpdyxn, dpdzx1, dpdzxn, source=dpdyx1)
    allocate(dpdxy1(xsize(1),xsize(3)))
    call var_zero(dpdxy1)
    allocate(dpdxyn, dpdzy1, dpdzyn, source=dpdxy1)
    allocate(dpdxz1(xsize(1),xsize(2)))
    call var_zero(dpdxz1)
    allocate(dpdxzn, dpdyz1, dpdyzn, source=dpdxz1)

    !Y PENCILS
    call alloc_y(ux2)
    call var_zero(ux2)
    allocate(uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, tj2, source=ux2)
    allocate(pp2(ph3%yst(1):ph3%yen(1),nym,ysize(3)))
    call var_zero(pp2)
    allocate(pgz2, source=pp2)
    allocate(ppi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    call var_zero(ppi2)
    allocate(pgy2, pgzi2, source=ppi2)
    allocate(duxdxp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    call var_zero(duxdxp2)
    allocate(uyp2, uzp2, source=duxdxp2)
    allocate(upi2(ph1%yst(1):ph1%yen(1),nym,ysize(3)))
    call var_zero(upi2)
    allocate(duydypi2, source=upi2)

    !Z PENCILS
    call alloc_z(ux3)
    call var_zero(ux3)
    allocate(uy3, uz3, ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, source=ux3)
    allocate(ppi3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    call var_zero(ppi3)
    allocate(pgz3, source=ppi3)

    allocate(duxydxyp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    call var_zero(duxydxyp3)
    allocate(uzp3, source=duxydxyp3)

    allocate(pp3(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzm, npress))
    call var_zero(pp3)

    call alloc_z(dv3, ph1, .true.)
    call var_zero(dv3)
    allocate(po3, source=dv3)

    !module waves
    allocate(zkz(nz/2+1))
    call var_zero(zkz)
    allocate(zk2, ezs, source=zkz)

    allocate(yky(ny))
    call var_zero(yky)
    allocate(yk2, eys, source=yky)

    allocate(xkx(nx))
    call var_zero(xkx)
    allocate(xk2, exs, source=xkx)

    !module mesh
    allocate(ppy(ny))
    call var_zero(ppy)
    allocate(pp2y, pp4y, ppyi, pp2yi, pp4yi, source=ppy)

    allocate(xp(nx))
    call var_zero(xp)
    allocate(xpi, source=xp)

    allocate(yp, ypi, del, source=ppy)

    allocate(zp(nz))
    call var_zero(zp)
    allocate(zpi, source=zp)

    allocate(yeta(ny))
    call var_zero(yeta)
    allocate(yetai, source=yeta)

    ! x-position
    do i=1,nx
      xp(i)=real(i-1,mytype)*dx
      xpi(i)=(real(i,mytype)-half)*dx
    enddo
    ! y-position
    if (istret.eq.0) then
       do j=1,ny
          yp(j)=real(j-1,mytype)*dy
          ypi(j)=(real(j,mytype)-half)*dy
       enddo
    endif
    ! z-position
    do k=1,nz
       zp(k)=real(k-1,mytype)*dz
       zpi(k)=(real(k,mytype)-half)*dz
    enddo
    !
    call var_zero(adt)
    call var_zero(bdt)
    call var_zero(cdt)
    call var_zero(gdt)

    if (itimescheme.eq.1) then ! Euler

       iadvance_time=1

       adt(1)=one*dt
       bdt(1)=zero
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 1
    endif
    allocate(dux1(xsize(1),xsize(2),xsize(3),ntime))
    call var_zero(dux1)
    allocate(duy1, duz1, source=dux1)

    call alloc_z(divu3, opt_global=.true.) !global indices
    divu3=zero

#ifdef DEBUG
    if (nrank ==  0) write(*,*) '# var_init done'
#endif

  end subroutine var_init

  !
  ! Zero a 1D array in parallel (ensure any first touch initialisation is performed)
  !
  subroutine var1_zero(v)

    use param, only : zero

    implicit none
    
    real(mytype), dimension(:), intent(inout) :: v

    integer :: i
    integer :: ni

    ni = size(v, 1)

    do concurrent (i = 1:ni)
       v(i) = zero
    end do
    
  end subroutine var1_zero
  subroutine cvar1_zero(v)

    use param, only : zero

    implicit none

    complex(mytype), dimension(:), intent(inout) :: v

    integer :: i
    integer :: ni

    ni = size(v, 1)

    do concurrent (i = 1:ni)
       v(i) = cmplx(zero, zero, kind=mytype)
    end do

  end subroutine cvar1_zero
  !
  ! Zero a 2D array in parallel (ensure any first touch initialisation is performed)
  !
  subroutine var2_zero(v)

    use param, only : zero

    implicit none
    
    real(mytype), dimension(:,:), intent(inout) :: v

    integer :: i, j
    integer :: ni, nj

    nj = size(v, 2)
    ni = size(v, 1)

    do concurrent (j = 1:nj, i = 1:ni)
       v(i, j) = zero
    end do
    
  end subroutine var2_zero
  !
  ! Zero a 3D array in parallel (ensure any first touch initialisation is performed)
  !
  subroutine var3_zero(v)

    use param, only : zero

    implicit none
    
    real(mytype), dimension(:,:,:), intent(inout) :: v

    integer :: i, j, k
    integer :: ni, nj, nk

    nk = size(v, 3)
    nj = size(v, 2)
    ni = size(v, 1)

    do concurrent (k = 1:nk, j = 1:nj, i = 1:ni)
       v(i, j, k) = zero
    end do
    
  end subroutine var3_zero
  !
  ! Zero a 4D array in parallel (ensure any first touch initialisation is performed)
  !
  subroutine var4_zero(v)

    use param, only : zero

    implicit none

    real(mytype), dimension(:,:,:,:), intent(inout) :: v

    integer :: i, j, k, l
    integer :: ni, nj, nk, nl

    nl = size(v, 4)
    nk = size(v, 3)
    nj = size(v, 2)
    ni = size(v, 1)

    do concurrent (l = 1:nl, k = 1:nk, j = 1:nj, i = 1:ni)
       v(i, j, k, l) = zero
    end do

  end subroutine var4_zero

  !
  ! Free memory
  !
  subroutine var_finalize()

    use variables

    implicit none

    !X PENCILS
    deallocate(ux1)
    deallocate(uy1)
    deallocate(uz1)
    deallocate(px1)
    deallocate(py1)
    deallocate(pz1)

    deallocate(ta1)
    deallocate(tb1)
    deallocate(tc1)
    deallocate(td1)
    deallocate(te1)
    deallocate(tf1)
    deallocate(tg1)
    deallocate(th1)
    deallocate(ti1)

    deallocate(pp1)
    deallocate(pgy1)
    deallocate(pgz1)

    !pre_correc 2d array
    deallocate(dpdyx1,dpdyxn)
    deallocate(dpdzx1,dpdzxn)
    deallocate(dpdxy1,dpdxyn)
    deallocate(dpdzy1,dpdzyn)
    deallocate(dpdxz1,dpdxzn)
    deallocate(dpdyz1,dpdyzn)

    !Y PENCILS
    deallocate(ux2)
    deallocate(uy2)
    deallocate(uz2)
    deallocate(ta2)
    deallocate(tb2)
    deallocate(tc2)
    deallocate(td2)
    deallocate(te2)
    deallocate(tf2)
    deallocate(tg2)
    deallocate(th2)
    deallocate(ti2)
    deallocate(tj2)
    deallocate(pgz2)
    deallocate(pp2)
    deallocate(ppi2)
    deallocate(pgy2)
    deallocate(pgzi2)
    deallocate(duxdxp2)
    deallocate(uyp2)
    deallocate(uzp2)
    deallocate(upi2)
    deallocate(duydypi2)

    !Z PENCILS
    deallocate(ux3)
    deallocate(uy3)
    deallocate(uz3)
    deallocate(ta3)
    deallocate(tb3)
    deallocate(tc3)
    deallocate(td3)
    deallocate(te3)
    deallocate(tf3)
    deallocate(tg3)
    deallocate(th3)
    deallocate(ti3)
    deallocate(pgz3)
    deallocate(ppi3)

    deallocate(duxydxyp3)
    deallocate(uzp3)

    deallocate(pp3)

    deallocate(dv3)
    deallocate(po3)

    !module waves
    deallocate(zkz)
    deallocate(zk2)
    deallocate(ezs)

    deallocate(yky)
    deallocate(yk2)
    deallocate(eys)

    deallocate(xkx)
    deallocate(xk2)
    deallocate(exs)

    !module mesh
    deallocate(ppy)
    deallocate(pp2y)
    deallocate(pp4y)

    deallocate(ppyi)
    deallocate(pp2yi)
    deallocate(pp4yi)

    deallocate(xp)
    deallocate(xpi)

    deallocate(yp)
    deallocate(ypi)
    deallocate(del)

    deallocate(zp)
    deallocate(zpi)

    deallocate(yeta)
    deallocate(yetai)

    deallocate(dux1)
    deallocate(duy1)
    deallocate(duz1)

    deallocate(divu3)

  end subroutine var_finalize

end module var
