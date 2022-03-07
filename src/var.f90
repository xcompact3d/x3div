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
    call alloc_x(ux1, opt_global=.true.) !global indices
    ux1 = zero
    allocate(uy1, uz1, px1, py1, pz1, source=ux1)

    call alloc_x(ta1)
    ta1 = zero
    allocate(tb1, tc1, td1, te1, tf1, tg1, th1, ti1, source=ta1)

    allocate(pp1(nxm,xsize(2),xsize(3)))
    pp1 = zero
    allocate(pgy1, source=pp1)
    allocate(pgz1, source=pp1)

    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)))
    dpdyx1=zero
    allocate(dpdyxn, dpdzx1, dpdzxn, source=dpdyx1)
    allocate(dpdxy1(xsize(1),xsize(3)))
    dpdxy1=zero
    allocate(dpdxyn, dpdzy1, dpdzyn, source=dpdxy1)
    allocate(dpdxz1(xsize(1),xsize(2)))
    dpdxz1=zero
    allocate(dpdxzn, dpdyz1, dpdyzn, source=dpdxz1)

    !Y PENCILS
    call alloc_y(ux2)
    ux2=zero
    allocate(uy2, uz2, ta2, tb2, tc2, td2, te2, tf2, tg2, th2, ti2, tj2, source=ux2)
    allocate(pp2(ph3%yst(1):ph3%yen(1),nym,ysize(3)))
    pp2=zero
    allocate(pgz2, source=pp2)
    allocate(ppi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    ppi2=zero
    allocate(pgy2, pgzi2, source=ppi2)
    allocate(duxdxp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    duxdxp2=zero
    allocate(uyp2, uzp2, source=duxdxp2)
    allocate(upi2(ph1%yst(1):ph1%yen(1),nym,ysize(3)))
    upi2=zero
    allocate(duydypi2, source=upi2)

    !Z PENCILS
    call alloc_z(ux3)
    ux3=zero
    allocate(uy3, uz3, ta3, tb3, tc3, td3, te3, tf3, tg3, th3, ti3, source=ux3)
    allocate(ppi3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    ppi3=zero
    allocate(pgz3, source=ppi3)

    allocate(duxydxyp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    duxydxyp3=zero
    allocate(uzp3, source=duxydxyp3)

    allocate(pp3(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzm, npress))
    pp3=zero

    call alloc_z(dv3, ph1, .true.)
    dv3=zero
    allocate(po3, source=dv3)


    !module derivative
    allocate(sx(xsize(2),xsize(3)))
    sx=zero
    allocate(vx, source=sx)

    allocate(sy(ysize(1),ysize(3)))
    sy=zero
    allocate(vy, source=sy)

    allocate(sz(zsize(1),zsize(2)))
    sz=zero
    allocate(vz, source=sz)

    !module waves
    allocate(zkz(nz/2+1))
    zkz=zero
    allocate(zk2, ezs, source=zkz)

    allocate(yky(ny))
    yky=zero
    allocate(yk2, eys, source=yky)

    allocate(xkx(nx))
    xkx=zero
    allocate(xk2, exs, source=xkx)

    !module mesh
    allocate(ppy(ny))
    ppy=zero
    allocate(pp2y, pp4y, ppyi, pp2yi, pp4yi, source=ppy)

    allocate(xp(nx))
    xp=zero
    allocate(xpi, source=xp)

    allocate(yp, ypi, del, source=ppy)

    allocate(zp(nz))
    zp=zero
    allocate(zpi, source=zp)

    allocate(yeta(ny))
    yeta=zero
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
    adt=zero
    bdt=zero
    cdt=zero
    gdt=zero

    if (itimescheme.eq.1) then ! Euler

       iadvance_time=1

       adt(1)=one*dt
       bdt(1)=zero
       gdt(1)=adt(1)+bdt(1)
       gdt(3)=gdt(1)

       ntime = 1
    endif
    allocate(dux1(xsize(1),xsize(2),xsize(3),ntime))
    dux1=zero
    allocate(duy1, duz1, source=dux1)

    call alloc_z(divu3, opt_global=.true.) !global indices
    divu3=zero

#ifdef DEBUG
    if (nrank ==  0) write(*,*) '# var_init done'
#endif

  end subroutine var_init

  !
  ! Free memory
  !
  subroutine var_finalize()

    use variables
    use param

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

    !module derivative
    deallocate(sx)
    deallocate(vx)

    deallocate(sy)
    deallocate(vy)

    deallocate(sz)
    deallocate(vz)

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
