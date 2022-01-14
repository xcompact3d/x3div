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
  !use decomp_2d, only : DECOMP_INFO
  !use decomp_2d
  !USE variables
  !USE param

  implicit none
  
  ! define all major arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3, po3, dv3
  !double precision, save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3, po3, dv3
  real(mytype), save, allocatable, dimension(:,:,:,:) :: pp3
  real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
  real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
  real(mytype), save, allocatable, dimension(:,:,:) :: divu3
  real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
  real(mytype), save, allocatable, dimension(:,:,:,:) :: dux1,duy1,duz1  ! Output of convdiff


  ! define all work arrays here
  real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
       te1,tf1,tg1,th1,ti1,di1
  real(mytype), save, allocatable, dimension(:,:,:) :: pp1,pgy1,pgz1
  real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
       te2,tf2,tg2,th2,ti2,tj2,di2
  real(mytype), save, allocatable, dimension(:,:,:) :: pp2,ppi2,pgy2,pgz2,pgzi2,dip2,dipp2,duxdxp2,uyp2,uzp2,upi2,duydypi2
  real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
       te3,tf3,tg3,th3,ti3,di3
  real(mytype), save, allocatable, dimension(:,:,:) :: pgz3,ppi3,dip3,dipp3,duxydxyp3,uzp3

  integer, save :: nxmsize, nymsize, nzmsize

contains

  subroutine init_variables

    !use decomp_2d  
    use variables
    use param
    use decomp_2d, only : DECOMP_INFO
    use decomp_2d , only : decomp_info_init 
    use decomp_2d , only : alloc_x, alloc_y, alloc_z 
    use decomp_2d , only : xsize, ysize, zsize, ph1, ph3
    use decomp_2d , only : nrank
    
    implicit none

    TYPE(DECOMP_INFO), save :: ph! decomposition object

    integer :: i, j , k

#ifdef DEBG
    if (nrank == 0) write(*,*) '# Init_variables start'
#endif

    if (nrank == 0) write(*,*) '# Initializing variables...'

    if (nclx) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)
    !xsize(i), ysize(i), zsize(i), i=1,2,3 - sizes of the sub-domains held by the current process. The first letter refers to the pencil orientation and the three 1D array elements contain the sub-domain sizes in X, Y and Z directions, respectively. In a 2D pencil decomposition, there is always one dimension which completely resides in local memory. So by definition xsize(1)==nx_global, ysize(2)==ny_global and zsize(3)==nz_global.

    !xstart(i), ystart(i), zstart(i), xend(i), yend(i), zend(i), i=1,2,3 - the starting and ending indices for each sub-domain, as in the global coordinate system. Obviously, it can be seen that xsize(i)=xend(i)-xstart(i)+1. It may be convenient for certain applications to use global coordinate (for example when extracting a 2D plane from a 3D domain, it is easier to know which process owns the plane if global index is used).


    !X PENCILS
    call alloc_x(ux1, opt_global=.true.) !global indices
    ux1 = zero
    call alloc_x(uy1, opt_global=.true.) !global indices
    uy1 = zero
    call alloc_x(uz1, opt_global=.true.) !global indices
    uz1 = zero
    call alloc_x(px1, opt_global=.true.) !global indices
    px1 = zero
    call alloc_x(py1, opt_global=.true.) !global indices
    py1 = zero
    call alloc_x(pz1, opt_global=.true.) !global indices
    pz1 = zero


    call alloc_x(ta1)
    ta1 = zero
    call alloc_x(tb1)
    tb1 = zero
    call alloc_x(tc1)
    tc1 = zero
    call alloc_x(td1)
    td1 = zero
    call alloc_x(te1)
    te1 = zero
    call alloc_x(tf1)
    tf1 = zero
    call alloc_x(tg1)
    tg1 = zero
    call alloc_x(th1)
    th1 = zero
    call alloc_x(ti1)
    ti1 = zero
    call alloc_x(di1)
    di1 = zero

    allocate(pp1(nxmsize,xsize(2),xsize(3)))
    pp1 = zero
    allocate(pgy1(nxmsize,xsize(2),xsize(3)))
    pgy1 = zero
    allocate(pgz1(nxmsize,xsize(2),xsize(3)))
    pgz1 = zero

    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    dpdyx1=zero
    dpdyxn=zero
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    dpdzx1=zero
    dpdzxn=zero
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    dpdxy1=zero
    dpdxyn=zero
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    dpdzy1=zero
    dpdzyn=zero
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    dpdxz1=zero
    dpdxzn=zero
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))
    dpdyz1=zero
    dpdyzn=zero

    !Y PENCILS
    call alloc_y(ux2)
    ux2=zero
    call alloc_y(uy2)
    uy2=zero
    call alloc_y(uz2)
    uz2=zero
    call alloc_y(ta2)
    ta2=zero
    call alloc_y(tb2)
    tb2=zero
    call alloc_y(tc2)
    tc2=zero
    call alloc_y(td2)
    td2=zero
    call alloc_y(te2)
    te2=zero
    call alloc_y(tf2)
    tf2=zero
    call alloc_y(tg2)
    tg2=zero
    call alloc_y(th2)
    th2=zero
    call alloc_y(ti2)
    ti2=zero
    call alloc_y(tj2)
    tj2=zero
    call alloc_y(di2)
    di2=zero
    allocate(pgz2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    pgz2=zero
    allocate(pp2(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)))
    pp2=zero
    allocate(dip2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    dip2=zero
    allocate(ppi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    ppi2=zero
    allocate(pgy2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    pgy2=zero
    allocate(pgzi2(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)))
    pgzi2=zero
    allocate(duxdxp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    duxdxp2=zero
    allocate(uyp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    uyp2=zero
    allocate(uzp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    uzp2=zero
    allocate(dipp2(ph1%yst(1):ph1%yen(1),ysize(2),ysize(3)))
    dipp2=zero
    allocate(upi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))
    upi2=zero
    allocate(duydypi2(ph1%yst(1):ph1%yen(1),nymsize,ysize(3)))
    duydypi2=zero

    !Z PENCILS
    call alloc_z(ux3)
    ux3=zero
    call alloc_z(uy3)
    uy3=zero
    call alloc_z(uz3)
    uz3=zero
    call alloc_z(ta3)
    ta3=zero
    call alloc_z(tb3)
    tb3=zero
    call alloc_z(tc3)
    tc3=zero
    call alloc_z(td3)
    td3=zero
    call alloc_z(te3)
    te3=zero
    call alloc_z(tf3)
    tf3=zero
    call alloc_z(tg3)
    tg3=zero
    call alloc_z(th3)
    th3=zero
    call alloc_z(ti3)
    ti3=zero
    call alloc_z(di3)
    di3=zero
    allocate(pgz3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    pgz3=zero
    allocate(ppi3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    ppi3=zero
    allocate(dip3(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)))
    dip3=zero

    allocate(duxydxyp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    duxydxyp3=zero
    allocate(uzp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    uzp3=zero
    allocate(dipp3(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),zsize(3)))
    dipp3=zero


    allocate(pp3(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress))
    pp3=zero

    call alloc_z(dv3,ph,.true.)
    dv3=zero
    call alloc_z(po3,ph,.true.)
    po3=zero


    !module derivative
    allocate(ffx(nx))
    ffx=zero
    allocate(sfx(nx))
    sfx=zero
    allocate(fsx(nx))
    fsx=zero
    allocate(fwx(nx))
    fwx=zero
    allocate(ssx(nx))
    ssx=zero
    allocate(swx(nx))
    swx=zero

    allocate(ffxp(nx))
    ffxp=zero
    allocate(sfxp(nx))
    sfxp=zero
    allocate(fsxp(nx))
    fsxp=zero
    allocate(fwxp(nx))
    fwxp=zero
    allocate(ssxp(nx))
    ssxp=zero
    allocate(swxp(nx))
    swxp=zero

    allocate(ffy(ny))
    ffy=zero
    allocate(sfy(ny))
    sfy=zero
    allocate(fsy(ny))
    fsy=zero
    allocate(fwy(ny))
    fwy=zero
    allocate(ssy(ny))
    ssy=zero
    allocate(swy(ny))
    swy=zero

    allocate(ffyp(ny))
    ffyp=zero
    allocate(sfyp(ny))
    sfyp=zero
    allocate(fsyp(ny))
    fsyp=zero
    allocate(fwyp(ny))
    fwyp=zero
    allocate(ssyp(ny))
    ssyp=zero
    allocate(swyp(ny))
    swyp=zero

    allocate(ffz(nz))
    ffz=zero
    allocate(sfz(nz))
    sfz=zero
    allocate(fsz(nz))
    fsz=zero
    allocate(fwz(nz))
    fwz=zero
    allocate(ssz(nz))
    ssz=zero
    allocate(swz(nz))
    swz=zero

    allocate(ffzp(nz))
    ffzp=zero
    allocate(sfzp(nz))
    sfzp=zero
    allocate(fszp(nz))
    fszp=zero
    allocate(fwzp(nz))
    fwzp=zero
    allocate(sszp(nz))
    sszp=zero
    allocate(swzp(nz))
    swzp=zero

    allocate(ffxS(nx))
    ffxS=zero
    allocate(sfxS(nx))
    sfxS=zero
    allocate(fsxS(nx))
    fsxS=zero
    allocate(fwxS(nx))
    fwxS=zero
    allocate(ssxS(nx))
    ssxS=zero
    allocate(swxS(nx))
    swxS=zero

    allocate(ffxpS(nx))
    ffxpS=zero
    allocate(sfxpS(nx))
    sfxpS=zero
    allocate(fsxpS(nx))
    fsxpS=zero
    allocate(fwxpS(nx))
    fwxpS=zero
    allocate(ssxpS(nx))
    ssxpS=zero
    allocate(swxpS(nx))
    swxpS=zero

    allocate(ffyS(ny))
    ffyS=zero
    allocate(sfyS(ny))
    sfyS=zero
    allocate(fsyS(ny))
    fsyS=zero
    allocate(fwyS(ny))
    fwyS=zero
    allocate(ssyS(ny))
    ssyS=zero
    allocate(swyS(ny))
    swyS=zero

    allocate(ffypS(ny))
    ffypS=zero
    allocate(sfypS(ny))
    sfypS=zero
    allocate(fsypS(ny))
    fsypS=zero
    allocate(fwypS(ny))
    fwypS=zero
    allocate(ssypS(ny))
    ssypS=zero
    allocate(swypS(ny))
    swypS=zero

    allocate(ffzS(nz))
    ffzS=zero
    allocate(sfzS(nz))
    sfzS=zero
    allocate(fszS(nz))
    fszS=zero
    allocate(fwzS(nz))
    fwzS=zero
    allocate(sszS(nz))
    sszS=zero
    allocate(swzS(nz))
    swzS=zero

    allocate(ffzpS(nz))
    ffzpS=zero
    allocate(sfzpS(nz))
    sfzpS=zero
    allocate(fszpS(nz))
    fszpS=zero
    allocate(fwzpS(nz))
    fwzpS=zero
    allocate(sszpS(nz))
    sszpS=zero
    allocate(swzpS(nz))
    swzpS=zero

    allocate(sx(xsize(2),xsize(3)))
    sx=zero
    allocate(vx(xsize(2),xsize(3)))
    vx=zero

    allocate(sy(ysize(1),ysize(3)))
    sy=zero
    allocate(vy(ysize(1),ysize(3)))
    vy=zero

    allocate(sz(zsize(1),zsize(2)))
    sz=zero
    allocate(vz(zsize(1),zsize(2)))
    vz=zero


    !module derpres
    allocate(cfx6(nxm))
    cfx6=zero
    allocate(ccx6(nxm))
    ccx6=zero
    allocate(cbx6(nxm))
    cbx6=zero
    allocate(cfxp6(nxm))
    cfxp6=zero
    allocate(ciwxp6(nxm))
    ciwxp6=zero
    allocate(csxp6(nxm))
    csxp6=zero
    allocate(cwxp6(nxm))
    cwxp6=zero
    allocate(csx6(nxm))
    csx6=zero
    allocate(cwx6(nxm))
    cwx6=zero
    allocate(cifx6(nxm))
    cifx6=zero
    allocate(cicx6(nxm))
    cicx6=zero
    allocate(cisx6(nxm))
    cisx6=zero

    allocate(cibx6(nxm))
    cibx6=zero
    allocate(cifxp6(nxm))
    cifxp6=zero
    allocate(cisxp6(nxm))
    cisxp6=zero
    allocate(ciwx6(nxm))
    ciwx6=zero

    allocate(cfi6(nx))
    cfi6=zero
    allocate(cci6(nx))
    cci6=zero
    allocate(cbi6(nx))
    cbi6=zero
    allocate(cfip6(nx))
    cfip6=zero
    allocate(csip6(nx))
    csip6=zero
    allocate(cwip6(nx))
    cwip6=zero
    allocate(csi6(nx))
    csi6=zero
    allocate(cwi6(nx))
    cwi6=zero
    allocate(cifi6(nx))
    cifi6=zero
    allocate(cici6(nx))
    cici6=zero
    allocate(cibi6(nx))
    cibi6=zero
    allocate(cifip6(nx))
    cifip6=zero

    allocate(cisip6(nx))
    cisip6=zero
    allocate(ciwip6(nx))
    ciwip6=zero
    allocate(cisi6(nx))
    cisi6=zero
    allocate(ciwi6(nx))
    ciwi6=zero

    allocate(cfy6(nym))
    cfy6=zero
    allocate(ccy6(nym))
    ccy6=zero
    allocate(cby6(nym))
    cby6=zero
    allocate(cfyp6(nym))
    cfyp6=zero
    allocate(csyp6(nym))
    csyp6=zero
    allocate(cwyp6(nym))
    cwyp6=zero
    allocate(csy6(nym))
    csy6=zero

    allocate(cwy6(nym))
    cwy6=zero
    allocate(cify6(nym))
    cify6=zero
    allocate(cicy6(nym))
    cicy6=zero
    allocate(ciby6(nym))
    ciby6=zero
    allocate(cifyp6(nym))
    cifyp6=zero
    allocate(cisyp6(nym))
    cisyp6=zero
    allocate(ciwyp6(nym))
    ciwyp6=zero
    allocate(cisy6(nym))
    cisy6=zero
    allocate(ciwy6(nym))
    ciwy6=zero

    allocate(cfi6y(ny))
    cfi6y=zero
    allocate(cci6y(ny))
    cci6y=zero
    allocate(cbi6y(ny))
    cbi6y=zero
    allocate(cfip6y(ny))
    cfip6y=zero
    allocate(csip6y(ny))
    csip6y=zero
    allocate(cwip6y(ny))
    cwip6y=zero
    allocate(csi6y(ny))
    csi6y=zero
    allocate(cwi6y(ny))
    cwi6y=zero
    allocate(cifi6y(ny))
    cifi6y=zero
    allocate(cici6y(ny))
    cici6y=zero

    allocate(cibi6y(ny))
    cibi6y=zero
    allocate(cifip6y(ny))
    cifip6y=zero
    allocate(cisip6y(ny))
    cisip6y=zero
    allocate(ciwip6y(ny))
    ciwip6y=zero
    allocate(cisi6y(ny))
    cisi6y=zero
    allocate(ciwi6y(ny))
    ciwi6y=zero

    allocate(cfz6(nzm))
    cfz6=zero
    allocate(ccz6(nzm))
    ccz6=zero
    allocate(cbz6(nzm))
    cbz6=zero
    allocate(cfzp6(nzm))
    cfzp6=zero
    allocate(cszp6(nzm))
    cszp6=zero
    allocate(cwzp6(nzm))
    cwzp6=zero
    allocate(csz6(nzm))
    csz6=zero

    allocate(cwz6(nzm))
    cwz6=zero
    allocate(cifz6(nzm))
    cifz6=zero
    allocate(cicz6(nzm))
    cicz6=zero
    allocate(cibz6(nzm))
    cibz6=zero
    allocate(cifzp6(nzm))
    cifzp6=zero
    allocate(ciszp6(nzm))
    ciszp6=zero
    allocate(ciwzp6(nzm))
    ciwzp6=zero
    allocate(cisz6(nzm))
    cisz6=zero
    allocate(ciwz6(nzm))
    ciwz6=zero

    allocate(cfi6z(nz))
    cfi6z=zero
    allocate(cci6z(nz))
    cci6z=zero
    allocate(cbi6z(nz))
    cbi6z=zero
    allocate(cfip6z(nz))
    cfip6z=zero
    allocate(csip6z(nz))
    csip6z=zero
    allocate(cwip6z(nz))
    cwip6z=zero
    allocate(csi6z(nz))
    csi6z=zero
    allocate(cwi6z(nz))
    cwi6z=zero
    allocate(cifi6z(nz))
    cifi6z=zero
    allocate(cici6z(nz))
    cici6z=zero

    allocate(cibi6z(nz))
    cibi6z=zero
    allocate(cifip6z(nz))
    cifip6z=zero
    allocate(cisip6z(nz))
    cisip6z=zero
    allocate(ciwip6z(nz))
    ciwip6z=zero
    allocate(cisi6z(nz))
    cisi6z=zero
    allocate(ciwi6z(nz))
    ciwi6z=zero

    !module waves
    allocate(zkz(nz/2+1))
    zkz=zero
    allocate(zk2(nz/2+1))
    zk2=zero
    allocate(ezs(nz/2+1))
    ezs=zero

    allocate(yky(ny))
    yky=zero
    allocate(yk2(ny))
    yk2=zero
    allocate(eys(ny))
    eys=zero

    allocate(xkx(nx))
    xkx=zero
    allocate(xk2(nx))
    xk2=zero
    allocate(exs(nx))
    exs=zero

    !module mesh
    allocate(ppy(ny))
    ppy=zero
    allocate(pp2y(ny))
    pp2y=zero
    allocate(pp4y(ny))
    pp4y=zero

    allocate(ppyi(ny))
    ppyi=zero
    allocate(pp2yi(ny))
    pp2yi=zero
    allocate(pp4yi(ny))
    pp4yi=zero

    allocate(xp(nx))
    xp=zero
    allocate(xpi(nx))
    xpi=zero

    allocate(yp(ny))
    yp=zero
    allocate(ypi(ny))
    ypi=zero
    allocate(del(ny))
    del=zero

    allocate(zp(nz))
    zp=zero
    allocate(zpi(nz))
    zpi=zero

    allocate(yeta(ny))
    yeta=zero
    allocate(yetai(ny))
    yetai=zero

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
    allocate(duy1(xsize(1),xsize(2),xsize(3),ntime))
    duy1=zero
    allocate(duz1(xsize(1),xsize(2),xsize(3),ntime))
    duz1=zero

    call alloc_z(divu3, opt_global=.true.) !global indices
    divu3=zero



#ifdef DEBG
    if (nrank ==  0) write(*,*) '# init_variables done'
#endif
    return
  end subroutine init_variables
end module var
