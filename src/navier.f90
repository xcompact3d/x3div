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
module navier

  implicit none

  private

  public :: solve_poisson, divergence
  public :: cor_vel
  public :: gradp

contains
  !############################################################################
  !!  SUBROUTINE: solve_poisson
  !!      AUTHOR: Paul Bartholomew
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
  !############################################################################
  SUBROUTINE solve_poisson(pp3, px1, py1, pz1, ux1, uy1, uz1)

    use x3dprecision, only : mytype
    USE decomp_2d, ONLY : xsize, zsize, ph1
    USE decomp_2d_poisson, ONLY : poisson
    USE var, ONLY : nzmsize
    USE param, ONLY : npress
    use nvtx

    implicit none

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzmsize, npress) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    INTEGER :: nlock

    nlock = 1 !! Corresponds to computing div(u*)

    call nvtxStartRange("divergence")
    CALL divergence(pp3(:,:,:,1),ux1,uy1,uz1,nlock)
    call nvtxEndRange
    !
    call nvtxStartRange("poisson")
    CALL poisson(pp3(:,:,:,1))
    call nvtxEndRange
    !
    call nvtxStartRange("gradp")
    CALL gradp(px1,py1,pz1,pp3(:,:,:,1))
    call nvtxEndRange

  END SUBROUTINE solve_poisson
  !############################################################################
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !written by SL 2018
  !############################################################################
  subroutine cor_vel (ux,uy,uz,px,py,pz)

    use x3dprecision, only : mytype
    use decomp_2d, only : xsize
    USE variables
    USE param

    implicit none

    integer :: i,j,k

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: px,py,pz

    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      ux(i,j,k)=ux(i,j,k)-px(i,j,k)
      uy(i,j,k)=uy(i,j,k)-py(i,j,k)
      uz(i,j,k)=uz(i,j,k)-pz(i,j,k)
    enddo

    return
  end subroutine cor_vel
  !############################################################################
  !subroutine DIVERGENCe
  !Calculation of div u* for nlock=1 and of div u^{n+1} for nlock=2
  ! input : ux1,uy1,uz1,ep1 (on velocity mesh)
  ! output : pp3 (on pressure mesh)
  !written by SL 2018
  !############################################################################
  subroutine divergence (pp3,ux1,uy1,uz1,nlock)

    use x3dprecision, only : mytype, real_type
    use param
    use decomp_2d, only : nrank, ph1, ph3, ph4, nproc
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : nx_global, ny_global, nz_global
    use decomp_2d, only : transpose_x_to_y, &
                          transpose_y_to_z, &
                          transpose_z_to_y, &
                          transpose_y_to_x
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, di1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, dipp2, &
         duxydxyp3, uzp3, po3, dipp3, nxmsize, nymsize, nzmsize
    USE MPI

    implicit none


    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzmsize) :: pp3

    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,tmax1,tmoy1

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzmsize

    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))    
       ta1(i,j,k) = ux1(i,j,k)
       tb1(i,j,k) = uy1(i,j,k)
       tc1(i,j,k) = uz1(i,j,k)
   enddo

    !WORK X-PENCILS

    call derxvp(pp1,ta1,di1,sx,cfx6,csx6,cwx6,xsize(1),nxmsize,xsize(2),xsize(3),0)

    call interxvp(pgy1,tb1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)
    call interxvp(pgz1,tc1,di1,sx,cifxp6,cisxp6,ciwxp6,xsize(1),nxmsize,xsize(2),xsize(3),1)

    call transpose_x_to_y(pp1,duxdxp2,ph4)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph4)
    call transpose_x_to_y(pgz1,uzp2,ph4)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)
    call deryvp(duydypi2,uyp2,dipp2,sy,cfy6,csy6,cwy6,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),0)

    !! Compute sum dudx + dvdy !ph1%yst(1):ph1%yen(1),nymsize,ysize(3)
    do concurrent (k=1:ysize(3), j=1:nymsize, i=ph1%yst(1):ph1%yen(1))
      duydypi2(i,j,k) = duydypi2(i,j,k) + upi2(i,j,k)
    enddo

    call interyvp(upi2,uzp2,dipp2,sy,cifyp6,cisyp6,ciwyp6,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nymsize,ysize(3),1)

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,dipp3,sz,cifzp6,ciszp6,ciwzp6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,1)
    call derzvp(po3,uzp3,dipp3,sz,cfz6,csz6,cwz6,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzmsize,0)

    !! Compute sum dudx + dvdy + dwdz
    do concurrent (k=1:nzmsize, j=ph1%zst(2):ph1%zen(2), i=ph1%zst(1):ph1%zen(1))
       pp3(i,j,k) = pp3(i,j,k) + po3(i,j,k)
    enddo

    if (nlock==2) then
       do concurrent (k=1:nzmsize, j=ph1%zst(2):ph1%zen(2),i=ph1%zst(1):ph1%zen(1))
          pp3(i,j,k)=pp3(i,j,k)-pp3(ph1%zst(1),ph1%zst(2),nzmsize)
       enddo
    endif

    tmax=-1609._mytype
    tmoy=zero
    do concurrent (k=1:nzmsize, j=ph1%zst(2):ph1%zen(2),i=ph1%zst(1):ph1%zen(1))
       if (pp3(i,j,k).gt.tmax) tmax=pp3(i,j,k)
       tmoy=tmoy+abs(pp3(i,j,k))
    enddo
    tmoy=tmoy/nvect3


    if (test_mode) then
       call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
       call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
       if ((nrank==0).and.(nlock.gt.0)) then
          if (nlock==2) then
             print *,'DIV U  max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
          else
             print *,'DIV U* max mean=',real(tmax1,4),real(tmoy1/real(nproc),4)
          endif
       endif
    end if

    return
  end subroutine divergence
  !############################################################################
  !subroutine GRADP
  !Computation of the pressure gradient from the pressure mesh to the
  !velocity mesh
  !Saving pressure gradients on boundaries for correct imposition of
  !BCs on u* via the fractional step methodi (it is not possible to
  !impose BC after correction by pressure gradient otherwise lost of
  !incompressibility--> BCs are imposed on u*
  !
  ! input: pp3 - pressure field (on pressure mesh)
  ! output: px1, py1, pz1 - pressure gradients (on velocity mesh)
  !written by SL 2018
  !############################################################################
  subroutine gradp(px1,py1,pz1,pp3)

    use x3dprecision, only: mytype 
    USE param
    USE decomp_2d, only: xsize, ysize, zsize, ph2, ph3
    use decomp_2d, only : transpose_x_to_y, &
                          transpose_y_to_z, &
                          transpose_z_to_y, &
                          transpose_y_to_x
    use decomp_2d, only: xstart, xend, ystart, yend, zstart, zend
    USE variables
    USE var, only: pp1,pgy1,pgz1,di1,pp2,ppi2,pgy2,pgz2,pgzi2,dip2,&
         pgz3,ppi3,dip3,nxmsize,nymsize,nzmsize

    implicit none

    integer :: i,j,k

    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1

    !WORK Z-PENCILS
    call interzpv(ppi3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
    call derzpv(pgz3,pp3,dip3,sz,cfip6z,csip6z,cwip6z,cfz6,csz6,cwz6,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)

    !WORK Y-PENCILS
    call transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
    call transpose_z_to_y(ppi3,pp2,ph3)

    call interypv(ppi2,pp2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call derypv(pgy2,pp2,dip2,sy,cfip6y,csip6y,cwip6y,cfy6,csy6,cwy6,ppy,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
    call interypv(pgzi2,pgz2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
         (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)

    !WORK X-PENCILS

    call transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call transpose_y_to_x(pgy2,pgy1,ph2)
    call transpose_y_to_x(pgzi2,pgz1,ph2)

    call derxpv(px1,pp1,di1,sx,cfip6,csip6,cwip6,cfx6,csx6,cwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(py1,pgy1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)
    call interxpv(pz1,pgz1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
         nxmsize,xsize(1),xsize(2),xsize(3),1)

    !we are in X pencils:
    if (nclx1.eq.2) then
       do concurrent (k=1:xsize(3),j=1:xsize(2))
          dpdyx1(j,k)=py1(1,j,k)/gdt(itr)
          dpdzx1(j,k)=pz1(1,j,k)/gdt(itr)
       enddo
    endif
    if (nclxn.eq.2) then
       do concurrent (k=1:xsize(3),j=1:xsize(2))
          dpdyxn(j,k)=py1(nx,j,k)/gdt(itr)
          dpdzxn(j,k)=pz1(nx,j,k)/gdt(itr)
       enddo
    endif

    if (ncly1.eq.2) then
       if (xstart(2)==1) then
          do concurrent (k=1:xsize(3),i=1:xsize(1))
             dpdxy1(i,k)=px1(i,1,k)/gdt(itr)
             dpdzy1(i,k)=pz1(i,1,k)/gdt(itr)
          enddo
       endif
    endif
    if (nclyn.eq.2) then
       if (xend(2)==ny) then
          do concurrent (k=1:xsize(3),i=1:xsize(1))
             dpdxyn(i,k)=px1(i,xsize(2),k)/gdt(itr)
             dpdzyn(i,k)=pz1(i,xsize(2),k)/gdt(itr)
          enddo
       endif
    endif

    if (nclz1.eq.2) then
       if (xstart(3)==1) then
          do concurrent (j=1:xsize(2),i=1:xsize(1))
             dpdxz1(i,j)=py1(i,j,1)/gdt(itr)
             dpdyz1(i,j)=pz1(i,j,1)/gdt(itr)
          enddo
       endif
    endif
    if (nclzn.eq.2) then
       if (xend(3)==nz) then
          do concurrent (j=1:xsize(2),i=1:xsize(1))
             dpdxzn(i,j)=py1(i,j,xsize(3))/gdt(itr)
             dpdyzn(i,j)=pz1(i,j,xsize(3))/gdt(itr)
          enddo
       endif
    endif

    return
  end subroutine gradp
  !############################################################################
  !############################################################################
endmodule navier
