!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module navier

  implicit none

  private

  public :: solve_poisson, divergence
  public :: cor_vel
  public :: gradp

contains
  !############################################################################
  !!  SUBROUTINE: solve_poisson
  !! DESCRIPTION: Takes the intermediate momentum field as input,
  !!              computes div and solves pressure-Poisson equation.
  !############################################################################
  SUBROUTINE solve_poisson(pp3, px1, py1, pz1, ux1, uy1, uz1)

    use decomp_2d, only : mytype
    USE decomp_2d, ONLY : xsize, zsize, ph1
    USE decomp_2d_poisson, ONLY : poisson
    USE variables, ONLY : nzm
    USE param, ONLY : npress

    implicit none

    !! Inputs
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)), INTENT(IN) :: ux1, uy1, uz1

    !! Outputs
    REAL(mytype), DIMENSION(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzm, npress) :: pp3
    REAL(mytype), DIMENSION(xsize(1), xsize(2), xsize(3)) :: px1, py1, pz1

    !! Locals
    INTEGER :: nlock

    nlock = 1 !! Corresponds to computing div(u*)

    CALL divergence(pp3(:,:,:,1),ux1,uy1,uz1,nlock)
    !
    CALL poisson(pp3(:,:,:,1))
    !
    CALL gradp(px1,py1,pz1,pp3(:,:,:,1))

  END SUBROUTINE solve_poisson
  !############################################################################
  !subroutine COR_VEL
  !Correction of u* by the pressure gradient to get a divergence free
  !field
  ! input : px,py,pz
  ! output : ux,uy,uz
  !############################################################################
  subroutine cor_vel (ux,uy,uz,px,py,pz)

    use decomp_2d, only : mytype
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
  !############################################################################
  subroutine divergence (pp3,ux1,uy1,uz1,nlock)

    use x3d_operator_1d
    use x3d_staggered
    use decomp_2d, only : mytype, real_type, decomp_2d_warning
    use param
    use decomp_2d, only : nrank, ph1, ph2, ph3, nproc
    use decomp_2d, only : xsize, ysize, zsize
    use decomp_2d, only : nx_global, ny_global, nz_global
    use decomp_2d, only : transpose_x_to_y, &
                          transpose_y_to_z, &
                          transpose_z_to_y, &
                          transpose_y_to_x
    USE variables
    USE var, ONLY: ta1, tb1, tc1, pp1, pgy1, pgz1, &
         duxdxp2, uyp2, uzp2, duydypi2, upi2, ta2, &
         duxydxyp3, uzp3, po3
    USE MPI

    implicit none


    !X PENCILS NX NY NZ  -->NXM NY NZ
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)),intent(in) :: ux1,uy1,uz1
    !Z PENCILS NXM NYM NZ  -->NXM NYM NZM
    real(mytype),dimension(ph1%zst(1):ph1%zen(1),ph1%zst(2):ph1%zen(2),nzm) :: pp3

    integer :: nvect3,i,j,k,nlock
    integer :: code
    real(mytype) :: tmax,tmoy,tmax1,tmoy1

    nvect3=(ph1%zen(1)-ph1%zst(1)+1)*(ph1%zen(2)-ph1%zst(2)+1)*nzm

    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))    
       ta1(i,j,k) = ux1(i,j,k)
       tb1(i,j,k) = uy1(i,j,k)
       tc1(i,j,k) = uz1(i,j,k)
   enddo

    !WORK X-PENCILS

    
    call derxvp(pp1,ta1,sx,x3d_op_derxvp,xsize(1),nxm,xsize(2),xsize(3))

    call interxvp(pgy1,tb1,sx,x3d_op_intxvp,xsize(1),nxm,xsize(2),xsize(3))
    call interxvp(pgz1,tc1,sx,x3d_op_intxvp,xsize(1),nxm,xsize(2),xsize(3))

    call transpose_x_to_y(pp1,duxdxp2,ph2)!->NXM NY NZ
    call transpose_x_to_y(pgy1,uyp2,ph2)
    call transpose_x_to_y(pgz1,uzp2,ph2)

    !WORK Y-PENCILS
    call interyvp(upi2,duxdxp2,sy,x3d_op_intyvp,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nym,ysize(3))
    call deryvp(duydypi2,uyp2,sy,x3d_op_deryvp,ppyi,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nym,ysize(3))

    !! Compute sum dudx + dvdy !ph1%yst(1):ph1%yen(1),nym,ysize(3)
    do concurrent (k=1:ysize(3), j=1:nym, i=ph1%yst(1):ph1%yen(1))
      duydypi2(i,j,k) = duydypi2(i,j,k) + upi2(i,j,k)
    enddo

    call interyvp(upi2,uzp2,sy,x3d_op_intyvp,(ph1%yen(1)-ph1%yst(1)+1),ysize(2),nym,ysize(3))

    call transpose_y_to_z(duydypi2,duxydxyp3,ph3)!->NXM NYM NZ
    call transpose_y_to_z(upi2,uzp3,ph3)

    !WORK Z-PENCILS
    call interzvp(pp3,duxydxyp3,sz,x3d_op_intzvp,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzm)
    call derzvp(po3,uzp3,sz,x3d_op_derzvp,(ph1%zen(1)-ph1%zst(1)+1),&
         (ph1%zen(2)-ph1%zst(2)+1),zsize(3),nzm)

    !! Compute sum dudx + dvdy + dwdz
    do concurrent (k=1:nzm, j=ph1%zst(2):ph1%zen(2), i=ph1%zst(1):ph1%zen(1))
       pp3(i,j,k) = pp3(i,j,k) + po3(i,j,k)
    enddo

    tmax = maxval(abs(pp3))
    tmoy = sum(abs(pp3)) / nvect3

    call MPI_REDUCE(tmax,tmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_REDUCE")
    call MPI_REDUCE(tmoy,tmoy1,1,real_type,MPI_SUM,0,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_warning(__FILE__, __LINE__, code, "MPI_REDUCE")

    if ((nrank==0).and.(nlock.gt.0)) then
       if (nlock==2) then
          print *,'DIV U  max mean=', tmax1, tmoy1
       else
          print *,'DIV U* max mean=', tmax1, tmoy1
       endif
    endif

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
  !############################################################################
  subroutine gradp(px1,py1,pz1,pp3)

    use x3d_operator_1d
    use x3d_staggered
    use x3d_transpose
    USE param
    USE decomp_2d, only: mytype, xsize, ysize, zsize, ph2, ph3
    use decomp_2d, only: xstart, xend, ystart, yend, zstart, zend
    USE variables
    USE var, only: pp1,pgy1,pgz1,pp2,ppi2,pgy2,pgz2,pgzi2,&
         pgz3,ppi3

    implicit none

    integer :: i,j,k

    real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzm) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1

    !WORK Z-PENCILS
    call interzpv(ppi3,pp3,sz,x3d_op_intzpv,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzm,zsize(3))
    call derzpv(pgz3,pp3,sz,x3d_op_derzpv,&
         (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzm,zsize(3))

    !WORK Y-PENCILS
    call x3d_transpose_z_to_y(pgz3,pgz2,ph3) !nxm nym nz
    call x3d_transpose_z_to_y(ppi3,pp2,ph3)

    call interypv(ppi2,pp2,sy,x3d_op_intypv,&
         (ph3%yen(1)-ph3%yst(1)+1),nym,ysize(2),ysize(3))
    call derypv(pgy2,pp2,sy,x3d_op_derypv,ppy,&
         (ph3%yen(1)-ph3%yst(1)+1),nym,ysize(2),ysize(3))
    call interypv(pgzi2,pgz2,sy,x3d_op_intypv,&
         (ph3%yen(1)-ph3%yst(1)+1),nym,ysize(2),ysize(3))

    !WORK X-PENCILS

    call x3d_transpose_y_to_x(ppi2,pp1,ph2) !nxm ny nz
    call x3d_transpose_y_to_x(pgy2,pgy1,ph2)
    call x3d_transpose_y_to_x(pgzi2,pgz1,ph2)

    call derxpv(px1,pp1,sx,x3d_op_derxpv,nxm,xsize(1),xsize(2),xsize(3))
    call interxpv(py1,pgy1,sx,x3d_op_intxpv,nxm,xsize(1),xsize(2),xsize(3))
    call interxpv(pz1,pgz1,sx,x3d_op_intxpv,nxm,xsize(1),xsize(2),xsize(3))

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
