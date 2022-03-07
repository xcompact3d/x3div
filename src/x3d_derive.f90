!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_derive

  use decomp_2d, only : mytype
  use x3d_operator_1d, only : x3doperator1d
  use param
  use thomas

  implicit none

  ! Make everything public unless declared private
  public

  ABSTRACT INTERFACE
     SUBROUTINE DERIVATIVE_X(t,u,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_X
     SUBROUTINE DERIVATIVE_Y(t,u,x3dop,pp,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(in), dimension(ny):: pp
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_Y
     SUBROUTINE DERIVATIVE_YY(t,u,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_YY
     SUBROUTINE DERIVATIVE_Z(t,u,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_Z
  END INTERFACE

  PROCEDURE (DERIVATIVE_X), POINTER :: derx=>null(), derxx=>null(), &
                                       derxS=>null(), derxxS=>null()
  PROCEDURE (DERIVATIVE_Y), POINTER :: dery=>null(), deryS=>null()
  PROCEDURE (DERIVATIVE_YY), POINTER :: deryy=>null(), deryyS=>null()
  PROCEDURE (DERIVATIVE_Z), POINTER :: derz=>null(), derzz=>null(), &
                                       derzS=>null(), derzzS=>null()

contains

  !
  ! Associate pointers with subroutines
  !
  subroutine x3d_derive_init()

    implicit none

    ! Velocity
    ! First derivative
    if (nclx1.eq.0.and.nclxn.eq.0) derx => derx_00
    if (nclx1.eq.1.and.nclxn.eq.1) derx => derx_11
    if (nclx1.eq.1.and.nclxn.eq.2) derx => derx_12
    if (nclx1.eq.2.and.nclxn.eq.1) derx => derx_21
    if (nclx1.eq.2.and.nclxn.eq.2) derx => derx_22
    !
    if (ncly1.eq.0.and.nclyn.eq.0) dery => dery_00
    if (ncly1.eq.1.and.nclyn.eq.1) dery => dery_11
    if (ncly1.eq.1.and.nclyn.eq.2) dery => dery_12
    if (ncly1.eq.2.and.nclyn.eq.1) dery => dery_21
    if (ncly1.eq.2.and.nclyn.eq.2) dery => dery_22
    !
    if (nclz1.eq.0.and.nclzn.eq.0) derz => derz_00
    if (nclz1.eq.1.and.nclzn.eq.1) derz => derz_11
    if (nclz1.eq.1.and.nclzn.eq.2) derz => derz_12
    if (nclz1.eq.2.and.nclzn.eq.1) derz => derz_21
    if (nclz1.eq.2.and.nclzn.eq.2) derz => derz_22
    ! Second derivative
    if (nclx1.eq.0.and.nclxn.eq.0) derxx => derxx_00
    if (nclx1.eq.1.and.nclxn.eq.1) derxx => derxx_11
    if (nclx1.eq.1.and.nclxn.eq.2) derxx => derxx_12
    if (nclx1.eq.2.and.nclxn.eq.1) derxx => derxx_21
    if (nclx1.eq.2.and.nclxn.eq.2) derxx => derxx_22
    !
    if (ncly1.eq.0.and.nclyn.eq.0) deryy => deryy_00
    if (ncly1.eq.1.and.nclyn.eq.1) deryy => deryy_11
    if (ncly1.eq.1.and.nclyn.eq.2) deryy => deryy_12
    if (ncly1.eq.2.and.nclyn.eq.1) deryy => deryy_21
    if (ncly1.eq.2.and.nclyn.eq.2) deryy => deryy_22
    !
    if (nclz1.eq.0.and.nclzn.eq.0) derzz => derzz_00
    if (nclz1.eq.1.and.nclzn.eq.1) derzz => derzz_11
    if (nclz1.eq.1.and.nclzn.eq.2) derzz => derzz_12
    if (nclz1.eq.2.and.nclzn.eq.1) derzz => derzz_21
    if (nclz1.eq.2.and.nclzn.eq.2) derzz => derzz_22

    ! Scalars
    if (iscalar.ne.0) then
      ! First derivative
      if (nclxS1.eq.0.and.nclxSn.eq.0) derxS => derx_00
      if (nclxS1.eq.1.and.nclxSn.eq.1) derxS => derx_11
      if (nclxS1.eq.1.and.nclxSn.eq.2) derxS => derx_12
      if (nclxS1.eq.2.and.nclxSn.eq.1) derxS => derx_21
      if (nclxS1.eq.2.and.nclxSn.eq.2) derxS => derx_22
      !
      if (nclyS1.eq.0.and.nclySn.eq.0) deryS => dery_00
      if (nclyS1.eq.1.and.nclySn.eq.1) deryS => dery_11
      if (nclyS1.eq.1.and.nclySn.eq.2) deryS => dery_12
      if (nclyS1.eq.2.and.nclySn.eq.1) deryS => dery_21
      if (nclyS1.eq.2.and.nclySn.eq.2) deryS => dery_22
      !
      if (nclzS1.eq.0.and.nclzSn.eq.0) derzS => derz_00
      if (nclzS1.eq.1.and.nclzSn.eq.1) derzS => derz_11
      if (nclzS1.eq.1.and.nclzSn.eq.2) derzS => derz_12
      if (nclzS1.eq.2.and.nclzSn.eq.1) derzS => derz_21
      if (nclzS1.eq.2.and.nclzSn.eq.2) derzS => derz_22
      ! Second derivative
      if (nclxS1.eq.0.and.nclxSn.eq.0) derxxS => derxx_00
      if (nclxS1.eq.1.and.nclxSn.eq.1) derxxS => derxx_11
      if (nclxS1.eq.1.and.nclxSn.eq.2) derxxS => derxx_12
      if (nclxS1.eq.2.and.nclxSn.eq.1) derxxS => derxx_21
      if (nclxS1.eq.2.and.nclxSn.eq.2) derxxS => derxx_22
      !
      if (nclyS1.eq.0.and.nclySn.eq.0) deryyS => deryy_00
      if (nclyS1.eq.1.and.nclySn.eq.1) deryyS => deryy_11
      if (nclyS1.eq.1.and.nclySn.eq.2) deryyS => deryy_12
      if (nclyS1.eq.2.and.nclySn.eq.1) deryyS => deryy_21
      if (nclyS1.eq.2.and.nclySn.eq.2) deryyS => deryy_22
      !
      if (nclzS1.eq.0.and.nclzSn.eq.0) derzzS => derzz_00
      if (nclzS1.eq.1.and.nclzSn.eq.1) derzzS => derzz_11
      if (nclzS1.eq.1.and.nclzSn.eq.2) derzzS => derzz_12
      if (nclzS1.eq.2.and.nclzSn.eq.1) derzzS => derzz_21
      if (nclzS1.eq.2.and.nclzSn.eq.2) derzzS => derzz_22
    endif

  end subroutine x3d_derive_init

  !
  ! Associate pointers with subroutines
  !
  subroutine x3d_derive_finalize()

    implicit none

    ! Velocity
    ! First derivative
    nullify(derx)
    nullify(dery)
    nullify(derz)
    ! Second derivative
    nullify(derxx)
    nullify(deryy)
    nullify(derzz)

    ! Scalars
    if (iscalar.ne.0) then
      ! First derivative
      nullify(derxS)
      nullify(deryS)
      nullify(derzS)
      ! Second derivative
      nullify(derxxS)
      nullify(deryyS)
      nullify(derzzS)
    endif

  end subroutine x3d_derive_finalize

subroutine derx_00(tx,ux,x3dop,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  do concurrent (k=1:nz, j=1:ny)
     ! Compute r.h.s.
     tx(1,j,k) = afix*(ux(2,j,k)-ux(nx,j,k)) &
               + bfix*(ux(3,j,k)-ux(nx-1,j,k))
     tx(2,j,k) = afix*(ux(3,j,k)-ux(1,j,k)) &
               + bfix*(ux(4,j,k)-ux(nx,j,k))
     do concurrent (i=3:nx-2)
        tx(i,j,k) = afix*(ux(i+1,j,k)-ux(i-1,j,k)) &
                  + bfix*(ux(i+2,j,k)-ux(i-2,j,k))
     enddo
     tx(nx-1,j,k) = afix*(ux(nx,j,k)-ux(nx-2,j,k)) &
                  + bfix*(ux(1,j,k)-ux(nx-3,j,k))
     tx(nx,j,k) = afix*(ux(1,j,k)-ux(nx-1,j,k)) &
                + bfix*(ux(2,j,k)-ux(nx-2,j,k))

  enddo

  ! Solve tri-diagonal system
  call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derx_00

!********************************************************************
!
subroutine derx_ij(tx,ux,ff,fs,fw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(in), dimension(nx):: ff, fs, fw

  ! Local variables
  integer :: i, j, k

  do concurrent (k=1:nz, j=1:ny)
     ! Compute r.h.s.
     if (ncl1==1) then
        if (npaire==1) then
           tx(1,j,k) = zero
           tx(2,j,k) = afix*(ux(3,j,k)-ux(1,j,k)) &
                     + bfix*(ux(4,j,k)-ux(2,j,k))
        else
           tx(1,j,k) = afix*(ux(2,j,k)+ux(2,j,k)) &
                     + bfix*(ux(3,j,k)+ux(3,j,k))
           tx(2,j,k) = afix*(ux(3,j,k)-ux(1,j,k)) &
                     + bfix*(ux(4,j,k)+ux(2,j,k))
        endif
     else
        tx(1,j,k) = af1x*ux(1,j,k) + bf1x*ux(2,j,k) + cf1x*ux(3,j,k)
        tx(2,j,k) = af2x*(ux(3,j,k)-ux(1,j,k))
     endif
     do concurrent (i=3:nx-2)
        tx(i,j,k) = afix*(ux(i+1,j,k)-ux(i-1,j,k)) &
                  + bfix*(ux(i+2,j,k)-ux(i-2,j,k))
     enddo
     ! nx-1 <= i <= nx
     if (ncln==1) then
        if (npaire==1) then
           tx(nx-1,j,k) = afix*(ux(nx,j,k)-ux(nx-2,j,k)) &
                        + bfix*(ux(nx-1,j,k)-ux(nx-3,j,k))
           tx(nx,j,k) = zero
        else
           tx(nx-1,j,k) = afix*(ux(nx,j,k)-ux(nx-2,j,k)) &
                        + bfix*((-ux(nx-1,j,k))-ux(nx-3,j,k))
           tx(nx,j,k) = afix*((-ux(nx-1,j,k))-ux(nx-1,j,k)) &
                      + bfix*((-ux(nx-2,j,k))-ux(nx-2,j,k))
        endif
     else
        tx(nx-1,j,k) = afmx*(ux(nx,j,k)-ux(nx-2,j,k))
        tx(nx,j,k) = - afnx*ux(nx,j,k) - bfnx*ux(nx-1,j,k) - cfnx*ux(nx-2,j,k)
     endif
  enddo

  ! Solve tri-diagonal system
  call xthomas(tx, ff, fs, fw, nx, ny, nz)

end subroutine derx_ij

!********************************************************************
!
subroutine derx_11(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine derx_11

!********************************************************************
!
subroutine derx_12(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derx_12

!********************************************************************
!
subroutine derx_21(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derx_21

!********************************************************************
!
subroutine derx_22(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derx_22

!********************************************************************
!
subroutine dery_00(ty,uy,x3dop,ppy,nx,ny,nz)
  !
  !********************************************************************

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  ! Compute r.h.s.
  do concurrent (k=1:nz)
     do concurrent (i=1:nx)
        ty(i,1,k) = afjy*(uy(i,2,k)-uy(i,ny,k)) &
                  + bfjy*(uy(i,3,k)-uy(i,ny-1,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,2,k) = afjy*(uy(i,3,k)-uy(i,1,k)) &
                  + bfjy*(uy(i,4,k)-uy(i,ny,k))
     enddo
     do concurrent (j=3:ny-2, i=1:nx)
        ty(i,j,k) = afjy*(uy(i,j+1,k)-uy(i,j-1,k)) &
                  + bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-1,k) = afjy*(uy(i,ny,k)-uy(i,ny-2,k)) &
                     + bfjy*(uy(i,1,k)-uy(i,ny-3,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny,k) = afjy*(uy(i,1,k)-uy(i,ny-1,k)) &
                   + bfjy*(uy(i,2,k)-uy(i,ny-2,k))
     enddo
  enddo

  ! Solve tri-diagonal system
  call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  ! Apply stretching if needed
  if (istret /= 0) then
     do concurrent (k=1:nz, j=1:ny, i=1:nx)
        ty(i,j,k) = ty(i,j,k) * ppy(j)
     enddo
  endif

end subroutine dery_00

!********************************************************************
!
subroutine dery_ij(ty,uy,ff,fs,fw,ppy,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ff,fs,fw,ppy

  ! Local variables
  integer :: i, j, k

  do concurrent (k=1:nz)

     ! Compute r.h.s.
     if (ncl1==1) then
        if (npaire==1) then
           do concurrent (i=1:nx)
              ty(i,1,k) = zero
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = afjy*(uy(i,3,k)-uy(i,1,k)) &
                        + bfjy*(uy(i,4,k)-uy(i,2,k))
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,1,k) = afjy*(uy(i,2,k)+uy(i,2,k)) &
                        + bfjy*(uy(i,3,k)+uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = afjy*(uy(i,3,k)-uy(i,1,k)) &
                        + bfjy*(uy(i,4,k)+uy(i,2,k))
           enddo
        endif
     else
        do concurrent (i=1:nx)
           ty(i,1,k) = af1y*uy(i,1,k)+bf1y*uy(i,2,k)+cf1y*uy(i,3,k)
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = af2y*(uy(i,3,k)-uy(i,1,k))
        enddo
     endif
     do concurrent (j=3:ny-2, i=1:nx)
        ty(i,j,k) = afjy*(uy(i,j+1,k)-uy(i,j-1,k)) &
                  + bfjy*(uy(i,j+2,k)-uy(i,j-2,k))
     enddo
     if (ncln==1) then
        if (npaire==1) then
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = afjy*(uy(i,ny,k)-uy(i,ny-2,k)) &
                           + bfjy*(uy(i,ny-1,k)-uy(i,ny-3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny,k) = zero
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = afjy*(uy(i,ny,k)-uy(i,ny-2,k)) &
                           + bfjy*((-uy(i,ny-1,k))-uy(i,ny-3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny,k) = afjy*((-uy(i,ny-1,k))-uy(i,ny-1,k)) &
                         + bfjy*((-uy(i,ny-2,k))-uy(i,ny-2,k))
           enddo
        endif
     else
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = afmy*(uy(i,ny,k)-uy(i,ny-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny,k) = -afny*uy(i,ny,k)-bfny*uy(i,ny-1,k)-cfny*uy(i,ny-2,k)
        enddo
     endif
  enddo

  ! Solve tri-diagonal system
  call ythomas(ty, ff, fs, fw, nx, ny, nz)

  ! Apply stretching if needed
  if (istret /= 0) then
     do concurrent (k=1:nz, j=1:ny, i=1:nx)
        ty(i,j,k) = ty(i,j,k) * ppy(j)
     enddo
  endif

end subroutine dery_ij

!********************************************************************
!
subroutine dery_11(ty,uy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,1,1)

end subroutine dery_11

!********************************************************************
!
subroutine dery_12(ty,uy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,1,2)

end subroutine dery_12

!********************************************************************
!
subroutine dery_21(ty,uy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,2,1)

end subroutine dery_21

!********************************************************************
!
subroutine dery_22(ty,uy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,2,2)

end subroutine dery_22

!********************************************************************
!
subroutine derz_00(tz,uz,x3dop,nx,ny,nz)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = zero
     enddo
     return
  endif

  ! Compute r.h.s.
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,1) = afkz*(uz(i,j,2)-uz(i,j,nz  )) &
               + bfkz*(uz(i,j,3)-uz(i,j,nz-1))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,2) = afkz*(uz(i,j,3)-uz(i,j,1 )) &
               + bfkz*(uz(i,j,4)-uz(i,j,nz))
  enddo
  do concurrent (k=3:nz-2, j=1:ny, i=1:nx)
     tz(i,j,k) = afkz*(uz(i,j,k+1)-uz(i,j,k-1)) &
               + bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,nz-1) = afkz*(uz(i,j,nz)-uz(i,j,nz-2)) &
                  + bfkz*(uz(i,j,1 )-uz(i,j,nz-3))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,nz) = afkz*(uz(i,j,1)-uz(i,j,nz-1)) &
                + bfkz*(uz(i,j,2)-uz(i,j,nz-2))
  enddo

  ! Solve tri-diagonal system
  call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derz_00

!********************************************************************
!
subroutine derz_ij(tz,uz,ff,fs,fw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(in), dimension(nz) :: ff,fs,fw

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = zero
     enddo
     return
  endif

  ! Compute r.h.s.
  if (ncl1==1) then
     if (npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = zero
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = afkz*(uz(i,j,3)-uz(i,j,1)) &
                     + bfkz*(uz(i,j,4)-uz(i,j,2))
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = afkz*(uz(i,j,2)+uz(i,j,2)) &
                     + bfkz*(uz(i,j,3)+uz(i,j,3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = afkz*(uz(i,j,3)-uz(i,j,1)) &
                     + bfkz*(uz(i,j,4)+uz(i,j,2))
        enddo
     endif
  else
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = af1z*uz(i,j,1) + bf1z*uz(i,j,2) &
                  + cf1z*uz(i,j,3)
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = af2z*(uz(i,j,3)-uz(i,j,1))
     enddo
  endif
  do concurrent (k=3:nz-2, j=1:ny, i=1:nx)
     tz(i,j,k) = afkz*(uz(i,j,k+1)-uz(i,j,k-1)) &
               + bfkz*(uz(i,j,k+2)-uz(i,j,k-2))
  enddo
  if (ncln==1) then
     if (npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = afkz*(uz(i,j,nz  )-uz(i,j,nz-2)) &
                        + bfkz*(uz(i,j,nz-1)-uz(i,j,nz-3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz) = zero
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = afkz*( uz(i,j,nz  )-uz(i,j,nz-2)) &
                        + bfkz*(-uz(i,j,nz-1)-uz(i,j,nz-3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz) = afkz*(-uz(i,j,nz-1)-uz(i,j,nz-1)) &
                      + bfkz*(-uz(i,j,nz-2)-uz(i,j,nz-2))
        enddo
     endif
  else
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = afmz*(uz(i,j,nz)-uz(i,j,nz-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz) = - afnz*uz(i,j,nz) - bfnz*uz(i,j,nz-1) &
                     - cfnz*uz(i,j,nz-2)
     enddo
  endif

  ! Solve tri-diagonal system
  call zthomas(tz, ff, fs, fw, nx, ny, nz)

end subroutine derz_ij

!********************************************************************
!
subroutine derz_11(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine derz_11

!********************************************************************
!
subroutine derz_12(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derz_12

!********************************************************************
!
subroutine derz_21(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derz_21

!********************************************************************
!
subroutine derz_22(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derz_22

!********************************************************************
!
subroutine derxx_00(tx,ux,x3dop,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  ! Compute r.h.s.
  do concurrent (k=1:nz, j=1:ny)
     tx(1,j,k) = asix*(ux(2,j,k)-ux(1   ,j,k) &
                      -ux(1,j,k)+ux(nx  ,j,k)) &
               + bsix*(ux(3,j,k)-ux(1   ,j,k) &
                      -ux(1,j,k)+ux(nx-1,j,k)) &
               + csix*(ux(4,j,k)-ux(1   ,j,k) &
                      -ux(1,j,k)+ux(nx-2,j,k)) &
               + dsix*(ux(5,j,k)-ux(1   ,j,k) &
                      -ux(1,j,k)+ux(nx-3,j,k))
     tx(2,j,k) = asix*(ux(3,j,k)-ux(2   ,j,k) &
                      -ux(2,j,k)+ux(1   ,j,k)) &
               + bsix*(ux(4,j,k)-ux(2   ,j,k) &
                      -ux(2,j,k)+ux(nx  ,j,k)) &
               + csix*(ux(5,j,k)-ux(2   ,j,k) &
                      -ux(2,j,k)+ux(nx-1,j,k)) &
               + dsix*(ux(6,j,k)-ux(2   ,j,k) &
                      -ux(2,j,k)+ux(nx-2,j,k))
     tx(3,j,k) = asix*(ux(4,j,k)-ux(3 ,j,k) &
                      -ux(3,j,k)+ux(2 ,j,k)) &
               + bsix*(ux(5,j,k)-ux(3 ,j,k) &
                      -ux(3,j,k)+ux(1 ,j,k)) &
               + csix*(ux(6,j,k)-ux(3 ,j,k) &
                      -ux(3,j,k)+ux(nx,j,k)) &
               + dsix*(ux(7,j,k)-ux(3 ,j,k) &
                      -ux(3,j,k)+ux(nx-1,j,k))
     tx(4,j,k) = asix*(ux(5,j,k)-ux(4 ,j,k) &
                      -ux(4,j,k)+ux(3 ,j,k)) &
               + bsix*(ux(6,j,k)-ux(4 ,j,k) &
                      -ux(4,j,k)+ux(2,j,k)) &
               + csix*(ux(7,j,k)-ux(4 ,j,k) &
                      -ux(4,j,k)+ux(1,j,k)) &
               + dsix*(ux(8,j,k)-ux(4 ,j,k) &
                      -ux(4,j,k)+ux(nx,j,k))
     do concurrent (i=5:nx-4)
        tx(i,j,k) = asix*(ux(i+1,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-1,j,k)) &
                  + bsix*(ux(i+2,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-2,j,k)) &
                  + csix*(ux(i+3,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-3,j,k)) &
                  + dsix*(ux(i+4,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-4,j,k))
     enddo
     tx(nx-3,j,k) = asix*(ux(nx-2,j,k)-ux(nx-3,j,k) &
                         -ux(nx-3,j,k)+ux(nx-4,j,k)) &
                  + bsix*(ux(nx-1,j,k)-ux(nx-3,j,k) &
                         -ux(nx-3,j,k)+ux(nx-5,j,k)) &
                  + csix*(ux(nx  ,j,k)-ux(nx-3,j,k) &
                         -ux(nx-3,j,k)+ux(nx-6,j,k)) &
                  + dsix*(ux(1   ,j,k)-ux(nx-3,j,k) &
                         -ux(nx-3,j,k)+ux(nx-7,j,k))
     tx(nx-2,j,k) = asix*(ux(nx-1,j,k)-ux(nx-2,j,k) &
                         -ux(nx-2,j,k)+ux(nx-3,j,k)) &
                  + bsix*(ux(nx  ,j,k)-ux(nx-2,j,k) &
                         -ux(nx-2,j,k)+ux(nx-4,j,k)) &
                  + csix*(ux(1   ,j,k)-ux(nx-2,j,k) &
                         -ux(nx-2,j,k)+ux(nx-5,j,k)) &
                  + dsix*(ux(2   ,j,k)-ux(nx-2,j,k) &
                         -ux(nx-2,j,k)+ux(nx-6,j,k))
     tx(nx-1,j,k) = asix*(ux(nx  ,j,k)-ux(nx-1,j,k) &
                         -ux(nx-1,j,k)+ux(nx-2,j,k)) &
                  + bsix*(ux(1   ,j,k)-ux(nx-1,j,k) &
                         -ux(nx-1,j,k)+ux(nx-3,j,k)) &
                  + csix*(ux(2   ,j,k)-ux(nx-1,j,k) &
                         -ux(nx-1,j,k)+ux(nx-4,j,k)) &
                  + dsix*(ux(3   ,j,k)-ux(nx-1,j,k) &
                         -ux(nx-1,j,k)+ux(nx-5,j,k))
     tx(nx  ,j,k) = asix*(ux(1 ,j,k)-ux(nx  ,j,k) &
                         -ux(nx,j,k)+ux(nx-1,j,k)) &
                  + bsix*(ux(2 ,j,k)-ux(nx  ,j,k) &
                         -ux(nx,j,k)+ux(nx-2,j,k)) &
                  + csix*(ux(3 ,j,k)-ux(nx  ,j,k) &
                         -ux(nx,j,k)+ux(nx-3,j,k)) &
                  + dsix*(ux(4 ,j,k)-ux(nx  ,j,k) &
                         -ux(nx,j,k)+ux(nx-4,j,k))
  enddo

  ! Solve tri-diagonal system
  call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derxx_00

!********************************************************************
!
subroutine derxx_ij(tx,ux,sf,ss,sw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(in), dimension(nx):: sf,ss,sw

  ! Local variables
  integer :: i, j, k

  do concurrent (k=1:nz, j=1:ny)

     ! Compute r.h.s.
     if (ncl1==1) then
        if (npaire==1) then
           tx(1,j,k) = asix*(ux(2,j,k)-ux(1,j,k) &
                            -ux(1,j,k)+ux(2,j,k)) &
                     + bsix*(ux(3,j,k)-ux(1,j,k) &
                            -ux(1,j,k)+ux(3,j,k)) &
                     + csix*(ux(4,j,k)-ux(1,j,k) &
                            -ux(1,j,k)+ux(4,j,k)) &
                     + dsix*(ux(5,j,k)-ux(1,j,k) &
                            -ux(1,j,k)+ux(5,j,k))
           tx(2,j,k) = asix*(ux(3,j,k)-ux(2,j,k) &
                            -ux(2,j,k)+ux(1,j,k)) &
                     + bsix*(ux(4,j,k)-ux(2,j,k) &
                            -ux(2,j,k)+ux(2,j,k)) &
                     + csix*(ux(5,j,k)-ux(2,j,k) &
                            -ux(2,j,k)+ux(3,j,k)) &
                     + dsix*(ux(6,j,k)-ux(2,j,k) &
                            -ux(2,j,k)+ux(4,j,k))
           tx(3,j,k) = asix*(ux(4,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(2,j,k)) &
                     + bsix*(ux(5,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(1,j,k)) &
                     + csix*(ux(6,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(2,j,k)) &
                     + dsix*(ux(7,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(3,j,k))
           tx(4,j,k) = asix*(ux(5,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(3,j,k)) &
                     + bsix*(ux(6,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(2,j,k)) &
                     + csix*(ux(7,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(1,j,k)) &
                     + dsix*(ux(8,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(2,j,k))
        else
           tx(1,j,k) = zero
           tx(2,j,k) = asix*(ux(3,j,k)-ux(2,j,k) &
                            -ux(2,j,k)+ux(1,j,k)) &
                     + bsix*(ux(4,j,k)-ux(2,j,k) &
                            -ux(2,j,k)-ux(2,j,k)) &
                     + csix*(ux(5,j,k)-ux(2,j,k) &
                            -ux(2,j,k)-ux(3,j,k)) &
                     + dsix*(ux(6,j,k)-ux(2,j,k) &
                            -ux(2,j,k)-ux(4,j,k))
           tx(3,j,k) = asix*(ux(4,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(2,j,k)) &
                     + bsix*(ux(5,j,k)-ux(3,j,k) &
                            -ux(3,j,k)+ux(1,j,k)) &
                     + csix*(ux(6,j,k)-ux(3,j,k) &
                            -ux(3,j,k)-ux(2,j,k)) &
                     + dsix*(ux(7,j,k)-ux(3,j,k) &
                            -ux(3,j,k)-ux(3,j,k))
           tx(4,j,k) = asix*(ux(5,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(3,j,k)) &
                     + bsix*(ux(6,j,k)-ux(4,j,k) &
                            -ux(4,j,k)+ux(2,j,k)) &
                     + csix*(ux(7,j,k)-ux(4,j,k) &
                            -ux(4,j,k)-ux(1,j,k)) &
                     + dsix*(ux(8,j,k)-ux(4,j,k) &
                            -ux(4,j,k)-ux(2,j,k))
        endif
     else
        tx(1,j,k) = as1x*ux(1,j,k) + bs1x*ux(2,j,k) &
                  + cs1x*ux(3,j,k) + ds1x*ux(4,j,k)
        tx(2,j,k) = as2x*(ux(3,j,k)-ux(2,j,k) &
                         -ux(2,j,k)+ux(1,j,k))
        tx(3,j,k) = as3x*(ux(4,j,k)-ux(3,j,k) &
                        -ux(3,j,k)+ux(2,j,k)) &
                  + bs3x*(ux(5,j,k)-ux(3,j,k) &
                         -ux(3,j,k)+ux(1,j,k))
        tx(4,j,k) = as4x*(ux(5,j,k)-ux(4,j,k) &
                         -ux(4,j,k)+ux(3,j,k)) &
                  + bs4x*(ux(6,j,k)-ux(4,j,k) &
                         -ux(4,j,k)+ux(2,j,k)) &
                  + cs4x*(ux(7,j,k)-ux(4,j,k) &
                         -ux(4,j,k)+ux(1,j,k))
     endif
     do concurrent (i=5:nx-4)
        tx(i,j,k) = asix*(ux(i+1,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-1,j,k)) &
                  + bsix*(ux(i+2,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-2,j,k)) &
                  + csix*(ux(i+3,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-3,j,k)) &
                  + dsix*(ux(i+4,j,k)-ux(i  ,j,k) &
                         -ux(i  ,j,k)+ux(i-4,j,k))
     enddo
     if (ncln == 1) then
        if (npaire==1) then
           tx(nx-3,j,k) = asix*(ux(nx-2,j,k)-ux(nx-3,j,k) &
                               -ux(nx-3,j,k)+ux(nx-4,j,k)) &
                        + bsix*(ux(nx-1,j,k)-ux(nx-3,j,k) &
                               -ux(nx-3,j,k)+ux(nx-5,j,k)) &
                        + csix*(ux(nx  ,j,k)-ux(nx-3,j,k) &
                               -ux(nx-3,j,k)+ux(nx-6,j,k)) &
                        + dsix*(ux(nx-1,j,k)-ux(nx-3,j,k) &
                               -ux(nx-3,j,k)+ux(nx-7,j,k))
           tx(nx-2,j,k) = asix*(ux(nx-1,j,k)-ux(nx-2,j,k) &
                               -ux(nx-2,j,k)+ux(nx-3,j,k)) &
                        + bsix*(ux(nx  ,j,k)-ux(nx-2,j,k) &
                               -ux(nx-2,j,k)+ux(nx-4,j,k)) &
                        + csix*(ux(nx-1,j,k)-ux(nx-2,j,k) &
                               -ux(nx-2,j,k)+ux(nx-5,j,k)) &
                        + dsix*(ux(nx-2,j,k)-ux(nx-2,j,k) &
                               -ux(nx-2,j,k)+ux(nx-6,j,k))
           tx(nx-1,j,k) = asix*(ux(nx  ,j,k)-ux(nx-1,j,k) &
                               -ux(nx-1,j,k)+ux(nx-2,j,k)) &
                        + bsix*(ux(nx-1,j,k)-ux(nx-1,j,k) &
                               -ux(nx-1,j,k)+ux(nx-3,j,k)) &
                        + csix*(ux(nx-2,j,k)-ux(nx-1,j,k) &
                               -ux(nx-1,j,k)+ux(nx-4,j,k)) &
                        + dsix*(ux(nx-3,j,k)-ux(nx-1,j,k) &
                               -ux(nx-1,j,k)+ux(nx-5,j,k))
           tx(nx  ,j,k) = asix*(ux(nx-1,j,k)-ux(nx  ,j,k) &
                               -ux(nx  ,j,k)+ux(nx-1,j,k)) &
                        + bsix*(ux(nx-2,j,k)-ux(nx  ,j,k) &
                               -ux(nx  ,j,k)+ux(nx-2,j,k)) &
                        + csix*(ux(nx-3,j,k)-ux(nx  ,j,k) &
                               -ux(nx  ,j,k)+ux(nx-3,j,k)) &
                        + dsix*(ux(nx-4,j,k)-ux(nx  ,j,k) &
                               -ux(nx  ,j,k)+ux(nx-4,j,k))
        else
           tx(nx-3,j,k) = asix*( ux(nx-2,j,k)-ux(nx-3,j,k) &
                                -ux(nx-3,j,k)+ux(nx-4,j,k)) &
                        + bsix*( ux(nx-1,j,k)-ux(nx-3,j,k) &
                                -ux(nx-3,j,k)+ux(nx-5,j,k)) &
                        + csix*(-ux(nx  ,j,k)-ux(nx-3,j,k) &
                                -ux(nx-3,j,k)+ux(nx-6,j,k)) &
                        + dsix*(-ux(nx-1,j,k)-ux(nx-3,j,k) &
                                -ux(nx-3,j,k)+ux(nx-7,j,k))
           tx(nx-2,j,k) = asix*( ux(nx-1,j,k)-ux(nx-2,j,k) &
                                -ux(nx-2,j,k)+ux(nx-3,j,k)) &
                        + bsix*( ux(nx  ,j,k)-ux(nx-2,j,k) &
                                -ux(nx-2,j,k)+ux(nx-4,j,k)) &
                        + csix*(-ux(nx-1,j,k)-ux(nx-2,j,k) &
                                -ux(nx-2,j,k)+ux(nx-5,j,k)) &
                        + dsix*(-ux(nx-2,j,k)-ux(nx-2,j,k) &
                                -ux(nx-2,j,k)+ux(nx-6,j,k))
           tx(nx-1,j,k) = asix*( ux(nx  ,j,k)-ux(nx-1,j,k) &
                                -ux(nx-1,j,k)+ux(nx-2,j,k)) &
                        + bsix*(-ux(nx-1,j,k)-ux(nx-1,j,k) &
                                -ux(nx-1,j,k)+ux(nx-3,j,k)) &
                        + csix*(-ux(nx-2,j,k)-ux(nx-1,j,k) &
                                -ux(nx-1,j,k)+ux(nx-4,j,k)) &
                        + dsix*(-ux(nx-3,j,k)-ux(nx-1,j,k) &
                                -ux(nx-1,j,k)+ux(nx-5,j,k))
           tx(nx  ,j,k) = zero
        endif
     else
        tx(nx-3,j,k) = asttx*(ux(nx-2,j,k)-ux(nx-3,j,k) &
                             -ux(nx-3,j,k)+ux(nx-4,j,k)) &
                     + bsttx*(ux(nx-1,j,k)-ux(nx-3,j,k) &
                             -ux(nx-3,j,k)+ux(nx-5,j,k)) &
                     + csttx*(ux(nx,j,k)-ux(nx-3,j,k) &
                             -ux(nx-3,j,k)+ux(nx-6,j,k))
        tx(nx-2,j,k) = astx*(ux(nx-1,j,k)-ux(nx-2,j,k) &
                            -ux(nx-2,j,k)+ux(nx-3,j,k)) &
                     + bstx*(ux(nx  ,j,k)-ux(nx-2,j,k) &
                            -ux(nx-2,j,k)+ux(nx-4,j,k))
        tx(nx-1,j,k) = asmx*(ux(nx  ,j,k)-ux(nx-1,j,k) &
                            -ux(nx-1,j,k)+ux(nx-2,j,k))
        tx(nx  ,j,k) = asnx*ux(nx  ,j,k) + bsnx*ux(nx-1,j,k) &
                     + csnx*ux(nx-2,j,k) + dsnx*ux(nx-3,j,k)
     endif
  enddo

  ! Solve tri-diagonal system
  call xthomas(tx, sf, ss, sw, nx, ny, nz)

end subroutine derxx_ij

!********************************************************************
!
subroutine derxx_11(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derxx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine derxx_11

!********************************************************************
!
subroutine derxx_12(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derxx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derxx_12

!********************************************************************
!
subroutine derxx_21(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derxx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derxx_21

!********************************************************************
!
subroutine derxx_22(tx,ux,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  call derxx_ij(tx,ux,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derxx_22

!********************************************************************
!
subroutine deryy_00(ty,uy,x3dop,nx,ny,nz)

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  ! Compute r.h.s.
  do concurrent (k=1:nz)
     do concurrent (i=1:nx)
        ty(i,1,k) = asjy*(uy(i,2,k)-uy(i,1,k) &
                         -uy(i,1,k)+uy(i,ny,k)) &
                  + bsjy*(uy(i,3,k)-uy(i,1,k) &
                         -uy(i,1,k)+uy(i,ny-1,k)) &
                  + csjy*(uy(i,4,k)-uy(i,1,k) &
                         -uy(i,1,k)+uy(i,ny-2,k)) &
                  + dsjy*(uy(i,5,k)-uy(i,1,k) &
                         -uy(i,1,k)+uy(i,ny-3,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,2,k) = asjy*(uy(i,3,k)-uy(i,2,k) &
                         -uy(i,2,k)+uy(i,1,k)) &
                  + bsjy*(uy(i,4,k)-uy(i,2,k) &
                         -uy(i,2,k)+uy(i,ny,k)) &
                  + csjy*(uy(i,5,k)-uy(i,2,k) &
                         -uy(i,2,k)+uy(i,ny-1,k)) &
                  + dsjy*(uy(i,6,k)-uy(i,2,k) &
                         -uy(i,2,k)+uy(i,ny-2,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,3,k) = asjy*(uy(i,4,k)-uy(i,3,k) &
                         -uy(i,3,k)+uy(i,2,k)) &
                  + bsjy*(uy(i,5,k)-uy(i,3,k) &
                         -uy(i,3,k)+uy(i,1,k)) &
                  + csjy*(uy(i,6,k)-uy(i,3,k) &
                         -uy(i,3,k)+uy(i,ny,k)) &
                  + dsjy*(uy(i,7,k)-uy(i,3,k) &
                         -uy(i,3,k)+uy(i,ny-1,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,4,k) = asjy*(uy(i,5,k)-uy(i,4,k) &
                         -uy(i,4,k)+uy(i,3,k)) &
                  + bsjy*(uy(i,6,k)-uy(i,4,k) &
                         -uy(i,4,k)+uy(i,2,k)) &
                  + csjy*(uy(i,7,k)-uy(i,4,k) &
                         -uy(i,4,k)+uy(i,1,k)) &
                  + dsjy*(uy(i,8,k)-uy(i,4,k) &
                         -uy(i,4,k)+uy(i,ny,k))
     enddo
     do concurrent (j=5:ny-4, i=1:nx)
        ty(i,j,k) = asjy*(uy(i,j+1,k)-uy(i,j,k) &
                         -uy(i,j,k)+uy(i,j-1,k)) &
                  + bsjy*(uy(i,j+2,k)-uy(i,j,k) &
                         -uy(i,j,k)+uy(i,j-2,k)) &
                  + csjy*(uy(i,j+3,k)-uy(i,j,k) &
                         -uy(i,j,k)+uy(i,j-3,k)) &
                  + dsjy*(uy(i,j+4,k)-uy(i,j,k) &
                         -uy(i,j,k)+uy(i,j-4,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-3,k) = asjy*(uy(i,ny-2,k)-uy(i,ny-3,k) &
                            -uy(i,ny-3,k)+uy(i,ny-4,k)) &
                     + bsjy*(uy(i,ny-1,k)-uy(i,ny-3,k) &
                            -uy(i,ny-3,k)+uy(i,ny-5,k)) &
                     + csjy*(uy(i,ny  ,k)-uy(i,ny-3,k) &
                            -uy(i,ny-3,k)+uy(i,ny-6,k)) &
                     + dsjy*(uy(i,1   ,k)-uy(i,ny-3,k) &
                            -uy(i,ny-3,k)+uy(i,ny-7,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-2,k) = asjy*(uy(i,ny-1,k)-uy(i,ny-2,k) &
                            -uy(i,ny-2,k)+uy(i,ny-3,k)) &
                     + bsjy*(uy(i,ny  ,k)-uy(i,ny-2,k) &
                            -uy(i,ny-2,k)+uy(i,ny-4,k)) &
                     + csjy*(uy(i,1   ,k)-uy(i,ny-2,k) &
                            -uy(i,ny-2,k)+uy(i,ny-5,k)) &
                     + dsjy*(uy(i,2   ,k)-uy(i,ny-2,k) &
                            -uy(i,ny-2,k)+uy(i,ny-6,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny-1,k) = asjy*(uy(i,ny  ,k)-uy(i,ny-1,k) &
                            -uy(i,ny-1,k)+uy(i,ny-2,k)) &
                     + bsjy*(uy(i,1   ,k)-uy(i,ny-1,k) &
                            -uy(i,ny-1,k)+uy(i,ny-3,k)) &
                     + csjy*(uy(i,2   ,k)-uy(i,ny-1,k) &
                            -uy(i,ny-1,k)+uy(i,ny-4,k)) &
                     + dsjy*(uy(i,3   ,k)-uy(i,ny-1,k) &
                            -uy(i,ny-1,k)+uy(i,ny-5,k))
     enddo
     do concurrent (i=1:nx)
        ty(i,ny  ,k) = asjy*(uy(i,1 ,k)-uy(i,ny  ,k) &
                            -uy(i,ny,k)+uy(i,ny-1,k)) &
                     + bsjy*(uy(i,2 ,k)-uy(i,ny  ,k) &
                            -uy(i,ny,k)+uy(i,ny-2,k)) &
                     + csjy*(uy(i,3 ,k)-uy(i,ny  ,k) &
                            -uy(i,ny,k)+uy(i,ny-3,k)) &
                     + dsjy*(uy(i,4 ,k)-uy(i,ny  ,k) &
                            -uy(i,ny,k)+uy(i,ny-4,k))
     enddo
  enddo

  ! Solve tri-diagonal system
  call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine deryy_00

!********************************************************************
!
subroutine deryy_ij(ty,uy,sf,ss,sw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: sf,ss,sw

  ! Local variables
  integer :: i, j, k

  ! Compute r.h.s.
  do concurrent (k=1:nz)
     if (ncl1==1) then
        if (npaire==1) then
           do concurrent (i=1:nx)
              ty(i,1,k) = asjy*(uy(i,2,k)-uy(i,1,k) &
                               -uy(i,1,k)+uy(i,2,k)) &
                        + bsjy*(uy(i,3,k)-uy(i,1,k) &
                               -uy(i,1,k)+uy(i,3,k)) &
                        + csjy*(uy(i,4,k)-uy(i,1,k) &
                               -uy(i,1,k)+uy(i,4,k)) &
                        + dsjy*(uy(i,5,k)-uy(i,1,k) &
                               -uy(i,1,k)+uy(i,5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = asjy*(uy(i,3,k)-uy(i,2,k) &
                               -uy(i,2,k)+uy(i,1,k)) &
                        + bsjy*(uy(i,4,k)-uy(i,2,k) &
                               -uy(i,2,k)+uy(i,2,k)) &
                        + csjy*(uy(i,5,k)-uy(i,2,k) &
                               -uy(i,2,k)+uy(i,3,k)) &
                        + dsjy*(uy(i,6,k)-uy(i,2,k) &
                               -uy(i,2,k)+uy(i,4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = asjy*(uy(i,4,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,2,k)) &
                        + bsjy*(uy(i,5,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,1,k)) &
                        + csjy*(uy(i,6,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,2,k)) &
                        + dsjy*(uy(i,7,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,4,k) = asjy*(uy(i,5,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,3,k)) &
                        + bsjy*(uy(i,6,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,2,k)) &
                        + csjy*(uy(i,7,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,1,k)) &
                        + dsjy*(uy(i,8,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,2,k))
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,1,k) = zero
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = asjy*(uy(i,3,k)-uy(i,2,k) &
                               -uy(i,2,k)+uy(i,1,k)) &
                        + bsjy*(uy(i,4,k)-uy(i,2,k) &
                               -uy(i,2,k)-uy(i,2,k)) &
                        + csjy*(uy(i,5,k)-uy(i,2,k) &
                               -uy(i,2,k)-uy(i,3,k)) &
                        + dsjy*(uy(i,6,k)-uy(i,2,k) &
                               -uy(i,2,k)-uy(i,4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = asjy*(uy(i,4,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,2,k)) &
                        + bsjy*(uy(i,5,k)-uy(i,3,k) &
                               -uy(i,3,k)+uy(i,1,k)) &
                        + csjy*(uy(i,6,k)-uy(i,3,k) &
                               -uy(i,3,k)-uy(i,2,k)) &
                        + dsjy*(uy(i,7,k)-uy(i,3,k) &
                               -uy(i,3,k)-uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,4,k) = asjy*(uy(i,5,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,3,k)) &
                        + bsjy*(uy(i,6,k)-uy(i,4,k) &
                               -uy(i,4,k)+uy(i,2,k)) &
                        + csjy*(uy(i,7,k)-uy(i,4,k) &
                               -uy(i,4,k)-uy(i,1,k)) &
                        + dsjy*(uy(i,8,k)-uy(i,4,k) &
                               -uy(i,4,k)-uy(i,2,k))
           enddo
        endif
     else
        do concurrent (i=1:nx)
           ty(i,1,k) = as1y*uy(i,1,k) + bs1y*uy(i,2,k) &
                     + cs1y*uy(i,3,k) + ds1y*uy(i,4,k)
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = as2y*(uy(i,3,k)-uy(i,2,k) &
                            -uy(i,2,k)+uy(i,1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,3,k) = as3y*(uy(i,4,k)-uy(i,3,k) &
                            -uy(i,3,k)+uy(i,2,k)) &
                     + bs3y*(uy(i,5,k)-uy(i,3,k) &
                            -uy(i,3,k)+uy(i,1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,4,k) = as4y*(uy(i,5,k)-uy(i,4,k) &
                            -uy(i,4  ,k)+uy(i,3,k)) &
                     + bs4y*(uy(i,6,k)-uy(i,4  ,k) &
                            -uy(i,4  ,k)+uy(i,2,k)) &
                     + cs4y*(uy(i,7,k)-uy(i,4  ,k) &
                            -uy(i,4  ,k)+uy(i,1,k))
        enddo
     endif
     do concurrent (j=5:ny-4, i=1:nx)
        ty(i,j,k) = asjy*(uy(i,j+1,k)-uy(i,j  ,k) &
                         -uy(i,j  ,k)+uy(i,j-1,k)) &
                  + bsjy*(uy(i,j+2,k)-uy(i,j  ,k) &
                         -uy(i,j  ,k)+uy(i,j-2,k)) &
                  + csjy*(uy(i,j+3,k)-uy(i,j  ,k) &
                         -uy(i,j  ,k)+uy(i,j-3,k)) &
                  + dsjy*(uy(i,j+4,k)-uy(i,j  ,k) &
                         -uy(i,j  ,k)+uy(i,j-4,k))
     enddo
     if (ncln==1) then
        if (npaire==1) then
           do concurrent (i=1:nx)
              ty(i,ny-3,k) = asjy*(uy(i,ny-2,k)-uy(i,ny-3,k) &
                                  -uy(i,ny-3,k)+uy(i,ny-4,k)) &
                           + bsjy*(uy(i,ny-1,k)-uy(i,ny-3,k) &
                                  -uy(i,ny-3,k)+uy(i,ny-5,k)) &
                           + csjy*(uy(i,ny  ,k)-uy(i,ny-3,k) &
                                  -uy(i,ny-3,k)+uy(i,ny-6,k)) &
                           + dsjy*(uy(i,ny-1,k)-uy(i,ny-3,k) &
                                  -uy(i,ny-3,k)+uy(i,ny-7,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-2,k) = asjy*(uy(i,ny-1,k)-uy(i,ny-2,k) &
                                  -uy(i,ny-2,k)+uy(i,ny-3,k)) &
                           + bsjy*(uy(i,ny  ,k)-uy(i,ny-2,k) &
                                  -uy(i,ny-2,k)+uy(i,ny-4,k)) &
                           + csjy*(uy(i,ny-1,k)-uy(i,ny-2,k) &
                                  -uy(i,ny-2,k)+uy(i,ny-5,k)) &
                           + dsjy*(uy(i,ny-2,k)-uy(i,ny-2,k) &
                                  -uy(i,ny-2,k)+uy(i,ny-6,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = asjy*(uy(i,ny  ,k)-uy(i,ny-1,k) &
                                  -uy(i,ny-1,k)+uy(i,ny-2,k)) &
                           + bsjy*(uy(i,ny-1,k)-uy(i,ny-1,k) &
                                  -uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + csjy*(uy(i,ny-2,k)-uy(i,ny-1,k) &
                                  -uy(i,ny-1,k)+uy(i,ny-4,k)) &
                           + dsjy*(uy(i,ny-3,k)-uy(i,ny-1,k) &
                                  -uy(i,ny-1,k)+uy(i,ny-5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny  ,k) = asjy*(uy(i,ny-1,k)-uy(i,ny  ,k) &
                                  -uy(i,ny  ,k)+uy(i,ny-1,k)) &
                           + bsjy*(uy(i,ny-2,k)-uy(i,ny  ,k) &
                                  -uy(i,ny  ,k)+uy(i,ny-2,k)) &
                           + csjy*(uy(i,ny-3,k)-uy(i,ny  ,k) &
                                  -uy(i,ny  ,k)+uy(i,ny-3,k)) &
                           + dsjy*(uy(i,ny-4,k)-uy(i,ny  ,k) &
                                  -uy(i,ny  ,k)+uy(i,ny-4,k))
           enddo
        else
           do concurrent (i=1:nx)
              ty(i,ny-3,k) = asjy*( uy(i,ny-2,k)-uy(i,ny-3,k) &
                                   -uy(i,ny-3,k)+uy(i,ny-4,k)) &
                           + bsjy*( uy(i,ny-1,k)-uy(i,ny-3,k) &
                                   -uy(i,ny-3,k)+uy(i,ny-5,k)) &
                           + csjy*(-uy(i,ny ,k)-uy(i,ny-3,k) &
                                   -uy(i,ny-3,k)+uy(i,ny-6,k)) &
                           + dsjy*(-uy(i,ny-1,k)-uy(i,ny-3,k) &
                                   -uy(i,ny-3,k)+uy(i,ny-7,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-2,k) = asjy*( uy(i,ny-1,k)-uy(i,ny-2,k) &
                                   -uy(i,ny-2,k)+uy(i,ny-3,k)) &
                           + bsjy*( uy(i,ny  ,k)-uy(i,ny-2,k) &
                                   -uy(i,ny-2,k)+uy(i,ny-4,k)) &
                           + csjy*(-uy(i,ny-1,k)-uy(i,ny-2,k) &
                                   -uy(i,ny-2,k)+uy(i,ny-5,k)) &
                           + dsjy*(-uy(i,ny-2,k)-uy(i,ny-2,k) &
                                   -uy(i,ny-2,k)+uy(i,ny-6,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = asjy*( uy(i,ny  ,k)-uy(i,ny-1,k) &
                                   -uy(i,ny-1,k)+uy(i,ny-2,k)) &
                           + bsjy*(-uy(i,ny-1,k)-uy(i,ny-1,k) &
                                   -uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + csjy*(-uy(i,ny-2,k)-uy(i,ny-1,k) &
                                   -uy(i,ny-1,k)+uy(i,ny-4,k)) &
                           + dsjy*(-uy(i,ny-3,k)-uy(i,ny-1,k) &
                                   -uy(i,ny-1,k)+uy(i,ny-5,k))
              ty(i,ny  ,k) = zero
           enddo
        endif
     else
        do concurrent (i=1:nx)
           ty(i,ny-3,k) = astty*(uy(i,ny-2,k)-uy(i,ny-3  ,k) &
                                -uy(i,ny-3  ,k)+uy(i,ny-4,k)) &
                        + bstty*(uy(i,ny-1,k)-uy(i,ny-3  ,k) &
                                -uy(i,ny-3  ,k)+uy(i,ny-5,k)) &
                        + cstty*(uy(i,ny,k)-uy(i,ny-3  ,k) &
                                -uy(i,ny-3  ,k)+uy(i,ny-6,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-2,k) = asty*(uy(i,ny-1,k)-uy(i,ny-2,k) &
                               -uy(i,ny-2,k)+uy(i,ny-3,k)) &
                        + bsty*(uy(i,ny  ,k)-uy(i,ny-2,k) &
                               -uy(i,ny-2,k)+uy(i,ny-4,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = asmy*(uy(i,ny  ,k)-uy(i,ny-1,k) &
                               -uy(i,ny-1,k)+uy(i,ny-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny  ,k) = asny*uy(i,ny  ,k) + bsny*uy(i,ny-1,k) &
                        + csny*uy(i,ny-2,k) + dsny*uy(i,ny-3,k)
        enddo
     endif
  enddo

  ! Solve tri-diagonal system
  call ythomas(ty, sf, ss, sw, nx, ny, nz)

end subroutine deryy_ij

!********************************************************************
!
subroutine deryy_11(ty,uy,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  call deryy_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine deryy_11

!********************************************************************
!
subroutine deryy_12(ty,uy,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  call deryy_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine deryy_12

!********************************************************************
!
subroutine deryy_21(ty,uy,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  call deryy_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine deryy_21

!********************************************************************
!
subroutine deryy_22(ty,uy,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  call deryy_ij(ty,uy,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine deryy_22

!********************************************************************
!
subroutine derzz_00(tz,uz,x3dop,nx,ny,nz)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
    do concurrent(k=1:nz, j=1:ny, i=1:nx)
      tz(i,j,k) = zero
    enddo
    return
  endif

  ! Compute r.h.s.
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,1) = askz*(uz(i,j,2)-uz(i,j,1   ) &
                      -uz(i,j,1)+uz(i,j,nz  )) &
               + bskz*(uz(i,j,3)-uz(i,j,1   ) &
                      -uz(i,j,1)+uz(i,j,nz-1)) &
               + cskz*(uz(i,j,4)-uz(i,j,1   ) &
                      -uz(i,j,1)+uz(i,j,nz-2)) &
               + dskz*(uz(i,j,5)-uz(i,j,1   ) &
                      -uz(i,j,1)+uz(i,j,nz-3))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,2) = askz*(uz(i,j,3)-uz(i,j,2 ) &
                      -uz(i,j,2)+uz(i,j,1 )) &
               + bskz*(uz(i,j,4)-uz(i,j,2 ) &
                      -uz(i,j,2)+uz(i,j,nz)) &
               + cskz*(uz(i,j,5)-uz(i,j,2 ) &
                      -uz(i,j,2)+uz(i,j,nz-1)) &
               + dskz*(uz(i,j,6)-uz(i,j,2 ) &
                      -uz(i,j,2)+uz(i,j,nz-2))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,3) = askz*(uz(i,j,4)-uz(i,j,3 ) &
                      -uz(i,j,3)+uz(i,j,2 )) &
               + bskz*(uz(i,j,5)-uz(i,j,3 ) &
                      -uz(i,j,3)+uz(i,j,1 )) &
               + cskz*(uz(i,j,6)-uz(i,j,3 ) &
                      -uz(i,j,3)+uz(i,j,nz)) &
               + dskz*(uz(i,j,7)-uz(i,j,3 ) &
                      -uz(i,j,3)+uz(i,j,nz-1))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,4) = askz*(uz(i,j,5)-uz(i,j,4 ) &
                      -uz(i,j,4)+uz(i,j,3 )) &
               + bskz*(uz(i,j,6)-uz(i,j,4 ) &
                      -uz(i,j,4)+uz(i,j,2 )) &
               + cskz*(uz(i,j,7)-uz(i,j,4 ) &
                      -uz(i,j,4)+uz(i,j,1)) &
               + dskz*(uz(i,j,8)-uz(i,j,4 ) &
                      -uz(i,j,4)+uz(i,j,nz))
  enddo
  do concurrent (k=5:nz-4, j=1:ny, i=1:nx)
     tz(i,j,k) = askz*(uz(i,j,k+1)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-1)) &
               + bskz*(uz(i,j,k+2)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-2)) &
               + cskz*(uz(i,j,k+3)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-3)) &
               + dskz*(uz(i,j,k+4)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-4))
  enddo
  do concurrent (j=1:ny, i=1:nx)
     tz(i,j,nz-3) = askz*(uz(i,j,nz-2)-uz(i,j,nz-3) &
                         -uz(i,j,nz-3)+uz(i,j,nz-4)) &
                  + bskz*(uz(i,j,nz-1 )-uz(i,j,nz-3) &
                         -uz(i,j,nz-3)+uz(i,j,nz-5)) &
                  + cskz*(uz(i,j,nz  )-uz(i,j,nz-3) &
                         -uz(i,j,nz-3)+uz(i,j,nz-6)) &
                  + dskz*(uz(i,j,1   )-uz(i,j,nz-3) &
                         -uz(i,j,nz-3)+uz(i,j,nz-7))
     tz(i,j,nz-2) = askz*(uz(i,j,nz-1)-uz(i,j,nz-2) &
                         -uz(i,j,nz-2)+uz(i,j,nz-3)) &
                  + bskz*(uz(i,j,nz  )-uz(i,j,nz-2) &
                         -uz(i,j,nz-2)+uz(i,j,nz-4)) &
                  + cskz*(uz(i,j,1   )-uz(i,j,nz-2) &
                         -uz(i,j,nz-2)+uz(i,j,nz-5)) &
                  + dskz*(uz(i,j,2   )-uz(i,j,nz-2) &
                         -uz(i,j,nz-2)+uz(i,j,nz-6))
     tz(i,j,nz-1) = askz*(uz(i,j,nz  )-uz(i,j,nz-1) &
                         -uz(i,j,nz-1)+uz(i,j,nz-2)) &
                  + bskz*(uz(i,j,1   )-uz(i,j,nz-1) &
                         -uz(i,j,nz-1)+uz(i,j,nz-3)) &
                  + cskz*(uz(i,j,2   )-uz(i,j,nz-1) &
                         -uz(i,j,nz-1)+uz(i,j,nz-4)) &
                  + dskz*(uz(i,j,3   )-uz(i,j,nz-1) &
                         -uz(i,j,nz-1)+uz(i,j,nz-5))
     tz(i,j,nz  ) = askz*(uz(i,j,1 )-uz(i,j,nz  ) &
                         -uz(i,j,nz)+uz(i,j,nz-1)) &
                  + bskz*(uz(i,j,2 )-uz(i,j,nz  ) &
                         -uz(i,j,nz)+uz(i,j,nz-2)) &
                  + cskz*(uz(i,j,3 )-uz(i,j,nz  ) &
                         -uz(i,j,nz)+uz(i,j,nz-3)) &
                  + dskz*(uz(i,j,4 )-uz(i,j,nz  ) &
                         -uz(i,j,nz)+uz(i,j,nz-4))
  enddo

  ! Solve tri-diagonal system
  call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derzz_00

!********************************************************************
!
subroutine derzz_ij(tz,uz,sf,ss,sw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(in), dimension(nz) :: sf,ss,sw

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
    do concurrent(k=1:nz, j=1:ny, i=1:nx)
      tz(i,j,k) = zero
    enddo
  endif

  ! Compute r.h.s.
  if (ncl1==1) then
     if (npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = askz*(uz(i,j,2)-uz(i,j,1) &
                            -uz(i,j,1)+uz(i,j,2)) &
                     + bskz*(uz(i,j,3)-uz(i,j,1) &
                            -uz(i,j,1)+uz(i,j,3)) &
                     + cskz*(uz(i,j,4)-uz(i,j,1) &
                            -uz(i,j,1)+uz(i,j,4)) &
                     + dskz*(uz(i,j,5)-uz(i,j,1) &
                            -uz(i,j,1)+uz(i,j,5))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = askz*(uz(i,j,3)-uz(i,j,2) &
                            -uz(i,j,2)+uz(i,j,1)) &
                     + bskz*(uz(i,j,4)-uz(i,j,2) &
                            -uz(i,j,2)+uz(i,j,2)) &
                     + cskz*(uz(i,j,5)-uz(i,j,2) &
                            -uz(i,j,2)+uz(i,j,3)) &
                     + dskz*(uz(i,j,6)-uz(i,j,2) &
                            -uz(i,j,2)+uz(i,j,4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,3) = askz*(uz(i,j,4)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,2)) &
                     + bskz*(uz(i,j,5)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,1)) &
                     + cskz*(uz(i,j,6)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,2)) &
                     + dskz*(uz(i,j,7)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,4) = askz*(uz(i,j,5)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,3)) &
                     + bskz*(uz(i,j,6)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,2)) &
                     + cskz*(uz(i,j,7)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,1)) &
                     + dskz*(uz(i,j,8)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,2))
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = zero
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = askz*(uz(i,j,3)-uz(i,j,2) &
                            -uz(i,j,2)+uz(i,j,1)) &
                     + bskz*(uz(i,j,4)-uz(i,j,2) &
                            -uz(i,j,2)-uz(i,j,2)) &
                     + cskz*(uz(i,j,5)-uz(i,j,2) &
                            -uz(i,j,2)-uz(i,j,3)) &
                     + dskz*(uz(i,j,6)-uz(i,j,2) &
                            -uz(i,j,2)-uz(i,j,4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,3) = askz*(uz(i,j,4)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,2)) &
                     + bskz*(uz(i,j,5)-uz(i,j,3) &
                            -uz(i,j,3)+uz(i,j,1)) &
                     + cskz*(uz(i,j,6)-uz(i,j,3) &
                            -uz(i,j,3)-uz(i,j,2)) &
                     + dskz*(uz(i,j,7)-uz(i,j,3) &
                            -uz(i,j,3)-uz(i,j,3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,4) = askz*(uz(i,j,5)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,3)) &
                     + bskz*(uz(i,j,6)-uz(i,j,4) &
                            -uz(i,j,4)+uz(i,j,2)) &
                     + cskz*(uz(i,j,7)-uz(i,j,4) &
                            -uz(i,j,4)-uz(i,j,1)) &
                     + dskz*(uz(i,j,8)-uz(i,j,4) &
                            -uz(i,j,4)-uz(i,j,2))
        enddo
     endif
  else
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = as1z*uz(i,j,1) + bs1z*uz(i,j,2) &
                  + cs1z*uz(i,j,3) + ds1z*uz(i,j,4)
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = as2z*(uz(i,j,3)-uz(i,j,2) &
                         -uz(i,j,2)+uz(i,j,1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,3) = as3z*(uz(i,j,4)-uz(i,j,3) &
                         -uz(i,j,3)+uz(i,j,2)) &
                  + bs3z*(uz(i,j,5)-uz(i,j,3) &
                         -uz(i,j,3)+uz(i,j,1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,4) = as4z*(uz(i,j,5)-uz(i,j,4  ) &
                         -uz(i,j,4  )+uz(i,j,3)) &
                  + bs4z*(uz(i,j,6)-uz(i,j,4 ) &
                         -uz(i,j,4 )+uz(i,j,2)) &
                  + cs4z*(uz(i,j,7)-uz(i,j,4  ) &
                         -uz(i,j,4  )+uz(i,j,1))
     enddo
  endif
  do concurrent (k=5:nz-4, j=1:ny, i=1:nx)
     tz(i,j,k) = askz*(uz(i,j,k+1)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-1)) &
               + bskz*(uz(i,j,k+2)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-2)) &
               + cskz*(uz(i,j,k+3)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-3)) &
               + dskz*(uz(i,j,k+4)-uz(i,j,k  ) &
                      -uz(i,j,k  )+uz(i,j,k-4))
  enddo
  if (ncln==1) then
     if (npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-3) = askz*(uz(i,j,nz-2)-uz(i,j,nz-3) &
                               -uz(i,j,nz-3)+uz(i,j,nz-4)) &
                        + bskz*(uz(i,j,nz-1)-uz(i,j,nz-3) &
                               -uz(i,j,nz-3)+uz(i,j,nz-5)) &
                        + cskz*(uz(i,j,nz  )-uz(i,j,nz-3) &
                               -uz(i,j,nz-3)+uz(i,j,nz-6)) &
                        + dskz*(uz(i,j,nz-1)-uz(i,j,nz-3) &
                               -uz(i,j,nz-3)+uz(i,j,nz-7))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-2) = askz*(uz(i,j,nz-1)-uz(i,j,nz-2) &
                               -uz(i,j,nz-2)+uz(i,j,nz-3)) &
                        + bskz*(uz(i,j,nz  )-uz(i,j,nz-2) &
                               -uz(i,j,nz-2)+uz(i,j,nz-4)) &
                        + cskz*(uz(i,j,nz-1)-uz(i,j,nz-2) &
                               -uz(i,j,nz-2)+uz(i,j,nz-5)) &
                        + dskz*(uz(i,j,nz-2)-uz(i,j,nz-2) &
                               -uz(i,j,nz-2)+uz(i,j,nz-6))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = askz*(uz(i,j,nz  )-uz(i,j,nz-1) &
                               -uz(i,j,nz-1)+uz(i,j,nz-2)) &
                        + bskz*(uz(i,j,nz-1)-uz(i,j,nz-1) &
                               -uz(i,j,nz-1)+uz(i,j,nz-3)) &
                        + cskz*(uz(i,j,nz-2)-uz(i,j,nz-1) &
                               -uz(i,j,nz-1)+uz(i,j,nz-4)) &
                        + dskz*(uz(i,j,nz-3)-uz(i,j,nz-1) &
                               -uz(i,j,nz-1)+uz(i,j,nz-5))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz  ) = askz*(uz(i,j,nz-1)-uz(i,j,nz  ) &
                               -uz(i,j,nz  )+uz(i,j,nz-1)) &
                        + bskz*(uz(i,j,nz-2)-uz(i,j,nz  ) &
                               -uz(i,j,nz  )+uz(i,j,nz-2)) &
                        + cskz*(uz(i,j,nz-3)-uz(i,j,nz  ) &
                               -uz(i,j,nz  )+uz(i,j,nz-3)) &
                        + dskz*(uz(i,j,nz-4)-uz(i,j,nz  ) &
                               -uz(i,j,nz  )+uz(i,j,nz-4))
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-3) = askz*( uz(i,j,nz-2)-uz(i,j,nz-3) &
                                -uz(i,j,nz-3)+uz(i,j,nz-4)) &
                        + bskz*( uz(i,j,nz-1)-uz(i,j,nz-3) &
                                -uz(i,j,nz-3)+uz(i,j,nz-5)) &
                        + cskz*(-uz(i,j,nz  )-uz(i,j,nz-3) &
                                -uz(i,j,nz-3)+uz(i,j,nz-6)) &
                        + dskz*(-uz(i,j,nz-1)-uz(i,j,nz-3) &
                                -uz(i,j,nz-3)+uz(i,j,nz-7))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-2) = askz*( uz(i,j,nz-1)-uz(i,j,nz-2) &
                                -uz(i,j,nz-2)+uz(i,j,nz-3)) &
                        + bskz*( uz(i,j,nz  )-uz(i,j,nz-2) &
                                -uz(i,j,nz-2)+uz(i,j,nz-4)) &
                        + cskz*(-uz(i,j,nz-1)-uz(i,j,nz-2) &
                                -uz(i,j,nz-2)+uz(i,j,nz-5)) &
                        + dskz*(-uz(i,j,nz-2)-uz(i,j,nz-2) &
                                -uz(i,j,nz-2)+uz(i,j,nz-6))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = askz*( uz(i,j,nz  )-uz(i,j,nz-1) &
                                -uz(i,j,nz-1)+uz(i,j,nz-2)) &
                        + bskz*(-uz(i,j,nz-1)-uz(i,j,nz-1) &
                                -uz(i,j,nz-1)+uz(i,j,nz-3)) &
                        + cskz*(-uz(i,j,nz-2)-uz(i,j,nz-1) &
                                -uz(i,j,nz-1)+uz(i,j,nz-4)) &
                        + dskz*(-uz(i,j,nz-3)-uz(i,j,nz-1) &
                                -uz(i,j,nz-1)+uz(i,j,nz-5))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz  ) = zero
        enddo
     endif
  else
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-3) = asttz*(uz(i,j,nz-2)-uz(i,j,nz-3  ) &
                             -uz(i,j,nz-3  )+uz(i,j,nz-4)) &
                     + bsttz*(uz(i,j,nz-1)-uz(i,j,nz-3  ) &
                             -uz(i,j,nz-3  )+uz(i,j,nz-5)) &
                     + csttz*(uz(i,j,nz)-uz(i,j,nz-3  ) &
                             -uz(i,j,nz-3  )+uz(i,j,nz-6))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-2) = astz*(uz(i,j,nz-1)-uz(i,j,nz-2) &
                            -uz(i,j,nz-2)+uz(i,j,nz-3)) &
                     + bstz*(uz(i,j,nz  )-uz(i,j,nz-2) &
                            -uz(i,j,nz-2)+uz(i,j,nz-4))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = asmz*(uz(i,j,nz  )-uz(i,j,nz-1) &
                            -uz(i,j,nz-1)+uz(i,j,nz-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz  ) = asnz*uz(i,j,nz  ) + bsnz*uz(i,j,nz-1) &
                     + csnz*uz(i,j,nz-2) + dsnz*uz(i,j,nz-3)
     enddo
  endif

  ! Solve tri-diagonal system
  call zthomas(tz, sf, ss, sw, nx, ny, nz)

end subroutine derzz_ij

!********************************************************************
!
subroutine derzz_11(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derzz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine derzz_11

!********************************************************************
!
subroutine derzz_12(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derzz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derzz_12

!********************************************************************
!
subroutine derzz_21(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derzz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derzz_21

!********************************************************************
!
subroutine derzz_22(tz,uz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  call derzz_ij(tz,uz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derzz_22

end module x3d_derive
