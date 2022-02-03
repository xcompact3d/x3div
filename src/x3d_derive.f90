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
!##################################################################

module x3d_derive

  use decomp_2d, only : mytype
  use x3d_operator_1d, only : x3doperator1d
  use param
  use thomas
  
  implicit none

  ! Make everything public unless declared private
  public

  ABSTRACT INTERFACE
     SUBROUTINE DERIVATIVE_X(t,u,s,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(ny,nz):: s
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_X
     SUBROUTINE DERIVATIVE_Y(t,u,s,x3dop,pp,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(nx,nz):: s
       real(mytype), intent(in), dimension(ny):: pp
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_Y
     SUBROUTINE DERIVATIVE_YY(t,u,s,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(nx,nz):: s
       type(x3doperator1d), intent(in) :: x3dop
     END SUBROUTINE DERIVATIVE_YY
     SUBROUTINE DERIVATIVE_Z(t,u,s,x3dop,nx,ny,nz)
       use decomp_2d, only : mytype
       use x3d_operator_1d, only : x3doperator1d
       integer, intent(in) :: nx,ny,nz
       real(mytype), intent(out), dimension(nx,ny,nz) :: t
       real(mytype), intent(in), dimension(nx,ny,nz) :: u
       real(mytype), intent(out), dimension(nx,ny):: s
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

    use param, only : nclx1, ncly1, nclz1, nclxn, nclyn, nclzn

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

  end subroutine x3d_derive_finalize

subroutine derx_00(tx,ux,sx,x3dop,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
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
  call xthomas(tx, sx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derx_00

!********************************************************************
!
subroutine derx_ij(tx,ux,sx,ff,fs,fw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
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
subroutine derx_11(tx,ux,sx,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,sx,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)
  
end subroutine derx_11

!********************************************************************
!
subroutine derx_12(tx,ux,sx,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz            
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,sx,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derx_12

!********************************************************************
!
subroutine derx_21(tx,ux,sx,x3dop,nx,ny,nz) 

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,sx,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derx_21

!********************************************************************
!
subroutine derx_22(tx,ux,sx,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz):: sx
  type(x3doperator1d), intent(in) :: x3dop

  call derx_ij(tx,ux,sx,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derx_22

!********************************************************************
!
subroutine dery_00(ty,uy,sy,x3dop,ppy,nx,ny,nz) 
  !
  !********************************************************************

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
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
  call ythomas(ty, sy, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  ! Apply stretching if needed
  if (istret /= 0) then
     do concurrent (k=1:nz, j=1:ny, i=1:nx)
        ty(i,j,k) = ty(i,j,k) * ppy(j)
     enddo
  endif

end subroutine dery_00

!********************************************************************
!
subroutine dery_ij(ty,uy,sy,ff,fs,fw,ppy,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
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
subroutine dery_11(ty,uy,sy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,sy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,1,1)

end subroutine dery_11

!********************************************************************
!
subroutine dery_12(ty,uy,sy,x3dop,ppy,nx,ny,nz) 

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,sy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,1,2)

end subroutine dery_12

!********************************************************************
!
subroutine dery_21(ty,uy,sy,x3dop,ppy,nx,ny,nz) 

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,sy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,2,1)

end subroutine dery_21

!********************************************************************
!
subroutine dery_22(ty,uy,sy,x3dop,ppy,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz)  :: sy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  call dery_ij(ty,uy,sy,x3dop%f,x3dop%s,x3dop%w,ppy,nx,ny,nz,x3dop%npaire,2,2)

end subroutine dery_22

!********************************************************************
!
subroutine derz_00(tz,uz,sz,x3dop,nx,ny,nz)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx,ny,nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
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
  call zthomas(tz, sz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

end subroutine derz_00

!********************************************************************
!
subroutine derz_ij(tz,uz,sz,ff,fs,fw,nx,ny,nz,npaire,ncl1,ncln)

  use x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, npaire, ncl1, ncln
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
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
subroutine derz_11(tz,uz,sz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,sz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,1)

end subroutine derz_11

!********************************************************************
!
subroutine derz_12(tz,uz,sz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,sz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,1,2)

end subroutine derz_12

!********************************************************************
!
subroutine derz_21(tz,uz,sz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,sz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,1)

end subroutine derz_21

!********************************************************************
!
subroutine derz_22(tz,uz,sz,x3dop,nx,ny,nz)

  implicit none

  integer, intent(in) :: nx, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  call derz_ij(tz,uz,sz,x3dop%f,x3dop%s,x3dop%w,nx,ny,nz,x3dop%npaire,2,2)

end subroutine derz_22

end module x3d_derive
