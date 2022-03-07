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

module x3d_staggered

  use decomp_2d, only : mytype
  use x3d_operator_1d, only : x3doperator1d
  use param
  use thomas
  
  implicit none

  ! Make everything public unless declared private
  public


contains


!********************************************************************
!
subroutine derxvp(tx,ux,sx,x3dop,nx,nxm,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nxm,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: sx
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny)

        ! Compute r.h.s.
        tx(1,j,k) = acix6*(ux(2,j,k)-ux(1 ,j,k)) &
                  + bcix6*(ux(3,j,k)-ux(nx,j,k))
        tx(2,j,k) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                  + bcix6*(ux(4,j,k)-ux(1,j,k))
        do concurrent (i=3:nx-2)
           tx(i,j,k) = acix6*(ux(i+1,j,k)-ux(i  ,j,k)) &
                     + bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
        enddo
        tx(nx-1,j,k) = acix6*(ux(nx,j,k)-ux(nx-1,j,k)) &
                     + bcix6*(ux(1 ,j,k)-ux(nx-2,j,k))
        tx(nx  ,j,k) = acix6*(ux(1,j,k)-ux(nx  ,j,k)) &
                     + bcix6*(ux(2,j,k)-ux(nx-1,j,k))
     enddo

     ! Solve tri-diagonal system
     call xthomas(tx, sx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nxm = nx-1
     do concurrent (k=1:nz, j=1:ny)

        ! Compute r.h.s.
        if (x3dop%npaire==1) then
           tx(1,j,k) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                     + bcix6*(ux(3,j,k)-ux(2,j,k))
           tx(2,j,k) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                     + bcix6*(ux(4,j,k)-ux(1,j,k))
        else
           tx(1,j,k) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                     + bcix6*(ux(3,j,k)-two*ux(1,j,k)+ux(2,j,k))
           tx(2,j,k) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                     + bcix6*(ux(4,j,k)-ux(1,j,k))
        endif
        do concurrent (i=3:nxm-2)
           tx(i,j,k) = acix6*(ux(i+1,j,k)-ux(i  ,j,k)) &
                     + bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
        enddo
        if (x3dop%npaire==1) then
           tx(nxm-1,j,k) = acix6*(ux(nxm,j,k)-ux(nxm-1,j,k)) &
                         + bcix6*(ux(nx ,j,k)-ux(nxm-2,j,k))
           tx(nxm,j,k) = acix6*(ux(nx ,j,k)-ux(nxm  ,j,k)) &
                       + bcix6*(ux(nxm,j,k)-ux(nxm-1,j,k))
        else
           tx(nxm-1,j,k) = acix6*(ux(nxm,j,k)-ux(nxm-1,j,k)) &
                         + bcix6*(ux(nx ,j,k)-ux(nxm-2,j,k))
           tx(nxm,j,k) = acix6*(ux(nx,j,k)-ux(nxm,j,k)) &
                       + bcix6*(two*ux(nx,j,k)-ux(nxm,j,k)-ux(nxm-1,j,k))
        endif
     enddo

     ! Solve tri-diagonal system
     call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, nxm, ny, nz)

  endif

end subroutine derxvp

!********************************************************************
!
subroutine interxvp(tx,ux,sx,x3dop,nx,nxm,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nxm,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: sx
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny)

        ! Compute r.h.s.
        tx(1,j,k) = aicix6*(ux(2,j,k)+ux(1  ,j,k)) &
                  + bicix6*(ux(3,j,k)+ux(nx,j,k)) &
                  + cicix6*(ux(4,j,k)+ux(nx-1,j,k)) &
                  + dicix6*(ux(5,j,k)+ux(nx-2,j,k))
        tx(2,j,k) = aicix6*(ux(3,j,k)+ux(2 ,j,k)) &
                  + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                  + cicix6*(ux(5,j,k)+ux(nx,j,k)) &
                  + dicix6*(ux(6,j,k)+ux(nx-1,j,k))
        tx(3,j,k) = aicix6*(ux(4,j,k)+ux(3 ,j,k)) &
                  + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                  + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                  + dicix6*(ux(7,j,k)+ux(nx,j,k))
        do concurrent (i=4:nx-4)
           tx(i,j,k) = aicix6*(ux(i+1,j,k)+ux(i,j,k)) &
                     + bicix6*(ux(i+2,j,k)+ux(i-1,j,k)) &
                     + cicix6*(ux(i+3,j,k)+ux(i-2,j,k)) &
                     + dicix6*(ux(i+4,j,k)+ux(i-3,j,k))
        enddo
        tx(nx-3,j,k) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                     + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                     + cicix6*(ux(nx,j,k)+ux(nx-5,j,k)) &
                     + dicix6*(ux(1,j,k)+ux(nx-6,j,k))
        tx(nx-2,j,k) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                     + bicix6*(ux(nx ,j,k)+ux(nx-3,j,k)) &
                     + cicix6*(ux(1,j,k)+ux(nx-4,j,k)) &
                     + dicix6*(ux(2,j,k)+ux(nx-5,j,k))
        tx(nx-1,j,k) = aicix6*(ux(nx,j,k)+ux(nx-1,j,k)) &
                     + bicix6*(ux(1 ,j,k)+ux(nx-2,j,k)) &
                     + cicix6*(ux(2,j,k)+ux(nx-3,j,k)) &
                     + dicix6*(ux(3,j,k)+ux(nx-4,j,k))
        tx(nx  ,j,k) = aicix6*(ux(1,j,k)+ux(nx,j,k)) &
                     + bicix6*(ux(2,j,k)+ux(nx-1,j,k)) &
                     + cicix6*(ux(3,j,k)+ux(nx-2,j,k)) &
                     + dicix6*(ux(4,j,k)+ux(nx-3,j,k))
     enddo

     ! Solve tri-diagonal system
     call xthomas(tx, sx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny)

           ! Compute r.h.s.
           tx(1,j,k) = aicix6*(ux(2,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + dicix6*(ux(5,j,k)+ux(4,j,k))
           tx(2,j,k) = aicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + dicix6*(ux(6,j,k)+ux(3,j,k))
           tx(3,j,k) = aicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(7,j,k)+ux(2,j,k))
           do concurrent (i=4:nxm-3)
              tx(i,j,k) = aicix6*(ux(i+1,j,k)+ux(i,j,k)) &
                        + bicix6*(ux(i+2,j,k)+ux(i-1,j,k)) &
                        + cicix6*(ux(i+3,j,k)+ux(i-2,j,k)) &
                        + dicix6*(ux(i+4,j,k)+ux(i-3,j,k))
           enddo
           tx(nxm-2,j,k) = aicix6*(ux(nxm-1,j,k)+ux(nxm-2,j,k)) &
                         + bicix6*(ux(nxm,j,k)+ux(nxm-3,j,k)) &
                         + cicix6*(ux(nx,j,k)+ux(nxm-4,j,k)) &
                         + dicix6*(ux(nxm,j,k)+ux(nxm-5,j,k))
           tx(nxm-1,j,k) = aicix6*(ux(nxm,j,k)+ux(nxm-1,j,k)) &
                         + bicix6*(ux(nx,j,k)+ux(nxm-2,j,k)) &
                         + cicix6*(ux(nxm,j,k)+ux(nxm-3,j,k)) &
                         + dicix6*(ux(nxm-1,j,k)+ux(nxm-4,j,k))
           tx(nxm  ,j,k) = aicix6*(ux(nx,j,k)+ux(nxm,j,k)) &
                         + bicix6*(ux(nxm,j,k)+ux(nxm-1,j,k)) &
                         + cicix6*(ux(nxm-1,j,k)+ux(nxm-2,j,k)) &
                         + dicix6*(ux(nxm-2,j,k)+ux(nxm-3,j,k))
        enddo

        ! Solve tri-diagonal system
        call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, nxm, ny, nz)

     endif
  endif

end subroutine interxvp

!********************************************************************
!
subroutine derxpv(tx,ux,sx,x3dop,nxm,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nxm,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: sx
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny)

        ! Compute r.h.s.
        tx(1,j,k) = acix6*(ux(1,j,k)-ux(nx  ,j,k)) &
                  + bcix6*(ux(2,j,k)-ux(nx-1,j,k))
        tx(2,j,k) = acix6*(ux(2,j,k)-ux(1 ,j,k)) &
                  + bcix6*(ux(3,j,k)-ux(nx,j,k))
        do concurrent (i=3:nx-2)
           tx(i,j,k) = acix6*(ux(i,j,k)-ux(i-1,j,k)) &
                     + bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
        enddo
        tx(nx-1,j,k) = acix6*(ux(nx-1,j,k)-ux(nx-2,j,k)) &
                     + bcix6*(ux(nx ,j,k)-ux(nx-3,j,k))
        tx(nx  ,j,k) = acix6*(ux(nx,j,k)-ux(nx-1,j,k)) &
                     + bcix6*(ux(1,j,k)-ux(nx-2,j,k))
     enddo

     ! Solve tri-diagonal system
     call xthomas(tx, sx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny)

           ! Compute r.h.s.
           tx(1,j,k) = zero
           tx(2,j,k) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                     + bcix6*(ux(3,j,k)-ux(1,j,k))
           do concurrent (i=3:nx-2)
              tx(i,j,k) = acix6*(ux(i,j,k)-ux(i-1,j,k)) &
                        + bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
           enddo
           tx(nx-1,j,k) = acix6*(ux(nx-1,j,k)-ux(nx-2,j,k)) &
                        + bcix6*(ux(nx-1,j,k)-ux(nx-3,j,k))
           tx(nx,j,k) = zero
        enddo

        ! Solve tri-diagonal system
        call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

end subroutine derxpv

!********************************************************************
!
subroutine interxpv(tx,ux,sx,x3dop,nxm,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nxm,ny,nz) :: ux
  real(mytype), intent(out), dimension(ny,nz) :: sx
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny)

        ! Compute r.h.s.
        tx(1,j,k) = aicix6*(ux(1,j,k)+ux(nx  ,j,k)) &
                  + bicix6*(ux(2,j,k)+ux(nx-1,j,k)) &
                  + cicix6*(ux(3,j,k)+ux(nx-2,j,k)) &
                  + dicix6*(ux(4,j,k)+ux(nx-3,j,k))
        tx(2,j,k) = aicix6*(ux(2,j,k)+ux(1 ,j,k)) &
                  + bicix6*(ux(3,j,k)+ux(nx,j,k)) &
                  + cicix6*(ux(4,j,k)+ux(nx-1,j,k)) &
                  + dicix6*(ux(5,j,k)+ux(nx-2,j,k))
        tx(3,j,k) = aicix6*(ux(3,j,k)+ux(2 ,j,k)) &
                  + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                  + cicix6*(ux(5,j,k)+ux(nx,j,k)) &
                  + dicix6*(ux(6,j,k)+ux(nx-1,j,k))
        tx(4,j,k) = aicix6*(ux(4,j,k)+ux(3 ,j,k)) &
                  + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                  + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                  + dicix6*(ux(7,j,k)+ux(nx,j,k))
        do concurrent (i=5:nx-3)
           tx(i,j,k) = aicix6*(ux(i,j,k)+ux(i-1,j,k)) &
                     + bicix6*(ux(i+1,j,k)+ux(i-2,j,k)) &
                     + cicix6*(ux(i+2,j,k)+ux(i-3,j,k)) &
                     + dicix6*(ux(i+3,j,k)+ux(i-4,j,k))
        enddo
        tx(nx-2,j,k) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                     + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                     + cicix6*(ux(nx,j,k)+ux(nx-5,j,k)) &
                     + dicix6*(ux(1,j,k)+ux(nx-6,j,k))
        tx(nx-1,j,k) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                     + bicix6*(ux(nx ,j,k)+ux(nx-3,j,k)) &
                     + cicix6*(ux(1,j,k)+ux(nx-4,j,k)) &
                     + dicix6*(ux(2,j,k)+ux(nx-5,j,k))
        tx(nx  ,j,k) = aicix6*(ux(nx,j,k)+ux(nx-1,j,k)) &
                     + bicix6*(ux(1,j,k)+ux(nx-2,j,k)) &
                     + cicix6*(ux(2,j,k)+ux(nx-3,j,k)) &
                     + dicix6*(ux(3,j,k)+ux(nx-4,j,k))
     enddo

     ! Solve tri-diagonal system
     call xthomas(tx, sx, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny)

           ! Compute r.h.s.
           tx(1,j,k) = aicix6*(ux(1,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(2,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(3,j,k)+ux(3,j,k)) &
                     + dicix6*(ux(4,j,k)+ux(4,j,k))
           tx(2,j,k) = aicix6*(ux(2,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(3,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(4,j,k)+ux(2,j,k)) &
                     + dicix6*(ux(5,j,k)+ux(3,j,k))
           tx(3,j,k) = aicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(5,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(6,j,k)+ux(2,j,k))
           tx(4,j,k) = aicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(7,j,k)+ux(1,j,k))
           do concurrent (i=5:nx-4)
              tx(i,j,k) = aicix6*(ux(i,j,k)+ux(i-1,j,k)) &
                        + bicix6*(ux(i+1,j,k)+ux(i-2,j,k)) &
                        + cicix6*(ux(i+2,j,k)+ux(i-3,j,k)) &
                        + dicix6*(ux(i+3,j,k)+ux(i-4,j,k))
           enddo
           tx(nx-3,j,k) = aicix6*(ux(nx-3,j,k)+ux(nx-4,j,k)) &
                        + bicix6*(ux(nx-2,j,k)+ux(nx-5,j,k)) &
                        + cicix6*(ux(nx-1,j,k)+ux(nx-6,j,k)) &
                        + dicix6*(ux(nx-1,j,k)+ux(nx-7,j,k))
           tx(nx-2,j,k) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                        + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                        + cicix6*(ux(nx-1,j,k)+ux(nx-5,j,k)) &
                        + dicix6*(ux(nx-2,j,k)+ux(nx-6,j,k))
           tx(nx-1,j,k) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                        + bicix6*(ux(nx-1,j,k)+ux(nx-3,j,k)) &
                        + cicix6*(ux(nx-2,j,k)+ux(nx-4,j,k)) &
                        + dicix6*(ux(nx-3,j,k)+ux(nx-5,j,k))
           tx(nx  ,j,k) = aicix6*(ux(nx-1,j,k)+ux(nx-1,j,k)) &
                        + bicix6*(ux(nx-2,j,k)+ux(nx-2,j,k)) &
                        + cicix6*(ux(nx-3,j,k)+ux(nx-3,j,k)) &
                        + dicix6*(ux(nx-4,j,k)+ux(nx-4,j,k))
        enddo

        ! Solve tri-diagonal system
        call xthomas(tx, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

end subroutine interxpv

!********************************************************************
!
subroutine interyvp(ty,uy,sy,x3dop,nx,ny,nym,nz)

  USE x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,nym,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz) :: sy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz)

        ! Compute r.h.s.
        do concurrent (i=1:nx)
           ty(i,1,k) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                     + biciy6*(uy(i,3,k)+uy(i,ny,k)) &
                     + ciciy6*(uy(i,4,k)+uy(i,ny-1,k)) &
                     + diciy6*(uy(i,5,k)+uy(i,ny-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                     + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                     + ciciy6*(uy(i,5,k)+uy(i,ny,k)) &
                     + diciy6*(uy(i,6,k)+uy(i,ny-1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,3,k) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                     + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                     + diciy6*(uy(i,7,k)+uy(i,ny,k))
        enddo
        do concurrent (j=4:ny-4, i=1:nx)
           ty(i,j,k) = aiciy6*(uy(i,j+1,k)+uy(i,j,k)) &
                     + biciy6*(uy(i,j+2,k)+uy(i,j-1,k)) &
                     + ciciy6*(uy(i,j+3,k)+uy(i,j-2,k)) &
                     + diciy6*(uy(i,j+4,k)+uy(i,j-3,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-3,k) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                        + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                        + ciciy6*(uy(i,ny,k)+uy(i,ny-5,k)) &
                        + diciy6*(uy(i,1,k)+uy(i,ny-6,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-2,k) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                        + biciy6*(uy(i,ny,k)+uy(i,ny-3,k)) &
                        + ciciy6*(uy(i,1,k)+uy(i,ny-4,k)) &
                        + diciy6*(uy(i,2,k)+uy(i,ny-5,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = aiciy6*(uy(i,ny,k)+uy(i,ny-1,k)) &
                        + biciy6*(uy(i,1,k)+uy(i,ny-2,k)) &
                        + ciciy6*(uy(i,2,k)+uy(i,ny-3,k)) &
                        + diciy6*(uy(i,3,k)+uy(i,ny-4,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny  ,k) = aiciy6*(uy(i,1,k)+uy(i,ny,k)) &
                        + biciy6*(uy(i,2,k)+uy(i,ny-1,k)) &
                        + ciciy6*(uy(i,3,k)+uy(i,ny-2,k)) &
                        + diciy6*(uy(i,4,k)+uy(i,ny-3,k))
        enddo
     enddo

     ! Solve tri-diagonal system
     call ythomas(ty, sy, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz)

           ! Compute r.h.s.
           do concurrent (i=1:nx)
              ty(i,1,k) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                        + biciy6*(uy(i,3,k)+uy(i,2,k)) &
                        + ciciy6*(uy(i,4,k)+uy(i,3,k)) &
                        + diciy6*(uy(i,5,k)+uy(i,4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                        + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                        + ciciy6*(uy(i,5,k)+uy(i,2,k)) &
                        + diciy6*(uy(i,6,k)+uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                        + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                        + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                        + diciy6*(uy(i,7,k)+uy(i,2,k))
           enddo
           do concurrent (j=4:nym-3, i=1:nx)
              ty(i,j,k) = aiciy6*(uy(i,j+1,k)+uy(i,j,k)) &
                        + biciy6*(uy(i,j+2,k)+uy(i,j-1,k)) &
                        + ciciy6*(uy(i,j+3,k)+uy(i,j-2,k)) &
                        + diciy6*(uy(i,j+4,k)+uy(i,j-3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,nym-2,k) = aiciy6*(uy(i,nym-1,k)+uy(i,nym-2,k)) &
                            + biciy6*(uy(i,nym,k)+uy(i,nym-3,k)) &
                            + ciciy6*(uy(i,ny,k)+uy(i,nym-4,k)) &
                            + diciy6*(uy(i,nym,k)+uy(i,nym-5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,nym-1,k) = aiciy6*(uy(i,nym,k)+uy(i,nym-1,k)) &
                            + biciy6*(uy(i,ny,k)+uy(i,nym-2,k)) &
                            + ciciy6*(uy(i,nym,k)+uy(i,nym-3,k)) &
                            + diciy6*(uy(i,nym-1,k)+uy(i,nym-4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,nym  ,k) = aiciy6*(uy(i,ny,k)+uy(i,nym,k)) &
                            + biciy6*(uy(i,nym,k)+uy(i,nym-1,k)) &
                            + ciciy6*(uy(i,nym-1,k)+uy(i,nym-2,k)) &
                            + diciy6*(uy(i,nym-2,k)+uy(i,nym-3,k))
           enddo
        enddo

        ! Solve tri-diagonal system
        call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, nx, nym, nz)

     endif
  endif

end subroutine interyvp

!********************************************************************
!
subroutine deryvp(ty,uy,sy,x3dop,ppyi,nx,ny,nym,nz)

  USE x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,nym,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz) :: sy
  real(mytype), intent(in), dimension(nym) :: ppyi
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz)

        ! Compute r.h.s.
        do concurrent (i=1:nx)
           ty(i,1,k) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                     + bciy6*(uy(i,3,k)-uy(i,ny,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = aciy6*(uy(i,3,k)-uy(i,2,k)) &
                     + bciy6*(uy(i,4,k)-uy(i,1,k))
        enddo
        do concurrent (j=3:ny-2, i=1:nx)
           ty(i,j,k) = aciy6*(uy(i,j+1,k)-uy(i,j,k)) &
                     + bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = aciy6*(uy(i,ny,k)-uy(i,ny-1,k)) &
                        + bciy6*(uy(i,1,k)-uy(i,ny-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny  ,k) = aciy6*(uy(i,1,k)-uy(i,ny,k)) &
                        + bciy6*(uy(i,2,k)-uy(i,ny-1,k))
        enddo
     enddo

     ! Solve tri-diagonal system
     call ythomas(ty, sy, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nym = ny-1
     if (x3dop%npaire==0) then
        do concurrent (k=1:nz)

           ! Compute r.h.s.
           do concurrent (i=1:nx)
              ty(i,1,k) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                        + bciy6*(uy(i,3,k)-two*uy(i,1,k)+uy(i,2,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = aciy6*(uy(i,3,k)-uy(i,2,k)) &
                        + bciy6*(uy(i,4,k)-uy(i,1,k))
           enddo
           do concurrent (j=3:nym-2, i=1:nx)
              ty(i,j,k) = aciy6*(uy(i,j+1,k)-uy(i,j,k)) &
                        + bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,nym-1,k) = aciy6*(uy(i,nym,k)-uy(i,nym-1,k)) &
                            + bciy6*(uy(i,ny,k)-uy(i,nym-2,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,nym  ,k) = aciy6*(uy(i,ny,k)-uy(i,nym,k)) &
                            + bciy6*(two*uy(i,ny,k)-uy(i,nym,k)-uy(i,nym-1,k))
           enddo
        enddo

        ! Solve tri-diagonal system
        call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, nx, nym, nz)

     endif
  endif

  if (istret /= 0) then
     do concurrent (k=1:nz, j=1:nym, i=1:nx)
        ty(i,j,k) = ty(i,j,k) * ppyi(j)
     enddo
  endif

end subroutine deryvp

!********************************************************************
!
subroutine interypv(ty,uy,sy,x3dop,nx,nym,ny,nz)

  USE x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,nym,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz) :: sy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz)

        ! Compute r.h.s.
        do concurrent (i=1:nx)
           ty(i,1,k) = aiciy6*(uy(i,1,k)+uy(i,ny,k)) &
                     + biciy6*(uy(i,2,k)+uy(i,ny-1,k)) &
                     + ciciy6*(uy(i,3,k)+uy(i,ny-2,k)) &
                     + diciy6*(uy(i,4,k)+uy(i,ny-3,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                     + biciy6*(uy(i,3,k)+uy(i,ny,k)) &
                     + ciciy6*(uy(i,4,k)+uy(i,ny-1,k)) &
                     + diciy6*(uy(i,5,k)+uy(i,ny-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,3,k) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                     + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                     + ciciy6*(uy(i,5,k)+uy(i,ny,k)) &
                     + diciy6*(uy(i,6,k)+uy(i,ny-1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,4,k) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                     + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                     + diciy6*(uy(i,7,k)+uy(i,ny,k))
        enddo
        do concurrent (j=5:ny-3, i=1:nx)
           ty(i,j,k) = aiciy6*(uy(i,j,k)+uy(i,j-1,k)) &
                     + biciy6*(uy(i,j+1,k)+uy(i,j-2,k)) &
                     + ciciy6*(uy(i,j+2,k)+uy(i,j-3,k)) &
                     + diciy6*(uy(i,j+3,k)+uy(i,j-4,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-2,k) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                        + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                        + ciciy6*(uy(i,ny,k)+uy(i,ny-5,k)) &
                        + diciy6*(uy(i,1,k)+uy(i,ny-6,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                        + biciy6*(uy(i,ny,k)+uy(i,ny-3,k)) &
                        + ciciy6*(uy(i,1,k)+uy(i,ny-4,k)) &
                        + diciy6*(uy(i,2,k)+uy(i,ny-5,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny  ,k) = aiciy6*(uy(i,ny,k)+uy(i,ny-1,k)) &
                        + biciy6*(uy(i,1,k)+uy(i,ny-2,k)) &
                        + ciciy6*(uy(i,2,k)+uy(i,ny-3,k)) &
                        + diciy6*(uy(i,3,k)+uy(i,ny-4,k))
        enddo
     enddo

     ! Solve tri-diagonal system
     call ythomas(ty, sy, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz)

           ! Compute r.h.s.
           do concurrent (i=1:nx)
              ty(i,1,k) = aiciy6*(uy(i,1,k)+uy(i,1,k)) &
                        + biciy6*(uy(i,2,k)+uy(i,2,k)) &
                        + ciciy6*(uy(i,3,k)+uy(i,3,k)) &
                        + diciy6*(uy(i,4,k)+uy(i,4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                        + biciy6*(uy(i,3,k)+uy(i,1,k)) &
                        + ciciy6*(uy(i,4,k)+uy(i,2,k)) &
                        + diciy6*(uy(i,5,k)+uy(i,3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,3,k) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                        + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                        + ciciy6*(uy(i,5,k)+uy(i,1,k)) &
                        + diciy6*(uy(i,6,k)+uy(i,2,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,4,k) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                        + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                        + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                        + diciy6*(uy(i,7,k)+uy(i,1,k))
           enddo
           do concurrent (j=5:ny-4, i=1:nx)
              ty(i,j,k) = aiciy6*(uy(i,j,k)+uy(i,j-1,k)) &
                        + biciy6*(uy(i,j+1,k)+uy(i,j-2,k)) &
                        + ciciy6*(uy(i,j+2,k)+uy(i,j-3,k)) &
                        + diciy6*(uy(i,j+3,k)+uy(i,j-4,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-3,k) = aiciy6*(uy(i,ny-3,k)+uy(i,ny-4,k)) &
                           + biciy6*(uy(i,ny-2,k)+uy(i,ny-5,k)) &
                           + ciciy6*(uy(i,ny-1,k)+uy(i,ny-6,k)) &
                           + diciy6*(uy(i,ny-1,k)+uy(i,ny-7,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-2,k) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                           + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                           + ciciy6*(uy(i,ny-1,k)+uy(i,ny-5,k)) &
                           + diciy6*(uy(i,ny-2,k)+uy(i,ny-6,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                           + biciy6*(uy(i,ny-1,k)+uy(i,ny-3,k)) &
                           + ciciy6*(uy(i,ny-2,k)+uy(i,ny-4,k)) &
                           + diciy6*(uy(i,ny-3,k)+uy(i,ny-5,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny  ,k) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-1,k)) &
                           + biciy6*(uy(i,ny-2,k)+uy(i,ny-2,k)) &
                           + ciciy6*(uy(i,ny-3,k)+uy(i,ny-3,k)) &
                           + diciy6*(uy(i,ny-4,k)+uy(i,ny-4,k))
           enddo
        enddo

        ! Solve tri-diagonal system
        call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

end subroutine interypv

!********************************************************************
!
subroutine derypv(ty,uy,sy,x3dop,ppy,nx,nym,ny,nz)

  USE x3d_operator_y_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,nym,nz) :: uy
  real(mytype), intent(out), dimension(nx,nz) :: sy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz)

        ! Compute r.h.s.
        do concurrent (i=1:nx)
           ty(i,1,k) = aciy6*(uy(i,1,k)-uy(i,ny,k)) &
                     + bciy6*(uy(i,2,k)-uy(i,ny-1,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,2,k) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                     + bciy6*(uy(i,3,k)-uy(i,ny,k))
        enddo
        do concurrent (j=3:ny-2, i=1:nx)
           ty(i,j,k) = aciy6*(uy(i,j,k)-uy(i,j-1,k)) &
                     + bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny-1,k) = aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k)) &
                        + bciy6*(uy(i,ny,k)-uy(i,ny-3,k))
        enddo
        do concurrent (i=1:nx)
           ty(i,ny,k) = aciy6*(uy(i,ny,k)-uy(i,ny-1,k)) &
                      + bciy6*(uy(i,1,k)-uy(i,ny-2,k))
        enddo
     enddo

     ! Solve tri-diagonal system
     call ythomas(ty, sy, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz)

           ! Compute r.h.s.
           do concurrent (i=1:nx)
              ty(i,1,k) = zero
           enddo
           do concurrent (i=1:nx)
              ty(i,2,k) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                        + bciy6*(uy(i,3,k)-uy(i,1,k))
           enddo
           do concurrent (j=3:ny-2, i=1:nx)
              ty(i,j,k) = aciy6*(uy(i,j,k)-uy(i,j-1,k)) &
                        + bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny-1,k) = aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k)) &
                           + bciy6*(uy(i,ny-1,k)-uy(i,ny-3,k))
           enddo
           do concurrent (i=1:nx)
              ty(i,ny,k) = zero
           enddo
        enddo

        ! Solve tri-diagonal system
        call ythomas(ty, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

  if (istret /= 0) then
     do concurrent (k=1:nz, j=1:ny, i=1:nx)
        ty(i,j,k) = ty(i,j,k) * ppy(j)
     enddo
  endif

end subroutine derypv

!********************************************************************
!
subroutine derzvp(tz,uz,sz,x3dop,nx,ny,nz,nzm)

  USE x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nzm) :: tz
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

  if (nclz) then
     ! nzm = nz

     ! Compute r.h.s.
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                  + bciz6*(uz(i,j,3)-uz(i,j,nz))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = aciz6*(uz(i,j,3)-uz(i,j,2)) &
                  + bciz6*(uz(i,j,4)-uz(i,j,1))
     enddo
     do concurrent (k=3:nz-2, j=1:ny, i=1:nx)
        tz(i,j,k) = aciz6*(uz(i,j,k+1)-uz(i,j,k)) &
                  + bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                     + bciz6*(uz(i,j,1)-uz(i,j,nz-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz  ) = aciz6*(uz(i,j,1)-uz(i,j,nz)) &
                     + bciz6*(uz(i,j,2)-uz(i,j,nz-1))
     enddo

     ! Solve tri-diagonal system
     call zthomas(tz, sz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nzm = nz-1

     ! Compute r.h.s.
     if (x3dop%npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-uz(i,j,2))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = aciz6*(uz(i,j,3)-uz(i,j,2))&
                     + bciz6*(uz(i,j,4)-uz(i,j,1))
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-two*uz(i,j,1)+uz(i,j,2))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = aciz6*(uz(i,j,3)-uz(i,j,2)) &
                     + bciz6*(uz(i,j,4)-uz(i,j,1))
        enddo
     endif
     do concurrent (k=3:nzm-2, j=1:ny, i=1:nx)
        tz(i,j,k) = aciz6*(uz(i,j,k+1)-uz(i,j,k)) &
                  + bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
     enddo
     if (x3dop%npaire==1) then
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm-1) = aciz6*(uz(i,j,nzm)-uz(i,j,nzm-1)) &
                         + bciz6*(uz(i,j,nz)-uz(i,j,nzm-2))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm  ) = aciz6*(uz(i,j,nz)-uz(i,j,nzm)) &
                         + bciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))
        enddo
     else
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                         + bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm  ) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                         + bciz6*(two*uz(i,j,nz)-uz(i,j,nz-1)-uz(i,j,nz-2))
        enddo
     endif

     ! Solve tri-diagonal system
     call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, nx, ny, nzm)

  endif

end subroutine derzvp

!********************************************************************
!
subroutine interzvp(tz,uz,sz,x3dop,nx,ny,nz,nzm)

  USE x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nzm) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = uz(i,j,k)
     enddo
     return
  endif

  if (nclz) then
     ! nzm = nz

     ! Compute r.h.s.
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                  + biciz6*(uz(i,j,3)+uz(i,j,nz)) &
                  + ciciz6*(uz(i,j,4)+uz(i,j,nz-1)) &
                  + diciz6*(uz(i,j,5)+uz(i,j,nz-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                  + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                  + ciciz6*(uz(i,j,5)+uz(i,j,nz)) &
                  + diciz6*(uz(i,j,6)+uz(i,j,nz-1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,3) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                  + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                  + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                  + diciz6*(uz(i,j,7)+uz(i,j,nz))
     enddo
     do concurrent (k=4:nz-4, j=1:ny, i=1:nx)
        tz(i,j,k) = aiciz6*(uz(i,j,k+1)+uz(i,j,k)) &
                  + biciz6*(uz(i,j,k+2)+uz(i,j,k-1)) &
                  + ciciz6*(uz(i,j,k+3)+uz(i,j,k-2)) &
                  + diciz6*(uz(i,j,k+4)+uz(i,j,k-3))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-3) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                     + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                     + ciciz6*(uz(i,j,nz)+uz(i,j,nz-5)) &
                     + diciz6*(uz(i,j,1)+uz(i,j,nz-6))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-2) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                     + biciz6*(uz(i,j,nz)+uz(i,j,nz-3)) &
                     + ciciz6*(uz(i,j,1)+uz(i,j,nz-4)) &
                     + diciz6*(uz(i,j,2)+uz(i,j,nz-5))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = aiciz6*(uz(i,j,nz)+uz(i,j,nz-1)) &
                     + biciz6*(uz(i,j,1)+uz(i,j,nz-2)) &
                     + ciciz6*(uz(i,j,2)+uz(i,j,nz-3)) &
                     + diciz6*(uz(i,j,3)+uz(i,j,nz-4))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz  ) = aiciz6*(uz(i,j,1)+uz(i,j,nz)) &
                     + biciz6*(uz(i,j,2)+uz(i,j,nz-1)) &
                     + ciciz6*(uz(i,j,3)+uz(i,j,nz-2)) &
                     + diciz6*(uz(i,j,4)+uz(i,j,nz-3))
     enddo

     ! Solve tri-diagonal system
     call zthomas(tz, sz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then

        ! Compute r.h.s.
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + diciz6*(uz(i,j,5)+uz(i,j,4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                     + ciciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + diciz6*(uz(i,j,6)+uz(i,j,3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,3) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,7)+uz(i,j,2))
        enddo
        do concurrent (k=4:nzm-3, j=1:ny, i=1:nx)
           tz(i,j,k) = aiciz6*(uz(i,j,k+1)+uz(i,j,k)) &
                     + biciz6*(uz(i,j,k+2)+uz(i,j,k-1)) &
                     + ciciz6*(uz(i,j,k+3)+uz(i,j,k-2)) &
                     + diciz6*(uz(i,j,k+4)+uz(i,j,k-3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm-2) = aiciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-2)) &
                         + biciz6*(uz(i,j,nzm)+uz(i,j,nzm-3)) &
                         + ciciz6*(uz(i,j,nz)+uz(i,j,nzm-4)) &
                         + diciz6*(uz(i,j,nzm)+uz(i,j,nzm-5))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm-1) = aiciz6*(uz(i,j,nzm)+uz(i,j,nzm-1)) &
                         + biciz6*(uz(i,j,nz)+uz(i,j,nzm-2)) &
                         + ciciz6*(uz(i,j,nzm)+uz(i,j,nzm-3)) &
                         + diciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nzm) = aiciz6*(uz(i,j,nz)+uz(i,j,nzm)) &
                       + biciz6*(uz(i,j,nzm)+uz(i,j,nzm-1)) &
                       + ciciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-2)) &
                       + diciz6*(uz(i,j,nzm-2)+uz(i,j,nzm-3))
        enddo

        ! Solve tri-diagonal system
        call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, nx, ny, nzm)

     endif
  endif

end subroutine interzvp

!********************************************************************
!
subroutine derzpv(tz,uz,sz,x3dop,nx,ny,nzm,nz)

  USE x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, nzm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nzm) :: uz
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

  if (nclz) then
     ! nzm = nz

     ! Compute r.h.s.
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = aciz6*(uz(i,j,1)-uz(i,j,nz)) &
                  + bciz6*(uz(i,j,2)-uz(i,j,nz-1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                  + bciz6*(uz(i,j,3)-uz(i,j,nz))
     enddo
     do concurrent (k=3:nz-2, j=1:ny, i=1:nx)
        tz(i,j,k) = aciz6*(uz(i,j,k)-uz(i,j,k-1)) &
                  + bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                     + bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                   + bciz6*(uz(i,j,1)-uz(i,j,nz-2))
     enddo

     ! Solve tri-diagonal system
     call zthomas(tz, sz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then

        ! Compute r.h.s.
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = zero
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-uz(i,j,1))
        enddo
        do concurrent (k=3:nz-2, j=1:ny, i=1:nx)
           tz(i,j,k) = aciz6*(uz(i,j,k)-uz(i,j,k-1)) &
                     + bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                        + bciz6*(uz(i,j,nz-1)-uz(i,j,nz-3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz) = zero
        enddo

        ! Solve tri-diagonal system
        call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

end subroutine derzpv

!********************************************************************
!
subroutine interzpv(tz,uz,sz,x3dop,nx,ny,nzm,nz)

  USE x3d_operator_z_data

  implicit none

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nzm) :: uz
  real(mytype), intent(out), dimension(nx,ny) :: sz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = uz(i,j,k)
     enddo
     return
  endif

  if (nclz) then
     ! nzm = nz

     ! Compute r.h.s.
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,1) = aiciz6*(uz(i,j,1)+uz(i,j,nz)) &
                  + biciz6*(uz(i,j,2)+uz(i,j,nz-1)) &
                  + ciciz6*(uz(i,j,3)+uz(i,j,nz-2)) &
                  + diciz6*(uz(i,j,4)+uz(i,j,nz-3))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,2) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                  + biciz6*(uz(i,j,3)+uz(i,j,nz)) &
                  + ciciz6*(uz(i,j,4)+uz(i,j,nz-1)) &
                  + diciz6*(uz(i,j,5)+uz(i,j,nz-2))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,3) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                  + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                  + ciciz6*(uz(i,j,5)+uz(i,j,nz)) &
                  + diciz6*(uz(i,j,6)+uz(i,j,nz-1))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,4) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                  + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                  + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                  + diciz6*(uz(i,j,7)+uz(i,j,nz))
     enddo
     do concurrent (k=5:nz-3, j=1:ny, i=1:nx)
        tz(i,j,k) = aiciz6*(uz(i,j,k)+uz(i,j,k-1)) &
                  + biciz6*(uz(i,j,k+1)+uz(i,j,k-2)) &
                  + ciciz6*(uz(i,j,k+2)+uz(i,j,k-3)) &
                  + diciz6*(uz(i,j,k+3)+uz(i,j,k-4))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-2) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                     + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                     + ciciz6*(uz(i,j,nz)+uz(i,j,nz-5)) &
                     + diciz6*(uz(i,j,1)+uz(i,j,nz-6))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz-1) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                     + biciz6*(uz(i,j,nz)+uz(i,j,nz-3)) &
                     + ciciz6*(uz(i,j,1)+uz(i,j,nz-4)) &
                     + diciz6*(uz(i,j,2)+uz(i,j,nz-5))
     enddo
     do concurrent (j=1:ny, i=1:nx)
        tz(i,j,nz  ) = aiciz6*(uz(i,j,nz)+uz(i,j,nz-1)) &
                     + biciz6*(uz(i,j,1)+uz(i,j,nz-2)) &
                     + ciciz6*(uz(i,j,2)+uz(i,j,nz-3)) &
                     + diciz6*(uz(i,j,3)+uz(i,j,nz-4))
     enddo

     ! Solve tri-diagonal system
     call zthomas(tz, sz, x3dop%f, x3dop%s, x3dop%w, x3dop%periodic, x3dop%alfa, nx, ny, nz)

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then

        ! Compute r.h.s.
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,1) = aiciz6*(uz(i,j,1)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,2)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,3)+uz(i,j,3)) &
                     + diciz6*(uz(i,j,4)+uz(i,j,4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,2) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,3)+uz(i,j,1))&
                     + ciciz6*(uz(i,j,4)+uz(i,j,2))&
                     + diciz6*(uz(i,j,5)+uz(i,j,3))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,3) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                     + ciciz6*(uz(i,j,5)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,6)+uz(i,j,2))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,4) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,7)+uz(i,j,1))
        enddo
        do concurrent (k=5:nz-4, j=1:ny, i=1:nx)
           tz(i,j,k) = aiciz6*(uz(i,j,k)+uz(i,j,k-1)) &
                     + biciz6*(uz(i,j,k+1)+uz(i,j,k-2)) &
                     + ciciz6*(uz(i,j,k+2)+uz(i,j,k-3)) &
                     + diciz6*(uz(i,j,k+3)+uz(i,j,k-4))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-3) = aiciz6*(uz(i,j,nz-3)+uz(i,j,nz-4)) &
                        + biciz6*(uz(i,j,nz-2)+uz(i,j,nz-5)) &
                        + ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-6)) &
                        + diciz6*(uz(i,j,nz-1)+uz(i,j,nz-7))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-2) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                        + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                        + ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-5)) &
                        + diciz6*(uz(i,j,nz-2)+uz(i,j,nz-6))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz-1) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                        + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-3)) &
                        + ciciz6*(uz(i,j,nz-2)+uz(i,j,nz-4)) &
                        + diciz6*(uz(i,j,nz-3)+uz(i,j,nz-5))
        enddo
        do concurrent (j=1:ny, i=1:nx)
           tz(i,j,nz  ) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-1)) &
                        + biciz6*(uz(i,j,nz-2)+uz(i,j,nz-2)) &
                        + ciciz6*(uz(i,j,nz-3)+uz(i,j,nz-3)) &
                        + diciz6*(uz(i,j,nz-4)+uz(i,j,nz-4))
        enddo

        ! Solve tri-diagonal system
        call zthomas(tz, x3dop%f, x3dop%s, x3dop%w, nx, ny, nz)

     endif
  endif

end subroutine interzpv

end module x3d_staggered
