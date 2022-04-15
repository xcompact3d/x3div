!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

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
subroutine derxvp(tx,ux,x3dop,nx,nxm,ny,nz)

   use x3d_operator_x_data

   implicit none

   !$acc routine(thomas1d_0, thomas1d_12) seq

   ! Arguments
   integer, intent(in) :: nx, nxm, ny, nz
   real(mytype), intent(out), dimension(nxm,ny,nz) :: tx
   real(mytype), intent(in), dimension(nx,ny,nz) :: ux
   type(x3doperator1d), intent(in) :: x3dop

   ! Local variables
   integer :: i, j, k
   real(mytype), dimension(nxm) :: buffer, ff, ss, ww, pp

   do concurrent (i=1:nxm)
      ff(i) = x3dop%f(i)
      ss(i) = x3dop%s(i)
      ww(i) = x3dop%w(i)
      if (allocated(x3dop%periodic)) pp(i) = x3dop%periodic(i)
   end do

   if (nclx) then
      ! nxm = nx
      do concurrent (k=1:nz, j=1:ny) local(buffer)
         ! Compute r.h.s.
         buffer(1) = acix6*(ux(2,j,k)-ux(1 ,j,k)) &
                   + bcix6*(ux(3,j,k)-ux(nx,j,k))
         buffer(2) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                   + bcix6*(ux(4,j,k)-ux(1,j,k))
         do concurrent (i=3:nx-2)
            buffer(i) = acix6*(ux(i+1,j,k)-ux(i  ,j,k)) &
                      + bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         buffer(nx-1) = acix6*(ux(nx,j,k)-ux(nx-1,j,k)) &
                      + bcix6*(ux(1 ,j,k)-ux(nx-2,j,k))
         buffer(nx  ) = acix6*(ux(1,j,k)-ux(nx  ,j,k)) &
                      + bcix6*(ux(2,j,k)-ux(nx-1,j,k))

         ! Solve tri-diagonal system
         call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nx)
         do concurrent (i=1:nx)
            tx(i,j,k) = buffer(i)
         enddo
      enddo

   else
      ! nxm = nx-1
      do concurrent (k=1:nz, j=1:ny)

         ! Compute r.h.s.
         if (x3dop%npaire==1) then
            buffer(1) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                      + bcix6*(ux(3,j,k)-ux(2,j,k))
            buffer(2) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                      + bcix6*(ux(4,j,k)-ux(1,j,k))
         else
            buffer(1) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                      + bcix6*(ux(3,j,k)-two*ux(1,j,k)+ux(2,j,k))
            buffer(2) = acix6*(ux(3,j,k)-ux(2,j,k)) &
                      + bcix6*(ux(4,j,k)-ux(1,j,k))
         endif
         do concurrent (i=3:nxm-2)
            buffer(i) = acix6*(ux(i+1,j,k)-ux(i  ,j,k)) &
                      + bcix6*(ux(i+2,j,k)-ux(i-1,j,k))
         enddo
         if (x3dop%npaire==1) then
            buffer(nxm-1) = acix6*(ux(nxm,j,k)-ux(nxm-1,j,k)) &
                          + bcix6*(ux(nx ,j,k)-ux(nxm-2,j,k))
            buffer(nxm) = acix6*(ux(nx ,j,k)-ux(nxm  ,j,k)) &
                         + bcix6*(ux(nx ,j,k)-ux(nxm-2,j,k))
            buffer(nxm) = acix6*(ux(nx,j,k)-ux(nxm,j,k)) &
                        + bcix6*(two*ux(nx,j,k)-ux(nxm,j,k)-ux(nxm-1,j,k))
         endif

         ! Solve tri-diagonal system
         call thomas1d(buffer, ff, ss, ww, nxm)
         do concurrent (i=1:nxm)
            tx(i,j,k) = buffer(i)
         enddo
      enddo
   endif

end subroutine derxvp

!********************************************************************
!
subroutine interxvp(tx,ux,x3dop,nx,nxm,ny,nz)

  use x3d_operator_x_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nxm,ny,nz) :: tx
  real(mytype), intent(in), dimension(nx,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nxm) :: buffer, ff, ss, ww, pp

  do concurrent (i=1:nxm)
      ff(i) = x3dop%f(i)
      ss(i) = x3dop%s(i)
      ww(i) = x3dop%w(i)
      if (allocated(x3dop%periodic)) pp(i) = x3dop%periodic(i)
   end do

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aicix6*(ux(2,j,k)+ux(1  ,j,k)) &
                  + bicix6*(ux(3,j,k)+ux(nx,j,k)) &
                  + cicix6*(ux(4,j,k)+ux(nx-1,j,k)) &
                  + dicix6*(ux(5,j,k)+ux(nx-2,j,k))
        buffer(2) = aicix6*(ux(3,j,k)+ux(2 ,j,k)) &
                  + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                  + cicix6*(ux(5,j,k)+ux(nx,j,k)) &
                  + dicix6*(ux(6,j,k)+ux(nx-1,j,k))
        buffer(3) = aicix6*(ux(4,j,k)+ux(3 ,j,k)) &
                  + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                  + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                  + dicix6*(ux(7,j,k)+ux(nx,j,k))
        do concurrent (i=4:nx-4)
           buffer(i) = aicix6*(ux(i+1,j,k)+ux(i,j,k)) &
                     + bicix6*(ux(i+2,j,k)+ux(i-1,j,k)) &
                     + cicix6*(ux(i+3,j,k)+ux(i-2,j,k)) &
                     + dicix6*(ux(i+4,j,k)+ux(i-3,j,k))
        enddo
        buffer(nx-3) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                     + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                     + cicix6*(ux(nx,j,k)+ux(nx-5,j,k)) &
                     + dicix6*(ux(1,j,k)+ux(nx-6,j,k))
        buffer(nx-2) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                     + bicix6*(ux(nx ,j,k)+ux(nx-3,j,k)) &
                     + cicix6*(ux(1,j,k)+ux(nx-4,j,k)) &
                     + dicix6*(ux(2,j,k)+ux(nx-5,j,k))
        buffer(nx-1) = aicix6*(ux(nx,j,k)+ux(nx-1,j,k)) &
                     + bicix6*(ux(1 ,j,k)+ux(nx-2,j,k)) &
                     + cicix6*(ux(2,j,k)+ux(nx-3,j,k)) &
                     + dicix6*(ux(3,j,k)+ux(nx-4,j,k))
        buffer(nx  ) = aicix6*(ux(1,j,k)+ux(nx,j,k)) &
                     + bicix6*(ux(2,j,k)+ux(nx-1,j,k)) &
                     + cicix6*(ux(3,j,k)+ux(nx-2,j,k)) &
                     + dicix6*(ux(4,j,k)+ux(nx-3,j,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nx)
        do concurrent (i=1:nx)
           tx(i,j,k) = buffer(i)
        enddo
     enddo

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aicix6*(ux(2,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + dicix6*(ux(5,j,k)+ux(4,j,k))
           buffer(2) = aicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + dicix6*(ux(6,j,k)+ux(3,j,k))
           buffer(3) = aicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(7,j,k)+ux(2,j,k))
           do concurrent (i=4:nxm-3)
              buffer(i) = aicix6*(ux(i+1,j,k)+ux(i,j,k)) &
                        + bicix6*(ux(i+2,j,k)+ux(i-1,j,k)) &
                        + cicix6*(ux(i+3,j,k)+ux(i-2,j,k)) &
                        + dicix6*(ux(i+4,j,k)+ux(i-3,j,k))
           enddo
           buffer(nxm-2) = aicix6*(ux(nxm-1,j,k)+ux(nxm-2,j,k)) &
                         + bicix6*(ux(nxm,j,k)+ux(nxm-3,j,k)) &
                         + cicix6*(ux(nx,j,k)+ux(nxm-4,j,k)) &
                         + dicix6*(ux(nxm,j,k)+ux(nxm-5,j,k))
           buffer(nxm-1) = aicix6*(ux(nxm,j,k)+ux(nxm-1,j,k)) &
                         + bicix6*(ux(nx,j,k)+ux(nxm-2,j,k)) &
                         + cicix6*(ux(nxm,j,k)+ux(nxm-3,j,k)) &
                         + dicix6*(ux(nxm-1,j,k)+ux(nxm-4,j,k))
           buffer(nxm  ) = aicix6*(ux(nx,j,k)+ux(nxm,j,k)) &
                         + bicix6*(ux(nxm,j,k)+ux(nxm-1,j,k)) &
                         + cicix6*(ux(nxm-1,j,k)+ux(nxm-2,j,k)) &
                         + dicix6*(ux(nxm-2,j,k)+ux(nxm-3,j,k))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nxm)
           do concurrent (i=1:nxm)
              tx(i,j,k) = buffer(i)
           enddo
        enddo

     endif
  endif

end subroutine interxvp

!********************************************************************
!
subroutine derxpv(tx,ux,x3dop,nxm,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nxm,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nx) :: buffer, ff, ss, ww, pp

  do concurrent (i=1:nx)
     ff(i) = x3dop%f(i)
     ss(i) = x3dop%s(i)
     ww(i) = x3dop%w(i)
     if (allocated(x3dop%periodic)) pp(i) = x3dop%periodic(i)
  end do

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny) local(buffer)

        ! Compute r.h.s.
        buffer(1) = acix6*(ux(1,j,k)-ux(nx  ,j,k)) &
                  + bcix6*(ux(2,j,k)-ux(nx-1,j,k))
        buffer(2) = acix6*(ux(2,j,k)-ux(1 ,j,k)) &
                  + bcix6*(ux(3,j,k)-ux(nx,j,k))
        do concurrent (i=3:nx-2)
           buffer(i) = acix6*(ux(i,j,k)-ux(i-1,j,k)) &
                     + bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
        enddo
        buffer(nx-1) = acix6*(ux(nx-1,j,k)-ux(nx-2,j,k)) &
                     + bcix6*(ux(nx ,j,k)-ux(nx-3,j,k))
        buffer(nx  ) = acix6*(ux(nx,j,k)-ux(nx-1,j,k)) &
                     + bcix6*(ux(1,j,k)-ux(nx-2,j,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nx)
        do concurrent (i=1:nx)
           tx(i,j,k) = buffer(i)
        enddo
     enddo

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny) local(buffer)

           ! Compute r.h.s.
           buffer(1) = zero
           buffer(2) = acix6*(ux(2,j,k)-ux(1,j,k)) &
                     + bcix6*(ux(3,j,k)-ux(1,j,k))
           do concurrent (i=3:nx-2)
              buffer(i) = acix6*(ux(i,j,k)-ux(i-1,j,k)) &
                        + bcix6*(ux(i+1,j,k)-ux(i-2,j,k))
           enddo
           buffer(nx-1) = acix6*(ux(nx-1,j,k)-ux(nx-2,j,k)) &
                        + bcix6*(ux(nx-1,j,k)-ux(nx-3,j,k))
           buffer(nx) = zero

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nx)
           do concurrent (i=1:nx)
              tx(i,j,k) = buffer(i)
           enddo
        enddo

     endif
  endif

end subroutine derxpv

!********************************************************************
!
subroutine interxpv(tx,ux,x3dop,nxm,nx,ny,nz)

  use x3d_operator_x_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, nxm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tx
  real(mytype), intent(in), dimension(nxm,ny,nz) :: ux
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nx) :: buffer, ff, ss, ww, pp

  do concurrent (i=1:nx)
      ff(i) = x3dop%f(i)
      ss(i) = x3dop%s(i)
      ww(i) = x3dop%w(i)
      if (allocated(x3dop%periodic)) pp(i) = x3dop%periodic(i)
   end do

  if (nclx) then
     ! nxm = nx
     do concurrent (k=1:nz, j=1:ny) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aicix6*(ux(1,j,k)+ux(nx  ,j,k)) &
                  + bicix6*(ux(2,j,k)+ux(nx-1,j,k)) &
                  + cicix6*(ux(3,j,k)+ux(nx-2,j,k)) &
                  + dicix6*(ux(4,j,k)+ux(nx-3,j,k))
        buffer(2) = aicix6*(ux(2,j,k)+ux(1 ,j,k)) &
                  + bicix6*(ux(3,j,k)+ux(nx,j,k)) &
                  + cicix6*(ux(4,j,k)+ux(nx-1,j,k)) &
                  + dicix6*(ux(5,j,k)+ux(nx-2,j,k))
        buffer(3) = aicix6*(ux(3,j,k)+ux(2 ,j,k)) &
                  + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                  + cicix6*(ux(5,j,k)+ux(nx,j,k)) &
                  + dicix6*(ux(6,j,k)+ux(nx-1,j,k))
        buffer(4) = aicix6*(ux(4,j,k)+ux(3 ,j,k)) &
                  + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                  + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                  + dicix6*(ux(7,j,k)+ux(nx,j,k))
        do concurrent (i=5:nx-3)
           buffer(i) = aicix6*(ux(i,j,k)+ux(i-1,j,k)) &
                     + bicix6*(ux(i+1,j,k)+ux(i-2,j,k)) &
                     + cicix6*(ux(i+2,j,k)+ux(i-3,j,k)) &
                     + dicix6*(ux(i+3,j,k)+ux(i-4,j,k))
        enddo
        buffer(nx-2) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                     + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                     + cicix6*(ux(nx,j,k)+ux(nx-5,j,k)) &
                     + dicix6*(ux(1,j,k)+ux(nx-6,j,k))
        buffer(nx-1) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                     + bicix6*(ux(nx ,j,k)+ux(nx-3,j,k)) &
                     + cicix6*(ux(1,j,k)+ux(nx-4,j,k)) &
                     + dicix6*(ux(2,j,k)+ux(nx-5,j,k))
        buffer(nx  ) = aicix6*(ux(nx,j,k)+ux(nx-1,j,k)) &
                     + bicix6*(ux(1,j,k)+ux(nx-2,j,k)) &
                     + cicix6*(ux(2,j,k)+ux(nx-3,j,k)) &
                     + dicix6*(ux(3,j,k)+ux(nx-4,j,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nx)
        do concurrent (i=1:nx)
           tx(i,j,k) = buffer(i)
        enddo
     enddo

  else
     ! nxm = nx-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, j=1:ny) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aicix6*(ux(1,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(2,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(3,j,k)+ux(3,j,k)) &
                     + dicix6*(ux(4,j,k)+ux(4,j,k))
           buffer(2) = aicix6*(ux(2,j,k)+ux(1,j,k)) &
                     + bicix6*(ux(3,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(4,j,k)+ux(2,j,k)) &
                     + dicix6*(ux(5,j,k)+ux(3,j,k))
           buffer(3) = aicix6*(ux(3,j,k)+ux(2,j,k)) &
                     + bicix6*(ux(4,j,k)+ux(1,j,k)) &
                     + cicix6*(ux(5,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(6,j,k)+ux(2,j,k))
           buffer(4) = aicix6*(ux(4,j,k)+ux(3,j,k)) &
                     + bicix6*(ux(5,j,k)+ux(2,j,k)) &
                     + cicix6*(ux(6,j,k)+ux(1,j,k)) &
                     + dicix6*(ux(7,j,k)+ux(1,j,k))
           do concurrent (i=5:nx-4)
              buffer(i) = aicix6*(ux(i,j,k)+ux(i-1,j,k)) &
                        + bicix6*(ux(i+1,j,k)+ux(i-2,j,k)) &
                        + cicix6*(ux(i+2,j,k)+ux(i-3,j,k)) &
                        + dicix6*(ux(i+3,j,k)+ux(i-4,j,k))
           enddo
           buffer(nx-3) = aicix6*(ux(nx-3,j,k)+ux(nx-4,j,k)) &
                        + bicix6*(ux(nx-2,j,k)+ux(nx-5,j,k)) &
                        + cicix6*(ux(nx-1,j,k)+ux(nx-6,j,k)) &
                        + dicix6*(ux(nx-1,j,k)+ux(nx-7,j,k))
           buffer(nx-2) = aicix6*(ux(nx-2,j,k)+ux(nx-3,j,k)) &
                        + bicix6*(ux(nx-1,j,k)+ux(nx-4,j,k)) &
                        + cicix6*(ux(nx-1,j,k)+ux(nx-5,j,k)) &
                        + dicix6*(ux(nx-2,j,k)+ux(nx-6,j,k))
           buffer(nx-1) = aicix6*(ux(nx-1,j,k)+ux(nx-2,j,k)) &
                        + bicix6*(ux(nx-1,j,k)+ux(nx-3,j,k)) &
                        + cicix6*(ux(nx-2,j,k)+ux(nx-4,j,k)) &
                        + dicix6*(ux(nx-3,j,k)+ux(nx-5,j,k))
           buffer(nx  ) = aicix6*(ux(nx-1,j,k)+ux(nx-1,j,k)) &
                        + bicix6*(ux(nx-2,j,k)+ux(nx-2,j,k)) &
                        + cicix6*(ux(nx-3,j,k)+ux(nx-3,j,k)) &
                        + dicix6*(ux(nx-4,j,k)+ux(nx-4,j,k))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nx)
           do concurrent (i=1:nx)
              tx(i,j,k) = buffer(i)
           enddo
        enddo

     endif
  endif

end subroutine interxpv

!********************************************************************
!
subroutine interyvp(ty,uy,x3dop,nx,ny,nym,nz)

  USE x3d_operator_y_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,nym,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nym) :: buffer, ff, ss, ww, pp

  do concurrent (j=1:nym)
      ff(j) = x3dop%f(j)
      ss(j) = x3dop%s(j)
      ww(j) = x3dop%w(j)
      if (allocated(x3dop%periodic)) pp(j) = x3dop%periodic(j)
   end do

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                  + biciy6*(uy(i,3,k)+uy(i,ny,k)) &
                  + ciciy6*(uy(i,4,k)+uy(i,ny-1,k)) &
                  + diciy6*(uy(i,5,k)+uy(i,ny-2,k))
        buffer(2) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                  + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                  + ciciy6*(uy(i,5,k)+uy(i,ny,k)) &
                  + diciy6*(uy(i,6,k)+uy(i,ny-1,k))
        buffer(3) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                  + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                  + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                  + diciy6*(uy(i,7,k)+uy(i,ny,k))
        do concurrent (j=4:ny-4)
           buffer(j) = aiciy6*(uy(i,j+1,k)+uy(i,j,k)) &
                     + biciy6*(uy(i,j+2,k)+uy(i,j-1,k)) &
                     + ciciy6*(uy(i,j+3,k)+uy(i,j-2,k)) &
                     + diciy6*(uy(i,j+4,k)+uy(i,j-3,k))
        enddo
        buffer(ny-3) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                     + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                     + ciciy6*(uy(i,ny,k)+uy(i,ny-5,k)) &
                     + diciy6*(uy(i,1,k)+uy(i,ny-6,k))
        buffer(ny-2) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                     + biciy6*(uy(i,ny,k)+uy(i,ny-3,k)) &
                     + ciciy6*(uy(i,1,k)+uy(i,ny-4,k)) &
                     + diciy6*(uy(i,2,k)+uy(i,ny-5,k))
        buffer(ny-1) = aiciy6*(uy(i,ny,k)+uy(i,ny-1,k)) &
                     + biciy6*(uy(i,1,k)+uy(i,ny-2,k)) &
                     + ciciy6*(uy(i,2,k)+uy(i,ny-3,k)) &
                     + diciy6*(uy(i,3,k)+uy(i,ny-4,k))
        buffer(ny  ) = aiciy6*(uy(i,1,k)+uy(i,ny,k)) &
                     + biciy6*(uy(i,2,k)+uy(i,ny-1,k)) &
                     + ciciy6*(uy(i,3,k)+uy(i,ny-2,k)) &
                     + diciy6*(uy(i,4,k)+uy(i,ny-3,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, ny)
        do concurrent (j=1:ny)
           ty(i,j,k) = buffer(j)
        enddo
     enddo

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                     + biciy6*(uy(i,3,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,4,k)+uy(i,3,k)) &
                     + diciy6*(uy(i,5,k)+uy(i,4,k))
           buffer(2) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                     + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                     + ciciy6*(uy(i,5,k)+uy(i,2,k)) &
                     + diciy6*(uy(i,6,k)+uy(i,3,k))
           buffer(3) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                     + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                     + diciy6*(uy(i,7,k)+uy(i,2,k))
           do concurrent (j=4:nym-3)
              buffer(j) = aiciy6*(uy(i,j+1,k)+uy(i,j,k)) &
                        + biciy6*(uy(i,j+2,k)+uy(i,j-1,k)) &
                        + ciciy6*(uy(i,j+3,k)+uy(i,j-2,k)) &
                        + diciy6*(uy(i,j+4,k)+uy(i,j-3,k))
           enddo
           buffer(nym-2) = aiciy6*(uy(i,nym-1,k)+uy(i,nym-2,k)) &
                         + biciy6*(uy(i,nym,k)+uy(i,nym-3,k)) &
                         + ciciy6*(uy(i,ny,k)+uy(i,nym-4,k)) &
                         + diciy6*(uy(i,nym,k)+uy(i,nym-5,k))
           buffer(nym-1) = aiciy6*(uy(i,nym,k)+uy(i,nym-1,k)) &
                         + biciy6*(uy(i,ny,k)+uy(i,nym-2,k)) &
                         + ciciy6*(uy(i,nym,k)+uy(i,nym-3,k)) &
                         + diciy6*(uy(i,nym-1,k)+uy(i,nym-4,k))
           buffer(nym  ) = aiciy6*(uy(i,ny,k)+uy(i,nym,k)) &
                         + biciy6*(uy(i,nym,k)+uy(i,nym-1,k)) &
                         + ciciy6*(uy(i,nym-1,k)+uy(i,nym-2,k)) &
                         + diciy6*(uy(i,nym-2,k)+uy(i,nym-3,k))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nym)
           do concurrent (j=1:nym)
              ty(i,j,k) = buffer(j)
           enddo
        enddo

     endif
  endif

end subroutine interyvp

!********************************************************************
!
subroutine deryvp(ty,uy,x3dop,ppyi,nx,ny,nym,nz)

  USE x3d_operator_y_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,nym,nz) :: ty
  real(mytype), intent(in), dimension(nx,ny,nz) :: uy
  real(mytype), intent(in), dimension(nym) :: ppyi
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nym) :: buffer, ff, ss, ww, pp

  do concurrent (j=1:nym)
     ff(j) = x3dop%f(j)
     ss(j) = x3dop%s(j)
     ww(j) = x3dop%w(j)
     if (allocated(x3dop%periodic)) pp(j) = x3dop%periodic(j)
  end do

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                  + bciy6*(uy(i,3,k)-uy(i,ny,k))
        buffer(2) = aciy6*(uy(i,3,k)-uy(i,2,k)) &
                  + bciy6*(uy(i,4,k)-uy(i,1,k))
        do concurrent (j=3:ny-2)
           buffer(j) = aciy6*(uy(i,j+1,k)-uy(i,j,k)) &
                     + bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
        enddo
        buffer(ny-1) = aciy6*(uy(i,ny,k)-uy(i,ny-1,k)) &
                     + bciy6*(uy(i,1,k)-uy(i,ny-2,k))
        buffer(ny  ) = aciy6*(uy(i,1,k)-uy(i,ny,k)) &
                     + bciy6*(uy(i,2,k)-uy(i,ny-1,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, ny)
        if (istret /= 0) then
           do concurrent (j=1:nym)
              buffer(j) = buffer(j) * ppyi(j)
           enddo
        endif
        do concurrent (j=1:ny)
           ty(i,j,k) = buffer(j)
        enddo
     enddo

  else
     ! nym = ny-1
     if (x3dop%npaire==0) then
        do concurrent (k=1:nz, i=1:nx) local(buffer)

           ! Compute r.h.s.
              buffer(1) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                        + bciy6*(uy(i,3,k)-two*uy(i,1,k)+uy(i,2,k))
              buffer(2) = aciy6*(uy(i,3,k)-uy(i,2,k)) &
                        + bciy6*(uy(i,4,k)-uy(i,1,k))
           do concurrent (j=3:nym-2)
              buffer(j) = aciy6*(uy(i,j+1,k)-uy(i,j,k)) &
                        + bciy6*(uy(i,j+2,k)-uy(i,j-1,k))
           enddo
              buffer(nym-1) = aciy6*(uy(i,nym,k)-uy(i,nym-1,k)) &
                            + bciy6*(uy(i,ny,k)-uy(i,nym-2,k))
              buffer(nym  ) = aciy6*(uy(i,ny,k)-uy(i,nym,k)) &
                            + bciy6*(two*uy(i,ny,k)-uy(i,nym,k)-uy(i,nym-1,k))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nym)
           if (istret /= 0) then
              do concurrent (j=1:nym)
                 buffer(j) = buffer(j) * ppyi(j)
              enddo
           endif
           do concurrent (j=1:nym)
              ty(i,j,k) = buffer(j)
           enddo
        enddo

     endif
  endif

end subroutine deryvp

!********************************************************************
!
subroutine interypv(ty,uy,x3dop,nx,nym,ny,nz)

  USE x3d_operator_y_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,nym,nz) :: uy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(ny) :: buffer, ff, ss, ww, pp

  do concurrent (j=1:ny)
     ff(j) = x3dop%f(j)
     ss(j) = x3dop%s(j)
     ww(j) = x3dop%w(j)
     if (allocated(x3dop%periodic)) pp(j) = x3dop%periodic(j)
  end do

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aiciy6*(uy(i,1,k)+uy(i,ny,k)) &
                  + biciy6*(uy(i,2,k)+uy(i,ny-1,k)) &
                  + ciciy6*(uy(i,3,k)+uy(i,ny-2,k)) &
                  + diciy6*(uy(i,4,k)+uy(i,ny-3,k))
        buffer(2) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                  + biciy6*(uy(i,3,k)+uy(i,ny,k)) &
                  + ciciy6*(uy(i,4,k)+uy(i,ny-1,k)) &
                  + diciy6*(uy(i,5,k)+uy(i,ny-2,k))
        buffer(3) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                  + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                  + ciciy6*(uy(i,5,k)+uy(i,ny,k)) &
                  + diciy6*(uy(i,6,k)+uy(i,ny-1,k))
        buffer(4) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                  + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                  + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                  + diciy6*(uy(i,7,k)+uy(i,ny,k))
        do concurrent (j=5:ny-3)
           buffer(j) = aiciy6*(uy(i,j,k)+uy(i,j-1,k)) &
                     + biciy6*(uy(i,j+1,k)+uy(i,j-2,k)) &
                     + ciciy6*(uy(i,j+2,k)+uy(i,j-3,k)) &
                     + diciy6*(uy(i,j+3,k)+uy(i,j-4,k))
        enddo
        buffer(ny-2) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                     + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                     + ciciy6*(uy(i,ny,k)+uy(i,ny-5,k)) &
                     + diciy6*(uy(i,1,k)+uy(i,ny-6,k))
        buffer(ny-1) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                     + biciy6*(uy(i,ny,k)+uy(i,ny-3,k)) &
                     + ciciy6*(uy(i,1,k)+uy(i,ny-4,k)) &
                     + diciy6*(uy(i,2,k)+uy(i,ny-5,k))
        buffer(ny  ) = aiciy6*(uy(i,ny,k)+uy(i,ny-1,k)) &
                     + biciy6*(uy(i,1,k)+uy(i,ny-2,k)) &
                     + ciciy6*(uy(i,2,k)+uy(i,ny-3,k)) &
                     + diciy6*(uy(i,3,k)+uy(i,ny-4,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, ny)
        do concurrent (j=1:ny)
           ty(i,j,k) = buffer(j)
        enddo
     enddo

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aiciy6*(uy(i,1,k)+uy(i,1,k)) &
                     + biciy6*(uy(i,2,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,3,k)+uy(i,3,k)) &
                     + diciy6*(uy(i,4,k)+uy(i,4,k))
           buffer(2) = aiciy6*(uy(i,2,k)+uy(i,1,k)) &
                     + biciy6*(uy(i,3,k)+uy(i,1,k)) &
                     + ciciy6*(uy(i,4,k)+uy(i,2,k)) &
                     + diciy6*(uy(i,5,k)+uy(i,3,k))
           buffer(3) = aiciy6*(uy(i,3,k)+uy(i,2,k)) &
                     + biciy6*(uy(i,4,k)+uy(i,1,k)) &
                     + ciciy6*(uy(i,5,k)+uy(i,1,k)) &
                     + diciy6*(uy(i,6,k)+uy(i,2,k))
           buffer(4) = aiciy6*(uy(i,4,k)+uy(i,3,k)) &
                     + biciy6*(uy(i,5,k)+uy(i,2,k)) &
                     + ciciy6*(uy(i,6,k)+uy(i,1,k)) &
                     + diciy6*(uy(i,7,k)+uy(i,1,k))
           do concurrent (j=5:ny-4)
              buffer(j) = aiciy6*(uy(i,j,k)+uy(i,j-1,k)) &
                        + biciy6*(uy(i,j+1,k)+uy(i,j-2,k)) &
                        + ciciy6*(uy(i,j+2,k)+uy(i,j-3,k)) &
                        + diciy6*(uy(i,j+3,k)+uy(i,j-4,k))
           enddo
           buffer(ny-3) = aiciy6*(uy(i,ny-3,k)+uy(i,ny-4,k)) &
                        + biciy6*(uy(i,ny-2,k)+uy(i,ny-5,k)) &
                        + ciciy6*(uy(i,ny-1,k)+uy(i,ny-6,k)) &
                        + diciy6*(uy(i,ny-1,k)+uy(i,ny-7,k))
           buffer(ny-2) = aiciy6*(uy(i,ny-2,k)+uy(i,ny-3,k)) &
                        + biciy6*(uy(i,ny-1,k)+uy(i,ny-4,k)) &
                        + ciciy6*(uy(i,ny-1,k)+uy(i,ny-5,k)) &
                        + diciy6*(uy(i,ny-2,k)+uy(i,ny-6,k))
           buffer(ny-1) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-2,k)) &
                        + biciy6*(uy(i,ny-1,k)+uy(i,ny-3,k)) &
                        + ciciy6*(uy(i,ny-2,k)+uy(i,ny-4,k)) &
                        + diciy6*(uy(i,ny-3,k)+uy(i,ny-5,k))
           buffer(ny  ) = aiciy6*(uy(i,ny-1,k)+uy(i,ny-1,k)) &
                        + biciy6*(uy(i,ny-2,k)+uy(i,ny-2,k)) &
                        + ciciy6*(uy(i,ny-3,k)+uy(i,ny-3,k)) &
                        + diciy6*(uy(i,ny-4,k)+uy(i,ny-4,k))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, ny)
           do concurrent (j=1:ny)
              ty(i,j,k) = buffer(j)
           enddo
        enddo

     endif
  endif

end subroutine interypv

!********************************************************************
!
subroutine derypv(ty,uy,x3dop,ppy,nx,nym,ny,nz)

  USE x3d_operator_y_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nym, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: ty
  real(mytype), intent(in), dimension(nx,nym,nz) :: uy
  real(mytype), intent(in), dimension(ny) :: ppy
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(ny) :: buffer, ff, ss, ww, pp

  do concurrent (j=1:ny)
     ff(j) = x3dop%f(j)
     ss(j) = x3dop%s(j)
     ww(j) = x3dop%w(j)
     if (allocated(x3dop%periodic)) pp(j) = x3dop%periodic(j)
  end do

  if (ncly) then
     ! nym = ny
     do concurrent (k=1:nz, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aciy6*(uy(i,1,k)-uy(i,ny,k)) &
                  + bciy6*(uy(i,2,k)-uy(i,ny-1,k))
        buffer(2) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                  + bciy6*(uy(i,3,k)-uy(i,ny,k))
        do concurrent (j=3:ny-2)
           buffer(j) = aciy6*(uy(i,j,k)-uy(i,j-1,k)) &
                     + bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
        enddo
        buffer(ny-1) = aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k)) &
                     + bciy6*(uy(i,ny,k)-uy(i,ny-3,k))
        buffer(ny) = aciy6*(uy(i,ny,k)-uy(i,ny-1,k)) &
                   + bciy6*(uy(i,1,k)-uy(i,ny-2,k))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, ny)
        if (istret /= 0) then
           do concurrent (j=1:ny)
              buffer(j) = buffer(j) * ppy(j)
           enddo
        endif
        do concurrent (j=1:ny)
           ty(i,j,k) = buffer(j)
        enddo
     enddo

  else
     ! nym = ny-1
     if (x3dop%npaire==1) then
        do concurrent (k=1:nz, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = zero
           buffer(2) = aciy6*(uy(i,2,k)-uy(i,1,k)) &
                     + bciy6*(uy(i,3,k)-uy(i,1,k))
           do concurrent (j=3:ny-2)
              buffer(j) = aciy6*(uy(i,j,k)-uy(i,j-1,k)) &
                        + bciy6*(uy(i,j+1,k)-uy(i,j-2,k))
           enddo
           buffer(ny-1) = aciy6*(uy(i,ny-1,k)-uy(i,ny-2,k)) &
                        + bciy6*(uy(i,ny-1,k)-uy(i,ny-3,k))
           buffer(ny) = zero

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, ny)
           if (istret /= 0) then
              do concurrent (j=1:ny)
                 buffer(j) = buffer(j) * ppy(j)
              enddo
           endif
           do concurrent (j=1:ny)
              ty(i,j,k) = buffer(j)
           enddo
        enddo

     endif
  endif

end subroutine derypv

!********************************************************************
!
subroutine derzvp(tz,uz,x3dop,nx,ny,nz,nzm)

  USE x3d_operator_z_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nzm) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nzm) :: buffer, ff, ss, ww, pp

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = zero
     enddo
     return
  endif

  do concurrent (k=1:nzm)
     ff(k) = x3dop%f(k)
     ss(k) = x3dop%s(k)
     ww(k) = x3dop%w(k)
     if (allocated(x3dop%periodic)) pp(k) = x3dop%periodic(k)
  end do

  if (nclz) then
     ! nzm = nz
     do concurrent (j=1:ny, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                  + bciz6*(uz(i,j,3)-uz(i,j,nz))
        buffer(2) = aciz6*(uz(i,j,3)-uz(i,j,2)) &
                  + bciz6*(uz(i,j,4)-uz(i,j,1))
        do concurrent (k=3:nz-2)
           buffer(k) = aciz6*(uz(i,j,k+1)-uz(i,j,k)) &
                     + bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
        enddo
        buffer(nz-1) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                     + bciz6*(uz(i,j,1)-uz(i,j,nz-2))
        buffer(nz  ) = aciz6*(uz(i,j,1)-uz(i,j,nz)) &
                     + bciz6*(uz(i,j,2)-uz(i,j,nz-1))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nz)
        do concurrent (k=1:nz)
           tz(i,j,k) = buffer(k)
        enddo
     enddo

  else
     ! nzm = nz-1
     do concurrent (j=1:ny, i=1:nx) local(buffer)

        ! Compute r.h.s.
        if (x3dop%npaire==1) then
           buffer(1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-uz(i,j,2))
           buffer(2) = aciz6*(uz(i,j,3)-uz(i,j,2))&
                     + bciz6*(uz(i,j,4)-uz(i,j,1))
        else
           buffer(1) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-two*uz(i,j,1)+uz(i,j,2))
           buffer(2) = aciz6*(uz(i,j,3)-uz(i,j,2)) &
                     + bciz6*(uz(i,j,4)-uz(i,j,1))
        endif
        do concurrent (k=3:nzm-2)
           buffer(k) = aciz6*(uz(i,j,k+1)-uz(i,j,k)) &
                     + bciz6*(uz(i,j,k+2)-uz(i,j,k-1))
        enddo
        if (x3dop%npaire==1) then
           buffer(nzm-1) = aciz6*(uz(i,j,nzm)-uz(i,j,nzm-1)) &
                         + bciz6*(uz(i,j,nz)-uz(i,j,nzm-2))
           buffer(nzm  ) = aciz6*(uz(i,j,nz)-uz(i,j,nzm)) &
                         + bciz6*(uz(i,j,nzm)-uz(i,j,nzm-1))
        else
           buffer(nzm-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                         + bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
           buffer(nzm  ) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                         + bciz6*(two*uz(i,j,nz)-uz(i,j,nz-1)-uz(i,j,nz-2))
        endif

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, nzm)
        do concurrent (k=1:nzm)
           tz(i,j,k) = buffer(k)
        enddo
     enddo

  endif

end subroutine derzvp

!********************************************************************
!
subroutine interzvp(tz,uz,x3dop,nx,ny,nz,nzm)

  USE x3d_operator_z_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nzm) :: tz
  real(mytype), intent(in), dimension(nx,ny,nz) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nzm) :: buffer, ff, ss, ww, pp

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = uz(i,j,k)
     enddo
     return
  endif

  do concurrent (k=1:nzm)
     ff(k) = x3dop%f(k)
     ss(k) = x3dop%s(k)
     ww(k) = x3dop%w(k)
     if (allocated(x3dop%periodic)) pp(k) = x3dop%periodic(k)
  end do

  if (nclz) then
     ! nzm = nz
     do concurrent (j=1:ny, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                  + biciz6*(uz(i,j,3)+uz(i,j,nz)) &
                  + ciciz6*(uz(i,j,4)+uz(i,j,nz-1)) &
                  + diciz6*(uz(i,j,5)+uz(i,j,nz-2))
        buffer(2) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                  + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                  + ciciz6*(uz(i,j,5)+uz(i,j,nz)) &
                  + diciz6*(uz(i,j,6)+uz(i,j,nz-1))
        buffer(3) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                  + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                  + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                  + diciz6*(uz(i,j,7)+uz(i,j,nz))
        do concurrent (k=4:nz-4)
           buffer(k) = aiciz6*(uz(i,j,k+1)+uz(i,j,k)) &
                     + biciz6*(uz(i,j,k+2)+uz(i,j,k-1)) &
                     + ciciz6*(uz(i,j,k+3)+uz(i,j,k-2)) &
                     + diciz6*(uz(i,j,k+4)+uz(i,j,k-3))
        enddo
        buffer(nz-3) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                     + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                     + ciciz6*(uz(i,j,nz)+uz(i,j,nz-5)) &
                     + diciz6*(uz(i,j,1)+uz(i,j,nz-6))
        buffer(nz-2) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                     + biciz6*(uz(i,j,nz)+uz(i,j,nz-3)) &
                     + ciciz6*(uz(i,j,1)+uz(i,j,nz-4)) &
                     + diciz6*(uz(i,j,2)+uz(i,j,nz-5))
        buffer(nz-1) = aiciz6*(uz(i,j,nz)+uz(i,j,nz-1)) &
                     + biciz6*(uz(i,j,1)+uz(i,j,nz-2)) &
                     + ciciz6*(uz(i,j,2)+uz(i,j,nz-3)) &
                     + diciz6*(uz(i,j,3)+uz(i,j,nz-4))
        buffer(nz  ) = aiciz6*(uz(i,j,1)+uz(i,j,nz)) &
                     + biciz6*(uz(i,j,2)+uz(i,j,nz-1)) &
                     + ciciz6*(uz(i,j,3)+uz(i,j,nz-2)) &
                     + diciz6*(uz(i,j,4)+uz(i,j,nz-3))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nz)
        do concurrent (k=1:nz)
           tz(i,j,k) = buffer(k)
        enddo
     enddo

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then
        do concurrent (j=1:ny, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + diciz6*(uz(i,j,5)+uz(i,j,4))
           buffer(2) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                     + ciciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + diciz6*(uz(i,j,6)+uz(i,j,3))
           buffer(3) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,7)+uz(i,j,2))
           do concurrent (k=4:nzm-3)
              buffer(k) = aiciz6*(uz(i,j,k+1)+uz(i,j,k)) &
                        + biciz6*(uz(i,j,k+2)+uz(i,j,k-1)) &
                        + ciciz6*(uz(i,j,k+3)+uz(i,j,k-2)) &
                        + diciz6*(uz(i,j,k+4)+uz(i,j,k-3))
           enddo
           buffer(nzm-2) = aiciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-2)) &
                         + biciz6*(uz(i,j,nzm)+uz(i,j,nzm-3)) &
                         + ciciz6*(uz(i,j,nz)+uz(i,j,nzm-4)) &
                         + diciz6*(uz(i,j,nzm)+uz(i,j,nzm-5))
           buffer(nzm-1) = aiciz6*(uz(i,j,nzm)+uz(i,j,nzm-1)) &
                         + biciz6*(uz(i,j,nz)+uz(i,j,nzm-2)) &
                         + ciciz6*(uz(i,j,nzm)+uz(i,j,nzm-3)) &
                         + diciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-4))
           buffer(nzm) = aiciz6*(uz(i,j,nz)+uz(i,j,nzm)) &
                       + biciz6*(uz(i,j,nzm)+uz(i,j,nzm-1)) &
                       + ciciz6*(uz(i,j,nzm-1)+uz(i,j,nzm-2)) &
                       + diciz6*(uz(i,j,nzm-2)+uz(i,j,nzm-3))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nzm)
           do concurrent (k=1:nzm)
              tz(i,j,k) = buffer(k)
           enddo
        enddo

     endif
  endif

end subroutine interzvp

!********************************************************************
!
subroutine derzpv(tz,uz,x3dop,nx,ny,nzm,nz)

  USE x3d_operator_z_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, nzm, ny, nz
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nzm) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nz) :: buffer, ff, ss, ww, pp

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = zero
     enddo
     return
  endif

  do concurrent (k=1:nz)
     ff(k) = x3dop%f(k)
     ss(k) = x3dop%s(k)
     ww(k) = x3dop%w(k)
     if (allocated(x3dop%periodic)) pp(k) = x3dop%periodic(k)
  end do

  if (nclz) then
     ! nzm = nz
     do concurrent (j=1:ny, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aciz6*(uz(i,j,1)-uz(i,j,nz)) &
                  + bciz6*(uz(i,j,2)-uz(i,j,nz-1))
        buffer(2) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                  + bciz6*(uz(i,j,3)-uz(i,j,nz))
        do concurrent (k=3:nz-2)
           buffer(k) = aciz6*(uz(i,j,k)-uz(i,j,k-1)) &
                     + bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
        enddo
        buffer(nz-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                     + bciz6*(uz(i,j,nz)-uz(i,j,nz-3))
        buffer(nz) = aciz6*(uz(i,j,nz)-uz(i,j,nz-1)) &
                   + bciz6*(uz(i,j,1)-uz(i,j,nz-2))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nz)
        do concurrent (k=1:nz)
           tz(i,j,k) = buffer(k)
        enddo
     enddo

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then
        do concurrent (j=1:ny, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = zero
           buffer(2) = aciz6*(uz(i,j,2)-uz(i,j,1)) &
                     + bciz6*(uz(i,j,3)-uz(i,j,1))
           do concurrent (k=3:nz-2)
              buffer(k) = aciz6*(uz(i,j,k)-uz(i,j,k-1)) &
                        + bciz6*(uz(i,j,k+1)-uz(i,j,k-2))
           enddo
           buffer(nz-1) = aciz6*(uz(i,j,nz-1)-uz(i,j,nz-2)) &
                        + bciz6*(uz(i,j,nz-1)-uz(i,j,nz-3))
           buffer(nz) = zero

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nz)
           do concurrent (k=1:nz)
              tz(i,j,k) = buffer(k)
           enddo
        enddo

     endif
  endif

end subroutine derzpv

!********************************************************************
!
subroutine interzpv(tz,uz,x3dop,nx,ny,nzm,nz)

  USE x3d_operator_z_data

  implicit none

  !$acc routine(thomas1d_0, thomas1d_12) seq

  ! Arguments
  integer, intent(in) :: nx, ny, nz, nzm
  real(mytype), intent(out), dimension(nx,ny,nz) :: tz
  real(mytype), intent(in), dimension(nx,ny,nzm) :: uz
  type(x3doperator1d), intent(in) :: x3dop

  ! Local variables
  integer :: i, j, k
  real(mytype), dimension(nz) :: buffer, ff, ss, ww, pp

  if (nz==1) then
     do concurrent(k=1:nz, j=1:ny, i=1:nx)
        tz(i,j,k) = uz(i,j,k)
     enddo
     return
  endif

  do concurrent (k=1:nz)
     ff(k) = x3dop%f(k)
     ss(k) = x3dop%s(k)
     ww(k) = x3dop%w(k)
     if (allocated(x3dop%periodic)) pp(k) = x3dop%periodic(k)
  end do

  if (nclz) then
     ! nzm = nz
     do concurrent (j=1:ny, i=1:nx) local(buffer)

        ! Compute r.h.s.
        buffer(1) = aiciz6*(uz(i,j,1)+uz(i,j,nz)) &
                  + biciz6*(uz(i,j,2)+uz(i,j,nz-1)) &
                  + ciciz6*(uz(i,j,3)+uz(i,j,nz-2)) &
                  + diciz6*(uz(i,j,4)+uz(i,j,nz-3))
        buffer(2) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                  + biciz6*(uz(i,j,3)+uz(i,j,nz)) &
                  + ciciz6*(uz(i,j,4)+uz(i,j,nz-1)) &
                  + diciz6*(uz(i,j,5)+uz(i,j,nz-2))
        buffer(3) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                  + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                  + ciciz6*(uz(i,j,5)+uz(i,j,nz)) &
                  + diciz6*(uz(i,j,6)+uz(i,j,nz-1))
        buffer(4) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                  + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                  + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                  + diciz6*(uz(i,j,7)+uz(i,j,nz))
        do concurrent (k=5:nz-3)
           buffer(k) = aiciz6*(uz(i,j,k)+uz(i,j,k-1)) &
                     + biciz6*(uz(i,j,k+1)+uz(i,j,k-2)) &
                     + ciciz6*(uz(i,j,k+2)+uz(i,j,k-3)) &
                     + diciz6*(uz(i,j,k+3)+uz(i,j,k-4))
        enddo
        buffer(nz-2) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                     + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                     + ciciz6*(uz(i,j,nz)+uz(i,j,nz-5)) &
                     + diciz6*(uz(i,j,1)+uz(i,j,nz-6))
        buffer(nz-1) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                     + biciz6*(uz(i,j,nz)+uz(i,j,nz-3)) &
                     + ciciz6*(uz(i,j,1)+uz(i,j,nz-4)) &
                     + diciz6*(uz(i,j,2)+uz(i,j,nz-5))
        buffer(nz  ) = aiciz6*(uz(i,j,nz)+uz(i,j,nz-1)) &
                     + biciz6*(uz(i,j,1)+uz(i,j,nz-2)) &
                     + ciciz6*(uz(i,j,2)+uz(i,j,nz-3)) &
                     + diciz6*(uz(i,j,3)+uz(i,j,nz-4))

        ! Solve tri-diagonal system
        call thomas1d(buffer, ff, ss, ww, pp, x3dop%alfa, nz)
        do concurrent (k=1:nz)
           tz(i,j,k) = buffer(k)
        enddo
     enddo

  else
     ! nzm = nz-1
     if (x3dop%npaire==1) then
        do concurrent (j=1:ny, i=1:nx) local(buffer)

           ! Compute r.h.s.
           buffer(1) = aiciz6*(uz(i,j,1)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,2)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,3)+uz(i,j,3)) &
                     + diciz6*(uz(i,j,4)+uz(i,j,4))
           buffer(2) = aiciz6*(uz(i,j,2)+uz(i,j,1)) &
                     + biciz6*(uz(i,j,3)+uz(i,j,1))&
                     + ciciz6*(uz(i,j,4)+uz(i,j,2))&
                     + diciz6*(uz(i,j,5)+uz(i,j,3))
           buffer(3) = aiciz6*(uz(i,j,3)+uz(i,j,2)) &
                     + biciz6*(uz(i,j,4)+uz(i,j,1)) &
                     + ciciz6*(uz(i,j,5)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,6)+uz(i,j,2))
           buffer(4) = aiciz6*(uz(i,j,4)+uz(i,j,3)) &
                     + biciz6*(uz(i,j,5)+uz(i,j,2)) &
                     + ciciz6*(uz(i,j,6)+uz(i,j,1)) &
                     + diciz6*(uz(i,j,7)+uz(i,j,1))
           do concurrent (k=5:nz-4)
              buffer(k) = aiciz6*(uz(i,j,k)+uz(i,j,k-1)) &
                        + biciz6*(uz(i,j,k+1)+uz(i,j,k-2)) &
                        + ciciz6*(uz(i,j,k+2)+uz(i,j,k-3)) &
                        + diciz6*(uz(i,j,k+3)+uz(i,j,k-4))
           enddo
           buffer(nz-3) = aiciz6*(uz(i,j,nz-3)+uz(i,j,nz-4)) &
                        + biciz6*(uz(i,j,nz-2)+uz(i,j,nz-5)) &
                        + ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-6)) &
                        + diciz6*(uz(i,j,nz-1)+uz(i,j,nz-7))
           buffer(nz-2) = aiciz6*(uz(i,j,nz-2)+uz(i,j,nz-3)) &
                        + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-4)) &
                        + ciciz6*(uz(i,j,nz-1)+uz(i,j,nz-5)) &
                        + diciz6*(uz(i,j,nz-2)+uz(i,j,nz-6))
           buffer(nz-1) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-2)) &
                        + biciz6*(uz(i,j,nz-1)+uz(i,j,nz-3)) &
                        + ciciz6*(uz(i,j,nz-2)+uz(i,j,nz-4)) &
                        + diciz6*(uz(i,j,nz-3)+uz(i,j,nz-5))
           buffer(nz  ) = aiciz6*(uz(i,j,nz-1)+uz(i,j,nz-1)) &
                        + biciz6*(uz(i,j,nz-2)+uz(i,j,nz-2)) &
                        + ciciz6*(uz(i,j,nz-3)+uz(i,j,nz-3)) &
                        + diciz6*(uz(i,j,nz-4)+uz(i,j,nz-4))

           ! Solve tri-diagonal system
           call thomas1d(buffer, ff, ss, ww, nz)
           do concurrent (k=1:nz)
              tz(i,j,k) = buffer(k)
           enddo
        enddo

     endif
  endif

end subroutine interzpv

end module x3d_staggered
