!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module time_integrators

  implicit none

  private
  public :: int_time

contains

  subroutine intt(var1,dvar1)

    use MPI
    use param    , only : itr, ntime, gdt
    use variables
    use decomp_2d, only : mytype, xsize

    implicit none

    !! INPUT / OUTPUT
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

    !! LOCAL
    integer :: i,j,k


    ! for the moment we just use euler
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      var1(i,j,k)=gdt(itr)*dvar1(i,j,k,1)+var1(i,j,k)
    enddo
    
    !if (iimplicit.ge.1) then
    !   !>>> (semi)implicit Y diffusion

    !   if (present(isc)) then
    !      is = isc
    !   else
    !      is = 0
    !   endif
    !   if (present(npaire).and.present(forcing1)) then
    !      call inttimp(var1, dvar1, npaire=npaire, isc=is, forcing1=forcing1)
    !   else if (present(npaire)) then
    !      call inttimp(var1, dvar1, npaire=npaire, isc=is)
    !   else
    !      if (nrank  == 0) write(*,*) "Error in intt call."
    !      call MPI_ABORT(MPI_COMM_WORLD,code,ierror); stop
    !   endif

    !elseif (itimescheme.eq.1) then
    !   !>>> Euler

    !   var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
    !elseif(itimescheme.eq.2) then
    !   !>>> Adam-Bashforth second order (AB2)

    !   ! Do first time step with Euler
    !   if(itime.eq.1.and.irestart.eq.0) then
    !      var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
    !   else
    !      var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
    !   endif
    !   dvar1(:,:,:,2)=dvar1(:,:,:,1)
    !elseif(itimescheme.eq.3) then
    !   !>>> Adams-Bashforth third order (AB3)

    !   ! Do first time step with Euler
    !   if(itime.eq.1.and.irestart.eq.0) then
    !      var1(:,:,:)=dt*dvar1(:,:,:,1)+var1(:,:,:)
    !   elseif(itime.eq.2.and.irestart.eq.0) then
    !      ! Do second time step with AB2
    !      var1(:,:,:)=onepfive*dt*dvar1(:,:,:,1)-half*dt*dvar1(:,:,:,2)+var1(:,:,:)
    !      dvar1(:,:,:,3)=dvar1(:,:,:,2)
    !   else
    !      ! Finally using AB3
    !      var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+cdt(itr)*dvar1(:,:,:,3)+var1(:,:,:)
    !      dvar1(:,:,:,3)=dvar1(:,:,:,2)
    !   endif
    !   dvar1(:,:,:,2)=dvar1(:,:,:,1)
    !elseif(itimescheme.eq.4) then
    !   !>>> Adams-Bashforth fourth order (AB4)

    !   if (nrank==0) then
    !      write(*,*) "AB4 not implemented!"
    !      stop
    !   endif

    !   !>>> Runge-Kutta (low storage) RK3
    !elseif(itimescheme.eq.5) then
    !   if(itr.eq.1) then
    !      var1(:,:,:)=gdt(itr)*dvar1(:,:,:,1)+var1(:,:,:)
    !   else
    !      var1(:,:,:)=adt(itr)*dvar1(:,:,:,1)+bdt(itr)*dvar1(:,:,:,2)+var1(:,:,:)
    !   endif
    !   dvar1(:,:,:,2)=dvar1(:,:,:,1)
    !   !>>> Runge-Kutta (low storage) RK4
    !elseif(itimescheme.eq.6) then

    !   if (nrank==0) then
    !      write(*,*) "RK4 not implemented!"
    !      STOP
    !   endif

    !else

    !   if (nrank==0) then
    !      write(*,*) "Unrecognised itimescheme: ", itimescheme
    !      STOP
    !   endif

    !endif


    return

  end subroutine intt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time
  !! DESCRIPTION: 
  !!      INPUTS: 
  !!     OUTPUTS:
  !!       NOTES: 
  !!      AUTHOR:  
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine int_time(ux1, uy1, uz1, dux1, duy1, duz1)

    use decomp_2d, only : mytype, xsize, nrank
    use param, only : ntime

    implicit none

    !! input/output
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    !! output
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)) :: ux1, uy1, uz1

    !! LOCAL
#ifdef DEBUG
    real(mytype) avg_param
    if (nrank .eq. 0) write(*,*)'## Init int_time'
#endif

    call int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)


  ENDSUBROUTINE int_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !!  SUBROUTINE: int_time_momentum
  !! DESCRIPTION: Integrates the momentum equations in time by calling time
  !!              integrator.
  !!      INPUTS: dux1, duy1, duz1 - the RHS(s) of the momentum equations
  !!     OUTPUTS: ux1,   uy1,  uz1 - the intermediate momentum state.
  !!       NOTES: This is integrating the MOMENTUM in time (!= velocity)
  !!      AUTHOR: Paul Bartholomew
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine int_time_momentum(ux1, uy1, uz1, dux1, duy1, duz1)

    use param
    use variables
    use decomp_2d

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

    call intt(ux1, dux1)
    call intt(uy1, duy1)
    call intt(uz1, duz1)

  endsubroutine int_time_momentum

end module time_integrators
