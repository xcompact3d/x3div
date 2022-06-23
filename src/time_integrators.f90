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
    use nvtx

    implicit none

    !! INPUT / OUTPUT
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: var1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dvar1

    !! LOCAL
    integer :: i,j,k
    integer :: imax,jmax,kmax
    imax = xsize(1)
    jmax = xsize(2)
    kmax = xsize(3)


    call nvtxStartRange("time loop orig")
    ! We have only Euler
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      var1(i,j,k)=gdt(itr)*dvar1(i,j,k,1)+var1(i,j,k)
    enddo
    call nvtxEndRange
    
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
