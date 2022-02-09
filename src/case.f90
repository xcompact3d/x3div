!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module case

  use param
  use decomp_2d, only : mytype
  use variables

  implicit none

  private ! All functions/subroutines private by default
  public :: init

contains
  !##################################################################
  subroutine init (ux1, uy1, uz1, dux1, duy1, duz1, &
       pp3, px1, py1, pz1)

    use mom, only : vel
    use decomp_2d, only : xsize, ph1

    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1
    real(mytype),dimension(ph1%zst(1):ph1%zen(1), ph1%zst(2):ph1%zen(2), nzm, npress) :: pp3
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: px1, py1, pz1

    integer :: it, is
    integer :: i, j, k 

    !! Zero out the pressure field
    do concurrent (k=1:nzm, j=ph1%zst(2):ph1%zen(2), i=ph1%zst(1):ph1%zen(1))
      pp3(i,j,k,1) = zero
    enddo
    
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1)) 
      px1(i,j,k) = zero
      py1(i,j,k) = zero
      pz1(i,j,k) = zero
    enddo

    call vel(ux1, uy1, uz1)
    
    !! Setup old arrays
    do it = 1, ntime
      do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1)) 
        dux1(i,j,k,it)=ux1(i,j,k)
        duy1(i,j,k,it)=uy1(i,j,k)
        duz1(i,j,k,it)=uz1(i,j,k)
      enddo
    enddo


  end subroutine init
end module case

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! case.f90 ends here
