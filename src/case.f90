!################################################################################
!This file is part of Xcompact3d
!
!Xcompact3d
!Copyright (c) 2012-2022, Xcompact3d
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!BSD 3-Clause License
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.!
!
!3. Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
