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
module tools

  implicit none

  private

  public :: error_l1_l2_linf

  interface error_l1_l2_linf
     module procedure error_l1_l2_linf_xsize
     module procedure error_l1_l2_linf_generic
  end interface error_l1_l2_linf

contains
  !##################################################################
  !##################################################################
  subroutine error_L1_L2_Linf_xsize(err, l1, l2, linf)

    USE variables   , only : nx, ny, nz
    use decomp_2d, only : mytype
    use decomp_2d   , only : xsize

    implicit none
      
    real(mytype),intent(in),dimension(xsize(1),xsize(2),xsize(3)) :: err
    real(mytype),intent(out) :: l1, l2, linf

    call error_L1_L2_Linf_generic(err, l1, l2, linf, xsize(1), xsize(2), xsize(3), nx*ny*nz)

  end subroutine error_L1_L2_Linf_xsize
  !##################################################################

  !##################################################################
  subroutine error_L1_L2_Linf_generic(err, l1, l2, linf, n1, n2, n3, ntot)

    ! Compute L1, L2 and Linf norm of given 3D array
    USE param, only : zero
    use decomp_2d, only : mytype, real_type, decomp_2d_abort
    USE MPI
    
    implicit none
      
    real(mytype),intent(in),dimension(n1,n2,n3) :: err
    real(mytype),intent(out) :: l1, l2, linf
    integer, intent(in) :: n1, n2, n3, ntot

    integer :: i, j, k, code
    real(mytype) :: l1l2(2)

    l1 = zero
    l2 = zero
    linf = zero

    do k = 1,n3
      do j = 1,n2
        do i = 1,n1
          l1 = l1 + abs(err(i,j,k))
          l2 = l2 + err(i,j,k)*err(i,j,k)
          linf = max(linf, abs(err(i,j,k)))
        enddo
      enddo
    enddo

    ! Parallel, MPI_SUM
    l1l2 = (/l1, l2/)
    code = 0
    call MPI_ALLREDUCE(MPI_IN_PLACE,l1l2,2,real_type,MPI_SUM,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")
    l1 = l1l2(1) / ntot
    l2 = sqrt(l1l2(2) / ntot)
    ! Parallel, MPI_MAX
    call MPI_ALLREDUCE(MPI_IN_PLACE,linf,1,real_type,MPI_MAX,MPI_COMM_WORLD,code)
    if (code /= 0) call decomp_2d_abort(__FILE__, __LINE__, code, "MPI_ALLREDUCE")

  end subroutine error_L1_L2_Linf_generic
end module tools
