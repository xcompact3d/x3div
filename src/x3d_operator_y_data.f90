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

module x3d_operator_y_data

   use decomp_2d, only: mytype
   use param, only: zero

   implicit none

   ! Make everything public unless declared private
   public

   ! This is needed for r.h.s. estimation
   real(mytype) :: alcaiy6, aciy6, bciy6
   real(mytype) :: ailcaiy6, aiciy6, biciy6, ciciy6, diciy6
   real(mytype) :: alfa1y, af1y, bf1y, cf1y, df1y, alfa2y, af2y, alfany, afny, bfny
   real(mytype) :: cfny, dfny, alfamy, afmy, alfajy, afjy, bfjy, alsa1y, as1y, bs1y
   real(mytype) :: cs1y, ds1y, alsa2y, as2y, alsany, asny, bsny, csny, dsny, alsamy
   real(mytype) :: asmy, alsa3y, as3y, bs3y, alsaty, asty, bsty
   !O6SVV
   real(mytype) :: alsa4y, as4y, bs4y, cs4y
   real(mytype) :: alsatty, astty, bstty, cstty
   real(mytype) :: alsajy, asjy, bsjy, csjy, dsjy
   real(mytype) :: alsajyt, asjyt, bsjyt, csjyt, dsjyt
   ! Filter
   real(mytype) :: fial1y, fia1y, fib1y, fic1y, fid1y, fie1y, fif1y ! Coefficients for filter at boundary point 1
   real(mytype) :: fial2y, fia2y, fib2y, fic2y, fid2y, fie2y, fif2y ! Coefficients for filter at boundary point 2
   real(mytype) :: fial3y, fia3y, fib3y, fic3y, fid3y, fie3y, fif3y ! Coefficients for filter at boundary point 3
   real(mytype) :: fialjy, fiajy, fibjy, ficjy, fidjy               ! Coefficient for filter at interior points
   real(mytype) :: fialny, fiany, fibny, ficny, fidny, fieny, fifny ! Coefficient for filter at boundary point n
   real(mytype) :: fialmy, fiamy, fibmy, ficmy, fidmy, fiemy, fifmy ! Coefficient for filter at boundary point m=n-1
   real(mytype) :: fialpy, fiapy, fibpy, ficpy, fidpy, fiepy, fifpy ! Coefficient for filter at boundary point p=n-2

   ! This is needed when solving the linear system
   real(mytype), target, allocatable, dimension(:) :: ffy, sfy, fsy, fwy, ssy, swy
   real(mytype), target, allocatable, dimension(:) :: ffyp, sfyp, fsyp, fwyp, ssyp, swyp

   real(mytype), allocatable, dimension(:) :: fiffy, fifcy, fifby, fisfy, fiscy, fisby, fifsy, fifwy, fissy, fiswy
   real(mytype), allocatable, dimension(:) :: fiffyp, fifsyp, fifwyp, fisfyp, fissyp, fiswyp

   real(mytype), allocatable, dimension(:) :: ffyS, sfyS, fsyS, fwyS, ssyS, swyS
   real(mytype), allocatable, dimension(:) :: ffypS, sfypS, fsypS, fwypS, ssypS, swypS

   real(mytype), target, allocatable, dimension(:) :: cfy6, ccy6, cby6, cfyp6, csyp6, cwyp6, csy6
   real(mytype), target, allocatable, dimension(:) :: cwy6, cify6, cicy6, ciby6, cifyp6, cisyp6, &
                                                      ciwyp6, cisy6, ciwy6
   real(mytype), target, allocatable, dimension(:) :: cfi6y, cci6y, cbi6y, cfip6y, csip6y, cwip6y, &
                                                      csi6y, cwi6y, cifi6y, cici6y
   real(mytype), target, allocatable, dimension(:) :: cibi6y, cifip6y, cisip6y, ciwip6y, cisi6y, ciwi6y

contains

   !
   ! Allocate memory and prepare arrays
   !
   subroutine x3d_operator_y_data_init()

      use variables, only: ny, nym

      implicit none

      allocate (ffy(ny))
      ffy = zero
      allocate (sfy(ny))
      sfy = zero
      allocate (fsy(ny))
      fsy = zero
      allocate (fwy(ny))
      fwy = zero
      allocate (ssy(ny))
      ssy = zero
      allocate (swy(ny))
      swy = zero

      allocate (ffyp(ny))
      ffyp = zero
      allocate (sfyp(ny))
      sfyp = zero
      allocate (fsyp(ny))
      fsyp = zero
      allocate (fwyp(ny))
      fwyp = zero
      allocate (ssyp(ny))
      ssyp = zero
      allocate (swyp(ny))
      swyp = zero

      allocate (ffyS(ny))
      ffyS = zero
      allocate (sfyS(ny))
      sfyS = zero
      allocate (fsyS(ny))
      fsyS = zero
      allocate (fwyS(ny))
      fwyS = zero
      allocate (ssyS(ny))
      ssyS = zero
      allocate (swyS(ny))
      swyS = zero

      allocate (ffypS(ny))
      ffypS = zero
      allocate (sfypS(ny))
      sfypS = zero
      allocate (fsypS(ny))
      fsypS = zero
      allocate (fwypS(ny))
      fwypS = zero
      allocate (ssypS(ny))
      ssypS = zero
      allocate (swypS(ny))
      swypS = zero

      allocate (cfy6(nym))
      cfy6 = zero
      allocate (ccy6(nym))
      ccy6 = zero
      allocate (cby6(nym))
      cby6 = zero
      allocate (cfyp6(nym))
      cfyp6 = zero
      allocate (csyp6(nym))
      csyp6 = zero
      allocate (cwyp6(nym))
      cwyp6 = zero
      allocate (csy6(nym))
      csy6 = zero

      allocate (cwy6(nym))
      cwy6 = zero
      allocate (cify6(nym))
      cify6 = zero
      allocate (cicy6(nym))
      cicy6 = zero
      allocate (ciby6(nym))
      ciby6 = zero
      allocate (cifyp6(nym))
      cifyp6 = zero
      allocate (cisyp6(nym))
      cisyp6 = zero
      allocate (ciwyp6(nym))
      ciwyp6 = zero
      allocate (cisy6(nym))
      cisy6 = zero
      allocate (ciwy6(nym))
      ciwy6 = zero

      allocate (cfi6y(ny))
      cfi6y = zero
      allocate (cci6y(ny))
      cci6y = zero
      allocate (cbi6y(ny))
      cbi6y = zero
      allocate (cfip6y(ny))
      cfip6y = zero
      allocate (csip6y(ny))
      csip6y = zero
      allocate (cwip6y(ny))
      cwip6y = zero
      allocate (csi6y(ny))
      csi6y = zero
      allocate (cwi6y(ny))
      cwi6y = zero
      allocate (cifi6y(ny))
      cifi6y = zero
      allocate (cici6y(ny))
      cici6y = zero

      allocate (cibi6y(ny))
      cibi6y = zero
      allocate (cifip6y(ny))
      cifip6y = zero
      allocate (cisip6y(ny))
      cisip6y = zero
      allocate (ciwip6y(ny))
      ciwip6y = zero
      allocate (cisi6y(ny))
      cisi6y = zero
      allocate (ciwi6y(ny))
      ciwi6y = zero

   end subroutine x3d_operator_y_data_init

   !
   ! Free memory
   !
   subroutine x3d_operator_y_data_finalize()

      implicit none

      deallocate (ffy)
      deallocate (sfy)
      deallocate (fsy)
      deallocate (fwy)
      deallocate (ssy)
      deallocate (swy)

      deallocate (ffyp)
      deallocate (sfyp)
      deallocate (fsyp)
      deallocate (fwyp)
      deallocate (ssyp)
      deallocate (swyp)

      deallocate (ffyS)
      deallocate (sfyS)
      deallocate (fsyS)
      deallocate (fwyS)
      deallocate (ssyS)
      deallocate (swyS)

      deallocate (ffypS)
      deallocate (sfypS)
      deallocate (fsypS)
      deallocate (fwypS)
      deallocate (ssypS)
      deallocate (swypS)

      deallocate (cfy6)
      deallocate (ccy6)
      deallocate (cby6)
      deallocate (cfyp6)
      deallocate (csyp6)
      deallocate (cwyp6)
      deallocate (csy6)

      deallocate (cwy6)
      deallocate (cify6)
      deallocate (cicy6)
      deallocate (ciby6)
      deallocate (cifyp6)
      deallocate (cisyp6)
      deallocate (ciwyp6)
      deallocate (cisy6)
      deallocate (ciwy6)

      deallocate (cfi6y)
      deallocate (cci6y)
      deallocate (cbi6y)
      deallocate (cfip6y)
      deallocate (csip6y)
      deallocate (cwip6y)
      deallocate (csi6y)
      deallocate (cwi6y)
      deallocate (cifi6y)
      deallocate (cici6y)

      deallocate (cibi6y)
      deallocate (cifip6y)
      deallocate (cisip6y)
      deallocate (ciwip6y)
      deallocate (cisi6y)
      deallocate (ciwi6y)

   end subroutine x3d_operator_y_data_finalize

end module x3d_operator_y_data
