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

module x3d_operator_z_data

   use decomp_2d, only: mytype
   use param, only: zero

   implicit none

   ! Make everything public unless declared private
   public

   ! This is needed for r.h.s. estimation
   real(mytype) :: alcaiz6, aciz6, bciz6
   real(mytype) :: ailcaiz6, aiciz6, biciz6, ciciz6, diciz6
   real(mytype) :: alfa1z, af1z, bf1z, cf1z, df1z, alfa2z, af2z, alfanz, afnz, bfnz
   real(mytype) :: cfnz, dfnz, alfamz, afmz, alfakz, afkz, bfkz, alsa1z, as1z, bs1z
   real(mytype) :: cs1z, ds1z, alsa2z, as2z, alsanz, asnz, bsnz, csnz, dsnz, alsamz
   real(mytype) :: asmz, alsa3z, as3z, bs3z, alsatz, astz, bstz
   !O6SVV
   real(mytype) :: alsa4z, as4z, bs4z, cs4z
   real(mytype) :: alsattz, asttz, bsttz, csttz
   real(mytype) :: alsakz, askz, bskz, cskz, dskz
   real(mytype) :: alsakzt, askzt, bskzt, cskzt, dskzt
   ! Filter
   real(mytype) :: fial1z, fia1z, fib1z, fic1z, fid1z, fie1z, fif1z ! Coefficients for filter at boundary point 1
   real(mytype) :: fial2z, fia2z, fib2z, fic2z, fid2z, fie2z, fif2z ! Coefficients for filter at boundary point 2
   real(mytype) :: fial3z, fia3z, fib3z, fic3z, fid3z, fie3z, fif3z ! Coefficients for filter at boundary point 3
   real(mytype) :: fialkz, fiakz, fibkz, fickz, fidkz               ! Coefficient for filter at interior points
   real(mytype) :: fialnz, fianz, fibnz, ficnz, fidnz, fienz, fifnz ! Coefficient for filter at boundary point n
   real(mytype) :: fialmz, fiamz, fibmz, ficmz, fidmz, fiemz, fifmz ! Coefficient for filter at boundary point m=n-1
   real(mytype) :: fialpz, fiapz, fibpz, ficpz, fidpz, fiepz, fifpz ! Coefficient for filter at boundary point p=n-2

   ! This is needed when solving the linear system
   real(mytype), target, allocatable, dimension(:) :: ffz, sfz, fsz, fwz, ssz, swz
   real(mytype), target, allocatable, dimension(:) :: ffzp, sfzp, fszp, fwzp, sszp, swzp

   real(mytype), allocatable, dimension(:) :: fiffz, fifcz, fifbz, fisfz, fiscz, fisbz, fifsz, fifwz, fissz, fiswz
   real(mytype), allocatable, dimension(:) :: fiffzp, fifszp, fifwzp, fisfzp, fisszp, fiswzp

   real(mytype), allocatable, dimension(:) :: ffzS, sfzS, fszS, fwzS, sszS, swzS
   real(mytype), allocatable, dimension(:) :: ffzpS, sfzpS, fszpS, fwzpS, sszpS, swzpS

   real(mytype), target, allocatable, dimension(:) :: cfz6, ccz6, cbz6, cfzp6, cszp6, cwzp6, csz6
   real(mytype), target, allocatable, dimension(:) :: cwz6, cifz6, cicz6, cibz6, cifzp6, ciszp6, &
                                                      ciwzp6, cisz6, ciwz6
   real(mytype), target, allocatable, dimension(:) :: cfi6z, cci6z, cbi6z, cfip6z, csip6z, cwip6z, &
                                                      csi6z, cwi6z, cifi6z, cici6z
   real(mytype), target, allocatable, dimension(:) :: cibi6z, cifip6z, cisip6z, ciwip6z, cisi6z, ciwi6z

contains

   !
   ! Allocate memory and prepare arrays
   !
   subroutine x3d_operator_z_data_init()

      use variables, only: nz, nzm

      implicit none

      allocate (ffz(nz))
      ffz = zero
      allocate (sfz(nz))
      sfz = zero
      allocate (fsz(nz))
      fsz = zero
      allocate (fwz(nz))
      fwz = zero
      allocate (ssz(nz))
      ssz = zero
      allocate (swz(nz))
      swz = zero

      allocate (ffzp(nz))
      ffzp = zero
      allocate (sfzp(nz))
      sfzp = zero
      allocate (fszp(nz))
      fszp = zero
      allocate (fwzp(nz))
      fwzp = zero
      allocate (sszp(nz))
      sszp = zero
      allocate (swzp(nz))
      swzp = zero

      allocate (ffzS(nz))
      ffzS = zero
      allocate (sfzS(nz))
      sfzS = zero
      allocate (fszS(nz))
      fszS = zero
      allocate (fwzS(nz))
      fwzS = zero
      allocate (sszS(nz))
      sszS = zero
      allocate (swzS(nz))
      swzS = zero

      allocate (ffzpS(nz))
      ffzpS = zero
      allocate (sfzpS(nz))
      sfzpS = zero
      allocate (fszpS(nz))
      fszpS = zero
      allocate (fwzpS(nz))
      fwzpS = zero
      allocate (sszpS(nz))
      sszpS = zero
      allocate (swzpS(nz))
      swzpS = zero

      allocate (cfz6(nzm))
      cfz6 = zero
      allocate (ccz6(nzm))
      ccz6 = zero
      allocate (cbz6(nzm))
      cbz6 = zero
      allocate (cfzp6(nzm))
      cfzp6 = zero
      allocate (cszp6(nzm))
      cszp6 = zero
      allocate (cwzp6(nzm))
      cwzp6 = zero
      allocate (csz6(nzm))
      csz6 = zero

      allocate (cwz6(nzm))
      cwz6 = zero
      allocate (cifz6(nzm))
      cifz6 = zero
      allocate (cicz6(nzm))
      cicz6 = zero
      allocate (cibz6(nzm))
      cibz6 = zero
      allocate (cifzp6(nzm))
      cifzp6 = zero
      allocate (ciszp6(nzm))
      ciszp6 = zero
      allocate (ciwzp6(nzm))
      ciwzp6 = zero
      allocate (cisz6(nzm))
      cisz6 = zero
      allocate (ciwz6(nzm))
      ciwz6 = zero

      allocate (cfi6z(nz))
      cfi6z = zero
      allocate (cci6z(nz))
      cci6z = zero
      allocate (cbi6z(nz))
      cbi6z = zero
      allocate (cfip6z(nz))
      cfip6z = zero
      allocate (csip6z(nz))
      csip6z = zero
      allocate (cwip6z(nz))
      cwip6z = zero
      allocate (csi6z(nz))
      csi6z = zero
      allocate (cwi6z(nz))
      cwi6z = zero
      allocate (cifi6z(nz))
      cifi6z = zero
      allocate (cici6z(nz))
      cici6z = zero

      allocate (cibi6z(nz))
      cibi6z = zero
      allocate (cifip6z(nz))
      cifip6z = zero
      allocate (cisip6z(nz))
      cisip6z = zero
      allocate (ciwip6z(nz))
      ciwip6z = zero
      allocate (cisi6z(nz))
      cisi6z = zero
      allocate (ciwi6z(nz))
      ciwi6z = zero

   end subroutine x3d_operator_z_data_init

   !
   ! Free memory
   !
   subroutine x3d_operator_z_data_finalize()

      implicit none

      deallocate (ffz)
      deallocate (sfz)
      deallocate (fsz)
      deallocate (fwz)
      deallocate (ssz)
      deallocate (swz)

      deallocate (ffzp)
      deallocate (sfzp)
      deallocate (fszp)
      deallocate (fwzp)
      deallocate (sszp)
      deallocate (swzp)

      deallocate (ffzS)
      deallocate (sfzS)
      deallocate (fszS)
      deallocate (fwzS)
      deallocate (sszS)
      deallocate (swzS)

      deallocate (ffzpS)
      deallocate (sfzpS)
      deallocate (fszpS)
      deallocate (fwzpS)
      deallocate (sszpS)
      deallocate (swzpS)

      deallocate (cfz6)
      deallocate (ccz6)
      deallocate (cbz6)
      deallocate (cfzp6)
      deallocate (cszp6)
      deallocate (cwzp6)
      deallocate (csz6)

      deallocate (cwz6)
      deallocate (cifz6)
      deallocate (cicz6)
      deallocate (cibz6)
      deallocate (cifzp6)
      deallocate (ciszp6)
      deallocate (ciwzp6)
      deallocate (cisz6)
      deallocate (ciwz6)

      deallocate (cfi6z)
      deallocate (cci6z)
      deallocate (cbi6z)
      deallocate (cfip6z)
      deallocate (csip6z)
      deallocate (cwip6z)
      deallocate (csi6z)
      deallocate (cwi6z)
      deallocate (cifi6z)
      deallocate (cici6z)

      deallocate (cibi6z)
      deallocate (cifip6z)
      deallocate (cisip6z)
      deallocate (ciwip6z)
      deallocate (cisi6z)
      deallocate (ciwi6z)

   end subroutine x3d_operator_z_data_finalize

end module x3d_operator_z_data
