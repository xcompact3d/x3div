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

module x3d_operator_x_data

   use decomp_2d, only: mytype
   use param, only: zero

   implicit none

   ! Make everything public unless declared private
   public

   ! This is needed for r.h.s. estimation
   real(mytype) :: alcaix6, acix6, bcix6
   real(mytype) :: ailcaix6, aicix6, bicix6, cicix6, dicix6
   real(mytype) :: alfa1x, af1x, bf1x, cf1x, df1x, alfa2x, af2x, alfanx, afnx, bfnx
   real(mytype) :: cfnx, dfnx, alfamx, afmx, alfaix, afix, bfix, alsa1x, as1x, bs1x
   real(mytype) :: cs1x, ds1x, alsa2x, as2x, alsanx, asnx, bsnx, csnx, dsnx, alsamx
   real(mytype) :: asmx, alsa3x, as3x, bs3x, alsatx, astx, bstx
   !O6SVV
   real(mytype) :: alsa4x, as4x, bs4x, cs4x
   real(mytype) :: alsattx, asttx, bsttx, csttx
   real(mytype) :: alsaix, asix, bsix, csix, dsix
   real(mytype) :: alsaixt, asixt, bsixt, csixt, dsixt
   ! Filter
   real(mytype) :: fial1x, fia1x, fib1x, fic1x, fid1x, fie1x, fif1x  ! Coefficients for filter at boundary point 1
   real(mytype) :: fial2x, fia2x, fib2x, fic2x, fid2x, fie2x, fif2x  ! Coefficients for filter at boundary point 2
   real(mytype) :: fial3x, fia3x, fib3x, fic3x, fid3x, fie3x, fif3x  ! Coefficients for filter at boundary point 3
   real(mytype) :: fialix, fiaix, fibix, ficix, fidix                ! Coefficient for filter at interior points
   real(mytype) :: fialnx, fianx, fibnx, ficnx, fidnx, fienx, fifnx  ! Coefficient for filter at boundary point n
   real(mytype) :: fialmx, fiamx, fibmx, ficmx, fidmx, fiemx, fifmx  ! Coefficient for filter at boundary point m=n-1
   real(mytype) :: fialpx, fiapx, fibpx, ficpx, fidpx, fiepx, fifpx  ! Coefficient for filter at boundary point p=n-2

   ! This is needed when solving the linear system
   real(mytype), target, allocatable, dimension(:) :: ffx, sfx, fsx, fwx, ssx, swx
   real(mytype), target, allocatable, dimension(:) :: ffxp, sfxp, fsxp, fwxp, ssxp, swxp
   !
   real(mytype), allocatable, dimension(:) :: fiffx, fifcx, fifbx, fisfx, fiscx, fisbx, fifsx, fifwx, fissx, fiswx
   real(mytype), allocatable, dimension(:) :: fiffxp, fifsxp, fifwxp, fisfxp, fissxp, fiswxp
   !
   real(mytype), allocatable, dimension(:) :: ffxS, sfxS, fsxS, fwxS, ssxS, swxS
   real(mytype), allocatable, dimension(:) :: ffxpS, sfxpS, fsxpS, fwxpS, ssxpS, swxpS
   !
   real(mytype), target, allocatable, dimension(:) :: cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                                                      cwxp6, csx6, cwx6, cifx6, cicx6, cisx6
   real(mytype), target, allocatable, dimension(:) :: cibx6, cifxp6, cisxp6, ciwx6
   real(mytype), target, allocatable, dimension(:) :: cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                                                      cwi6, cifi6, cici6, cibi6, cifip6
   real(mytype), target, allocatable, dimension(:) :: cisip6, ciwip6, cisi6, ciwi6

contains

   !
   ! Allocate memory
   !
   subroutine x3d_operator_x_data_init()

      use variables, only: nx, nxm

      implicit none

      allocate (ffx(nx))
      ffx = zero
      allocate (sfx, source=ffx)
      allocate (fsx, source=ffx)
      allocate (fwx, source=ffx)
      allocate (ssx, source=ffx)
      allocate (swx, source=ffx)

      allocate (ffxp, source=ffx)
      allocate (sfxp, source=ffx)
      allocate (fsxp, source=ffx)
      allocate (fwxp, source=ffx)
      allocate (ssxp, source=ffx)
      allocate (swxp, source=ffx)

      allocate (ffxS, source=ffx)
      allocate (sfxS, source=ffx)
      allocate (fsxS, source=ffx)
      allocate (fwxS, source=ffx)
      allocate (ssxS, source=ffx)
      allocate (swxS, source=ffx)

      allocate (ffxpS, source=ffx)
      allocate (sfxpS, source=ffx)
      allocate (fsxpS, source=ffx)
      allocate (fwxpS, source=ffx)
      allocate (ssxpS, source=ffx)
      allocate (swxpS, source=ffx)

      allocate (cfx6(nxm))
      cfx6 = zero
      allocate (ccx6, source=cfx6)
      allocate (cbx6, source=cfx6)
      allocate (cfxp6, source=cfx6)
      allocate (ciwxp6, source=cfx6)
      allocate (csxp6, source=cfx6)
      allocate (cwxp6, source=cfx6)
      allocate (csx6, source=cfx6)
      allocate (cwx6, source=cfx6)
      allocate (cifx6, source=cfx6)
      allocate (cicx6, source=cfx6)
      allocate (cisx6, source=cfx6)

      allocate (cibx6, source=cfx6)
      allocate (cifxp6, source=cfx6)
      allocate (cisxp6, source=cfx6)
      allocate (ciwx6, source=cfx6)

      allocate (cfi6, source=ffx)
      allocate (cci6, source=ffx)
      allocate (cbi6, source=ffx)
      allocate (cfip6, source=ffx)
      allocate (csip6, source=ffx)
      allocate (cwip6, source=ffx)
      allocate (csi6, source=ffx)
      allocate (cwi6, source=ffx)
      allocate (cifi6, source=ffx)
      allocate (cici6, source=ffx)
      allocate (cibi6, source=ffx)
      allocate (cifip6, source=ffx)

      allocate (cisip6, source=ffx)
      allocate (ciwip6, source=ffx)
      allocate (cisi6, source=ffx)
      allocate (ciwi6, source=ffx)

   end subroutine x3d_operator_x_data_init

   !
   ! Free memory
   !
   subroutine x3d_operator_x_data_finalize()

      implicit none

      deallocate (ffx)
      deallocate (sfx)
      deallocate (fsx)
      deallocate (fwx)
      deallocate (ssx)
      deallocate (swx)

      deallocate (ffxp)
      deallocate (sfxp)
      deallocate (fsxp)
      deallocate (fwxp)
      deallocate (ssxp)
      deallocate (swxp)

      deallocate (ffxS)
      deallocate (sfxS)
      deallocate (fsxS)
      deallocate (fwxS)
      deallocate (ssxS)
      deallocate (swxS)

      deallocate (ffxpS)
      deallocate (sfxpS)
      deallocate (fsxpS)
      deallocate (fwxpS)
      deallocate (ssxpS)
      deallocate (swxpS)

      deallocate (cfx6)
      deallocate (ccx6)
      deallocate (cbx6)
      deallocate (cfxp6)
      deallocate (ciwxp6)
      deallocate (csxp6)
      deallocate (cwxp6)
      deallocate (csx6)
      deallocate (cwx6)
      deallocate (cifx6)
      deallocate (cicx6)
      deallocate (cisx6)

      deallocate (cibx6)
      deallocate (cifxp6)
      deallocate (cisxp6)
      deallocate (ciwx6)

      deallocate (cfi6)
      deallocate (cci6)
      deallocate (cbi6)
      deallocate (cfip6)
      deallocate (csip6)
      deallocate (cwip6)
      deallocate (csi6)
      deallocate (cwi6)
      deallocate (cifi6)
      deallocate (cici6)
      deallocate (cibi6)
      deallocate (cifip6)

      deallocate (cisip6)
      deallocate (ciwip6)
      deallocate (cisi6)
      deallocate (ciwi6)

   end subroutine x3d_operator_x_data_finalize

end module x3d_operator_x_data
