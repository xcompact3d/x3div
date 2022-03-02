!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_operator_z_data

   use decomp_2d, only: mytype

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

   real(mytype), target, allocatable, dimension(:) :: ffzS, sfzS, fszS, fwzS, sszS, swzS
   real(mytype), target, allocatable, dimension(:) :: ffzpS, sfzpS, fszpS, fwzpS, sszpS, swzpS

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
   subroutine x3d_operator_z_data_init(nz, nzm)

      implicit none

      integer, intent(in) :: nz, nzm

      if (nz==1) return

      allocate (ffz(nz))
      ffz = 0._mytype
      allocate (sfz, source=ffz)
      allocate (fsz, source=ffz)
      allocate (fwz, source=ffz)
      allocate (ssz, source=ffz)
      allocate (swz, source=ffz)

      allocate (ffzp, source=ffz)
      allocate (sfzp, source=ffz)
      allocate (fszp, source=ffz)
      allocate (fwzp, source=ffz)
      allocate (sszp, source=ffz)
      allocate (swzp, source=ffz)

      allocate (ffzS, source=ffz)
      allocate (sfzS, source=ffz)
      allocate (fszS, source=ffz)
      allocate (fwzS, source=ffz)
      allocate (sszS, source=ffz)
      allocate (swzS, source=ffz)

      allocate (ffzpS, source=ffz)
      allocate (sfzpS, source=ffz)
      allocate (fszpS, source=ffz)
      allocate (fwzpS, source=ffz)
      allocate (sszpS, source=ffz)
      allocate (swzpS, source=ffz)

      allocate (cfz6(nzm))
      cfz6 = 0._mytype
      allocate (ccz6, source=cfz6)
      allocate (cbz6, source=cfz6)
      allocate (cfzp6, source=cfz6)
      allocate (cszp6, source=cfz6)
      allocate (cwzp6, source=cfz6)
      allocate (csz6, source=cfz6)

      allocate (cwz6, source=cfz6)
      allocate (cifz6, source=cfz6)
      allocate (cicz6, source=cfz6)
      allocate (cibz6, source=cfz6)
      allocate (cifzp6, source=cfz6)
      allocate (ciszp6, source=cfz6)
      allocate (ciwzp6, source=cfz6)
      allocate (cisz6, source=cfz6)
      allocate (ciwz6, source=cfz6)

      allocate (cfi6z, source=ffz)
      allocate (cci6z, source=ffz)
      allocate (cbi6z, source=ffz)
      allocate (cfip6z, source=ffz)
      allocate (csip6z, source=ffz)
      allocate (cwip6z, source=ffz)
      allocate (csi6z, source=ffz)
      allocate (cwi6z, source=ffz)
      allocate (cifi6z, source=ffz)
      allocate (cici6z, source=ffz)

      allocate (cibi6z, source=ffz)
      allocate (cifip6z, source=ffz)
      allocate (cisip6z, source=ffz)
      allocate (ciwip6z, source=ffz)
      allocate (cisi6z, source=ffz)
      allocate (ciwi6z, source=ffz)

   end subroutine x3d_operator_z_data_init

   !
   ! Free memory
   !
   subroutine x3d_operator_z_data_finalize()

      implicit none

      if (.not.allocated(ffz)) return

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
