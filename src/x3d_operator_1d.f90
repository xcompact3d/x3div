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
!##################################################################

module x3d_operator_1d

  use decomp_2d, only: mytype, nrank

  implicit none

  type, public :: x3doperator1d
    ! Size
    integer :: n = 0
    ! Odd or even
    integer :: npaire
    ! Extra-diag coefficient
    real(mytype) :: alfa
    ! Arrays needed by the Thomas solver
    ! See subroutine prepare(b,c,f,s,w,n) in schemes.f90
    real(mytype), dimension(:), pointer :: f=>null(), s=>null(), w=>null()
    ! Array needed by the optimized Thomas solver
    real(mytype), dimension(:), allocatable :: periodic

    contains

      procedure, public :: init
      procedure, public :: finalize

  end type x3doperator1d

  ! First derivative on the velocity grid
  type(x3doperator1d), save, public :: x3d_op_derx, x3d_op_derxp
  type(x3doperator1d), save, public :: x3d_op_dery, x3d_op_deryp
  type(x3doperator1d), save, public :: x3d_op_derz, x3d_op_derzp

  ! First derivative from velocity grid => pressure grid
  type(x3doperator1d), save, public :: x3d_op_derxvp, x3d_op_deryvp, x3d_op_derzvp

  ! First derivative from pressure grid => velocity grid
  type(x3doperator1d), save, public :: x3d_op_derxpv, x3d_op_derypv, x3d_op_derzpv

  ! Interpolation from velocity grid => pressure grid
  type(x3doperator1d), save, public :: x3d_op_intxvp, x3d_op_intyvp, x3d_op_intzvp

  ! Interpolation from pressure grid => velocity grid
  type(x3doperator1d), save, public :: x3d_op_intxpv, x3d_op_intypv, x3d_op_intzpv

  ! Make everything private unless declared public
  private
  public :: x3d_operator_1d_init, x3d_operator_1d_finalize


contains

  !
  ! Associate pointers and pre-compute some arrays for periodic cases
  !
  subroutine x3d_operator_1d_init()

    use variables
    use param, only : nclx, ncly, nclz
    use x3d_operator_x_data
    use x3d_operator_y_data
    use x3d_operator_z_data
    use param, only: dx, dx2, nclx, nclx1, nclxn, &
                     dy, dy2, ncly, ncly1, nclyn, &
                     dz, dz2, nclz, nclz1, nclzn

    implicit none

#ifdef DEBUG
      if (nrank == 0) write (*, *) '# x3d_operator_1d_init start'
#endif

      call first_derivative(alfa1x, af1x, bf1x, cf1x, df1x, alfa2x, af2x, alfanx, afnx, bfnx, &
                            cfnx, dfnx, alfamx, afmx, alfaix, afix, bfix, &
                            ffx, fsx, fwx, ffxp, fsxp, fwxp, dx, nx, nclx1, nclxn)
      call first_derivative(alfa1y, af1y, bf1y, cf1y, df1y, alfa2y, af2y, alfany, afny, bfny, &
                            cfny, dfny, alfamy, afmy, alfajy, afjy, bfjy, &
                            ffy, fsy, fwy, ffyp, fsyp, fwyp, dy, ny, ncly1, nclyn)
      call first_derivative(alfa1z, af1z, bf1z, cf1z, df1z, alfa2z, af2z, alfanz, afnz, bfnz, &
                            cfnz, dfnz, alfamz, afmz, alfakz, afkz, bfkz, &
                            ffz, fsz, fwz, ffzp, fszp, fwzp, dz, nz, nclz1, nclzn)

      call second_derivative(alsa1x, as1x, bs1x, &
                             cs1x, ds1x, alsa2x, as2x, alsanx, asnx, bsnx, csnx, dsnx, alsamx, &
                             asmx, alsa3x, as3x, bs3x, alsatx, astx, bstx, &
                             alsa4x, as4x, bs4x, cs4x, &
                             alsattx, asttx, bsttx, csttx, &
                             alsaix, asix, bsix, csix, dsix, &
                             sfx, ssx, swx, sfxp, ssxp, swxp, dx2, nx, nclx1, nclxn)
      call second_derivative(alsa1y, as1y, bs1y, &
                             cs1y, ds1y, alsa2y, as2y, alsany, asny, bsny, csny, dsny, alsamy, &
                             asmy, alsa3y, as3y, bs3y, alsaty, asty, bsty, &
                             alsa4y, as4y, bs4y, cs4y, &
                             alsatty, astty, bstty, cstty, &
                             alsajy, asjy, bsjy, csjy, dsjy, &
                             sfy, ssy, swy, sfyp, ssyp, swyp, dy2, ny, ncly1, nclyn)
      call second_derivative(alsa1z, as1z, bs1z, &
                             cs1z, ds1z, alsa2z, as2z, alsanz, asnz, bsnz, csnz, dsnz, alsamz, &
                             asmz, alsa3z, as3z, bs3z, alsatz, astz, bstz, &
                             alsa4z, as4z, bs4z, cs4z, &
                             alsattz, asttz, bsttz, csttz, &
                             alsakz, askz, bskz, cskz, dskz, &
                             sfz, ssz, swz, sfzp, sszp, swzp, dz2, nz, nclz1, nclzn)

      call interpolation(dx, nxm, nx, nclx1, nclxn, &
                         alcaix6, acix6, bcix6, &
                         ailcaix6, aicix6, bicix6, cicix6, dicix6, &
                         cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                         cwxp6, csx6, cwx6, cifx6, cicx6, cisx6, &
                         cibx6, cifxp6, cisxp6, ciwx6, &
                         cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                         cwi6, cifi6, cici6, cibi6, cifip6, &
                         cisip6, ciwip6, cisi6, ciwi6)
      call interpolation(dy, nym, ny, ncly1, nclyn, &
                         alcaiy6, aciy6, bciy6, &
                         ailcaiy6, aiciy6, biciy6, ciciy6, diciy6, &
                         cfy6, ccy6, cby6, cfyp6, ciwyp6, csyp6, &
                         cwyp6, csy6, cwy6, cify6, cicy6, cisy6, &
                         ciby6, cifyp6, cisyp6, ciwy6, &
                         cfi6y, cci6y, cbi6y, cfip6y, csip6y, cwip6y, csi6y, &
                         cwi6y, cifi6y, cici6y, cibi6y, cifip6y, &
                         cisip6y, ciwip6y, cisi6y, ciwi6y)
      call interpolation(dz, nzm, nz, nclz1, nclzn, &
                         alcaiz6, aciz6, bciz6, &
                         ailcaiz6, aiciz6, biciz6, ciciz6, diciz6, &
                         cfz6, ccz6, cbz6, cfzp6, ciwzp6, cszp6, &
                         cwzp6, csz6, cwz6, cifz6, cicz6, cisz6, &
                         cibz6, cifzp6, ciszp6, ciwz6, &
                         cfi6z, cci6z, cbi6z, cfip6z, csip6z, cwip6z, csi6z, &
                         cwi6z, cifi6z, cici6z, cibi6z, cifip6z, &
                         cisip6z, ciwip6z, cisi6z, ciwi6z)

    ! derx operators for the velocity
    call init(x3d_op_derx,  ffx,  fsx,  fwx,  nx, nclx, alfaix, 0)
    call init(x3d_op_derxp, ffxp, fsxp, fwxp, nx, nclx, alfaix, 1)

    ! dery operators for the velocity
    call init(x3d_op_dery,  ffy,  fsy,  fwy,  ny, ncly, alfajy, 0)
    call init(x3d_op_deryp, ffyp, fsyp, fwyp, ny, ncly, alfajy, 1)

    ! derz operators for the velocity
    call init(x3d_op_derz,  ffz,  fsz,  fwz,  nz, nclz, alfakz, 0)
    call init(x3d_op_derzp, ffzp, fszp, fwzp, nz, nclz, alfakz, 1)

    ! Staggered derivative velocity => pressure
    call init(x3d_op_derxvp, cfx6, csx6, cwx6, nxm, nclx, alcaix6, 1)
    call init(x3d_op_deryvp, cfy6, csy6, cwy6, nym, ncly, alcaiy6, 1)
    call init(x3d_op_derzvp, cfz6, csz6, cwz6, nzm, nclz, alcaiz6, 1)

    ! Staggered derivative pressure => velocity
    if (nclx) then
      call init(x3d_op_derxpv,  cfx6,  csx6,  cwx6, nx, nclx, alcaix6, 0)
    else
      call init(x3d_op_derxpv, cfip6, csip6, cwip6, nx, nclx, alcaix6, 0)
    endif
    if (ncly) then
      call init(x3d_op_derypv,   cfy6,   csy6,   cwy6, ny, ncly, alcaiy6, 0)
    else
      call init(x3d_op_derypv, cfip6y, csip6y, cwip6y, ny, ncly, alcaiy6, 0)
    endif
    if (nclz) then
      call init(x3d_op_derzpv,   cfz6,   csz6,   cwz6, nz, nclz, alcaiz6, 0)
    else
      call init(x3d_op_derzpv, cfip6z, csip6z, cwip6z, nz, nclz, alcaiz6, 0)
    endif

    ! Interpolation velocity => pressure
    call init(x3d_op_intxvp, cifxp6, cisxp6, ciwxp6, nxm, nclx, ailcaix6, 1)
    call init(x3d_op_intyvp, cifyp6, cisyp6, ciwyp6, nym, ncly, ailcaiy6, 1)
    call init(x3d_op_intzvp, cifzp6, ciszp6, ciwzp6, nzm, nclz, ailcaiz6, 1)

    ! Interpolation pressure => velocity
    if (nclx) then
      call init(x3d_op_intxpv,  cifx6,  cisx6,  ciwx6, nx, nclx, ailcaix6, 1)
    else
      call init(x3d_op_intxpv, cifip6, cisip6, ciwip6, nx, nclx, ailcaix6, 1)
    endif
    if (ncly) then
      call init(x3d_op_intypv,   cify6,   cisy6,   ciwy6, ny, ncly, ailcaiy6, 1)
    else
      call init(x3d_op_intypv, cifip6y, cisip6y, ciwip6y, ny, ncly, ailcaiy6, 1)
    endif
    if (nclz) then
      call init(x3d_op_intzpv,   cifz6,   cisz6,   ciwz6, nz, nclz, ailcaiz6, 1)
    else
      call init(x3d_op_intzpv, cifip6z, cisip6z, ciwip6z, nz, nclz, ailcaiz6, 1)
    endif

  end subroutine x3d_operator_1d_init

  !
  ! Nullify pointers and free allocated memory
  !
  subroutine x3d_operator_1d_finalize()

    implicit none

    call finalize(x3d_op_derx)
    call finalize(x3d_op_derxp)

    call finalize(x3d_op_dery)
    call finalize(x3d_op_deryp)

    call finalize(x3d_op_derz)
    call finalize(x3d_op_derzp)

    call finalize(x3d_op_derxvp)
    call finalize(x3d_op_deryvp)
    call finalize(x3d_op_derzvp)

    call finalize(x3d_op_derxpv)
    call finalize(x3d_op_derypv)
    call finalize(x3d_op_derzpv)

    call finalize(x3d_op_intxvp)
    call finalize(x3d_op_intyvp)
    call finalize(x3d_op_intzvp)

    call finalize(x3d_op_intxpv)
    call finalize(x3d_op_intypv)
    call finalize(x3d_op_intzpv)

  end subroutine x3d_operator_1d_finalize

  !
  ! Associate pointers with the given targets
  !
  subroutine init(x3dop, f, s, w, n, ncl, alfa, paire)

    use param, only : zero, one
    use thomas, only : thomas1d

    implicit none

    ! Arguments
    class(x3doperator1d) :: x3dop
    real(mytype), dimension(:), target, intent(in) :: f, s, w
    integer, intent(in) :: n, paire
    logical, intent(in) :: ncl
    real(mytype), intent(in) :: alfa

    ! Local variable
    integer :: i

    ! Nothing to do when n=1 (nz=1 for instance)
    if (n==1) return

    x3dop%n = n
    x3dop%npaire = paire
    x3dop%f => f
    x3dop%s => s
    x3dop%w => w
    x3dop%alfa = alfa
    if (ncl) then
      allocate(x3dop%periodic(n))
      x3dop%periodic = (/-one, (zero, i=2, n-1), alfa/)
      call thomas1d(x3dop%periodic, f, s, w, n)
    endif

  end subroutine init

  !
  ! Nullify pointer and free allocated memory of the given operator
  !
  subroutine finalize(x3dop)

    implicit none

    class(x3doperator1d) :: x3dop

    x3dop%n = 0
    x3dop%npaire = 0
    nullify(x3dop%f)
    nullify(x3dop%s)
    nullify(x3dop%w)
    if (allocated(x3dop%periodic)) deallocate(x3dop%periodic)

  end subroutine finalize

   !
   ! Prepare Thomas algorithm for a tri-diagonal matrix M
   ! See https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Method
   ! Here :
   !     ( c(1) f(1) 0    0    0 ...
   ! M = ( b(1) c(2) f(2) 0    0 ...
   !     ( 0    b(2) c(3) f(3) 0 ...
   !
   ! b, in, lower diagonal
   ! c, in, diagonal coefficient
   ! f, in, upper diagonal
   ! s, out, vector for the forward step of the Thomas algo.
   ! w, out, vector for the backward step of the Thomas algo.
   ! n, in, size of the problem
   !
   subroutine prepare(b, c, f, s, w, n)

      use decomp_2d, only: mytype
      use param, only: one

      implicit none

      integer, intent(in) :: n
      real(mytype), dimension(n), intent(in) :: b, c, f
      real(mytype), dimension(n), intent(out) :: s, w
      integer :: i

      do i = 1, n
         w(i) = c(i)
      end do
      do i = 2, n
         s(i) = b(i - 1)/w(i - 1)
         w(i) = w(i) - f(i - 1)*s(i)
      end do
      do i = 1, n
         w(i) = one/w(i)
      end do

   end subroutine prepare

   subroutine first_derivative(alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                               cfn, dfn, alfam, afm, alfai, afi, bfi, &
                               ff, fs, fw, ffp, fsp, fwp, d, n, ncl1, ncln)

      use decomp_2d, only: mytype, decomp_2d_abort
      use param

      implicit none

      real(mytype), intent(in) :: d
      integer, intent(in) :: n, ncl1, ncln
      real(mytype), dimension(n), intent(out) :: ff, fs, fw, ffp, fsp, fwp
      real(mytype), intent(out) :: alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                                   cfn, dfn, alfam, afm, alfai, afi, bfi
      integer :: i
      real(mytype), dimension(n) :: fb, fc

      if (n==1) return

      ff = zero; fs = zero; fw = zero; ffp = zero; fsp = zero; fwp = zero
      fb = zero; fc = zero

      if (ifirstder == 1) then    ! Second-order central
         alfai = zero
         afi = one/(two*d)
         bfi = zero
      elseif (ifirstder == 2) then ! Fourth-order central
         call decomp_2d_abort(__FILE__, __LINE__, ifirstder, "Set of coefficients not ready yet")
      elseif (ifirstder == 3) then ! Fourth-order compact
         call decomp_2d_abort(__FILE__, __LINE__, ifirstder, "Set of coefficients not ready yet")
      elseif (ifirstder == 4) then ! Sixth-order compact
         alfai = one/three
         afi = (seven/nine)/d
         bfi = (one/36._mytype)/d
      else
         call decomp_2d_abort(__FILE__, __LINE__, ifirstder, "This is not an option. Please use ifirstder=1,2,3,4")
      end if

      if (ifirstder == 1) then
         alfa1 = zero
         af1 = zero
         bf1 = zero
         cf1 = zero
         df1 = zero
         alfa2 = zero
         af2 = zero

         alfam = zero
         afm = zero
         alfan = zero
         afn = zero
         bfn = zero
         cfn = zero
         dfn = zero
      else
         alfa1 = two
         af1 = -(five/two)/d
         bf1 = (two)/d
         cf1 = (half)/d
         df1 = zero

         alfa2 = one/four
         af2 = (three/four)/d
         alfan = two
         afn = -(five/two)/d
         bfn = (two)/d
         cfn = (half)/d
         dfn = zero
         alfam = one/four
         afm = (three/four)/d
      end if

      if (ncl1 .eq. 0) then !Periodic
         ff(1) = alfai
         ff(2) = alfai
         fc(1) = two
         fc(2) = one
         fb(1) = alfai
         fb(2) = alfai
      elseif (ncl1 .eq. 1) then !Free-slip
         ff(1) = alfai + alfai
         ff(2) = alfai
         fc(1) = one
         fc(2) = one
         fb(1) = alfai
         fb(2) = alfai
      elseif (ncl1 .eq. 2) then !Dirichlet
         ff(1) = alfa1
         ff(2) = alfa2
         fc(1) = one
         fc(2) = one
         fb(1) = alfa2
         fb(2) = alfai
      end if
      if (ncln .eq. 0) then !Periodic
         ff(n - 2) = alfai
         ff(n - 1) = alfai
         ff(n) = zero
         fc(n - 2) = one
         fc(n - 1) = one
         fc(n) = one + alfai*alfai
         fb(n - 2) = alfai
         fb(n - 1) = alfai
         fb(n) = zero
      elseif (ncln .eq. 1) then !Free-slip
         ff(n - 2) = alfai
         ff(n - 1) = alfai
         ff(n) = zero
         fc(n - 2) = one
         fc(n - 1) = one
         fc(n) = one
         fb(n - 2) = alfai
         fb(n - 1) = alfai + alfai
         fb(n) = zero
      elseif (ncln .eq. 2) then !Dirichlet
         ff(n - 2) = alfai
         ff(n - 1) = alfam
         ff(n) = zero
         fc(n - 2) = one
         fc(n - 1) = one
         fc(n) = one
         fb(n - 2) = alfam
         fb(n - 1) = alfan
         fb(n) = zero
      end if
      do i = 3, n - 3
         ff(i) = alfai
         fc(i) = one
         fb(i) = alfai
      end do

      do i = 1, n
         ffp(i) = ff(i)
      end do

      call prepare(fb, fc, ff, fs, fw, n)

      if (ncl1 .eq. 1) then
         ffp(1) = zero
      end if
      if (ncln .eq. 1) then
         fb(n - 1) = zero
      end if

      call prepare(fb, fc, ffp, fsp, fwp, n)

   end subroutine first_derivative

   subroutine second_derivative(alsa1, as1, bs1, &
                                cs1, ds1, alsa2, as2, alsan, asn, bsn, csn, dsn, alsam, &
                                asm, alsa3, as3, bs3, alsat, ast, bst, &
                                alsa4, as4, bs4, cs4, &
                                alsatt, astt, bstt, cstt, &
                                alsai, asi, bsi, csi, dsi, &
                                sf, ss, sw, sfp, ssp, swp, d2, n, ncl1, ncln)

      use decomp_2d, only: mytype, nrank, decomp_2d_abort
      use x3dprecision, only: pi, twopi
      use param
      use variables, only: nu0nu, cnu

      implicit none

      real(mytype), intent(in) :: d2
      integer, intent(in) :: n, ncl1, ncln
      real(mytype), dimension(n), intent(out) :: sf, ss, sw, sfp, ssp, swp
      real(mytype), intent(out) :: alsa1, as1, bs1, &
                                   cs1, ds1, alsa2, as2, alsan, asn, bsn, csn, dsn, alsam, &
                                   asm, alsa3, as3, bs3, alsat, ast, bst, &
                                   alsa4, as4, bs4, cs4, &
                                   alsatt, astt, bstt, cstt, &
                                   alsai, asi, bsi, csi, dsi
      integer :: i
      real(mytype), dimension(n) :: sb, sc
      real(mytype) :: dpis3, kppkc, kppkm, xnpi2, xmpi2, den

      if (n==1) return

      sf = zero; ss = zero; sw = zero; sfp = zero; ssp = zero; swp = zero

      ! Define coefficients based on the desired formal accuracy of the numerical schemes
      if (isecondder == 1) then    ! Second-order central
         alsai = zero
         asi = one/d2  !((six-nine*alsai)/four)/d2
         bsi = zero !((-three+twentyfour*alsai)/five)/(four*d2)
         csi = zero !((two-eleven*alsai)/twenty)/(nine*d2)
         dsi = zero

         alsa4 = alsai
         as4 = asi
         bs4 = bsi
         cs4 = csi

         alsatt = alsai
         astt = asi
         bstt = bsi
         cstt = csi
      elseif (isecondder == 2) then ! Fourth-order central
         call decomp_2d_abort(__FILE__, __LINE__, isecondder, "Set of coefficients not ready yet")
      elseif (isecondder == 3) then ! Fourth-order compact
         call decomp_2d_abort(__FILE__, __LINE__, isecondder, "Set of coefficients not ready yet")
      elseif (isecondder == 4) then ! Sixth-order compact Lele style (no extra dissipation)
         alsai = two/11._mytype
         asi = (12._mytype/11._mytype)/d2
         bsi = (three/44.0_mytype)/d2
         csi = zero
         dsi = zero

         alsa4 = alsai
         as4 = asi
         bs4 = bsi
         cs4 = csi

         alsatt = alsai
         astt = asi
         bstt = bsi
         cstt = csi
      elseif (isecondder == 5) then ! Sixth-order Hyperviscous operator
         if (nrank == 0) write (*, *) 'Using the hyperviscous operator with (nu_0/nu,c_nu) = ', '(', nu0nu, ',', cnu, ')'
         dpis3 = twopi/three
         kppkc = pi*pi*(one + nu0nu)
         kppkm = dpis3*dpis3*(one + cnu*nu0nu) !exp(-((pi-dpis3)/(zpthree*pi-dpis3))**two)/xxnu+dpis3*dpis3
         xnpi2 = kppkc
         xmpi2 = kppkm

         den = 405._mytype*xnpi2 - 640._mytype*xmpi2 + 144._mytype

         alsai = half - (320._mytype*xmpi2 - 1296._mytype)/den
         asi = -(4329._mytype*xnpi2/eight - 32._mytype*xmpi2 - 140._mytype*xnpi2*xmpi2 + 286._mytype)/den/d2
         bsi = (2115._mytype*xnpi2 - 1792._mytype*xmpi2 - 280._mytype*xnpi2*xmpi2 + 1328._mytype)/den/(four*d2)
         csi = -(7695._mytype*xnpi2/eight + 288._mytype*xmpi2 - 180._mytype*xnpi2*xmpi2 - 2574._mytype)/den/(nine*d2)
         dsi = (198._mytype*xnpi2 + 128._mytype*xmpi2 - 40._mytype*xnpi2*xmpi2 - 736._mytype)/den/(four**2*d2)
      else
         if (nrank == 0) then
            write (*, *) 'This is not an option.'
         end if
      end if

      ! Defined for the bounadies when dirichlet conditions are used
      alsa1 = 11._mytype
      as1 = (13._mytype)/d2
      bs1 = -(27._mytype)/d2
      cs1 = (15._mytype)/d2
      ds1 = -(one)/d2

      if (isecondder == 1) then
         alsa2 = zero
         as2 = one/d2
      else
         alsa2 = 0.1_mytype
         as2 = (six/five)/d2
      end if

      alsa3 = two/11._mytype
      as3 = (12._mytype/11._mytype)/d2
      bs3 = (three/44._mytype)/d2

      alsa4 = two/11._mytype
      as4 = (12._mytype/11._mytype)/d2
      bs4 = (three/44._mytype)/d2
      cs4 = zero

      alsan = 11._mytype
      asn = (13._mytype)/d2
      bsn = -(27._mytype)/d2
      csn = (15._mytype)/d2
      dsn = -(one)/d2

      if (isecondder == 1) then
         alsam = zero
         asm = one/d2
      else
         alsam = 0.1_mytype
         asm = (six/five)/d2
      end if

      alsat = two/11._mytype
      ast = (12._mytype/11._mytype)/d2
      bst = (three/44._mytype)/d2

      alsatt = two/11._mytype
      astt = (12._mytype/11._mytype)/d2
      bstt = (three/44._mytype)/d2
      cstt = zero

      if (ncl1 .eq. 0) then !Periodic
         sf(1) = alsai
         sf(2) = alsai
         sf(3) = alsai
         sf(4) = alsai
         sc(1) = two
         sc(2) = one
         sc(3) = one
         sc(4) = one
         sb(1) = alsai
         sb(2) = alsai
         sb(3) = alsai
         sb(4) = alsai
      elseif (ncl1 .eq. 1) then !Free-slip
         sf(1) = alsai + alsai
         sf(2) = alsai
         sf(3) = alsai
         sf(4) = alsai
         sc(1) = one
         sc(2) = one
         sc(3) = one
         sc(4) = one
         sb(1) = alsai
         sb(2) = alsai
         sb(3) = alsai
         sb(4) = alsai
      elseif (ncl1 .eq. 2) then !Dirichlet
         sf(1) = alsa1
         sf(2) = alsa2
         sf(3) = alsa3
         sf(4) = alsa4
         sc(1) = one
         sc(2) = one
         sc(3) = one
         sc(4) = one
         sb(1) = alsa2
         sb(2) = alsa3
         sb(3) = alsa4
         sb(4) = alsai
      end if
      if (ncln .eq. 0) then !Periodic
         sf(n - 4) = alsai
         sf(n - 3) = alsai
         sf(n - 2) = alsai
         sf(n - 1) = alsai
         sf(n) = zero
         sc(n - 4) = one
         sc(n - 3) = one
         sc(n - 2) = one
         sc(n - 1) = one
         sc(n) = one + alsai*alsai
         sb(n - 4) = alsai
         sb(n - 3) = alsai
         sb(n - 2) = alsai
         sb(n - 1) = alsai
         sb(n) = zero
      elseif (ncln .eq. 1) then !Free-slip
         sf(n - 4) = alsai
         sf(n - 3) = alsai
         sf(n - 2) = alsai
         sf(n - 1) = alsai
         sf(n) = zero
         sc(n - 4) = one
         sc(n - 3) = one
         sc(n - 2) = one
         sc(n - 1) = one
         sc(n) = one
         sb(n - 4) = alsai
         sb(n - 3) = alsai
         sb(n - 2) = alsai
         sb(n - 1) = alsai + alsai
         sb(n) = zero
      elseif (ncln .eq. 2) then !Dirichlet
         sf(n - 4) = alsai
         sf(n - 3) = alsatt
         sf(n - 2) = alsat
         sf(n - 1) = alsam
         sf(n) = zero
         sc(n - 4) = one
         sc(n - 3) = one
         sc(n - 2) = one
         sc(n - 1) = one
         sc(n) = one
         sb(n - 4) = alsatt
         sb(n - 3) = alsat
         sb(n - 2) = alsam
         sb(n - 1) = alsan
         sb(n) = zero
      end if
      do i = 5, n - 5
         sf(i) = alsai
         sc(i) = one
         sb(i) = alsai
      end do

      do i = 1, n
         sfp(i) = sf(i)
      end do

      if (ncl1 .eq. 1) then
         sf(1) = zero
      end if

      call prepare(sb, sc, sf, ss, sw, n)
      call prepare(sb, sc, sfp, ssp, swp, n)

      if (ncln .eq. 1) then
         sb(n - 1) = zero
         call prepare(sb, sc, sf, ss, sw, n)
      end if

   end subroutine second_derivative

   subroutine interpolation(dx, nxm, nx, nclx1, nclxn, &
                            alcaix6, acix6, bcix6, &
                            ailcaix6, aicix6, bicix6, cicix6, dicix6, &
                            cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                            cwxp6, csx6, cwx6, cifx6, cicx6, cisx6, &
                            cibx6, cifxp6, cisxp6, ciwx6, &
                            cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                            cwi6, cifi6, cici6, cibi6, cifip6, &
                            cisip6, ciwip6, cisi6, ciwi6)

      use decomp_2d, only: mytype
      use param, only: zero, half, one, two, three, four, nine, ten, ipinter, ifirstder

      implicit none

      real(mytype), intent(in) :: dx
      integer, intent(in) :: nxm, nx, nclx1, nclxn
      real(mytype) :: alcaix6, acix6, bcix6
      real(mytype) :: ailcaix6, aicix6, bicix6, cicix6, dicix6
      real(mytype), dimension(nxm) :: cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                                      cwxp6, csx6, cwx6, cifx6, cicx6, cisx6
      real(mytype), dimension(nxm) :: cibx6, cifxp6, cisxp6, ciwx6
      real(mytype), dimension(nx) :: cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                                     cwi6, cifi6, cici6, cibi6, cifip6
      real(mytype), dimension(nx) :: cisip6, ciwip6, cisi6, ciwi6

      integer :: i

      if (nx==1) return

      if (ifirstder == 1) then
         alcaix6 = zero
         acix6 = one/dx
         bcix6 = zero
      else
         alcaix6 = nine/62._mytype
         acix6 = (63._mytype/62._mytype)/dx
         bcix6 = (17._mytype/62._mytype)/three/dx
      end if

      cfx6(1) = alcaix6
      cfx6(2) = alcaix6
      cfx6(nxm - 2) = alcaix6
      cfx6(nxm - 1) = alcaix6
      cfx6(nxm) = zero
      if (nclx1 == 0) ccx6(1) = two
      if (nclx1 == 1) ccx6(1) = one + alcaix6
      if (nclx1 == 2) ccx6(1) = one + alcaix6
      ccx6(2) = one
      ccx6(nxm - 2) = one
      ccx6(nxm - 1) = one
      if (nclxn == 0) ccx6(nxm) = one + alcaix6*alcaix6
      if (nclxn == 1) ccx6(nxm) = one + alcaix6
      if (nclxn == 2) ccx6(nxm) = one + alcaix6
      cbx6(1) = alcaix6
      cbx6(2) = alcaix6
      cbx6(nxm - 2) = alcaix6
      cbx6(nxm - 1) = alcaix6
      cbx6(nxm) = zero
      do i = 3, nxm - 3
         cfx6(i) = alcaix6
         ccx6(i) = one
         cbx6(i) = alcaix6
      end do

      cfi6(1) = alcaix6 + alcaix6
      cfi6(2) = alcaix6
      cfi6(nx - 2) = alcaix6
      cfi6(nx - 1) = alcaix6
      cfi6(nx) = zero
      cci6(1) = one
      cci6(2) = one
      cci6(nx - 2) = one
      cci6(nx - 1) = one
      cci6(nx) = one
      cbi6(1) = alcaix6
      cbi6(2) = alcaix6
      cbi6(nx - 2) = alcaix6
      cbi6(nx - 1) = alcaix6 + alcaix6
      cbi6(nx) = zero
      do i = 3, nx - 3
         cfi6(i) = alcaix6
         cci6(i) = one
         cbi6(i) = alcaix6
      end do

      if (ifirstder == 1) then
         ailcaix6 = zero
         aicix6 = half
         bicix6 = zero
         cicix6 = zero
         dicix6 = zero
      else if (ipinter .eq. 1) then
         ailcaix6 = three/ten
         aicix6 = three/four
         bicix6 = one/(two*ten)
         cicix6 = zero
         dicix6 = zero
      else if (ipinter .eq. 2) then
         ailcaix6 = 0.461658_mytype

         dicix6 = 0.00293016_mytype
         aicix6 = one/64._mytype*(75._mytype + 70._mytype*ailcaix6 - 320._mytype*dicix6)
         bicix6 = one/128._mytype*(126._mytype*ailcaix6 - 25._mytype + 1152._mytype*dicix6)
         cicix6 = one/128._mytype*(-ten*ailcaix6 + three - 640._mytype*dicix6)

         aicix6 = aicix6/two
         bicix6 = bicix6/two
         cicix6 = cicix6/two
         dicix6 = dicix6/two
      else if (ipinter .eq. 3) then
         ailcaix6 = 0.49_mytype
         aicix6 = one/128._mytype*(75._mytype + 70._mytype*ailcaix6)
         bicix6 = one/256._mytype*(126._mytype*ailcaix6 - 25._mytype)
         cicix6 = one/256._mytype*(-ten*ailcaix6 + three)
         dicix6 = zero
      end if

      cifx6(1) = ailcaix6
      cifx6(2) = ailcaix6
      cifx6(nxm - 2) = ailcaix6
      cifx6(nxm - 1) = ailcaix6
      cifx6(nxm) = zero
      if (nclx1 == 0) cicx6(1) = two
      if (nclx1 == 1) cicx6(1) = one + ailcaix6
      if (nclx1 == 2) cicx6(1) = one + ailcaix6
      cicx6(2) = one
      cicx6(nxm - 2) = one
      cicx6(nxm - 1) = one
      if (nclxn == 0) cicx6(nxm) = one + ailcaix6*ailcaix6
      if (nclxn == 1) cicx6(nxm) = one + ailcaix6
      if (nclxn == 2) cicx6(nxm) = one + ailcaix6
      cibx6(1) = ailcaix6
      cibx6(2) = ailcaix6
      cibx6(nxm - 2) = ailcaix6
      cibx6(nxm - 1) = ailcaix6
      cibx6(nxm) = zero
      do i = 3, nxm - 3
         cifx6(i) = ailcaix6
         cicx6(i) = one
         cibx6(i) = ailcaix6
      end do
      cifi6(1) = ailcaix6 + ailcaix6
      cifi6(2) = ailcaix6
      cifi6(nx - 2) = ailcaix6
      cifi6(nx - 1) = ailcaix6
      cifi6(nx) = zero
      cici6(1) = one
      cici6(2) = one
      cici6(nx - 2) = one
      cici6(nx - 1) = one
      cici6(nx) = one
      cibi6(1) = ailcaix6
      cibi6(2) = ailcaix6
      cibi6(nx - 2) = ailcaix6
      cibi6(nx - 1) = ailcaix6 + ailcaix6
      cibi6(nx) = zero
      do i = 3, nx - 3
         cifi6(i) = ailcaix6
         cici6(i) = one
         cibi6(i) = ailcaix6
      end do

      do i = 1, nxm
         cfxp6(i) = cfx6(i)
         cifxp6(i) = cifx6(i)
      end do
      do i = 1, nx
         cifip6(i) = cifi6(i)
         cfip6(i) = cfi6(i)
      end do
      cfxp6(1) = zero
      cfip6(1) = zero
      call prepare(cbx6, ccx6, cfx6, csx6, cwx6, nxm)
      call prepare(cbx6, ccx6, cfxp6, csxp6, cwxp6, nxm)
      call prepare(cibx6, cicx6, cifx6, cisx6, ciwx6, nxm)
      call prepare(cibx6, cicx6, cifxp6, cisxp6, ciwxp6, nxm)
      call prepare(cbi6, cci6, cfi6, csi6, cwi6, nx)
      call prepare(cbi6, cci6, cfip6, csip6, cwip6, nx)
      call prepare(cibi6, cici6, cifi6, cisi6, ciwi6, nx)
      call prepare(cibi6, cici6, cifip6, cisip6, ciwip6, nx)
      if (nclxn .eq. 1) then
         cbx6(nxm - 1) = zero
         cibx6(nxm) = zero
         cbi6(nx - 1) = zero
         cibi6(nx) = zero
         call prepare(cbx6, ccx6, cfxp6, csxp6, cwxp6, nxm)
         call prepare(cibx6, cicx6, cifxp6, cisxp6, ciwxp6, nxm)
         call prepare(cbi6, cci6, cfip6, csip6, cwip6, nx)
         call prepare(cibi6, cici6, cifip6, cisip6, ciwip6, nx)
      end if
      if (nclxn .eq. 2) then
         cbx6(nxm - 1) = zero
         cibx6(nxm) = zero
         cbi6(nx - 1) = zero
         cibi6(nx) = zero
         call prepare(cbx6, ccx6, cfxp6, csxp6, cwxp6, nxm)
         call prepare(cibx6, cicx6, cifxp6, cisxp6, ciwxp6, nxm)
         call prepare(cbi6, cci6, cfip6, csip6, cwip6, nx)
         call prepare(cibi6, cici6, cifip6, cisip6, ciwip6, nx)
      end if

   end subroutine interpolation

end module x3d_operator_1d
