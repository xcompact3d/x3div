!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module x3d_operator_1d

  use decomp_2d_constants, only: mytype
  use decomp_2d_mpi, only: nrank

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

      ! First derivative
      ! Velocity
      call first_deriv_exp_(alfa1x, af1x, bf1x, cf1x, df1x, alfa2x, af2x, alfanx, afnx, bfnx, &
                            cfnx, dfnx, alfamx, afmx, alfaix, afix, bfix, dx, nx, nclx1, nclxn)
      call first_deriv_imp_(alfa1x, af1x, bf1x, cf1x, df1x, alfa2x, af2x, alfanx, afnx, bfnx, &
                            cfnx, dfnx, alfamx, afmx, alfaix, afix, bfix, &
                            ffx, fsx, fwx, ffxp, fsxp, fwxp, dx, nx, nclx1, nclxn)
      call first_deriv_exp_(alfa1y, af1y, bf1y, cf1y, df1y, alfa2y, af2y, alfany, afny, bfny, &
                            cfny, dfny, alfamy, afmy, alfajy, afjy, bfjy, dy, ny, ncly1, nclyn)
      call first_deriv_imp_(alfa1y, af1y, bf1y, cf1y, df1y, alfa2y, af2y, alfany, afny, bfny, &
                            cfny, dfny, alfamy, afmy, alfajy, afjy, bfjy, &
                            ffy, fsy, fwy, ffyp, fsyp, fwyp, dy, ny, ncly1, nclyn)
      call first_deriv_exp_(alfa1z, af1z, bf1z, cf1z, df1z, alfa2z, af2z, alfanz, afnz, bfnz, &
                            cfnz, dfnz, alfamz, afmz, alfakz, afkz, bfkz, dz, nz, nclz1, nclzn)
      if (nz /= 1) then
         call first_deriv_imp_(alfa1z, af1z, bf1z, cf1z, df1z, alfa2z, af2z, alfanz, afnz, bfnz, &
                               cfnz, dfnz, alfamz, afmz, alfakz, afkz, bfkz, &
                               ffz, fsz, fwz, ffzp, fszp, fwzp, dz, nz, nclz1, nclzn)
      endif

      call interpol_exp_(dx, nxm, nx, nclx1, nclxn, &
                         alcaix6, acix6, bcix6, &
                         ailcaix6, aicix6, bicix6, cicix6, dicix6)
      call interpol_imp_(dx, nxm, nx, nclx1, nclxn, &
                         alcaix6, acix6, bcix6, &
                         ailcaix6, aicix6, bicix6, cicix6, dicix6, &
                         cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                         cwxp6, csx6, cwx6, cifx6, cicx6, cisx6, &
                         cibx6, cifxp6, cisxp6, ciwx6, &
                         cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                         cwi6, cifi6, cici6, cibi6, cifip6, &
                         cisip6, ciwip6, cisi6, ciwi6)
      call interpol_exp_(dy, nym, ny, ncly1, nclyn, &
                         alcaiy6, aciy6, bciy6, &
                         ailcaiy6, aiciy6, biciy6, ciciy6, diciy6)
      call interpol_imp_(dy, nym, ny, ncly1, nclyn, &
                         alcaiy6, aciy6, bciy6, &
                         ailcaiy6, aiciy6, biciy6, ciciy6, diciy6, &
                         cfy6, ccy6, cby6, cfyp6, ciwyp6, csyp6, &
                         cwyp6, csy6, cwy6, cify6, cicy6, cisy6, &
                         ciby6, cifyp6, cisyp6, ciwy6, &
                         cfi6y, cci6y, cbi6y, cfip6y, csip6y, cwip6y, csi6y, &
                         cwi6y, cifi6y, cici6y, cibi6y, cifip6y, &
                         cisip6y, ciwip6y, cisi6y, ciwi6y)
      call interpol_exp_(dz, nzm, nz, nclz1, nclzn, &
                         alcaiz6, aciz6, bciz6, &
                         ailcaiz6, aiciz6, biciz6, ciciz6, diciz6)
      if (nz /= 1) then
         call interpol_imp_(dz, nzm, nz, nclz1, nclzn, &
                            alcaiz6, aciz6, bciz6, &
                            ailcaiz6, aiciz6, biciz6, ciciz6, diciz6, &
                            cfz6, ccz6, cbz6, cfzp6, ciwzp6, cszp6, &
                            cwzp6, csz6, cwz6, cifz6, cicz6, cisz6, &
                            cibz6, cifzp6, ciszp6, ciwz6, &
                            cfi6z, cci6z, cbi6z, cfip6z, csip6z, cwip6z, csi6z, &
                            cwi6z, cifi6z, cici6z, cibi6z, cifip6z, &
                            cisip6z, ciwip6z, cisi6z, ciwi6z)
     endif

    ! derx operators for the velocity and scalars
    call init(x3d_op_derx,  ffx,  fsx,  fwx,  nx, nclx, alfaix, 0)
    call init(x3d_op_derxp, ffxp, fsxp, fwxp, nx, nclx, alfaix, 1)

    ! dery operators for the velocity and scalars
    call init(x3d_op_dery,  ffy,  fsy,  fwy,  ny, ncly, alfajy, 0)
    call init(x3d_op_deryp, ffyp, fsyp, fwyp, ny, ncly, alfajy, 1)

    ! derz operators for the velocity and scalars
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

#ifdef DEBUG
      if (nrank == 0) write (*, *) '# x3d_operator_1d_init done'
#endif

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

   subroutine first_deriv_exp_(alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                               cfn, dfn, alfam, afm, alfai, afi, bfi, &
                               d, n, ncl1, ncln)

      use decomp_2d_mpi, only: mytype, decomp_2d_abort
      use param

      implicit none

      real(mytype), intent(in) :: d
      integer, intent(in) :: n, ncl1, ncln
      real(mytype), intent(out) :: alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                                   cfn, dfn, alfam, afm, alfai, afi, bfi
      integer :: i

      if (ifirstder == 1) then    ! Second-order central
         alfai = zero
         afi = one/(two*d)
         bfi = zero
      elseif (ifirstder == 2) then ! Fourth-order central
         call decomp_2d_abort(__FILE__, __LINE__, ifirstder, &
                 "Set of coefficients not ready yet")
      elseif (ifirstder == 3) then ! Fourth-order compact
         call decomp_2d_abort(__FILE__, __LINE__, ifirstder, &
                 "Set of coefficients not ready yet")
      elseif (ifirstder == 4) then ! Sixth-order compact
         alfai = one/three
         afi = (seven/nine)/d
         bfi = (one/36._mytype)/d
      else
         call decomp_2d_abort(__FILE__, __LINE__, &
                 ifirstder, "This is not an option. Please use ifirstder=1,2,3,4")
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

   end subroutine first_deriv_exp_

   subroutine first_deriv_imp_(alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                               cfn, dfn, alfam, afm, alfai, afi, bfi, &
                               ff, fs, fw, ffp, fsp, fwp, d, n, ncl1, ncln)

      use param

      implicit none

      real(mytype), intent(in) :: d
      integer, intent(in) :: n, ncl1, ncln
      real(mytype), dimension(n), intent(out) :: ff, fs, fw, ffp, fsp, fwp
      real(mytype), intent(in) :: alfa1, af1, bf1, cf1, df1, alfa2, af2, alfan, afn, bfn, &
                                   cfn, dfn, alfam, afm, alfai, afi, bfi
      integer :: i
      real(mytype), dimension(n) :: fb, fc

      ff = zero; fs = zero; fw = zero; ffp = zero; fsp = zero; fwp = zero
      fb = zero; fc = zero

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

   end subroutine first_deriv_imp_

   subroutine interpol_exp_(dx, nxm, nx, nclx1, nclxn, &
                            alcaix6, acix6, bcix6, &
                            ailcaix6, aicix6, bicix6, cicix6, dicix6)

      use param, only: zero, half, one, two, three, four, nine, ten, ipinter, ifirstder

      implicit none

      real(mytype), intent(in) :: dx
      integer, intent(in) :: nxm, nx, nclx1, nclxn
      real(mytype) :: alcaix6, acix6, bcix6
      real(mytype) :: ailcaix6, aicix6, bicix6, cicix6, dicix6

      integer :: i

      if (ifirstder == 1) then
         alcaix6 = zero
         acix6 = one/dx
         bcix6 = zero
      else
         alcaix6 = nine/62._mytype
         acix6 = (63._mytype/62._mytype)/dx
         bcix6 = (17._mytype/62._mytype)/three/dx
      end if

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

   end subroutine interpol_exp_

   subroutine interpol_imp_(dx, nxm, nx, nclx1, nclxn, &
                            alcaix6, acix6, bcix6, &
                            ailcaix6, aicix6, bicix6, cicix6, dicix6, &
                            cfx6, ccx6, cbx6, cfxp6, ciwxp6, csxp6, &
                            cwxp6, csx6, cwx6, cifx6, cicx6, cisx6, &
                            cibx6, cifxp6, cisxp6, ciwx6, &
                            cfi6, cci6, cbi6, cfip6, csip6, cwip6, csi6, &
                            cwi6, cifi6, cici6, cibi6, cifip6, &
                            cisip6, ciwip6, cisi6, ciwi6)

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

   end subroutine interpol_imp_

end module x3d_operator_1d
