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

  use decomp_2d, only: mytype

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

    implicit none

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
    use thomas, only : thomas_optim, thomas1d

    implicit none

    ! Arguments
    class(x3doperator1d) :: x3dop
    real(mytype), dimension(:), target, intent(in) :: f, s, w
    integer, intent(in) :: n, paire
    logical, intent(in) :: ncl
    real(mytype), intent(in) :: alfa

    ! Local variable
    integer :: i

    x3dop%n = n
    x3dop%npaire = paire
    x3dop%f => f
    x3dop%s => s
    x3dop%w => w
    x3dop%alfa = alfa
    if (thomas_optim.and.ncl) then
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

end module x3d_operator_1d
