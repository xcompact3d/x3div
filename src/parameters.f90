!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module parameters

  use decomp_2d_constants, only : mytype
  use decomp_2d_mpi, only : decomp_2d_abort
  use decomp_2d_mpi, only : nrank, nproc
  use param
  use variables

  private :: parameter_defaults

contains

!###########################################################################
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!
!###########################################################################
subroutine parameter()

  use x3d_precision, only : pi

  implicit none

  integer :: longueur ,impi,j, is, total
#ifdef DEBUG
  if (nrank == 0) write(*,*) '# parameter start'
#endif

  call parameter_defaults()

  if (nz==1) then
     nclz1 = 0
     nclzn = 0
     p_row = nproc
     p_col = 1
  endif

  !! Set Scalar BCs same as fluid (may be overridden) [DEFAULT]
  nclxS1 = nclx1; nclxSn = nclxn
  nclyS1 = ncly1; nclySn = nclyn
  nclzS1 = nclz1; nclzSn = nclzn
  nu0nu=four
  cnu=0.44_mytype

  if (nclx1.eq.0.and.nclxn.eq.0) then
     nclx=.true.
     nxm=nx
  else
     nclx=.false.
     nxm=nx-1
  endif
  if (ncly1.eq.0.and.nclyn.eq.0) then
     ncly=.true.
     nym=ny
  else
     ncly=.false.
     nym=ny-1
  endif
  if (nclz1.eq.0.and.nclzn.eq.0) then
     nclz=.true.
     nzm=nz
  else
     nclz=.false.
     nzm=nz-1
  endif

  dx=xlx/real(nxm,mytype)
  dy=yly/real(nym,mytype)
  dz=zlz/real(nzm,mytype)

  dx2 = dx * dx
  dy2 = dy * dy
  dz2 = dz * dz

  if (abs(re) > epsilon(re)) then
     xnu = one / re
  else
     xnu = zero
  endif

  ! Some safety check
  if (itimescheme > 1) call decomp_2d_abort(itimescheme, "itimescheme must be specified as 1-6")

end subroutine parameter

!###########################################################################
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!
!###########################################################################
subroutine parameter_defaults()

  use x3d_precision, only : twopi

  implicit none

  integer :: i

  ifirstder = 4
  isecondder = 4
  ro = 99999999._mytype
  angle = zero
  u1 = 2
  u2 = 1
  init_noise = zero
  inflow_noise = zero
  iin = 0
  itimescheme = 1
  iimplicit = 0
  istret = 0
  ipinter=3
  beta = 0
  iscalar = 0
  cont_phi = 0
  irestart = 0
  itime0 = 0
  t0 = zero
  dt = zero

  xlx = one; yly = one; zlz = one
  nclx1 = 0; nclxn = 0
  ncly1 = 0; nclyn = 0
  nclz1 = 0; nclzn = 0

  npress = 1 !! By default people only need one pressure field
  imodulo2 = 1


end subroutine parameter_defaults

!
! Log / output
!
subroutine listing()

  use MPI
  use iso_fortran_env

  implicit none

  if (nrank==0) then
     print *,'==========================================================='
     print *,'========================x3div=============================='
     print *,'==========================================================='
     print *,'==========================================================='
     print *,'===============mini app of Xcompact3D======================'
     print *,'==========================================================='
     print *,'==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
     write(*,"(' xnu                    : ',F17.8)") xnu
     print *,'==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     print *,'==========================================================='
     write(*,"(' Time step dt           : ',F17.8)") dt
     !
     if (itimescheme.eq.1) then
       write(*,"(' Temporal scheme        : ',A20)") "Forwards Euler"
     else
       print *,'Error: itimescheme must be specified as 1-6'
       stop
     endif
     !
     write(*,*)'==========================================================='
     write(*,"(' ifirst                 : ',I17)") ifirst
     write(*,"(' ilast                  : ',I17)") ilast
     write(*,*)'==========================================================='
     write(*,"(' Lx                     : ',F17.8)") xlx
     write(*,"(' Ly                     : ',F17.8)") yly
     write(*,"(' Lz                     : ',F17.8)") zlz
     write(*,"(' nx                     : ',I17)") nx
     write(*,"(' ny                     : ',I17)") ny
     write(*,"(' nz                     : ',I17)") nz
     write(*,*)'==========================================================='
     write(*,"(' nu0nu                  : ',F17.8)") nu0nu
     write(*,"(' cnu                    : ',F17.8)") cnu
     write(*,*)'==========================================================='
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
     write(*,*)'==========================================================='
     ! Show the compile flags detected and the version of the MPI library
#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     write(*,*)'Numerical precision: Double, saving in single'
#else
     write(*,*)'Numerical precision: Double'
#endif
#else
     write(*,*)'Numerical precision: Single'
#endif
     write(*,*)'Compiled with ', compiler_version()
     write(*,*)'Compiler options : ', compiler_options()
     write(*,'(" Version of the MPI library : ",I0,".",I0)') MPI_VERSION, MPI_SUBVERSION
#ifdef DEBUG
     write(*,*)'Compile flag DEBUG detected'
#endif
#ifdef SHM
     write(*,*)'Compile flag SHM detected'
#endif
#ifdef EVEN
     write(*,*)'Compile flag EVEN detected'
#endif
#ifdef OCC
     write(*,*)'Compile flag OCC detected'
#endif
#ifdef OVERWRITE
     write(*,*)'Compile flag OVERWRITE detected'
#endif
#ifdef HALO_DEBUG
     write(*,*)'Compile flag HALO_DEBUG detected'
#endif
#ifdef SHM_DEBUG
     write(*,*)'Compile flag SHM_DEBUG detected'
#endif
     write(*,*)'==========================================================='

  endif

end subroutine listing

end module parameters
