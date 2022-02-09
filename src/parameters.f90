!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

!###########################################################################
!
!  SUBROUTINE: parameter
! DESCRIPTION: Reads the input.i3d file and sets the parameters of the
!              simulation.
!
!###########################################################################
subroutine parameter()

  use iso_fortran_env

  use decomp_2d, only : mytype
  use x3d_precision, only : pi
  use param
  use variables
  use decomp_2d, only : nrank, nproc

  implicit none

  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total
#ifdef DEBUG
  if (nrank == 0) write(*,*) '# parameter start'
#endif

  if (nrank==0) then
     write(*,*) '==========================================================='
     write(*,*) '======================Xcompact3D==========================='
     write(*,*) '===Copyright (c) 2022 Éric Lamballais and Sylvain Laizet==='
     write(*,*) '==========================================================='
#if defined(VERSION)
     write(*,*)'Git version        : ', VERSION
#else
     write(*,*)'Git version        : unknown'
#endif
  endif

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

  ! cfl = 0.2 and u = 1
  dt = 0.2_mytype * dx
  if (abs(re) > epsilon(re)) then
     xnu = one / re
  else
     xnu = one / 1600._mytype
  endif

  anglex = sin(pi*angle/180._mytype)
  angley = cos(pi*angle/180._mytype)
  !###########################################################################
  ! Log-output
  !###########################################################################

  if (nrank==0) then
     write(*,"(' Reynolds number Re     : ',F17.3)") re
     write(*,"(' xnu                    : ',F17.8)") xnu
     print *,'==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     print *,'==========================================================='
     write(*,"(' Time step dt           : ',F17.8)") dt
     !
     if (itimescheme.eq.1) then
       !print *,'Temporal scheme        : Forwards Euler'
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
     ! Show the compile flags detected
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

#ifdef DEBUG
  if (nrank .eq. 0) write(*,*)'# parameter done'
#endif

  return
end subroutine parameter

!###########################################################################
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!
!###########################################################################
subroutine parameter_defaults()

  use param
  use variables
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

  xlx = twopi; yly = twopi; zlz = one
  nclx1 = 0; nclxn = 0
  ncly1 = 0; nclyn = 0
  nclz1 = 0; nclzn = 0
  
  npress = 1 !! By default people only need one pressure field
  imodulo2 = 1

end subroutine parameter_defaults
