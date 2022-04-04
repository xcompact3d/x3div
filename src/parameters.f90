!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

subroutine parameter()

  use iso_fortran_env

  use param
  use variables
  use complex_geometry
  use decomp_2d
  use ibm_param

  use var, only : dphi1

  implicit none

  real(mytype) :: theta, cfl,cf2
  integer :: longueur ,impi,j, is, total
#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter start'
#endif

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
  endif

  call parameter_defaults()
  
  !! Set Scalar BCs same as fluid (may be overridden) [DEFAULT]
  nclxS1 = nclx1; nclxSn = nclxn
  nclyS1 = ncly1; nclySn = nclyn
  nclzS1 = nclz1; nclzSn = nclzn
  if(ilesmod==0) then
     nu0nu=four
     cnu=0.44_mytype
  endif

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

  xnu=one/re
  !! Constant pressure gradient, re = Re_tau -> use to compute Re_centerline
  if (icpg.eq.one) then
    re_cent = (re/0.116)**(1.0/0.88)
    xnu = one/re_cent ! viscosity based on Re_cent to keep same scaling as CFR
    !
    fcpg = two/yly * (re/re_cent)**2
  end if

  if (ilmn) then
     if (ivarcoeff) then
        npress = 2 !! Need current pressure and previous iterate
     else
        npress = 1
     endif
  endif

  if (iimplicit.ne.0) then
     if ((itimescheme.eq.5).or.(itimescheme.eq.6)) then
        print *,'Error: implicit Y diffusion not yet compatible with RK time schemes'
        stop
     endif
     if (isecondder.eq.5) then
        print *, "Warning : support for implicit Y diffusion and isecondder=5 is experimental"
     endif
     if (iimplicit.eq.1) then
        xcst = dt * xnu
     else if (iimplicit.eq.2) then
        xcst = dt * xnu * half
     else
        print *, 'Error: wrong value for iimplicit ', iimplicit
        stop
     endif
     if (iscalar.eq.1) xcst_sc = xcst / sc
  endif

  if (itype.eq.itype_tbl.and.A_tr .gt. 0.0)  print *, "TBL tripping is active"

#ifdef DOUBLE_PREC
  anglex = dsin(pi*angle/180._mytype)
  angley = dcos(pi*angle/180._mytype)
#else
  anglex = sin(pi*angle/180._mytype)
  angley = cos(pi*angle/180._mytype)
#endif
  !###########################################################################
  ! Log-output
  !###########################################################################
  if (nrank==0) call system('mkdir data out probes 2> /dev/null')

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter input.i3d done'
#endif
  if (nrank==0) then
     print *,'==========================================================='
     if (itype.eq.itype_user) then
        print *,'User-defined simulation'
     elseif (itype.eq.itype_lockexch) then
        print *,'Simulating lock-exchange'
     elseif (itype.eq.itype_tgv) then
        print *,'Simulating TGV'
     elseif (itype.eq.itype_channel) then
        print *,'Simulating channel'
     elseif (itype.eq.itype_hill) then
        print *,'Simulating periodic hill'
     elseif (itype.eq.itype_cyl) then
        print *,'Simulating cylinder'
     elseif (itype.eq.itype_dbg) then
        print *,'Debug schemes'
     elseif (itype.eq.itype_mixlayer) then
        print *,'Mixing layer'
     elseif (itype.eq.itype_jet) then
        print *,'Jet'
     elseif (itype.eq.itype_tbl) then
        print *,'Turbulent boundary layer'
     else
        print *,'Unknown itype: ', itype
        stop
     endif
     print *,'==========================================================='
     if (itype.eq.itype_channel) then
       if (icpg.eq.zero) then
         print *,'Channel forcing with constant flow rate (CFR)'
         write(*,"(' Re_cl                  : ',F17.3)") re
       else if (icpg.eq.one) then
         print *,'Channel forcing with constant pressure gradient (CPG)'
         write(*,"(' Re_tau                 : ',F17.3)") re
         write(*,"(' Re_cl (estimated)      : ',F17.3)") re_cent
         write(*,"(' fcpg                   : ',F17.8)") fcpg
       end if
     else
       write(*,"(' Reynolds number Re     : ',F17.3)") re
     endif
     write(*,"(' xnu                    : ',F17.8)") xnu
     print *,'==========================================================='
     write(*,"(' p_row, p_col           : ',I9, I8)") p_row, p_col
     print *,'==========================================================='
     write(*,"(' Time step dt           : ',F17.8)") dt
     !
     if (itimescheme.eq.1) then
       !print *,'Temporal scheme        : Forwards Euler'
       write(*,"(' Temporal scheme        : ',A20)") "Forwards Euler"
     elseif (itimescheme.eq.2) then
       !print *,'Temporal scheme        : Adams-bashforth 2'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 2"
     elseif (itimescheme.eq.3) then
       !print *,'Temporal scheme        : Adams-bashforth 3'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 3"
     elseif (itimescheme.eq.4) then
       !print *,'Temporal scheme        : Adams-bashforth 4'
       write(*,"(' Temporal scheme        : ',A20)") "Adams-bashforth 4"
       print *,'Error: Adams-bashforth 4 not implemented!'
       stop
     elseif (itimescheme.eq.5) then
       !print *,'Temporal scheme        : Runge-kutta 3'
       write(*,"(' Temporal scheme        : ',A20)") "Runge-kutta 3"
     elseif (itimescheme.eq.6) then
       !print *,'Temporal scheme        : Runge-kutta 4'
       write(*,"(' Temporal scheme        : ',A20)") "Runge-kutta 4"
       print *,'Error: Runge-kutta 4 not implemented!'
       stop
     else
       print *,'Error: itimescheme must be specified as 1-6'
       stop
     endif
     !
     if (iimplicit.ne.0) then
       if (iimplicit.eq.1) then
         write(*,"('                          ',A40)") "With backward Euler for Y diffusion"
       else if (iimplicit.eq.2) then
         write(*,"('                          ',A40)") "With CN for Y diffusion"
       endif
     endif
     !
     if (ilesmod.ne.0) then
       print *,'                   : DNS'
     else
       if (jles.eq.1) then
          print *,'                   : Phys Smag'
       else if (jles.eq.2) then
          print *,'                   : Phys WALE'
       else if (jles.eq.3) then
          print *,'                   : Phys dyn. Smag'
       else if (jles.eq.4) then
          print *,'                   : iSVV'
       else
       endif
     endif
     print *,'==========================================================='
     write(*,"(' ifirst                 : ',I17)") ifirst
     write(*,"(' ilast                  : ',I17)") ilast
     print *,'==========================================================='
     write(*,"(' Lx                     : ',F17.8)") xlx
     write(*,"(' Ly                     : ',F17.8)") yly
     write(*,"(' Lz                     : ',F17.8)") zlz
     write(*,"(' nx                     : ',I17)") nx
     write(*,"(' ny                     : ',I17)") ny
     write(*,"(' nz                     : ',I17)") nz
     print *,'==========================================================='
     write(*,"(' istret                 : ',I17)") istret
     write(*,"(' beta                   : ',F17.8)") beta
     print *,'==========================================================='
     write(*,"(' nu0nu                  : ',F17.8)") nu0nu
     write(*,"(' cnu                    : ',F17.8)") cnu
     print *,'==========================================================='
     if (iscalar.eq.0) write(*,"(' Scalar                 : ',A17)") "off"
     if (iscalar.eq.1) write(*,"(' Scalar                 : ',A17)") "on"
     write(*,"(' numscalar              : ',I17)") numscalar
     if (iscalar.eq.1) then
       do is=1, numscalar
          write(*,"(' Schmidt number sc(',I2,')  : ',F17.8)") is, sc(is)
          write(*,"(' Richardson n.  ri(',I2,')  : ',F17.8)") is, ri(is)
          if (scalar_lbound(is).gt.-huge(one)) then
             write(*,"(' Lower bound      (',I2,')  : ',F17.8)") is, scalar_lbound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (scalar_ubound(is).lt.huge(one)) then
             write(*,"(' Upper bound      (',I2,')  : ',F17.8)") is, scalar_ubound(is)
          else
             ! This is the default option, no information printed in the listing
          endif
          if (iscalar.eq.1) then
             if (nclxS1.eq.1 .or. nclxSn.eq.1 .or. &
                 nclyS1.eq.1 .or. nclySn.eq.1 .or. &
                 nclzS1.eq.1 .or. nclzSn.eq.1) then
                if (sc_even(is)) then
                   ! This is the default option, no information printed in the listing
                else
                   write(*,"(' Scalar ',I2,' is odd')") is
                endif
             endif
             if (sc_skew(is)) then
                write(*,"(' Scalar ',I2,' with skew-symmetric convection')") is
             else
                ! This is the default option, no information printed in the listing
             endif
          endif
       end do
     endif
     print *,'==========================================================='
     write(*,"(' spinup_time            : ',I17)") spinup_time
     write(*,"(' wrotation              : ',F17.8)") wrotation
     print *,'==========================================================='
     if (iibm.eq.0) write(*,"(' Immersed boundary      : ',A17)") "off"
     if (iibm.gt.1) then
      write(*,"(' Immersed boundary      : ',A17)") "on"
      write(*,"(' iibm                   : ',I17)") iibm
     end if
     if (iibm.eq.1) print *,'Simple immersed boundary method'
     if (iibm.eq.2) then
       print *,'Lagrangian polynomial reconstruction'
       print *,'==========================================================='
       write(*,"(' npif                   : ',I17)") npif
       write(*,"(' izap                   : ',I17)") izap
       write(*,"(' nraf                   : ',I17)") nraf
       write(*,"(' nobjmax                : ',I17)") nobjmax
     end if
     print *,'==========================================================='
     write(*,"(' Boundary condition velocity field: ')")
     write(*,"(' nclx1, nclxn           : ',I15,',',I1 )") nclx1,nclxn
     write(*,"(' ncly1, nclyn           : ',I15,',',I1 )") ncly1,nclyn
     write(*,"(' nclz1, nclzn           : ',I15,',',I1 )") nclz1,nclzn
     print *,'==========================================================='
     if ((iscalar.eq.1).or.(ilmn)) then
       write(*,"(' Boundary condition scalar field: ')")
       write(*,"(' nclxS1, nclxSn         : ',I15,',',I1 )") nclxS1,nclxSn
       write(*,"(' nclyS1, nclySn         : ',I15,',',I1 )") nclyS1,nclySn
       write(*,"(' nclzS1, nclzSn         : ',I15,',',I1 )") nclzS1,nclzSn
       print *,'==========================================================='
     endif

#ifdef DOUBLE_PREC
#ifdef SAVE_SINGLE
     print *,'Numerical precision: Double, saving in single'
#else
     print *,'Numerical precision: Double'
#endif
#else
     print *,'Numerical precision: Single'
#endif
     print *,'==========================================================='
     write(*,"(' High and low speed : u1=',F6.2,' and u2=',F6.2)") u1,u2
     write(*,"(' Gravity vector     : (gx, gy, gz)=(',F15.8,',',F15.8,',',F15.8,')')") gravx, gravy, gravz
     if (ilmn) then
        print *, "LMN                : Enabled"
        if (ivarcoeff) then
           print *, "LMN-Poisson solver : Variable-coefficient"
        else
           print *, "LMN-Poisson solver : Constant-coefficient"
        endif
        if (ilmn_bound) then
           print *, "LMN boundedness    : Enforced"
        else
           print *, "LMN boundedness    : Not enforced"
        endif
        write(*,"(' dens1 and dens2    : ',F6.2' ',F6.2)") dens1, dens2
        write(*,"(' Prandtl number Re  : ',F15.8)") prandtl
     endif
     if (angle.ne.0.) write(*,"(' Solid rotation     : ',F6.2)") angle
     print *, ' '
     print *,'==========================================================='
  endif

#ifdef DEBG
  if (nrank .eq. 0) print *,'# parameter done'
#endif

  return
end subroutine parameter

!###########################################################################
!
!  SUBROUTINE: parameter_defaults
! DESCRIPTION: Sets the default simulation parameters.
!      AUTHOR: Paul Bartholomew <paul.bartholomew08@imperial.ac.uk>
!
!###########################################################################
subroutine parameter_defaults()

  use param
  use variables
  use decomp_2d
  use complex_geometry

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
  filepath = './data/'
  irestart = 0
  itime0 = 0
  t0 = zero
  datapath = './data/'

  xlx = one; yly = one; zlz = one
  nclx1 = 0; nclxn = 0
  ncly1 = 0; nclyn = 0
  nclz1 = 0; nclzn = 0
  
  !! IBM stuff
  nraf = 0
  nobjmax = 0

  itrip = 0
  wrotation = zero
  irotation = 0
  itest=1

  !! Gravity field
  gravx = zero
  gravy = zero
  gravz = zero

  !! LMN stuff
  ilmn = .FALSE.
  ilmn_bound = .TRUE.
  pressure0 = one
  prandtl = one
  dens1 = one
  dens2 = one
  ivarcoeff = .FALSE.
  npress = 1 !! By default people only need one pressure field
  ilmn_solve_temp = .FALSE.
  imultispecies = .FALSE.
  Fr = zero
  ibirman_eos = .FALSE.

  primary_species = -1

  !! Channel
  icpg = 0
  icfr = 1

  !! IO
  ivisu = 1
  ipost = 0
  iprocessing = huge(i)
  initstat = huge(i)

  save_ux = 0
  save_uy = 0
  save_uz = 0
  save_phi = 0
  save_uxm = 0
  save_uym = 0
  save_uzm = 0
  save_phim = 0
  save_w = 0
  save_w1 = 0
  save_w2 = 0
  save_w3 = 0
  save_qc = 0
  save_pc = 0
  save_V = 0
  save_dudx = 0
  save_dudy = 0
  save_dudz = 0
  save_dvdx = 0
  save_dvdy = 0
  save_dvdz = 0
  save_dwdx = 0
  save_dwdy = 0
  save_dwdz = 0
  save_dphidx = 0
  save_dphidy = 0
  save_dphidz = 0
  save_pre = 0
  save_prem = 0
  save_dmap = 0
  save_utmap = 0
  save_ibm = 0

  ipost = 0
  iibm=0
  npif = 2
  izap = 1

  imodulo2 = 1

  !! CASE specific variables
  tgv_twod = .FALSE.

  !! TRIPPING
  A_tr=0.0
  xs_tr_tbl=1.402033
  ys_tr_tbl=0.350508
  ts_tr_tbl=1.402033
  x0_tr_tbl=3.505082

  !! X3DIV
  test_mode = .FALSE.
  
end subroutine parameter_defaults
