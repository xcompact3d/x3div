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
module variables
  !USE param
  !USE var
  use decomp_2d, only : mytype

  ! Boundary conditions : ncl = 2 --> Dirichlet
  ! Boundary conditions : ncl = 1 --> Free-slip
  ! Boundary conditions : ncl = 0 --> Periodic
  ! l: power of 2,3,4,5 and 6
  ! if ncl = 1 or 2, --> n  = 2l+ 1
  !                  --> nm = n - 1
  !                  --> m  = n + 1
  ! If ncl = 0,      --> n  = 2*l
  !                  --> nm = n
  !                  --> m  = n + 2
  !nstat = size arrays for statistic collection
  !2-->every 2 mesh nodes
  !4-->every 4 mesh nodes
  !nvisu = size for visualization collection
  !nprobe =  size for probe collection (energy spectra)

  !Possible n points: 3 5 7 9 11 13 17 19 21 25 31 33 37 41 49 51 55 61 65 73 81 91 97 101 109 121 129 145 151 161 163 181 193 201 217 241 251 257 271 289 301 321 325 361 385 401 433 451 481 487 501 513 541 577 601 641 649 721 751 769 801 811 865 901 961 973 1001 1025 1081 1153 1201 1251 1281 1297 1351 1441 1459 1501 1537 1601 1621 1729 1801 1921 1945 2001 2049 2161 2251 2305 2401 2431 2501 2561 2593 2701 2881 2917 3001 3073 3201 3241 3457 3601 3751 3841 3889 4001 4051 4097 4321 4375 4501 4609 4801 4861 5001 5121 5185 5401 5761 5833 6001 6145 6251 6401 6481 6751 6913 7201 7291 7501 7681 7777 8001 8101 8193 8641 8749 9001 9217 9601 9721 enough

  integer :: nx,ny,nz,numscalar,p_row,p_col,nxm,nym,nzm,spinup_time
  integer :: nstat=1,nvisu=1,nprobe=1,nlength=1

  real(mytype),allocatable,dimension(:) :: sc,uset,cp,ri,group
  real(mytype) :: nu0nu, cnu

#ifndef DOUBLE_PREC
  integer,parameter :: prec = 4
#else
#ifdef SAVE_SINGLE
  integer,parameter :: prec = 4
#else
  integer,parameter :: prec = 8
#endif
#endif

  !module filter
  real(mytype),dimension(200) :: idata

  real(mytype),allocatable,dimension(:,:) :: fisx,fivx
  real(mytype),allocatable,dimension(:,:) :: fisy,fivy
  real(mytype),allocatable,dimension(:,:) :: fisz,fivz


  real(mytype), save, allocatable, dimension(:,:) :: sx,vx
  real(mytype), save, allocatable, dimension(:,:) :: sy,vy
  real(mytype), save, allocatable, dimension(:,:) :: sz,vz

  !! X3DIV
  logical :: test_mode

  !O6SVV
  real(mytype),allocatable,dimension(:) :: newsm,newtm,newsmt,newtmt
  !real(mytype),allocatable,dimension(:) :: newrm,ttm,newrmt,ttmt
  real(mytype),allocatable,dimension(:) :: newrm,newrmt

  !module pressure
  real(mytype), save, allocatable, dimension(:,:) :: dpdyx1,dpdyxn,dpdzx1,dpdzxn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxy1,dpdxyn,dpdzy1,dpdzyn
  real(mytype), save, allocatable, dimension(:,:) :: dpdxz1,dpdxzn,dpdyz1,dpdyzn

  !module waves
  complex(mytype),allocatable,dimension(:) :: zkz,zk2,ezs
  complex(mytype),allocatable,dimension(:) :: yky,yk2,eys
  complex(mytype),allocatable,dimension(:) :: xkx,xk2,exs

  !module mesh
  real(mytype),allocatable,dimension(:) :: ppy,pp2y,pp4y
  real(mytype),allocatable,dimension(:) :: ppyi,pp2yi,pp4yi
  real(mytype),allocatable,dimension(:) :: xp,xpi,yp,ypi,dyp,zp,zpi,del
  real(mytype),allocatable,dimension(:) :: yeta,yetai
  real(mytype) :: alpha,beta

end module variables
!############################################################################
!############################################################################
module param

  use decomp_2d, only : mytype

  integer :: nclx1,nclxn,ncly1,nclyn,nclz1,nclzn
  integer :: nclxS1,nclxSn,nclyS1,nclySn,nclzS1,nclzSn

  !logical variable for boundary condition that is true in periodic case
  !and false otherwise
  logical :: nclx,ncly,nclz

  integer :: itype
  integer, parameter :: &
       itype_user = 0, &
       itype_lockexch = 1, &
       itype_tgv = 2, &
       itype_channel = 3, &
       itype_hill = 4, &
       itype_cyl = 5, &
       itype_dbg = 6, &
       itype_mixlayer = 7, &
       itype_jet = 8, &
       itype_tbl = 9

  integer :: cont_phi,itr,itime,itest,iprocessing
  integer :: ifft,istret,iforc_entree,iturb
  integer :: iin,itimescheme,iimplicit,ifirst,ilast,iles
  integer :: ntime ! How many (sub)timestpeps do we need to store?
  integer :: icheckpoint,irestart,idebmod,ioutput,imodulo2,idemarre,icommence,irecord
  integer :: itime0
  integer :: iscalar,nxboite,istat,iread,iadvance_time,irotation,iibm
  integer :: npif,izap
  integer :: ivisu, ipost, initstat
  real(mytype) :: xlx,yly,zlz,dx,dy,dz,dx2,dy2,dz2,t,xxk1,xxk2,t0
  real(mytype) :: dt,re,xnu,init_noise,inflow_noise,u1,u2,angle,anglex,angley
  real(mytype) :: wrotation,ro
  real(mytype) :: dens1, dens2


  !! Numerics control
  integer :: ifirstder,isecondder,ipinter

  !! CFL_diffusion parameter
  real(mytype) :: cfl_diff_x,cfl_diff_y,cfl_diff_z,cfl_diff_sum

  !!
  real(mytype) :: xcst
  real(mytype), allocatable, dimension(:) :: xcst_sc
  real(mytype), allocatable, dimension(:,:) :: alpha_sc, beta_sc, g_sc
  real(mytype) :: g_bl_inf, f_bl_inf


  integer :: npress
  real(mytype), dimension(5) :: adt,bdt,cdt,ddt,gdt
  !numbers

  real(mytype),parameter :: half=0.5_mytype
  real(mytype),parameter :: twothird=2._mytype/3._mytype
  real(mytype),parameter :: zero=0._mytype
  real(mytype),parameter :: one=1._mytype
  real(mytype),parameter :: onepfive=1.5_mytype
  real(mytype),parameter :: two=2._mytype
  real(mytype),parameter :: twopfive=2.5_mytype
  real(mytype),parameter :: three=3._mytype
  real(mytype),parameter :: threepfive=3.5_mytype
  real(mytype),parameter :: four=4._mytype
  real(mytype),parameter :: five=5._mytype
  real(mytype),parameter :: six=6._mytype
  real(mytype),parameter :: seven=7._mytype
  real(mytype),parameter :: eight=8._mytype
  real(mytype),parameter :: nine=9._mytype
  real(mytype),parameter :: ten=10._mytype
  
  complex(mytype),parameter :: cx_one_one=cmplx(one, one, kind=mytype)

end module param
!############################################################################

