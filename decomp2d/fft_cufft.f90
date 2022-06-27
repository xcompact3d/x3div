!=======================================================================
! This is part of the 2DECOMP&FFT library
! 
! 2DECOMP&FFT is a software framework for general-purpose 2D (pencil) 
! decomposition. It also implements a highly scalable distributed
! three-dimensional Fast Fourier Transform (FFT).
!
! Copyright (C) 2009-2011 Ning Li, the Numerical Algorithms Group (NAG)
!
!=======================================================================

! This is the FFTW (version 3.x) implementation of the FFT library

module decomp_2d_fft

  use decomp_2d  ! 2D decomposition module
  use iso_c_binding
  use cudafor
  use cufft

  implicit none

  private        ! Make everything private unless declared public

  ! engine-specific global variables
  ! integer, save :: plan_type = FFTW_MEASURE

  ! FFTW plans
  ! j=1,2,3 corresponds to the 1D FFTs in X,Y,Z direction, respectively
  ! For c2c transforms: 
  !     use plan(-1,j) for  forward transform; 
  !     use plan( 1,j) for backward transform;
  ! For r2c/c2r transforms:
  !     use plan(0,j) for r2c transforms;
  !     use plan(2,j) for c2r transforms;
  integer*4, save :: plan(-1:2,3)
  complex*8, device, allocatable, dimension(:) :: cufft_workspace

  ! common code used for all engines, including global variables, 
  ! generic interface definitions and several subroutines
#include "fft_common.f90"

  ! Return a cuFFT plan for multiple 1D FFTs in X direction
  subroutine plan_1m_x(plan1, decomp, cufft_type, worksize)
    
    implicit none

    integer*4, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: cufft_type

    integer :: istat
    integer(int_ptr_kind()), intent(out) :: worksize
    integer, pointer :: null_fptr
    call c_f_pointer( c_null_ptr, null_fptr )
   
    istat = cufftCreate(plan1)
    istat = istat + cufftSetAutoAllocation(plan1,0)
    istat = istat + cufftMakePlanMany(plan1,1,                           &
                              decomp%xsz(1),null_fptr,1,         &
                              decomp%xsz(1),null_fptr,1,         &
                              decomp%xsz(1),cufft_type,          &
                              decomp%xsz(2)*decomp%xsz(3),worksize)
    if (istat /= 0) &
       write (*,*) "Cannot create plan for plan_1m_x FFTs batch (cufftPlanMany)"

    return
  end subroutine plan_1m_x

  ! Return a cuFFT plan for multiple 1D FFTs in Y direction
  subroutine plan_1m_y(plan1, decomp, cufft_type,worksize)

    implicit none

    integer*4, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: cufft_type

    ! Due to memory pattern of 3D arrays, 1D FFTs along Y have to be
    ! done one Z-plane at a time. So plan for 2D data sets here.
    integer :: istat
    integer(int_ptr_kind()), intent(out) :: worksize
    integer, pointer :: null_fptr
    call c_f_pointer( c_null_ptr, null_fptr )
   
    istat = cufftCreate(plan1)
    istat = istat + cufftSetAutoAllocation(plan1,0)
    istat = istat + cufftMakePlanMany(plan1,1,                    &
                              decomp%ysz(2),null_fptr,1,  &
                              decomp%ysz(2),null_fptr,1,  &
                              decomp%ysz(2),cufft_type,   &
                              decomp%ysz(1),worksize)
    if (istat /= 0) &
       write (*,*) "Cannot create plan for plan_1m_y FFTs batch (cufftPlanMany)"

    return 
  end subroutine plan_1m_y

  ! Return a cuFFT plan for multiple 1D FFTs in Z direction
  subroutine plan_1m_z(plan1, decomp, cufft_type,worksize)

    implicit none

    integer*4, intent(OUT) :: plan1
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: cufft_type

    integer :: istat
    integer(int_ptr_kind()), intent(out) :: worksize
    integer, pointer :: null_fptr
    call c_f_pointer( c_null_ptr, null_fptr )
   
    istat = cufftCreate(plan1)
    istat = istat + cufftSetAutoAllocation(plan1,0)
    istat = istat + cufftMakePlanMany(plan1,1,                      &
                              decomp%zsz(3),null_fptr,1,    &
                              decomp%zsz(3),null_fptr,1,    &
                              decomp%zsz(3),cufft_type,     &
                              decomp%zsz(1)*decomp%zsz(2),worksize)

    if (istat /= 0) &
       write (*,*) "Cannot create plan for plan_1m_z FFTs batch (cufftPlanMany)"


    return
  end subroutine plan_1m_z

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time initialisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine init_fft_engine

    implicit none

    !integer*4 :: cufft_ws, ws
    integer(int_ptr_kind()) :: cufft_ws, ws
    integer   :: i, j, istat

    if (nrank==0) then
       write(*,*) ' '
       write(*,*) '***** Using the cuFFT engine *****'
       write(*,*) ' '
    end if
    
     cufft_ws = 0 
#ifdef DOUBLE_PREC
    if (format == PHYSICAL_IN_X) then
       ! For C2C transforms
       call plan_1m_x(plan(-1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(-1,2), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(-1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan( 1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan( 1,2), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan( 1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       ! For R2C/C2R tranforms
       call plan_1m_x(plan(0,1), ph, CUFFT_D2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(0,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(0,3), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(2,3), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(2,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(2,1), sp, CUFFT_Z2D,ws)
       cufft_ws = max (cufft_ws,ws)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       write(*,*) 'Create the plans'
       call plan_1m_z(plan(-1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(-1,2), ph, CUFFT_Z2Z,ws) 
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(-1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan( 1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan( 1,2), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan( 1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)

       ! For R2C/C2R tranforms
       call plan_1m_z(plan(0,3), ph, CUFFT_D2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(0,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(0,1), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(2,1), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(2,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(2,3), sp, CUFFT_Z2D,ws)
       cufft_ws = max (cufft_ws,ws)

    end if
#else
    if (format == PHYSICAL_IN_X) then
       ! For C2C transforms
       call plan_1m_x(plan(-1,1), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(-1,2), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(-1,3), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan( 1,3), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan( 1,2), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan( 1,1), ph, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       ! For R2C/C2R tranforms
       call plan_1m_x(plan(0,1), ph, CUFFT_R2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(0,2), sp, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(0,3), sp, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(2,3), sp, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(2,2), sp, CUFFT_C2C,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(2,1), sp, CUFFT_C2R,ws)
       cufft_ws = max (cufft_ws,ws)

    else if (format == PHYSICAL_IN_Z) then

       ! For C2C transforms
       call plan_1m_z(plan(-1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(-1,2), ph, CUFFT_Z2Z,ws) 
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(-1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan( 1,1), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan( 1,2), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan( 1,3), ph, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)

       ! For R2C/C2R tranforms
       call plan_1m_z(plan(0,3), ph, CUFFT_D2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(0,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(0,1), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_x(plan(2,1), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_y(plan(2,2), sp, CUFFT_Z2Z,ws)
       cufft_ws = max (cufft_ws,ws)
       call plan_1m_z(plan(2,3), sp, CUFFT_Z2R,ws)
       cufft_ws = max (cufft_ws,ws)
    end if
#endif
    allocate(cufft_workspace(cufft_ws))
    do i=-1,2
      do j=1,3
         istat = cufftSetWorkArea(plan(i,j),cufft_workspace)
      enddo
    enddo

    return
  end subroutine init_fft_engine


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !  This routine performs one-time finalisations for the FFT engine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine finalize_fft_engine

    implicit none

    integer :: i,j, istat

    do j=1,3
       do i=-1,2
          istat = cufftDestroy(plan(i,j))
       end do
    end do

    return
  end subroutine finalize_fft_engine


  ! Following routines calculate multiple one-dimensional FFTs to form 
  ! the basis of three-dimensional FFTs.

  ! c2c transform, multiple 1D FFTs in x direction
  subroutine c2c_1m_x(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*4, intent(IN) :: plan1
    integer :: istat

#ifdef DOUBLE_PREC
    !$acc host_data use_device(inout)
    istat = cufftExecZ2Z(plan1, inout, inout,isign)
    !$acc end host_data
#else
    !$acc host_data use_device(inout)
    istat = cufftExecC2C(plan1, inout, inout,isign)
    !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing c2c_1m_x"

    return
  end subroutine c2c_1m_x


  ! c2c transform, multiple 1D FFTs in y direction
  subroutine c2c_1m_y(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*4, intent(IN) :: plan1
    integer :: istat

    integer :: k, s3

    ! transform on one Z-plane at a time
    s3 = size(inout,3)
    do k=1,s3
#ifdef DOUBLE_PREC
       !$acc host_data use_device(inout)
       istat = cufftExecZ2Z(plan1, inout(:,:,k), inout(:,:,k),isign)
       !$acc end host_data
#else
       !$acc host_data use_device(inout)
       istat = cufftExecC2C(plan1, inout(:,:,k), inout(:,:,k),isign)
       !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing c2c_1m_y istat ", istat, isign
    end do

    return
  end subroutine c2c_1m_y

  ! c2c transform, multiple 1D FFTs in z direction
  subroutine c2c_1m_z(inout, isign, plan1)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: inout
    integer, intent(IN) :: isign
    integer*4, intent(IN) :: plan1
    integer :: istat

#ifdef DOUBLE_PREC
    !$acc host_data use_device(inout)
    istat = cufftExecZ2Z(plan1, inout, inout,isign)
    !$acc end host_data
#else
    !$acc host_data use_device(inout)
    istat = cufftExecC2C(plan1, inout, inout,isign)
    !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing c2c_1m_z"

    return
  end subroutine c2c_1m_z

  ! r2c transform, multiple 1D FFTs in x direction
  subroutine r2c_1m_x(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    integer :: istat


#ifdef DOUBLE_PREC
    !$acc host_data use_device(input,output)
    istat =  cufftExecD2Z(plan(0,1), input, output)
    !$acc end host_data
#else
    !$acc host_data use_device(input,output)
    istat =  cufftExecR2C(plan(0,1), input, output)
    !$acc end host_data
#endif    
    if (istat /= 0) &
       write (*,*) "Error in executing r2c_1m_x istat ", istat

    return

  end subroutine r2c_1m_x

  ! r2c transform, multiple 1D FFTs in z direction
  subroutine r2c_1m_z(input, output)

    implicit none

    real(mytype), dimension(:,:,:), intent(IN)  ::  input
    complex(mytype), dimension(:,:,:), intent(OUT) :: output
    integer :: istat

#ifdef DOUBLE_PREC
    !$acc host_data use_device(input,output)
    istat = cufftExecD2Z(plan(0,3), input, output)
    !$acc end host_data
#else
    !$acc host_data use_device(input,output)
    istat = cufftExecR2C(plan(0,3), input, output)
    !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing r2c_1m_z istat ", istat

    return

  end subroutine r2c_1m_z

  ! c2r transform, multiple 1D FFTs in x direction
  subroutine c2r_1m_x(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN)  ::  input
    real(mytype), dimension(:,:,:), intent(OUT) :: output
    integer :: istat

#ifdef DOUBLE_PREC
    !$acc host_data use_device(input,output)
    istat = cufftExecZ2D(plan(2,1), input, output)
    !$acc end host_data
#else
    !$acc host_data use_device(input,output)
    istat = cufftExecC2R(plan(2,1), input, output)
    !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing c2r_1m_x"

    return

  end subroutine c2r_1m_x

  ! c2r transform, multiple 1D FFTs in z direction
  subroutine c2r_1m_z(input, output)

    implicit none

    complex(mytype), dimension(:,:,:), intent(IN) :: input
    real(mytype), dimension(:,:,:), intent(OUT) :: output
    integer :: istat

#ifdef DOUBLE_PREC
    !$acc host_data use_device(input,output)
    istat = cufftExecZ2D(plan(2,3), input, output)
    !$acc end host_data
#else
    !$acc host_data use_device(input,output)
    istat = cufftExecC2R(plan(2,3), input, output)
    !$acc end host_data
#endif
    if (istat /= 0) &
       write (*,*) "Error in executing c2r_1m_z"

    return

  end subroutine c2r_1m_z



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D FFT - complex to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2c(in, out, isign)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: in
    complex(mytype), dimension(:,:,:), intent(OUT) :: out
    integer, intent(IN) :: isign

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_FORWARD .OR.  &
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_BACKWARD) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in,isign,plan(isign,1))
#else
       allocate (wk1(ph%xsz(1),ph%xsz(2),ph%xsz(3)))
       wk1 = in
       call c2c_1m_x(wk1,isign,plan(isign,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====

       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in,wk2_c2c,ph)
#else
          call transpose_x_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else
#ifdef OVERWRITE
          call c2c_1m_y(in,isign,plan(isign,2))
#else
          call c2c_1m_y(wk1,isign,plan(isign,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_c2c,out,ph)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in,out,ph)
#else
          call transpose_y_to_z(wk1,out,ph)
#endif
       end if
       call c2c_1m_z(out,isign,plan(isign,3))

    else if (format==PHYSICAL_IN_X .AND. isign==DECOMP_2D_FFT_BACKWARD &
         .OR. & 
         format==PHYSICAL_IN_Z .AND. isign==DECOMP_2D_FFT_FORWARD) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in,isign,plan(isign,3))
#else
       allocate (wk1(ph%zsz(1),ph%zsz(2),ph%zsz(3)))
       wk1 = in
       call c2c_1m_z(wk1,isign,plan(isign,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_z_to_y(in,wk2_c2c,ph)
#else
          call transpose_z_to_y(wk1,wk2_c2c,ph)
#endif
          call c2c_1m_y(wk2_c2c,isign,plan(isign,2))
       else  ! out==wk2_c2c if 1D decomposition
#ifdef OVERWRITE
          call transpose_z_to_y(in,out,ph)
#else
          call transpose_z_to_y(wk1,out,ph)
#endif
          call c2c_1m_y(out,isign,plan(isign,2))
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_c2c,out,ph)
       end if
       call c2c_1m_x(out,isign,plan(isign,1))

    end if

#ifndef OVERWRITE
    deallocate (wk1)
#endif

    return
  end subroutine fft_3d_c2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D forward FFT - real to complex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_r2c(in_r, out_c)
    use nvtx
    implicit none

    real(mytype), dimension(:,:,:), intent(IN) :: in_r
    complex(mytype), dimension(:,:,:), intent(OUT) :: out_c

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in X =====
       call r2c_1m_x(in_r,wk13)

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call transpose_x_to_y(wk13,wk2_r2c,sp)
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
       else
          call c2c_1m_y(wk13,-1,plan(0,2))
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,out_c,sp)
       else
          call transpose_y_to_z(wk13,out_c,sp)
       end if
       call c2c_1m_z(out_c,-1,plan(0,3))

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in Z =====
       call nvtxStartRange("Z r2c_1m_z")
       call r2c_1m_z(in_r,wk13)
       call nvtxEndRange

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
          call nvtxStartRange("Z1 transpose_z_to_y")
          call transpose_z_to_y(wk13,wk2_r2c,sp)
          call nvtxEndRange
          call nvtxStartRange("Z1 c2c_1m_y")
          call c2c_1m_y(wk2_r2c,-1,plan(0,2))
          call nvtxEndRange
       else  ! out_c==wk2_r2c if 1D decomposition
          call nvtxStartRange("Z transpose_z_to_y")
          call transpose_z_to_y(wk13,out_c,sp)
          call nvtxEndRange
          call nvtxStartRange("Z c2c_1m_y")
          call c2c_1m_y(out_c,-1,plan(0,2))
          call nvtxEndRange
       end if

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call nvtxStartRange("Z1 transpose_y_to_x")
          call transpose_y_to_x(wk2_r2c,out_c,sp)
          call nvtxEndRange
       end if
       call nvtxStartRange("c2c_1m_x")
       call c2c_1m_x(out_c,-1,plan(0,1))
       call nvtxEndRange

    end if

    return
  end subroutine fft_3d_r2c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 3D inverse FFT - complex to real
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine fft_3d_c2r(in_c, out_r)

    implicit none

    complex(mytype), dimension(:,:,:), intent(INOUT) :: in_c
    real(mytype), dimension(:,:,:), intent(OUT) :: out_r

#ifndef OVERWRITE
    complex(mytype), allocatable, dimension(:,:,:) :: wk1
#endif

    if (format==PHYSICAL_IN_X) then

       ! ===== 1D FFTs in Z =====
#ifdef OVERWRITE
       call c2c_1m_z(in_c,1,plan(2,3))       
#else
       allocate (wk1(sp%zsz(1),sp%zsz(2),sp%zsz(3)))
       wk1 = in_c
       call c2c_1m_z(wk1,1,plan(2,3))
#endif

       ! ===== Swap Z --> Y; 1D FFTs in Y =====
#ifdef OVERWRITE
       call transpose_z_to_y(in_c,wk2_r2c,sp)
#else
       call transpose_z_to_y(wk1,wk2_r2c,sp)
#endif
       call c2c_1m_y(wk2_r2c,1,plan(2,2))

       ! ===== Swap Y --> X; 1D FFTs in X =====
       if (dims(1)>1) then
          call transpose_y_to_x(wk2_r2c,wk13,sp)
          call c2r_1m_x(wk13,out_r)
       else
          call c2r_1m_x(wk2_r2c,out_r)
       end if

    else if (format==PHYSICAL_IN_Z) then

       ! ===== 1D FFTs in X =====
#ifdef OVERWRITE
       call c2c_1m_x(in_c,1,plan(2,1))
#else
       allocate(wk1(sp%xsz(1),sp%xsz(2),sp%xsz(3)))
       wk1 = in_c
       call c2c_1m_x(wk1,1,plan(2,1))
#endif

       ! ===== Swap X --> Y; 1D FFTs in Y =====
       if (dims(1)>1) then
#ifdef OVERWRITE
          call transpose_x_to_y(in_c,wk2_r2c,sp)
#else
          call transpose_x_to_y(wk1,wk2_r2c,sp)
#endif
          call c2c_1m_y(wk2_r2c,1,plan(2,2))
       else  ! in_c==wk2_r2c if 1D decomposition
#ifdef OVERWRITE
          call c2c_1m_y(in_c,1,plan(2,2))
#else
          call c2c_1m_y(wk1,1,plan(2,2))
#endif
       end if

       ! ===== Swap Y --> Z; 1D FFTs in Z =====
       if (dims(1)>1) then
          call transpose_y_to_z(wk2_r2c,wk13,sp)
       else
#ifdef OVERWRITE
          call transpose_y_to_z(in_c,wk13,sp)
#else
          call transpose_y_to_z(wk1,wk13,sp)
#endif
       end if
       call c2r_1m_z(wk13,out_r)

    end if

#ifndef OVERWRITE
    deallocate (wk1)
#endif

    return
  end subroutine fft_3d_c2r


end module decomp_2d_fft
