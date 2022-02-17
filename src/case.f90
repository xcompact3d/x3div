!Copyright (c) 2012-2022, Xcompact3d
!This file is part of Xcompact3d (xcompact3d.com)
!SPDX-License-Identifier: BSD 3-Clause

module case

  use param, only : itype, itype_tgv2d
  use decomp_2d, only : mytype, xsize

  use tgv2d

  implicit none

  private ! All functions/subroutines private by default
  public :: case_boot, &
            case_listing, &
            case_init, &
            case_bc, &
            case_forcing, &
            case_visu, &
            case_postprocess, &
            case_finalize

contains

  !
  ! Read case-specific parameters in the input file
  ! Initialize case-specific IO
  ! Allocate memory
  !
  subroutine case_boot()

    implicit none

    if (itype == itype_tgv2d) call tgv2d_boot()

  end subroutine case_boot

  !
  ! Print case-specific parameters in the listing
  !
  subroutine case_listing()

    implicit none

    if (itype == itype_tgv2d) call tgv2d_listing()

  end subroutine case_listing

  !
  ! Case-specific initialization
  !
  subroutine case_init(ux1, uy1, uz1)

    implicit none

    ! Arguments
    real(mytype),intent(out),dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

    ! Case-specific init
    if (itype == itype_tgv2d) call tgv2d_init(ux1, uy1, uz1)

  end subroutine case_init

  !
  ! Case-specific boundary conditions
  !
  subroutine case_bc(ux1, uy1, uz1)

    implicit none

    ! Arguments
    real(mytype), intent(inout), dimension(xsize(1),xsize(2),xsize(3)) :: ux1, uy1, uz1

  end subroutine case_bc

  !
  ! Add case-specific forcing term in the momentum r.h.s.
  !
  subroutine case_forcing(dux1, duy1, duz1)

    use param, only : ntime

    implicit none

    ! Arguments
    real(mytype), intent(inout), dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1, duy1, duz1

  end subroutine case_forcing

  !
  ! Visualization
  ! This is called when itime % ioutput = 0
  !
  subroutine case_visu()

    implicit none

  end subroutine case_visu

  !
  ! Case-specific post-processing
  ! This is called at the end of each time step
  !
  subroutine case_postprocess(ux1, uy1, uz1, ndt)

    use param, only : ivisu, ioutput

    implicit none

    ! Arguments
    real(mytype), dimension(xsize(1),xsize(2),xsize(3)), intent(in) :: ux1, uy1, uz1
    integer, intent(in) :: ndt

    if (itype == itype_tgv2d) call tgv2d_postprocess(ux1, uy1, uz1, ndt)

    if ((ivisu /= 0).and.(ioutput /= 0)) then
      if (mod(ndt, ioutput) == 0) call case_visu()
    endif

  end subroutine case_postprocess

  !
  ! Finalize case-specific IO
  ! Free memory
  !
  subroutine case_finalize()

    implicit none

    if (itype == itype_tgv2d) call tgv2d_finalize()

  end subroutine case_finalize

end module case
