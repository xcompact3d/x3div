!################################################################################
!This file is part of Xcompact3d
!
!Xcompact3d
!Copyright (c) 2012-2022, Xcompact3d
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!BSD 3-Clause License
!
!Redistribution and use in source and binary forms, with or without
!modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this
!   list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice,
!   this list of conditions and the following disclaimer in the documentation
!   and/or other materials provided with the distribution.!
!
!3. Neither the name of the copyright holder nor the names of its
!   contributors may be used to endorse or promote products derived from
!   this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
!DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
!FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
!DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
!SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
!CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
!OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

module transeq

  private
  public :: calculate_transeq_rhs

contains
  !############################################################################
  !!  SUBROUTINE: calculate_transeq_rhs
  !! DESCRIPTION: Calculates the right hand sides of all transport
  !!              equations - momentum, scalar transport, etc.
  !############################################################################
  subroutine calculate_transeq_rhs(dux1,duy1,duz1,ux1,uy1,uz1)

    use decomp_2d, only : mytype
    use decomp_2d, only : xsize, zsize
    use param, only : ntime

    implicit none

    !! Inputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3)), intent(in) :: ux1, uy1, uz1

    !! Outputs
    real(mytype), dimension(xsize(1), xsize(2), xsize(3), ntime) :: dux1, duy1, duz1

    !! Momentum equations
    call momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1)

  end subroutine calculate_transeq_rhs
  !############################################################################
  !############################################################################
  !!
  !!  subroutine: momentum_rhs_eq
  !! DESCRIPTION: Calculation of convective and diffusion terms of momentum
  !!              equation
  !!
  !############################################################################
  !############################################################################
  subroutine momentum_rhs_eq(dux1,duy1,duz1,ux1,uy1,uz1)

    use param
    use variables
    use x3d_operator_1d
    use decomp_2d, only : mytype
    use x3d_transpose
    use x3d_derive
    use decomp_2d , only : xsize, ysize, zsize
    use var, only : ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1
    use var, only : ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2
    use var, only : ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3

    implicit none

    !! INPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1

    !! OUTPUTS
    real(mytype),dimension(xsize(1),xsize(2),xsize(3),ntime) :: dux1,duy1,duz1


    integer :: i,j,k,is

    !SKEW SYMMETRIC FORM
    !WORK X-PENCILS
    
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      ta1(i,j,k) = ux1(i,j,k) * ux1(i,j,k)
      tb1(i,j,k) = ux1(i,j,k) * uy1(i,j,k)
      tc1(i,j,k) = ux1(i,j,k) * uz1(i,j,k)
    enddo

    call derx (td1,ta1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (te1,tb1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tf1,tc1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (ta1,ux1,sx,x3d_op_derx, xsize(1),xsize(2),xsize(3))
    call derx (tb1,uy1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))
    call derx (tc1,uz1,sx,x3d_op_derxp,xsize(1),xsize(2),xsize(3))

    ! Convective terms of x-pencil are stored in tg1,th1,ti1
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      dux1(i,j,k,1) = -half*(td1(i,j,k) + ux1(i,j,k) * ta1(i,j,k))
      duy1(i,j,k,1) = -half*(te1(i,j,k) + ux1(i,j,k) * tb1(i,j,k))
      duz1(i,j,k,1) = -half*(tf1(i,j,k) + ux1(i,j,k) * tc1(i,j,k))
    enddo

    call x3d_transpose_x_to_y(ux1,ux2)
    call x3d_transpose_x_to_y(uy1,uy2)
    call x3d_transpose_x_to_y(uz1,uz2)

    !WORK Y-PENCILS
    
    do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
      td2(i,j,k) = ux2(i,j,k) * uy2(i,j,k)
      te2(i,j,k) = uy2(i,j,k) * uy2(i,j,k)
      tf2(i,j,k) = uz2(i,j,k) * uy2(i,j,k)
    enddo

    call dery (tg2,td2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (th2,te2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call dery (ti2,tf2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (td2,ux2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))
    call dery (te2,uy2,sy,x3d_op_dery ,ppy,ysize(1),ysize(2),ysize(3))
    call dery (tf2,uz2,sy,x3d_op_deryp,ppy,ysize(1),ysize(2),ysize(3))

    ! Convective terms of y-pencil in tg2,th2,ti2
    do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
      tg2(i,j,k) = -half*(tg2(i,j,k) + uy2(i,j,k) * td2(i,j,k))
      th2(i,j,k) = -half*(th2(i,j,k) + uy2(i,j,k) * te2(i,j,k))
      ti2(i,j,k) = -half*(ti2(i,j,k) + uy2(i,j,k) * tf2(i,j,k))
    enddo

    ! Add a part of the diffusive term if needed
    if (istret /= 0 .and. xnu /= 0) then
      do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
        tg2(i,j,k) = tg2(i,j,k) - pp4y(j)*td2(i,j,k)
        th2(i,j,k) = th2(i,j,k) - pp4y(j)*te2(i,j,k)
        ti2(i,j,k) = ti2(i,j,k) - pp4y(j)*tf2(i,j,k)
      enddo
    endif
    
    call x3d_transpose_y_to_z(ux2,ux3)
    call x3d_transpose_y_to_z(uy2,uy3)
    call x3d_transpose_y_to_z(uz2,uz3)

    !WORK Z-PENCILS
    do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
      td3(i,j,k) = ux3(i,j,k) * uz3(i,j,k)
      te3(i,j,k) = uy3(i,j,k) * uz3(i,j,k)
      tf3(i,j,k) = uz3(i,j,k) * uz3(i,j,k)
    enddo

    call derz (tg3,td3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (th3,te3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))
    call derz (ti3,tf3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (td3,ux3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (te3,uy3,sz,x3d_op_derzp,zsize(1),zsize(2),zsize(3))
    call derz (tf3,uz3,sz,x3d_op_derz ,zsize(1),zsize(2),zsize(3))

    ! Convective terms of z-pencil in ta3,tb3,tc3
    do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
      ta3(i,j,k) = -half*(tg3(i,j,k) + uz3(i,j,k) * td3(i,j,k))
      tb3(i,j,k) = -half*(th3(i,j,k) + uz3(i,j,k) * te3(i,j,k))
      tc3(i,j,k) = -half*(ti3(i,j,k) + uz3(i,j,k) * tf3(i,j,k))
    enddo 

    ! If needed, compute and add diffusion
    if (xnu /= zero) then

      ! Compute diffusion in td3, te3, tf3
      call derzz(td3,ux3,sz,x3d_op_derzzp,zsize(1),zsize(2),zsize(3))
      call derzz(te3,uy3,sz,x3d_op_derzzp,zsize(1),zsize(2),zsize(3))
      call derzz(tf3,uz3,sz,x3d_op_derzz ,zsize(1),zsize(2),zsize(3))

      ! Add convective and diffusive terms of z-pencil
      do concurrent (k=1:zsize(3), j=1:zsize(2), i=1:zsize(1))
         ta3(i,j,k) = ta3(i,j,k) + xnu*td3(i,j,k)
         tb3(i,j,k) = tb3(i,j,k) + xnu*te3(i,j,k)
         tc3(i,j,k) = tc3(i,j,k) + xnu*tf3(i,j,k)
      enddo

    endif

    ! Send z-rhs (ta3,tb3,tc3) to y-pencil (td2,te2,tf2)
    call x3d_transpose_z_to_y(ta3,td2)
    call x3d_transpose_z_to_y(tb3,te2)
    call x3d_transpose_z_to_y(tc3,tf2)

    ! If needed, compute and add diffusion
    if (xnu /= 0) then

      ! Compute diffusion in ta2, tb2 and tc2
      call deryy(ta2,ux2,sy,x3d_op_deryyp,ysize(1),ysize(2),ysize(3))
      call deryy(tb2,uy2,sy,x3d_op_deryy ,ysize(1),ysize(2),ysize(3))
      call deryy(tc2,uz2,sy,x3d_op_deryyp,ysize(1),ysize(2),ysize(3))

      ! Add convective and diffusive terms of y-pencil
      if (istret /= 0) then
        ! In this case, a part of the y-diffusive term was added before
        do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))                                         
          tg2(i,j,k) = tg2(i,j,k) + xnu*ta2(i,j,k)*pp2y(j)                                               
          th2(i,j,k) = th2(i,j,k) + xnu*tb2(i,j,k)*pp2y(j)                                               
          ti2(i,j,k) = ti2(i,j,k) + xnu*tc2(i,j,k)*pp2y(j)                                               
        enddo
      else
        do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
          tg2(i,j,k) = tg2(i,j,k) + xnu*ta2(i,j,k)
          th2(i,j,k) = th2(i,j,k) + xnu*tb2(i,j,k)
          ti2(i,j,k) = ti2(i,j,k) + xnu*tc2(i,j,k)
        enddo
      endif

    endif

    ! Combine y-rhs with z-rhs
    do concurrent (k=1:ysize(3), j=1:ysize(2), i=1:ysize(1))
      tg2(i,j,k) = tg2(i,j,k) + td2(i,j,k)
      th2(i,j,k) = th2(i,j,k) + te2(i,j,k)
      ti2(i,j,k) = ti2(i,j,k) + tf2(i,j,k)
    enddo

    !WORK X-PENCILS
    call x3d_transpose_y_to_x(tg2,ta1)
    call x3d_transpose_y_to_x(th2,tb1)
    call x3d_transpose_y_to_x(ti2,tc1) !diff+conv. terms

    ! If needed, compute and add diffusion
    if (xnu /= 0) then

      ! Compute diffusion in td1, te1, tf1
      call derxx(td1,ux1,sx,x3d_op_derxx ,xsize(1),xsize(2),xsize(3))
      call derxx(te1,uy1,sx,x3d_op_derxxp,xsize(1),xsize(2),xsize(3))
      call derxx(tf1,uz1,sx,x3d_op_derxxp,xsize(1),xsize(2),xsize(3))

      ! Add convective and diffusive terms of x-pencil
      do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
        dux1(i,j,k,1) = dux1(i,j,k,1) + xnu*td1(i,j,k)
        duy1(i,j,k,1) = duy1(i,j,k,1) + xnu*te1(i,j,k)
        duz1(i,j,k,1) = duz1(i,j,k,1) + xnu*tf1(i,j,k)
      enddo
    endif

    !FINAL SUM: y and z rhs + x rhs
    do concurrent (k=1:xsize(3), j=1:xsize(2), i=1:xsize(1))
      dux1(i,j,k,1) = dux1(i,j,k,1) + ta1(i,j,k)
      duy1(i,j,k,1) = duy1(i,j,k,1) + tb1(i,j,k)
      duz1(i,j,k,1) = duz1(i,j,k,1) + tc1(i,j,k)
    enddo

  end subroutine momentum_rhs_eq
  !############################################################################
  !############################################################################
end module transeq
