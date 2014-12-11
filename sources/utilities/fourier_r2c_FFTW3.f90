MODULE fourier_r2c
!
! This module contains all necessary features for the use of FFTW 3.3.4 in HOS-ocean
! with comprehensive interface and normalization
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Copyright (C) 2014 - LHEEA Lab., Ecole Centrale de Nantes, UMR CNRS 6598
!
!    This program is part of HOS-ocean
!
!    HOS-ocean is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
USE type
USE variables_3D
USE, INTRINSIC :: iso_c_binding
!
IMPLICIT NONE
!
include 'fftw3.f03'
 !
 ! FFTW module procedure for original and dealiased domains
 !
 INTERFACE SPACE_2_FOURIER
    MODULE PROCEDURE FFTW_R2C
 END INTERFACE SPACE_2_FOURIER
 !
 INTERFACE FOURIER_2_SPACE
    MODULE PROCEDURE FFTW_C2R
 END INTERFACE FOURIER_2_SPACE
 !
 INTERFACE SPACE_2_FOURIER_big
    MODULE PROCEDURE FFTW_R2C_big
 END INTERFACE SPACE_2_FOURIER_big
 !
 INTERFACE FOURIER_2_SPACE_big
    MODULE PROCEDURE FFTW_C2R_big
 END INTERFACE FOURIER_2_SPACE_big
!
type(C_PTR)                              :: plan_R2C, plan_C2R
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: s_2_f, f_2_s
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: in
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: out
!
type(C_PTR)                              :: plan_R2C_big, plan_C2R_big
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: s_2_f_big, f_2_s_big
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: in_big
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: out_big
!
REAL(RP)                      :: twoon1, twooNd1, twoon2, twooNd2, oneon1, oneon2, oneoNd1, oneoNd2
!
!
!
CONTAINS
!
!
!
SUBROUTINE Fourier_ini(library)
!
! Initializes the Fourier Transforms
! default library: FFTW-3.3 (library=3)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN), OPTIONAL :: library
! Local variables
INTEGER :: lib, flag
!
! Input argument
IF (PRESENT(library)) THEN
   lib = library
ELSE
   ! default library is FFTW3
   lib = 3
END IF
!
! Memory allocation
ALLOCATE(in(n1,n2), out(n1o2p1,n2), s_2_f(n1o2p1,n2), f_2_s(n1o2p1,n2))
ALLOCATE(in_big(Nd1,Nd2), out_big(Nd1o2p1,Nd2), s_2_f_big(Nd1o2p1,Nd2), f_2_s_big(Nd1o2p1,Nd2))
!
! Initializing the FFT library
SELECT CASE (lib)
   ! FIXME: remove this need (see RF_solution...)
   ! FIXME: everything is done with FFTW now
   CASE (1)
    print*,'Not taken in charge any more'
    STOP
   CASE (3)
      ! FFTW-3.3 library
      ! FFTW_ESTIMATE is deterministic (useful for debug/testing) but not optimal
      flag = FFTW_ESTIMATE !FFTW_PATIENT ! FFTW_MEASURE ! FFTW_EXHAUSTIVE
      IF (n2 == 1) THEN
        ! 1D FFTs
        CALL dfftw_plan_dft_r2c_1d(plan_R2C, n1,        in,        out,          flag)
        CALL dfftw_plan_dft_c2r_1d(plan_C2R, n1,        out,       in,           flag)
        !
        CALL dfftw_plan_dft_r2c_1d(plan_R2C_big, Nd1,   in_big,     out_big,     flag)
        CALL dfftw_plan_dft_c2r_1d(plan_C2R_big, Nd1,   out_big,     in_big,     flag)
      ELSE
        ! 2D FFTs
        CALL dfftw_plan_dft_r2c_2d(plan_R2C, n1, n2,  in,     out,     flag)
        CALL dfftw_plan_dft_c2r_2d(plan_C2R, n1, n2,  out,     in,     flag)
        !
        CALL dfftw_plan_dft_r2c_2d(plan_R2C_big, Nd1, Nd2,  in_big,     out_big,     flag)
        CALL dfftw_plan_dft_c2r_2d(plan_C2R_big, Nd1, Nd2,  out_big,     in_big,     flag)
      END IF
   CASE DEFAULT
      STOP 'Unknown DFT library in Fourier_ini'
END SELECT
!
! Evaluating conversion coefficients
! when computing a space to Fourier Transform
! on the free surface
s_2_f(1,1)           = 1.0_rp / REAL(n1*n2,RP)
s_2_f(2:n1o2p1,1)    = 2.0_rp / REAL(n1*n2,RP)
!
s_2_f(1,2:n2)        = 1.0_rp / REAL(n1*n2,RP)
s_2_f(2:n1o2p1,2:n2) = 2.0_rp / REAL(n1*n2,RP)
!
IF(iseven(n1)) THEN
    s_2_f(n1o2p1,1:n2) = 1.0_rp / REAL(n1*n2,RP)
ENDIF
!
! for 'dealiased' domain
!
s_2_f_big(1,1)             = 1.0_rp / REAL(Nd1*Nd2,RP)
s_2_f_big(2:Nd1o2p1,1)     = 2.0_rp / REAL(Nd1*Nd2,RP)
!
s_2_f_big(1,2:Nd2)         = 1.0_rp / REAL(Nd1*Nd2,RP)
s_2_f_big(2:Nd1o2p1,2:Nd2) = 2.0_rp / REAL(Nd1*Nd2,RP)
!
IF(iseven(Nd1)) THEN
    s_2_f_big(Nd1o2p1,1:Nd2) = 1.0_rp / REAL(Nd1*Nd2,RP)
ENDIF
!
!  when computing a Fourier to space Transform
!  on the free surface
f_2_s(1,1)           = 1.0_rp
f_2_s(2:n1o2p1,1)    = 0.5_rp
!
f_2_s(1,2:n2)        = 1.0_rp
f_2_s(2:n1o2p1,2:n2) = 0.5_rp
!
IF(iseven(n1)) THEN
    f_2_s(n1o2p1,1:n2) = 1.0_rp
ENDIF
!
! for 'dealiased' domain
!
f_2_s_big(1,1)             = 1.0_rp
f_2_s_big(2:Nd1o2p1,1)     = 0.5_rp
!
f_2_s_big(1,2:Nd2)         = 1.0_rp
f_2_s_big(2:Nd1o2p1,2:Nd2) = 0.5_rp
!
IF(iseven(Nd1)) THEN
    f_2_s_big(Nd1o2p1,1:Nd2) = 1.0_rp
ENDIF
!
END SUBROUTINE Fourier_ini
!
!
!
SUBROUTINE Fourier_end(library)
!
! Desallocate the Fourier Transforms
! default library: FFTW-3.3 (library=3)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN), OPTIONAL :: library
! Local variables
INTEGER :: lib
!
! Input argument
IF (PRESENT(library)) THEN
   lib = library
ELSE
   ! default library is FFTW3
   lib = 3
END IF
!
! End the FFT library
SELECT CASE (lib)
   ! FIXME: remove this need (see RF_solution...)
   CASE (1)
        ! Nothing to do
   CASE (3)
      ! FFTW-3.3 library
      !
      CALL dfftw_destroy_plan(plan_R2C)
      CALL dfftw_destroy_plan(plan_C2R)
      CALL dfftw_destroy_plan(plan_R2C_big)
      CALL dfftw_destroy_plan(plan_C2R_big)
      !
      CALL fftw_cleanup
      !
   CASE DEFAULT
      STOP 'Unknown DFT library in Fourier_ini'
END SELECT
!
! Memory allocation
DEALLOCATE(in, out, s_2_f, f_2_s)
DEALLOCATE(in_big, out_big, s_2_f_big, f_2_s_big)
!
END SUBROUTINE Fourier_end
!
! Define the effective transormations
!
SUBROUTINE FFTW_R2C(x,y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m1,m2), INTENT(IN)         :: x
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(OUT) :: y
!
in   = x(1:n1,1:n2)
!
call dfftw_execute_dft_r2c(plan_R2C,in,out)
!
y(1:n1o2p1,1:n2) = out * s_2_f
!
END SUBROUTINE FFTW_R2C
!
SUBROUTINE FFTW_C2R(y,x)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(m1,m2), INTENT(OUT)       :: x
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: y
!
out   = y(1:n1o2p1,1:n2) * f_2_s
!
call dfftw_execute_dft_c2r(plan_C2R,out,in)
!
x(1:n1,1:n2) = in
!
END SUBROUTINE FFTW_C2R
!
! Dealiased domains
!
SUBROUTINE FFTW_R2C_big(x,y)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(md1,md2), INTENT(IN)         :: x
COMPLEX(CP), DIMENSION(md1o2p1,md2), INTENT(OUT) :: y
!
in_big   = x(1:Nd1,1:Nd2)
!
call dfftw_execute_dft_r2c(plan_R2C_big,in_big,out_big)
!
y(1:Nd1o2p1,1:Nd2) = out_big * s_2_f_big
!
END SUBROUTINE FFTW_R2C_big
!
SUBROUTINE FFTW_C2R_big(y,x)
!
IMPLICIT NONE
!
! Input variables, output variable
REAL(RP), DIMENSION(md1,md2), INTENT(OUT)       :: x
COMPLEX(CP), DIMENSION(md1o2p1,md2), INTENT(IN) :: y
!
out_big   = y(1:Nd1o2p1,1:Nd2) * f_2_s_big
!
call dfftw_execute_dft_c2r(plan_C2R_big,out_big,in_big)
!
x(1:Nd1,1:Nd2) = in_big
!
END SUBROUTINE FFTW_C2R_big
!
!
!
END MODULE fourier_r2c
