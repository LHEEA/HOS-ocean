MODULE RF_solution
!
! This module computes the solution of Rienecker & Fenton (J.F.M. Vol.104, 1981)
! It uses a set of coefficients describing eta and phi obtained from a specific R&F program
! from LHEEA Lab., ECN (Pierre Ferrant & Félicien Bonnefoy)
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
!
IMPLICIT NONE
!
TYPE RF_data
   ! see the Rienecker and Fenton solver (eg wave_rf_hinf_N12.f90)
   CHARACTER(LEN=100)                	:: file_name
   INTEGER                           	:: N_eta, N_phi
   REAL(RP)                          	:: ka, lambda, k, T, C, CC, CPH, H, omega
   REAL(RP)                          	:: depth, Q, R
   LOGICAL                           	:: inf_depth
   REAL(RP), ALLOCATABLE, DIMENSION(:)  :: A, B
   !
   INTEGER                           	:: N_space          ! number of points
   REAL(RP), ALLOCATABLE, DIMENSION(:)  :: x, eta, phis, W  ! values
   ! building in space on N>=N_space, then moving to Fourier, keeping the N_space first
   ! modes and back to space on N_space
   REAL(RP), ALLOCATABLE, DIMENSION(:)  :: eta_dealiased, phis_dealiased, W_dealiased
!    ! building in space on N_space from the N_space first modes of the RF data
!    REAL(RP), POINTER, DIMENSION(:)   :: eta_dealiased2, phis_dealiased2, W_dealiased2
   !
   INTEGER                           	:: N_space_press
   REAL(RP), ALLOCATABLE, DIMENSION(:)  :: pressure
   REAL(RP), ALLOCATABLE, DIMENSION(:,:):: coord_press
END TYPE RF_data
!
!
!
CONTAINS
!
!
!
SUBROUTINE read_RF_data(file_name, RF_obj, inf_depth)
!
IMPLICIT NONE
!
! Input
CHARACTER(LEN=100), INTENT(IN) :: file_name
TYPE(RF_data), INTENT(OUT)     :: RF_obj
LOGICAL, INTENT(IN)            :: inf_depth
! Local
INTEGER                        :: i_loop, i_tmp
!
! Data file name from the Rienecker and Fenton solver (eg wave_rf_hinf_N12.f90)
RF_obj%file_name = TRIM(file_name)
! Whether the depth is finite or not
RF_obj%inf_depth = inf_depth
!
! Reading the data from the file
OPEN(UNIT=967,FILE=TRIM(RF_obj%file_name))
READ(967,*) RF_obj%lambda, RF_obj%H, RF_obj%k, RF_obj%T, RF_obj%C, RF_obj%CPH, RF_obj%CC, RF_obj%N_phi, RF_obj%N_eta
!
! Steepness
RF_obj%ka = RF_obj%k * RF_obj%H * 0.5_rp
!
! Frequency
RF_obj%omega = TWOPI / RF_obj%T
!
IF(ALLOCATED(RF_obj%A)) DEALLOCATE(RF_obj%A,RF_obj%B)
ALLOCATE(RF_obj%A(RF_obj%N_eta+1), RF_obj%B(RF_obj%N_phi+1))
!
! Coefficients for the potential phi
DO i_loop = 1,RF_obj%N_phi+1
   READ(967,*) i_tmp, RF_obj%B(i_loop)
END DO
! Coefficients for the free surface elevation eta
DO i_loop = 1,RF_obj%N_eta+1
   READ(967,*) i_tmp, RF_obj%A(i_loop)
END DO
!
CLOSE(967)
!
END SUBROUTINE read_RF_data
!
!
!
SUBROUTINE build_RF_reference(RF_obj, N_space, N_lambda)
!
use, intrinsic :: iso_c_binding
!
IMPLICIT NONE
!
! Input
TYPE(RF_data), INTENT(INOUT)   :: RF_obj
INTEGER, INTENT(IN)            :: N_space
INTEGER, INTENT(IN), OPTIONAL  :: N_lambda
! Local
INTEGER                             :: i_loop, j_loop, N_work
REAL(RP)                            :: delx, kx, keta, jk
REAL(RP), DIMENSION(:), ALLOCATABLE :: x_w, eta_w, phis_w, W_w
REAL(RP), DIMENSION(:), ALLOCATABLE :: TF, TF_w
! for FFTW transforms
REAL(RP), DIMENSION(:), ALLOCATABLE    :: in_RF, in_RF2, s_2_f_RF, f_2_s_RF2
COMPLEX(CP), DIMENSION(:), ALLOCATABLE :: out_RF, out_RF2
type(C_PTR) :: plan_R2C_RF, plan_C2R_RF2
!
! Requested number of points in space
RF_obj%N_space = N_space
!
ALLOCATE(RF_obj%x(RF_obj%N_space))
ALLOCATE(RF_obj%phis(RF_obj%N_space), RF_obj%eta(RF_obj%N_space), RF_obj%W(RF_obj%N_space))
!
IF (PRESENT(N_lambda)) THEN
   ! Building the solution on several wavelengths
   delx = N_lambda * RF_obj%lambda / RF_obj%N_space
ELSE
   ! Building the solution on one wavelength
   delx = RF_obj%lambda / RF_obj%N_space
END IF
! Position vector
DO i_loop = 1, RF_obj%N_space
   RF_obj%x(i_loop) = REAL(i_loop - 1, RP) * delx
END DO
!
! Direct solution on N_space points from the modes of phi and eta
! (aliased if N_space < N_phi or N_eta)
N_work = RF_obj%N_space
ALLOCATE(x_w(N_work))
ALLOCATE(phis_w(N_work), eta_w(N_work), W_w(N_work))
!
x_w = RF_obj%x
!
DO i_loop = 1, N_work
   ! temp
   kx = RF_obj%k * x_w(i_loop)
   ! Free surface elevation
   eta_w(i_loop) = RF_obj%A(1) * 0.5_rp
   DO j_loop = 1, RF_obj%N_eta
      eta_w(i_loop) = eta_w(i_loop) + RF_obj%A(j_loop+1) * COS(j_loop*kx)
   END DO
   ! temp
   keta = RF_obj%k * eta_w(i_loop)
   ! Free surface potential
   phis_w(i_loop) = (RF_obj%C + RF_obj%B(1)) * x_w(i_loop)
   IF (RF_obj%inf_depth) THEN
      DO j_loop = 1, RF_obj%N_phi
         phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
      END DO
   ELSE
      DO j_loop = 1, RF_obj%N_phi
        IF (j_loop*RF_obj%k < 500.0_rp) THEN
           phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) &
           	* COSH(j_loop*(keta+RF_obj%k))/COSH(j_loop*RF_obj%k)
        ELSE
           phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
        END IF
      END DO
   END IF
   ! Vertical velocity on the free surface
   W_w(i_loop) = 0.0_rp
   IF (RF_obj%inf_depth) THEN
      DO j_loop = 1, RF_obj%N_phi
         W_w(i_loop) = W_w(i_loop) + REAL(j_loop, RP) * RF_obj%k * RF_obj%B(j_loop+1) &
         	* SIN(REAL(j_loop, RP)*kx) * EXP(REAL(j_loop, RP)*keta)
      END DO
   ELSE
      DO j_loop = 1, RF_obj%N_phi
         jk = REAL(j_loop, RP) * RF_obj%k
         IF (jk < 500.0_rp) THEN
            W_w(i_loop) = W_w(i_loop) + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * SINH(REAL(j_loop, RP)*keta+jk)/COSH(jk)
         ELSE
            W_w(i_loop) = W_w(i_loop) + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(REAL(j_loop, RP)*keta)
         END IF
      END DO
   END IF
END DO
!
RF_obj%eta  = eta_w
RF_obj%phis = phis_w
RF_obj%W    = W_w

! Dealiased solution on N_space points from the modes of phi and eta
! building in space on N>=N_space, then moving to Fourier, keeping the N_space first
! modes and back to space on N_space
IF (RF_obj%N_eta > RF_obj%N_space/2) THEN
   N_work = 2 * 2**nextpow2(RF_obj%N_eta)
   !
   DEALLOCATE(x_w, phis_w, eta_w, W_w)
   !
   ALLOCATE(x_w(N_work))
   ALLOCATE(phis_w(N_work), eta_w(N_work), W_w(N_work))
    !
    ! Create plan and other stuff for FFTW tranforms
    !
    ALLOCATE(in_RF(N_work), out_RF(N_work/2+1), s_2_f_RF(N_work/2+1))
    ALLOCATE(in_RF2(N_space), out_RF2(N_space/2+1), f_2_s_RF2(N_space/2+1))

    CALL dfftw_plan_dft_r2c_1d(plan_R2C_RF,  N_work,        in_RF,        out_RF,          64) !0 for FFTW_MEASURE or 64 for FFTW_ESTIMATE
    CALL dfftw_plan_dft_c2r_1d(plan_C2R_RF2, N_space,       out_RF2,      in_RF2,          64) !0 for FFTW_MEASURE or 64 for FFTW_ESTIMATE

    s_2_f_RF(1)            = 1.0_rp / REAL(N_work,RP)
    s_2_f_RF(2:N_work/2+1) = 2.0_rp / REAL(N_work,RP)
    s_2_f_RF(N_work/2+1)   = 1.0_rp / REAL(N_work,RP) ! N_work is even (see definition previously)

    f_2_s_RF2(1)             = 1.0_rp
    f_2_s_RF2(2:N_space/2+1) = 0.5_rp
    IF(iseven(N_space)) f_2_s_RF2(N_space/2+1)   = 1.0_rp
   !
   IF (PRESENT(N_lambda)) THEN
      delx = N_lambda * RF_obj%lambda / N_work
   ELSE
      delx = RF_obj%lambda / N_work
   END IF
   DO i_loop = 1, N_work
      x_w(i_loop) = REAL(i_loop - 1, RP) * delx
   END DO
   !
   DO i_loop = 1, N_work
      ! temp
      kx = RF_obj%k * x_w(i_loop)
      ! Free surface elevation
      eta_w(i_loop) = RF_obj%A(1) * 0.5_rp
      DO j_loop = 1, RF_obj%N_eta
         eta_w(i_loop) = eta_w(i_loop) + RF_obj%A(j_loop+1) * COS(j_loop*kx)
      END DO
      ! temp
      keta = RF_obj%k * eta_w(i_loop)
      ! Free surface potential
      phis_w(i_loop) = (RF_obj%C + RF_obj%B(1)) * x_w(i_loop)
      IF (RF_obj%inf_depth) THEN
         DO j_loop = 1, RF_obj%N_phi
            phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
         END DO
      ELSE
         DO j_loop = 1, RF_obj%N_phi
            IF (j_loop*RF_obj%k < 500.0_rp) THEN
               phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) &
               	* COSH(j_loop*(keta+RF_obj%k))/COSH(j_loop*RF_obj%k)
            ELSE
               phis_w(i_loop) = phis_w(i_loop) + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
            END IF
         END DO
      END IF
      ! Vertical velocity on the free surface
      W_w(i_loop) = 0.0_rp
      IF (RF_obj%inf_depth) THEN
         DO j_loop = 1, RF_obj%N_phi
            W_w(i_loop) = W_w(i_loop) + REAL(j_loop, RP) * RF_obj%k * RF_obj%B(j_loop+1) &
            	* SIN(REAL(j_loop, RP)*kx) * EXP(REAL(j_loop, RP)*keta)
         END DO
      ELSE
         DO j_loop = 1, RF_obj%N_phi
            jk = REAL(j_loop, RP) * RF_obj%k
            IF (jk < 500.0_rp) THEN
               W_w(i_loop) = W_w(i_loop) + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * SINH(REAL(j_loop, RP)*keta+jk)/COSH(jk)
            ELSE
               W_w(i_loop) = W_w(i_loop) + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(REAL(j_loop, RP)*keta)
            END IF
         END DO
      END IF
   END DO
   !
   ALLOCATE(TF_w(N_work), TF(RF_obj%N_space))
    !
    in_RF   = eta_w(1:N_work)
    call dfftw_execute_dft_r2c(plan_R2C_RF,in_RF,out_RF)
    !
    out_RF2(1:N_space/2+1) = out_RF(1:N_space/2+1) * s_2_f_RF(1:N_space/2+1)
    IF (iseven(N_space)) THEN ! last mode of eta will be a cosine without corresponding sine for phis and W
    	out_RF2(N_space/2+1)   = ((0.0_rp, 0.0_rp))
    ENDIF
    !
    ! Back Transform
    !
    DEALLOCATE(eta_w)
    ALLOCATE(eta_w(RF_obj%N_space))
    !
    out_RF2 = out_RF2 * f_2_s_RF2
    call dfftw_execute_dft_c2r(plan_C2R_RF2,out_RF2,in_RF2)
    eta_w = in_RF2
    !
    !
    !
    in_RF   = phis_w(1:N_work)
    call dfftw_execute_dft_r2c(plan_R2C_RF,in_RF,out_RF)
    !
    out_RF2(1:N_space/2+1) = out_RF(1:N_space/2+1) * s_2_f_RF(1:N_space/2+1)
    IF (iseven(N_space)) THEN ! last mode of eta will be a cosine without corresponding sine for phis and W
    	out_RF2(N_space/2+1)   = ((0.0_rp, 0.0_rp))
    ENDIF
    !
    ! Back Transform
    !
    DEALLOCATE(phis_w)
    ALLOCATE(phis_w(RF_obj%N_space))
    !
    out_RF2 = out_RF2 * f_2_s_RF2
    call dfftw_execute_dft_c2r(plan_C2R_RF2,out_RF2,in_RF2)
    phis_w = in_RF2
    !
    !
    !
    in_RF   = W_w(1:N_work)
    call dfftw_execute_dft_r2c(plan_R2C_RF,in_RF,out_RF)
    !
    out_RF2(1:N_space/2+1) = out_RF(1:N_space/2+1) * s_2_f_RF(1:N_space/2+1)
    IF (iseven(N_space)) THEN ! last mode of eta will be a cosine without corresponding sine for phis and W
    	out_RF2(N_space/2+1)   = ((0.0_rp, 0.0_rp))
    ENDIF
    !
    ! Back Transform
    !
    DEALLOCATE(W_w)
    ALLOCATE(W_w(RF_obj%N_space))
    !
    out_RF2 = out_RF2 * f_2_s_RF2
    call dfftw_execute_dft_c2r(plan_C2R_RF2,out_RF2,in_RF2)
    W_w = in_RF2
    !
    ! Deallocate temporary variables for FFTW tranforms + destroy plans
    !
    CALL dfftw_destroy_plan(plan_R2C_RF)
    CALL dfftw_destroy_plan(plan_C2R_RF2)
    DEALLOCATE(in_RF, out_RF, in_RF2, out_RF2, s_2_f_RF, f_2_s_RF2, TF_w, TF)
    CALL dfftw_cleanup
END IF
! GD : modif feb. 2013 put this outside the IF
ALLOCATE(RF_obj%phis_dealiased(RF_obj%N_space), RF_obj%eta_dealiased(RF_obj%N_space), &
RF_obj%W_dealiased(RF_obj%N_space))
!
RF_obj%eta_dealiased  = eta_w
RF_obj%phis_dealiased = phis_w
RF_obj%W_dealiased    = W_w
!
! Deallocate
DEALLOCATE(x_w, phis_w, eta_w, W_w)
!
END SUBROUTINE build_RF_reference
!
!
!
SUBROUTINE build_RF_pressure(RF_obj, coord, time, height_threshold)
!
IMPLICIT NONE
!
! Input
TYPE(RF_data), INTENT(INOUT)   :: RF_obj
REAL(RP), DIMENSION(:,:)       :: coord
REAL(RP)                       :: time, height_threshold
! Local
INTEGER                        :: i_loop, j_loop, N_work
REAL(RP)                       :: kx, kz, eta_w
REAL(RP), DIMENSION(:), ALLOCATABLE    :: phit_w, phix_w, phiz_w
!
RF_obj%N_space_press = SIZE(coord,1)
ALLOCATE(RF_obj%pressure(RF_obj%N_space_press), RF_obj%coord_press(RF_obj%N_space_press, 2))
!
RF_obj%coord_press = coord
!
N_work = RF_obj%N_space_press
ALLOCATE(phit_w(N_work), phix_w(N_work), phiz_w(N_work))
!
DO i_loop = 1, N_work
   ! temp
   kx = RF_obj%k * RF_obj%coord_press(i_loop, 1) - RF_obj%omega * time
   kz = RF_obj%k * RF_obj%coord_press(i_loop, 2)
   ! Free surface elevation
   eta_w = RF_obj%A(1) * 0.5_rp
   DO j_loop = 1, RF_obj%N_eta
      eta_w = eta_w + RF_obj%A(j_loop+1) * COS(j_loop*kx)
   END DO
   IF ( ((RF_obj%coord_press(i_loop, 2)>=0.0_rp) .AND. (RF_obj%coord_press(i_loop, 2) <= &
         (1.0_rp + height_threshold) * eta_w)) .OR. &
         ((RF_obj%coord_press(i_loop, 2)<0.0_rp) .AND. (RF_obj%coord_press(i_loop, 2) <= &
         (1.0_rp - height_threshold) * eta_w)) ) THEN
      ! Derivatives of the potential
      phit_w(i_loop) = 0.0_rp
      phix_w(i_loop) = RF_obj%C + RF_obj%B(1)
      phiz_w(i_loop) = 0.0_rp
      IF (RF_obj%inf_depth) THEN
         DO j_loop = 1, RF_obj%N_phi
            phit_w(i_loop) = phit_w(i_loop) - j_loop * RF_obj%omega * RF_obj%B(j_loop+1) * COS(j_loop*kx) * EXP(j_loop*kz)
            phix_w(i_loop) = phix_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * COS(j_loop*kx) * EXP(j_loop*kz)
            phiz_w(i_loop) = phiz_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*kz)
         END DO
      ELSE
         DO j_loop = 1, RF_obj%N_phi
           IF (j_loop*RF_obj%k < 500.0_rp) THEN
              phit_w(i_loop) = phit_w(i_loop) - j_loop * RF_obj%omega * RF_obj%B(j_loop+1) * COS(j_loop*kx) &
              	* COSH(j_loop*(kz+RF_obj%k))/COSH(j_loop*RF_obj%k)
              phix_w(i_loop) = phix_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * COS(j_loop*kx) &
              	* COSH(j_loop*(kz+RF_obj%k))/COSH(j_loop*RF_obj%k)
              phiz_w(i_loop) = phiz_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * SIN(j_loop*kx) &
              	* SINH(j_loop*(kz+RF_obj%k))/COSH(j_loop*RF_obj%k)
           ELSE
               phit_w(i_loop) = phit_w(i_loop) - j_loop * RF_obj%omega * RF_obj%B(j_loop+1) * COS(j_loop*kx) * EXP(j_loop*kz)
               phix_w(i_loop) = phix_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * COS(j_loop*kx) * EXP(j_loop*kz)
               phiz_w(i_loop) = phiz_w(i_loop) + j_loop * RF_obj%k     * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*kz)
           END IF
         END DO
      END IF
      RF_obj%pressure(i_loop)  = - 1.0_rp * RF_obj%coord_press(i_loop, 2) - phit_w(i_loop) &
      	- 0.5_rp * (phix_w(i_loop) * phix_w(i_loop) + phiz_w(i_loop) * phiz_w(i_loop))
   ELSE
      phit_w(i_loop) = 0.0_rp
      phix_w(i_loop) = 0.0_rp
      phiz_w(i_loop) = 0.0_rp
      RF_obj%pressure(i_loop)  = 0.0_rp
   END IF
END DO
!
END SUBROUTINE build_RF_pressure
!
!
!
SUBROUTINE display_RF_data(RF_obj)
!
IMPLICIT NONE
!
! Input
TYPE(RF_data), INTENT(IN)     :: RF_obj
!
WRITE(*,'(2A)')        'File name:      ', TRIM(RF_obj%file_name)
WRITE(*,'(A,ES25.16)') 'Wave length:    ', RF_obj%lambda
WRITE(*,'(A,ES25.16)') 'Wave height:    ', RF_obj%H
WRITE(*,'(A,ES25.16)') 'Wave number:    ', RF_obj%k
WRITE(*,'(A,ES25.16)') 'Phase velocity: ', RF_obj%C
WRITE(*,'(A,ES25.16)') 'Current?:       ', RF_obj%CPH
WRITE(*,'(A,ES25.16)') 'Current?:       ', RF_obj%CC
WRITE(*,'(A,I3)')      'Modes for phi:  ', RF_obj%N_phi
WRITE(*,'(A,I3)')      'Modes for eta:  ', RF_obj%N_eta
WRITE(*,'(A,L1)')      'Infinite depth: ', RF_obj%inf_depth
!
IF (ALLOCATED(RF_obj%A)) THEN
   WRITE(*,'(A,ES25.16)') 'Amplitude 0 for eta: ', RF_obj%A(1)
   WRITE(*,'(A,ES25.16)') 'Amplitude 1 for eta: ', RF_obj%A(2)
   WRITE(*,'(A,ES25.16)') 'Amplitude 2 for eta: ', RF_obj%A(3)
END IF
!
IF (ALLOCATED(RF_obj%B)) THEN
   WRITE(*,'(A,ES25.16)') 'Amplitude 0 for phi: ', RF_obj%B(1)
   WRITE(*,'(A,ES25.16)') 'Amplitude 1 for phi: ', RF_obj%B(2)
   WRITE(*,'(A,ES25.16)') 'Amplitude 2 for phi: ', RF_obj%B(3)
END IF
!
WRITE(*,'(A,I3)')      'Points for eta, phis and W: ', RF_obj%N_space
!
END SUBROUTINE display_RF_data
!
!
!
FUNCTION nextpow2(N)
!
INTEGER :: N, nextpow2, P, twopP
!
P=0
twopP = 1
DO WHILE (twopP < N)
   P=P+1
   twopP = twopP * 2
END DO
nextpow2 = P
!
END FUNCTION nextpow2
!
!
!
END MODULE RF_solution
