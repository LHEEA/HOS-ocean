MODULE initial_condition
!
! This module contains all subroutines necessary for initialization of the computation
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
USE variables_3d
USE fourier_r2c
USE energy_calc
USE RF_solution
USE runge_kutta
USE Bivar
USE input_HOS
USE linear_wave
!
IMPLICIT NONE
!
TYPE(RF_data) :: RF_obj
INTEGER :: n_lambda_x, n_lambda_y
! FIXME : is it necessary to have these?
INTEGER, PARAMETER ::  ifreq = 40
INTEGER, PARAMETER ::  ithet = 36
REAL(RP), DIMENSION(ifreq) :: freq_ww
REAL(RP), DIMENSION(ithet) :: thet_ww
REAL(RP), DIMENSION(ifreq,ithet):: phi_ww
!
!
!
CONTAINS
!
!
!
SUBROUTINE initiate_parameters()
!
IMPLICIT NONE
!
CHARACTER(LEN=100) :: filename
REAL(RP)           :: kp_real
!
! Checking numerical parameters
CALL check_range(M,'M ')
CALL check_range(p1,'p1')
CALL check_range(p2,'p2')
!
IF (p1 > M) STOP 'initiate_parameters: p1 is higher than M'
IF (p2 > M) STOP 'initiate_parameters: p2 is higher than M'
!
! Setting numerical parameters
!
! Check for infinite depth
IF (depth < 0.0_rp) THEN
   depth = 1.0e15_rp
END IF
!
SELECT CASE (i_case)
   CASE (1,2,21)
      WRITE(*,'(A,I3)') 'initiate_parameters: No parameters to set with i_case=',i_case
   CASE(81,82,83,84,809)
      WRITE(*,'(A)') 'initiate_parameters: Fructus case'
      SELECT CASE (i_case)
         CASE (81)
            filename = 'waverf_L628_inf_ka01_N15_30.cof'
         CASE (82)
            filename = 'waverf_L628_inf_ka02_N20_40.cof'
         CASE (83)
            filename = 'waverf_L628_inf_ka03_N25_50.cof'
         CASE (84)
            filename = 'waverf_L628_inf_ka04_N50_100.cof'
         CASE (809)
            filename = 'waverf_L628_inf_ka009_N50_100.cof'
      END SELECT
      ! Depth is infinite in those cases
      depth = 1.0e15_rp
      !
      CALL read_RF_data(filename, RF_obj, .TRUE.)
      !
      n_lambda_x = FLOOR(xlen)
      n_lambda_y = FLOOR(ylen)
      !
      CALL build_RF_reference(RF_obj, MAX(n1 / n_lambda_x,1))
      !
      ! input is non-dimensional
      !
      xlen   = REAL(FLOOR(xlen), RP) * RF_obj%lambda
      ylen   = REAL(FLOOR(ylen), RP) * RF_obj%lambda
      !
      T_stop = T_stop * RF_obj%T
      f_out = f_out / RF_obj%T
      !
      ! transform to dimensional for main program
      !
      IF (ABS(depth-1.0e15_rp) <= epsilon(1.0e15_rp)) THEN ! Infinite depth
      	! Time-scale change only (as if depth=1)
      	T_stop = T_stop/sqrt(grav)
      	f_out  = f_out*sqrt(grav)
      ELSE
      	xlen = xlen*depth
      	ylen = ylen*depth
      	!
      	T_stop = T_stop*sqrt(depth/grav)
      	f_out  = f_out*sqrt(grav/depth)
      ENDIF
   CASE (3,31,32)
      n_lambda_x = FLOOR(xlen)
      n_lambda_y = FLOOR(ylen)
      !
      !xlen   = REAL(FLOOR(xlen), RP) * TWOPI/kp
      !ylen   = REAL(FLOOR(ylen), RP) * TWOPI/kp
      !
      !grav = g ! gravity defined in input file
      !
      ! Everything should be dimensional here
      !
	  E_cible = (0.25_rp*Hs_real)**2
	  !
	  IF(ABS(depth-1.0e15_rp) <= epsilon(1.0e15_rp)) THEN
	  	kp_real = (TWOPI/Tp_real)**2/grav
	  ELSE
	  	kp_real = wave_number_r(1/Tp_real,depth,grav,1.0E-15_rp)
	  ENDIF
	  !
	  xlen   = REAL(FLOOR(xlen), RP) * TWOPI/kp_real
      ylen   = REAL(FLOOR(ylen), RP) * TWOPI/kp_real
      !
      T_stop = T_stop*Tp_real
      f_out  = f_out/Tp_real
      Ta     = Ta*Tp_real
      !
      IF(i_case == 31) CALL read_irreg_f
CASE DEFAULT
   WRITE(*,'(A)') 'initiate_parameters: Unknow case'
END SELECT
!
END SUBROUTINE initiate_parameters
!
!
!
SUBROUTINE initiate(time, eta, phis, RK_param)
!
IMPLICIT NONE
!
REAL(RP), INTENT(INOUT)                 :: time
REAL(RP), DIMENSION(m1,m2), INTENT(OUT) :: eta, phis
TYPE(RK_parameters), INTENT(IN)         :: RK_param
!
INTEGER                         		:: i1, i2,i_initiate_NL
REAL(RP)                        		:: phase, ampli, sig
COMPLEX(CP), DIMENSION(m1o2p1,n2) 		:: da_eta
!
time = 0.0_rp
!
! Elevation
SELECT CASE (i_case)
   CASE (1,9)
      WRITE(*,'(A)') 'initiate: Starting from rest'
      ! Rest
      eta(1:n1,1:n2) = 0.0_rp
   CASE (2,21)
      WRITE(*,'(A)') 'initiate: Starting with a natural mode'
      ! Natural mode, either progressive (2) or stationary (21)
      ! Number of the mode (may be negative for progressive mode)
      i1    = 1
      i2    = 0
      IF (n2 == 1) i2 = 0 ! Making sure we're in 2D
      IF (i_case == 21) THEN ! Stationary case
         i1 = ABS(i1)
      END IF
      ampli = 0.1_rp/L !0.1_rp * xlen_star/ (i1* TWOPI)
      !
      phase = 0.0_rp
      sig   = SIGN(1.0_rp,REAL(i1,RP))
      !
      a_eta(1:n1o2p1,1:n2) = 0.0_rp
      ! First order
      IF (i2 < 0) THEN
         a_eta(ABS(i1)+1,n2-(ABS(i2)+1)+2) = ampli * EXP(i*phase)
      ELSE
         a_eta(ABS(i1)+1,ABS(i2)+1)        = ampli * EXP(i*phase)
      END IF
      ! FIXME : create a case for NL initialization in this case? (2nd, 3rd order)
      !
      CALL fourier_2_space(a_eta, eta)
   CASE (81,82,83,84,809)
      IF (RF_dealias == 1) THEN
         DO i2 = 1, n2
            DO i1 = 0, n_lambda_x-1
               eta(i1*RF_obj%N_space+1:(i1+1)*RF_obj%N_space,i2)  = RF_obj%eta_dealiased(:)
            END DO
         END DO
      ELSE IF (RF_dealias == 0) THEN
         DO i2 = 1, n2
            DO i1 = 0, n_lambda_x-1
               eta(i1*RF_obj%N_space+1:(i1+1)*RF_obj%N_space,i2)  = RF_obj%eta(:)
            END DO
         END DO
      END IF
      !
      IF (iseven(n1)) THEN ! last mode of eta will be a cosine without corresponding sine for phis
         CALL space_2_fourier(eta, a_eta)
         a_eta(n1o2p1,:) = 0.0_cp
         CALL fourier_2_space(a_eta, eta)
      END IF
      !
      IF (n_lambda_x /= 1) THEN
         ! zero forcing to remove BF instabilities (in case L_x=2lambda)
         CALL space_2_fourier(eta, a_eta)
         DO i1=1,n_lambda_x-1
            a_eta(1+i1:n1o2p1:n_lambda_x,:) = 0.0_cp
         END DO
         CALL fourier_2_space(a_eta, eta)
      END IF
   CASE (3,31,32) ! Irregular sea-state
      WRITE(*,*) 'irregular sea-state'
      IF(i_case==3) CALL initiate_irreg(RK_param)
      IF(i_case==31) CALL initiate_irreg_f
      IF(i_case==32) CALL initiate_irreg_i
      !
      CALL fourier_2_space(a_phis, phis)
      CALL fourier_2_space(a_eta, eta)
!      !
!      ! Test pi phase-shift
!      !
!      eta  = - eta
!      phis = - phis
!      CALL space_2_fourier(phis, a_phis)
!      CALL space_2_fourier(eta,a_eta)
CASE DEFAULT
   WRITE(*,'(A)') 'initiate: maybe not correctly initiated'
   eta(1:n1,1:n2) = 0.0_rp
   STOP
END SELECT
!
! Velocity
SELECT CASE (i_case)
   CASE (1,9)
      !  Rest
      phis(1:n1,1:n2) = 0.0_rp
   CASE (2)
      ! Progressive natural mode
      a_phis(1:n1o2p1,1:n2) = 0.0_rp
      ! First order
      IF (i2 < 0)  THEN
         a_phis(ABS(i1)+1,n2-(ABS(i2)+1)+2) = - sig * i * ampli * EXP(i*sig*phase) / omega_n2(i1+1,i2+1)
      ELSE
         a_phis(ABS(i1)+1,ABS(i2)+1)        = - sig * i * ampli * EXP(i*sig*phase) / omega_n2(i1+1,i2+1)
      END IF
      !
      CALL fourier_2_space(a_phis, phis)
      ! FIXME : specific case for NL initialization + CHECK if the following is working?
      ! Change the initialization of a_eta and a_phis
      !
	  i_initiate_NL = 0
      IF(i_initiate_NL == 1)THEN
	    CALL initiate_NL
      ENDIF
      !
   CASE (21)
      !  Stationary natural mode, no velocity
      phis(1:n1,1:n2) = 0.0_rp
   CASE (81,82,83,84,809)
      IF (RF_dealias == 1) THEN
         DO i2 = 1, n2
            DO i1 = 0, n_lambda_x-1
               phis(i1*RF_obj%N_space+1:(i1+1)*RF_obj%N_space,i2) = RF_obj%phis_dealiased(:)
            END DO
         END DO
      ELSE IF (RF_dealias == 0) THEN
         DO i2 = 1, n2
            DO i1 = 0, n_lambda_x-1
               phis(i1*RF_obj%N_space+1:(i1+1)*RF_obj%N_space,i2) = RF_obj%phis(:)
            END DO
         END DO
      END IF
      !
      IF (iseven(n1)) THEN ! useless as cosine amplitude is already zero for phis
         CALL space_2_fourier(phis, a_phis)
         a_phis(n1o2p1,:) = 0.0_cp
         CALL fourier_2_space(a_phis, phis)
      END IF
      !
      IF (n_lambda_x /= 1) THEN
         ! zero forcing to remove BF instabilities (in case L_x=2lambda)
         CALL space_2_fourier(phis, a_phis)
         DO i1=1,n_lambda_x-1
            a_phis(1+i1:n1o2p1:n_lambda_x,:) = 0.0_rp
         END DO
         CALL fourier_2_space(a_phis, phis)
      END IF
   CASE (3,31,32)
      ! Initialisation already done... Nothing to do
 CASE DEFAULT
      phis(1:n1,1:n2) = 0.0_rp
END SELECT
!
! t=0, evaluates the energy
SELECT CASE (i_case)
   CASE (1,9,2,21)
      ! Rest or maximum
      da_eta(1:n1o2p1,1:n2) = 0.0_cp
   CASE(81,82,83,84,809,3,31,32)
      CALL space_2_fourier(eta,  a_eta)
      CALL space_2_fourier(phis, a_phis)
      CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
CASE DEFAULT
      da_eta(1:n1o2p1,1:n2) = 0.0_cp
END SELECT
!
CALL space_2_fourier(eta,  a_eta)
CALL space_2_fourier(phis, a_phis)
E_o = calc_energy(a_eta, a_phis, da_eta)
!
! Display the error on vertical velocity at initial stage
SELECT CASE (i_case)
	CASE(81,82,83,84,809)
	    print*,'erreur W abs.', MAXVAL(ABS(RF_obj%W-phiz(:,1))),(MAXVAL(ABS(RF_obj%W)))
END SELECT
!
END SUBROUTINE initiate
!
!
!
SUBROUTINE check_range(n,name)
!
! This subroutine check the values of parameters p1/p2 or M
! FIXME: With FFTW, maybe not necessary: test the loss of efficiency for large prime numbers
!
IMPLICIT NONE
!
INTEGER :: n
CHARACTER(LEN=2) :: name
!
SELECT CASE (n)
   CASE (1,2,3,4,5,7,8,9,11,14,15,17,19,23,29)
      WRITE(*,'(A,A,A)') 'initiate_parameters: value of ',name,' in correct range'
   CASE (6,10,12,13,16,18,20,21,22,24,25,26,27,28)
      WRITE(*,'(A,A,A)') 'initiate_parameters: forbidden value of ',name,'. cf readme.txt for instructions'
      STOP
   CASE DEFAULT ! n above 30
      WRITE(*,'(A,A,A)') 'initiate_parameters: forbidden value of ',name,'. cf readme.txt for instructions'
      STOP
END SELECT
!
END SUBROUTINE check_range
!
!
!
SUBROUTINE initiate_irreg(RK_param)
!
!-------------------------------------------
!
! Initialisation of irregular multi-directional sea state (linear)
! 	- Frequency spectrum is a JONSWAP
! 	- Directional spreading is extracted from Dysthe's publication
!		- 1/beta*COS(2*pi*theta/(4*beta))**2
! FIXME: Test the implemented directional function extracted from WWIII
! FIXME: create external subroutines for frequency spectrum and directional functions
!
!-------------------------------------------
!
IMPLICIT NONE

REAL(RP) :: sigma, theta, E, Cj, pioxlen, pioylen
REAL(RP) :: test, rnd, angle, angle1, angle2
REAL(RP), DIMENSION(m1o2p1,n2) :: phi_JONSWAP
INTEGER :: iseed,i1,i2,i_initiate_NL
REAL(RP), DIMENSION(4)          :: energy
TYPE(RK_parameters), INTENT(IN)         :: RK_param
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: da_eta,a_eta_temp,a_phi_temp
REAL(RP),DIMENSION(m1o2p1,m2) :: DD_WW
REAL(RP),DIMENSION(ikp) :: ones_k
INTEGER ::  i_dir_JSWP,jj
REAL(RP) :: beta_min,beta_max,frr,fp_w,s_ww,Norm_DD
!
pioxlen = TWOPI/xlen_star
pioylen = TWOPI/ylen_star
! FIXME : add this choice of normalized dir function in input file
! ----- Normalised Dir Function
DD_WW = 0.0_rp
i_dir_JSWP = 0
IF (i_dir_JSWP == 1) THEN
    beta_min = 4.06_rp
    beta_max = -2.34_rp
    !
    ones_k = 1.0_rp
    fp_w = T_out / Tp_real !FIXME: check this definition involving T_out
    DO i1=1,n1o2p1
        DO i2=1,n2
            Norm_DD = 0.0_rp
            frr = MIN(2.5_rp,omega_n2(i1,i2)/TWOPI/fp_w)
            IF(frr .lt. 1.0_rp)THEN
                s_ww = 9.77_rp * frr ** beta_min
            ELSE
                s_ww = 9.77_rp * frr ** beta_max
            ENDIF

            IF(cos(theta_abs(i1,i2)/2.0_rp)**2 .GT. 1.E-20 )THEN
                IF(s_ww*LOG(cos(theta_abs(i1,i2)/2.0_rp)**2) .GT.-170.0_rp) THEN
                    DO jj=1,ithp
                        Norm_DD = Norm_DD + (cos(theta_base(jj)/2.0_rp)**2)**s_ww * dth
                    ENDDO
                    DD_WW(i1,i2) = 1.0_rp/Norm_DD * (cos(theta_abs(i1,i2)/2.0_rp)**2)**s_ww
                ENDIF
            ENDIF
        ENDDO
    ENDDO
ELSE
    DO i1=1,n1o2p1
        DO i2=1,n2
            DD_WW(i1,i2) =1.0_rp/beta*(cos(TWOPI*theta_abs(i1,i2)/(4*beta)))**2.0_rp
        ENDDO
    ENDDO
ENDIF
! ----- End dir function

iseed = 4*n1o2p1*n2_all
!
Cj = 3.279_rp*E_cible
E     = E_cible
test  = 0.0_rp
!
DO WHILE(ABS(test-E_cible)/E_cible.GT.0.001)
   !
   Cj = 3.279_rp*E
   write(*,*) 'E = ', E
   call srand(iseed)
   !
   rnd = RAND(0)
   angle1 = rnd*TWOPI
   !
   rnd = RAND(0)
   angle2 =  rnd*TWOPI
   a_eta(1,1)  = 0.0_rp
   a_phis(1,1) = 0.0_rp
   !
   DO i1 = 2, n1o2p1
      rnd = RAND(0)
	  angle1 = rnd*TWOPI
      !
      rnd = RAND(0)
	  angle2 =  rnd*TWOPI
      angle = 0.0_rp
      !
      if(omega_n2(i1,1).LT.1.0_rp) then !FIXME: it is assumed that omega_p=1 (peak)
         sigma = 0.07_rp
      else
         sigma = 0.09_rp
      endif
      theta = ATAN(ky_n2(1)/kx(i1))
      !
      ! Directional spectrum : phi(omega,theta) = psi(omega) x G(theta)
      ! Here G(theta) = 1/beta x cos(2 pi / (4 beta))
      !
      if(ABS(theta).LE.beta) then
         phi_JONSWAP(i1,1) = Cj*omega_n2(i1,1)**(-5.0_rp)*exp(-5.0_rp/(4.0_rp*(omega_n2(i1,1))**(4.0_rp)))*gamma** &
              (exp(-(omega_n2(i1,1)-1.0_rp)**2/(2.0_rp*sigma**2)))*DD_WW(i1,1)
      ELSE
         phi_JONSWAP(i1,1) = 0.0_rp
      ENDIF
	!%%%%%%%%%%%%%%% Takes 1D / 2D into account
	IF(n2 == 1) THEN
        a_eta(i1,1)  = (2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,1))*phi_JONSWAP(i1,1)*pioxlen*pioylen)**(0.5_rp) &
           *exp(i*(angle1+angle2+angle))
        a_phis(i1,1) = -1.0_rp*i/omega_n2(i1,1)*(2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,1))*phi_JONSWAP(i1,1) &
           *pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
	ELSE
	    a_eta(i1,1)  = (2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,1)**3.0_rp)*phi_JONSWAP(i1,1)*pioxlen*pioylen)**(0.5_rp) &
           *exp(i*(angle1+angle2+angle))
      	a_phis(i1,1) = -1.0_rp*i/omega_n2(i1,1)*(2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,1)**3.0_rp)*phi_JONSWAP(i1,1) &
           *pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
	END IF
   ENDDO
   !
   DO i2 = 2, n2
      rnd = RAND(0)
	  angle1 =  rnd*TWOPI
      !
      rnd = RAND(0)
	  angle2 =  rnd*TWOPI
      angle = 0.0_rp
      !
      if(omega_n2(1,i2).LT.1.0_rp) then !FIXME: it is assumed that omega_p=1 (peak)
         sigma = 0.07_rp
      else
         sigma = 0.09_rp
      endif
      theta = PIO2 !plus ou moins pi/2 <- ATAN(ky_n2(i2)/kx(1)), kx nul
      !
      ! Directional spectrum : phi(omega,theta) = psi(omega) x G(theta)
      ! Here G(theta) = 1/beta x cos(2 pi / (4 beta))
      !
      if(ABS(theta).LE.beta) then
         phi_JONSWAP(1,i2) = Cj*omega_n2(1,i2)**(-5.0_rp)*exp(-5.0_rp/(4.0_rp*(omega_n2(1,i2))**(4.0_rp)))*gamma** &
              (exp(-(omega_n2(1,i2)-1.0_rp)**2/(2.0_rp*sigma**2)))*DD_WW(1,i2)
      else
         phi_JONSWAP(1,i2) = 0.0_rp
      endif
      a_eta(1,i2)  = (2.0_rp*1.0_rp/(2.0_rp*omega_n2(1,i2)**3.0_rp)*phi_JONSWAP(1,i2)*pioxlen*pioylen)**(0.5_rp) &
           *exp(i*(angle1+angle2+angle))
      a_phis(1,i2) = -1.0_rp*i/omega_n2(1,i2)*(2.0_rp*1.0_rp/(2.0_rp*omega_n2(1,i2)**3.0_rp)*phi_JONSWAP(1,i2) &
           *pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
      DO i1 = 2, n1o2p1
         rnd = RAND(0)
		 angle1 =  rnd*TWOPI
         !
         rnd = RAND(0)
		 angle2 =  rnd*TWOPI
         angle = 0.0_rp
         !
         if(omega_n2(i1,i2).LT.1.0_rp) then !FIXME: it is assumed that omega_p=1 (peak)
            sigma = 0.07_rp
         else
            sigma = 0.09_rp
         endif
         theta = ATAN(ky_n2(i2)/kx(i1))
         !
         ! Directional spectrum : phi(omega,theta) = psi(omega) x G(theta)
         ! Here G(theta) = 1/beta x cos(2 pi / (4 beta))
         !
         if(ABS(theta).LE.beta) then
            phi_JONSWAP(i1,i2) = Cj*omega_n2(i1,i2)**(-5.0_rp)*exp(-5.0_rp/(4.0_rp*(omega_n2(i1,i2))**(4.0_rp)))*gamma** &
                 (exp(-(omega_n2(i1,i2)-1.0_rp)**2/(2.0_rp*sigma**2)))*DD_WW(i1,i2)
         else
            phi_JONSWAP(i1,i2) = 0.0_rp
         endif
         a_eta(i1,i2)  = (2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,i2)**3.0_rp)*phi_JONSWAP(i1,i2)*pioxlen*pioylen)**(0.5_rp) &
              *exp(i*(angle1+angle2+angle))
         a_phis(i1,i2) = -1.0_rp*i/omega_n2(i1,i2)*(2.0_rp*1.0_rp/(2.0_rp*omega_n2(i1,i2)**3.0_rp)*phi_JONSWAP(i1,i2) &
              *pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
      ENDDO
   ENDDO
   ! n1o2p1 is special for n1 even (only a real part)
   IF (iseven(n1)) THEN
    a_eta(n1o2p1,:)  = 0.0_rp !i*ABS(a_eta(n1o2p1,:)) !0.0_rp !REAL(a_eta(n1o2p1,:))
    a_phis(n1o2p1,:) = 0.0_rp !i*ABS(a_phis(n1o2p1,:)) !0.0_rp !REAL(a_phis(n1o2p1,:))
   ENDIF
   !
   ! Convergence on energy...
   !
	a_eta_temp = a_eta
	a_phi_temp = a_phis
    ! FIXME : put the choice in input file?
	i_initiate_NL = 0
    IF(i_initiate_NL == 1)THEN
	    CALL initiate_NL
	    WRITE(*,*)'***************'
	    WRITE(*,*)'	test =',test
	    CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phi_temp, a_eta_temp, da_eta)
	    energy = calc_energy(a_eta, a_phis, da_eta)
	    WRITE(*,*)'***************'
	    test = energy(4)
     ELSE IF(i_initiate_NL == 2)THEN
        ! FIXME : depend on wind input/dissip ?
        !CALL initiate_NL_o2
        print*, 'check this test case'
        stop
        WRITE(*,*)'***************'
        WRITE(*,*)'      Comparative computation of Modal Pressure - Done'
        WRITE(*,*)'***************'
        !stop
      ELSE
	    CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phi_temp, a_eta_temp, da_eta)
	    energy = calc_energy(a_eta, a_phis, da_eta)
	    test = energy(4)
      ENDIF
      write(*,*) 'E_current=', test, 'E cible =', E_cible
      E = E * E_cible /test
ENDDO
E_tot=E_cible
!
END SUBROUTINE initiate_irreg
!
!
!
SUBROUTINE initiate_irreg_f
!
!-------------------------------------------
!
! Initialisation of irregular multi-directional sea state (linear)
! 	- From spectrum file from WAVEWATCH III
!	- Bilinear or bicubic interpolation
!
!-------------------------------------------
!
IMPLICIT NONE

REAL(RP) :: theta, E, pioxlen, pioylen,d1,d2,d3,d4
REAL(RP) :: rnd, angle, angle1, angle2
REAL(RP), DIMENSION(m1o2p1,m2) :: phi_E
REAL(RP), DIMENSION(m1o2p1,m2,4) :: coef_ww2cart
INTEGER, DIMENSION(m1o2p1,m2) :: ind_o_ww,ind_t_ww
INTEGER :: iseed,i1,i2,i_int,I_cpt,ndp
REAL(RP), DIMENSION(ifreq) :: ones_f
REAL(RP), DIMENSION(ithet) :: ones_t
INTEGER, DIMENSION(31*ifreq*ithet+1) :: IWK
REAL(RP), DIMENSION(8*ifreq*ithet) :: WK
REAL(RP), DIMENSION(ifreq*ithet) :: xd,yd,zd
REAL(RP), DIMENSION(1) :: xi,yi,zi
!% I_int = 1 bilinear interpol between WW3 intput and HOS grid
!% I_int = 2 bivar interpol
I_int=2

pioxlen = TWOPI/xlen_star
pioylen = TWOPI/ylen_star


 a_eta(1,1)  = 0.0_rp
 a_phis(1,1) = 0.0_rp
 phi_E = 0.0_rp
 ones_t(:)=1.0
 ones_f(:)=1.0

IF(i_int == 1) THEN
write(*,*) '----------Bilinear Interpolation from WW3 Input-----------------'
DO i1 = 1, n1o2p1
    DO i2 = 1, n2
    IF(omega_n2(i1,i2) .ge. freq_ww(1) .and. omega_n2(i1,i2) .le. freq_ww(ifreq)) THEN
    ind_o_ww(i1,i2)=MINLOC(abs(freq_ww(:)-omega_n2(i1,i2)*ones_f(:)),1)

        IF(freq_ww(ind_o_ww(i1,i2)) .gt. omega_n2(i1,i2)) ind_o_ww(i1,i2)=ind_o_ww(i1,i2)-1
!write(*,*) 'ind_o_ww(i1,i2)=',ind_o_ww(i1,i2),freq_ww(ind_o_ww(i1,i2)),omega_n2(i1,i2)
!
	IF(i1 == 1) THEN
	theta = PIO2
	ELSE
        theta = ATAN(ky_n2(i2)/kx(i1))
	IF(theta .lt. 0.0_rp) theta = theta + TWOPI
	END IF
	ind_t_ww(i1,i2)=MINLOC(abs(theta*ones_t-thet_ww),1)
	IF(thet_ww(ind_t_ww(i1,i2)) .gt. theta .and. ind_t_ww(i1,i2).ne. 1) ind_t_ww(i1,i2)=ind_t_ww(i1,i2)-1
    !
        IF(ind_t_ww(i1,i2) /= ithet)THEN
            d1 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2))**2*cos(thet_ww(ind_t_ww(i1,i2))))**2+ &
	            (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2))**2*sin(thet_ww(ind_t_ww(i1,i2))))**2)
            d2 = 	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*cos(thet_ww(ind_t_ww(i1,i2))))**2+ &
	            (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*sin(thet_ww(ind_t_ww(i1,i2))))**2)
            d3 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2))**2*cos(thet_ww(ind_t_ww(i1,i2)+1)))**2+ &
	            (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2))**2*sin(thet_ww(ind_t_ww(i1,i2)+1)))**2)
            d4 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*cos(thet_ww(ind_t_ww(i1,i2)+1)))**2+ &
	            (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*sin(thet_ww(ind_t_ww(i1,i2)+1)))**2)
            !
            coef_ww2cart(i1,i2,1)=d2*d3*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
            coef_ww2cart(i1,i2,2)=d1*d3*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
            coef_ww2cart(i1,i2,3)=d1*d2*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
            coef_ww2cart(i1,i2,4)=d1*d2*d3/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
            !
            phi_E(i1,i2)= coef_ww2cart(i1,i2,1) * phi_ww(ind_o_ww(i1,i2),ind_t_ww(i1,i2)) + &
		            coef_ww2cart(i1,i2,2) * phi_ww(ind_o_ww(i1,i2)+1,ind_t_ww(i1,i2)) + &
		            coef_ww2cart(i1,i2,3) * phi_ww(ind_o_ww(i1,i2),ind_t_ww(i1,i2)+1) + &
		            coef_ww2cart(i1,i2,4) * phi_ww(ind_o_ww(i1,i2)+1,ind_t_ww(i1,i2)+1)
            !IF(i2==1) write(*,*) 'phi_E(i1,i2)=',phi_E(i1,i2),phi_ww(ind_o_ww(i1,i2),ind_t_ww(i1,i2)),phi_ww(ind_o_ww(i1,i2),ind_t_ww(i1,i2)+1)
    ELSE
        d1 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2))**2*cos(thet_ww(ind_t_ww(i1,i2))))**2+ &
	        (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2))**2*sin(thet_ww(ind_t_ww(i1,i2))))**2)
        d2 = 	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*cos(thet_ww(ind_t_ww(i1,i2))))**2+ &
	        (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*sin(thet_ww(ind_t_ww(i1,i2))))**2)
        d3 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2))**2*cos(thet_ww(1)))**2+ &
	        (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2))**2*sin(thet_ww(1)))**2)
        d4 =	SQRT((omega_n2(i1,i2)**2*cos(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*cos(thet_ww(1)))**2+ &
	        (omega_n2(i1,i2)**2*sin(theta)-freq_ww(ind_o_ww(i1,i2)+1)**2*sin(thet_ww(1)))**2)
        !
        coef_ww2cart(i1,i2,1)=d2*d3*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
        coef_ww2cart(i1,i2,2)=d1*d3*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
        coef_ww2cart(i1,i2,3)=d1*d2*d4/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
        coef_ww2cart(i1,i2,4)=d1*d2*d3/(d2*d3*d4+d1*d3*d4+d1*d2*d4+d1*d2*d3)
        !
        !
        phi_E(i1,i2)= coef_ww2cart(i1,i2,1) * phi_ww(ind_o_ww(i1,i2),ind_t_ww(i1,i2)) + &
		        coef_ww2cart(i1,i2,2) * phi_ww(ind_o_ww(i1,i2)+1,ind_t_ww(i1,i2)) + &
		        coef_ww2cart(i1,i2,3) * phi_ww(ind_o_ww(i1,i2),1) + &
		        coef_ww2cart(i1,i2,4) * phi_ww(ind_o_ww(i1,i2)+1,1)

	ENDIF
ELSE
    phi_E(i1,i2)=0.0_rp
ENDIF
ENDDO
ENDDO
ELSEIF(i_int == 2) THEN
    write(*,*) '-----------Bicubic Interpolation from WW3 Input---------'
    I_cpt=1
    NDP=ithet*ifreq
    xd=0.0
    yd=0.0
    zd=0.0
	DO i1=1,ifreq
		DO i2=1,ithet
	        xd((i1-1)*ithet + i2) = freq_ww(i1)**2*cos(thet_ww(i2))
	        yd((i1-1)*ithet + i2) = freq_ww(i1)**2*sin(thet_ww(i2))
	        zd((i1-1)*ithet + i2) = phi_ww(i1,i2)
		ENDDO
	ENDDO
	DO i1 = 1, n1o2p1
       	DO i2 = 1, n2
	    IF(omega_n2(i1,i2) .ge. freq_ww(1) .and. omega_n2(i1,i2) .le. freq_ww(ifreq)) THEN
	        xi=k_abs(i1,i2)*cos(theta_abs(i1,i2))
	        yi=k_abs(i1,i2)*sin(theta_abs(i1,i2))
		        CALL IDBVIP(I_cpt,ndp,REAL(xd,RP),REAL(yd,RP),REAL(zd,RP),1,REAL(xi,RP),REAL(yi,RP),REAL(zi,RP),IWK,REAL(WK,RP))
	        phi_E(i1,i2) = zi(1)
	        I_cpt = 2
	    ENDIF
		ENDDO
	ENDDO
    IF(ISEVEN(n2)) phi_E(:,n2o2p1) = 0.0_rp
END IF
!
iseed=2*n1o2p1*n2
call srand(iseed)
!
a_eta=0.0_cp
a_phis=0.0_cp
E=0.0_rp
DO i1 = 1, n1o2p1
    DO i2 = 1, n2
	    IF(ABS(phi_E(i1,i2)).gt.tiny)THEN
	        IF(i1 /= 1 .or. i2 /= 1)THEN
              rnd = RAND(0)
	          angle1 = rnd*TWOPI
              !
              rnd = RAND(0)
	          angle2 =  rnd*TWOPI
              angle = 0.0_rp
              a_eta(i1,i2)  = (phi_E(i1,i2)/abs(phi_E(i1,i2)))* (1.0_rp/omega_n2(i1,i2)**3 &
              	*abs(phi_E(i1,i2))*pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
              a_phis(i1,i2) = (phi_E(i1,i2)/abs(phi_E(i1,i2)))* (-1.0_rp*i/omega_n2(i1,i2)) *(1.0_rp/omega_n2(i1,i2)**3 &
              	*abs(phi_E(i1,i2))*pioxlen*pioylen)**(0.5_rp)*exp(i*(angle1+angle2+angle))
	            E = E + 1.0_rp/(2.0_rp*omega_n2(i1,i2)**3)*phi_E(i1,i2)*pioxlen*pioylen
	        ENDIF
	    ENDIF
    ENDDO
ENDDO
!
E_tot = E
write(*,*) 'E_tot_ini =',E * L_out **2 ,'Hs_ini =',4*SQRT(E) * L_out, 'Tp_ini =',Tp_real
!
END SUBROUTINE initiate_irreg_f
!
!
!
SUBROUTINE read_irreg_f
!
!----------------------------------------------------------------------
!
! Reading of the spectrum given in a WAVEWATCH III simulation
! FIXME : describe the subroutine
! FIXME : do we keep it? Useful? Bring back global variables here?
!
!----------------------------------------------------------------------
!
IMPLICIT NONE
!
REAL(RP) ::dthet,O_p_ww,k_p_ww,E_int_tot
INTEGER :: i1,i2,ii,jj,ind_o_max
REAL(RP), DIMENSION(ifreq) :: phi_ww_int,S_in_int,S_diss_int,E_int,S_nl_int
REAL(RP), DIMENSION(ifreq,ithet) :: S_diss_ww, S_in_ww,S_nl_ww
REAL(RP),DIMENSION(7) :: phi_temp
!
! FIXME : name of file in input file?
OPEN(1002,file='/home/perignon/PHD/WW3/ww3.1m5s_nl.src')
!
READ(1002,*)
DO ii=1,FLOOR(ifreq/8.0)
    READ(1002,*) freq_ww(8*(ii-1)+1:8*ii)
ENDDO
DO ii=1,FLOOR(ithet/7.0)
  READ(1002,*) thet_ww(7*(ii-1)+1:7*ii)
ENDDO
READ(1002,*) thet_ww(7*(ii-1)+1:ithet)
dthet=abs(thet_ww(2)-thet_ww(1))
!
READ(1002,*)
READ(1002,*)
DO i1=1,FLOOR(ithet * ifreq / 7.0)
    READ(1002,*) phi_temp(1:7)
    DO i2 = 1,7
        jj=FLOOR(((i1-1)*7+i2-1)/REAL(ifreq,RP))+1
        ii=(i1-1)*7+i2 - (jj-1)*ifreq
        phi_ww(ii,jj)=phi_temp(i2)
    ENDDO
ENDDO
!
READ(1002,*) phi_ww(ifreq-(ithet *ifreq-7*FLOOR(ithet * ifreq / 7.0))+1:ifreq,ithet)

DO i1=1,FLOOR(ithet * ifreq / 7.0)
    READ(1002,*) phi_temp(1:7)
    DO i2 = 1,7
        jj=FLOOR(((i1-1)*7+i2-1)/REAL(ifreq,RP))+1
        ii=(i1-1)*7+i2 - (jj-1)*ifreq
        S_in_ww(ii,jj)=phi_temp(i2)
    ENDDO
ENDDO
READ(1002,*) S_in_ww(ifreq-(ithet *ifreq-7*FLOOR(ithet * ifreq / 7.0))+1:ifreq,ithet)
!
DO i1=1,FLOOR(ithet * ifreq / 7.0)
    READ(1002,*) phi_temp(1:7)
    DO i2 = 1,7
        jj=FLOOR(((i1-1)*7+i2-1)/REAL(ifreq,RP))+1
        ii=(i1-1)*7+i2 - (jj-1)*ifreq
        S_nl_ww(ii,jj)=phi_temp(i2)
    ENDDO
ENDDO
READ(1002,*) S_nl_ww(ifreq-(ithet *ifreq-7*FLOOR(ithet * ifreq / 7.0))+1:ifreq,ithet)
!
DO i1=1,FLOOR(ithet * ifreq / 7.0)
    READ(1002,*) phi_temp(1:7)
    DO i2 = 1,7
        jj=FLOOR(((i1-1)*7+i2-1)/REAL(ifreq,RP))+1
        ii=(i1-1)*7+i2 - (jj-1)*ifreq
        S_diss_ww(ii,jj)=phi_temp(i2)
    ENDDO
ENDDO
READ(1002,*) S_diss_ww(ifreq-(ithet *ifreq-7*FLOOR(ithet * ifreq / 7.0))+1:ifreq,ithet)
CLOSE(1002)
!
S_in_int = 0.0_rp
S_diss_int = 0.0_rp
E_int = 0.0_rp
S_nl_int = 0.0_rp
E_int_tot = 0.0_rp
DO jj=1,ithet
    DO ii = 1,ifreq
        S_in_int(ii) = S_in_int(ii) + S_in_ww(ii,jj) * dthet
        S_diss_int(ii) = S_diss_int(ii) + S_diss_ww(ii,jj) * dthet
        E_int(ii) = E_int(ii) + phi_ww(ii,jj) * dthet
        S_nl_int(ii) = S_nl_int(ii) + S_nl_ww(ii,jj) * dthet
    ENDDO
ENDDO
DO ii = 1,ifreq-1
    E_int_tot = E_int_tot + 0.5_rp*(E_int(ii)+E_int(ii+1))*(freq_ww(ii+1)-freq_ww(ii))
ENDDO
!
133 FORMAT(4(ES15.8,X),ES15.8)
103 FORMAT(A,F9.2,A,I5,A,I5)
!
! Transform frequency -> omega & theta grid to same ref. as HOS
phi_ww=phi_ww / TWOPI
thet_ww=PIO2-thet_ww
DO i1=1,ithet
    IF(thet_ww(i1) .lt. 0.0_rp) thet_ww(i1) = thet_ww(i1) + TWOPI
ENDDO
freq_ww=freq_ww * TWOPI
!
! ************
! Adim the variables
! ************
DO ii=1,ifreq
    phi_ww_int(ii)=SUM(phi_ww(ii,:)) * dthet
ENDDO
!
ind_o_max=MAXLOC(phi_ww_int,1)
O_p_ww=freq_ww(ind_o_max)
k_p_ww = O_p_ww**2 / grav ! FIXME: check if it is working if grav.ne.9.81
Tp_real = TWOPI / O_p_ww
L_out=1/k_p_ww **2
T_out=1/O_p_ww
Hs_real = 4.0_rp * sqrt(E_int_tot)
!%%%%%%%%%%%%%%%
OPEN(130,file='./S_BAJ_int.dat',status='unknown')
CALL write_input(130)
WRITE(130,'(A)')'TITLE=" 1D WW3 Source Term"'		!
WRITE(130,'(A)') 'VARIABLES="fx","E_int","S_in_ww3_int","S_diss_int","S_nl_int"'

IF (tecplot == 11) THEN
    WRITE(130,103)'ZONE SOLUTIONTIME = ',0.0_rp,', I=',ifreq
ELSE
    WRITE(130,103)'ZONE T = "',0.0_rp,'", I=',ifreq
END IF
!
DO ii = 1, ifreq
    WRITE(130,133) freq_ww(ii)/ TWOPI,E_int(ii),S_in_int(ii),-S_diss_int(ii),S_nl_int(ii)
END DO
CLOSE(130)
!%%%%%%%%%%%%%%%
freq_ww =freq_ww / SQRT(grav*k_p_ww) ! FIXME: check if it is working if grav.ne.9.81
phi_ww = phi_ww * k_p_ww **2 * O_p_ww
!
!
END SUBROUTINE read_irreg_f
!
!
!
SUBROUTINE initiate_irreg_i
!
!-------------------------------------------
!
! Initialisation of irregular multi-directional sea state (linear)
! 	- From HOS-ocean simulation, file '3d_ini.dat'
!
!-------------------------------------------
!
IMPLICIT NONE
!
INTEGER :: i1,i2,ii,unit,m_i,n_i
REAL(RP) :: time_i
CHARACTER :: A1,A2
!
unit = 1000
!
open(1000,FILE='3d_ini.dat')
DO ii= 1,63
 CALL read_blank_line(unit)
ENDDO
READ(unit,1003) A1,time_i,A2,m_i,n_i
IF(m_i == m1 .and. n_i == n2)THEN
    DO i2 =  1, n2
        DO i1 = 1, m1
            READ(unit,1004) eta(i1,i2),phis(i1,i2)
	    ENDDO
    ENDDO
    eta = eta /L_out
    phis = phis / (L_out**2/T_out)
    CALL space_2_fourier(eta, a_eta)
    CALL space_2_fourier(phis, a_phis)
ELSE
    WRITE(*,*) 'INPUT ERROR - WRONG DIMENSIONS'
    STOP
ENDIF
!
1003 FORMAT(A,F9.2,A,I5,A,I5)
1004 FORMAT(2(ES12.5,X),ES12.5)
END SUBROUTINE initiate_irreg_i
!
!
!
SUBROUTINE initiate_NL
!-----------
! Inital condition following an improved second order formulation - Ref Dalzell (3D) + Febo Thesis and Dunkan&Drake (2D)
! -  eta = eta_1 + eta_stokes + eta_pos + eta_neg
! 	(-  k_nbecomes k_n' according to 3rd order interactions ) not yet
!-----------
INTEGER :: i1,j1,i2,j2
REAL(RP), DIMENSION(m1,m2) :: eta_1,eta_s,eta_pos,eta_neg,eta_t,phi_1,phi_pos,phi_neg,phi_d,eta_3
COMPLEX(RP), DIMENSION(m1o2p1,m2) :: a_eta_s, a_eta_pos, a_eta_neg,a_phi_1,a_phi_pos,a_phi_neg,a_phi_d,a_eta_3
REAL(RP) :: E_NL, eta_mult,A_pos,A_neg,B_pos,B_neg
!
E_NL=E_tot
!
eta_1=0.0_rp
eta_s=0.0_rp
eta_3=0.0_rp
eta_pos=0.0_rp
eta_t=0.0_rp
eta_neg=0.0_rp
B_pos=0.0_rp
B_neg=0.0_rp
eta_mult = 1.0_rp
a_eta_pos=0.0_rp
a_eta_neg=0.0_rp
a_eta_s=0.0_rp
a_eta_3=0.0_rp
a_phi_1=0.0_rp
a_phi_pos=0.0_rp
a_phi_neg=0.0_rp
!
! Eta Computing
!
 CALL fourier_2_space(a_eta,eta_1)
!
DO j1=1,n2
    DO i1=1,n1o2p1
    	IF(i1 .ne.1 .or. j1 .ne.1)THEN
            a_phi_1(i1,j1)= - i * a_eta(i1,j1) * g_star / (omega_n2(i1,j1))
        	a_phi_d(i1,j1) = k_abs(i1,j1) * a_phi_1(i1,j1)
        	IF(2*i1-1 .le. n1o2p1)THEN
        	    IF(j1 .le. n2o2p1 .and. 2*j1-1 .le. n2o2p1)THEN
    		        a_eta_s(2*i1-1,2*j1-1)= + 0.5_rp * a_eta(i1,j1) ** 2 * k_abs(i1,j1)
	            ENDIF
	        IF(j1 .gt. n2o2p1 .and. 2*(n2-j1+1) .lt. n2o2p1) a_eta_s(2*i1-1,2*j1-n2-1)= 0.5_rp * a_eta(i1,j1) ** 2 * k_abs(i1,j1)
	    ENDIF
        !
	    DO j2=1,n2o2p1
	        DO i2=1,n1o2p1
		        IF(k_abs(i2,j2).ge.k_abs(i1,j1) .and. (i1 .ne. i2 .or. j1 .ne. j2))THEN
                    !
                    A_pos = - ((omega_n2(i1,j1)*omega_n2(i2,j2)*(omega_n2(i1,j1)+omega_n2(i2,j2))* &
                    	(1.0_rp-cos(theta_abs(i1,j1)-theta_abs(i2,j2)))) &
			            / ((omega_n2(i1,j1)+omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2)) )
		            A_neg = + ((omega_n2(i1,j1)*omega_n2(i2,j2)*(omega_n2(i2,j2)-omega_n2(i1,j1)) &
		            	*(1.0_rp+cos(theta_abs(i1,j1)-theta_abs(i2,j2)))) &
			            / ((omega_n2(i1,j1)-omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2)) )
                    !
		            B_pos = 0.5_rp/ g_star * (omega_n2(i1,j1)**2 + omega_n2(i2,j2)**2) &
		                - 0.5_rp / g_star * omega_n2(i1,j1) * omega_n2(i2,j2) * (1 - cos(theta_abs(i1,j1)-theta_abs(i2,j2))) &
		                * ( ((omega_n2(i1,j1) + omega_n2(i2,j2))**2 + g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2))&
		                /   ((omega_n2(i1,j1) + omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2)))
                    !
		            B_neg = 0.5_rp/ g_star * (omega_n2(i1,j1)**2 + omega_n2(i2,j2)**2) &
		                + 0.5_rp / g_star * omega_n2(i1,j1) * omega_n2(i2,j2) * (1 + cos(theta_abs(i1,j1)-theta_abs(i2,j2))) &
		                * ( ((omega_n2(i1,j1) - omega_n2(i2,j2))**2 + g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2))&
		                /   ((omega_n2(i1,j1) - omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2)))
		            IF(i1+i2 .le. n1o2p1)THEN
		                IF(j1+j2 .le. n2o2p1)THEN
                 	        a_eta_pos(i1+i2-1,j1+j2-1) = a_eta_pos (i1+i2-1,j1+j2-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j1+j2-1) = a_phi_pos (i1+i2-1,j1+j2-1) - i* A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
                        IF(j1.gt.n2o2p1 .and. j1+j2 .gt. n2+1)THEN
                            a_eta_pos(i1+i2-1,j1+j2-n2-1) = a_eta_pos (i1+i2-1,j1+j2-n2-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j1+j2-n2-1) = a_phi_pos (i1+i2-1,j1+j2-n2-1) - i* A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
                        IF(j1.gt.n2o2p1 .and. j1+j2 .le. n2+1)THEN
                            a_eta_pos(i1+i2-1,j1+j2-1) = a_eta_pos (i1+i2-1,j1+j2-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j1+j2-1) = a_phi_pos (i1+i2-1,j1+j2-1) - i* A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
		            ENDIF
		            IF(i2-i1 .ge. 0)THEN
		                IF(j2-j1 .ge. 0)THEN
                            a_eta_neg (i2-i1+1,j2-j1+1 ) = a_eta_neg (i2-i1+1,j2-j1+1 ) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+1 ) = a_phi_neg (i2-i1+1,j2-j1+1 ) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                        IF(j1 .le. n2o2p1 .and. j2-j1 .lt. 0)THEN
                            a_eta_neg (i2-i1+1,j2-j1+n2+1 ) = a_eta_neg (i2-i1+1,j2-j1+n2+1 ) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+n2+1 ) = a_phi_neg (i2-i1+1,j2-j1+n2+1 ) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                        IF(j1 .gt. n2o2p1 .and. j2-j1+n2+1 .le. n2o2p1)THEN
                            a_eta_neg (i2-i1+1,j2-j1+n2+1 ) = a_eta_neg (i2-i1+1,j2-j1+n2+1 ) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+n2+1 ) = a_phi_neg (i2-i1+1,j2-j1+n2+1 ) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                    ENDIF
		        ENDIF
	        ENDDO
	    ENDDO
	    DO j2=n2o2p1+1,n2
	        DO i2=1,n1o2p1
		        IF(k_abs(i2,j2) .ge. k_abs(i1,j1).and. (i1 .ne. i2 .or. j1 .ne. j2))THEN
                    !
                    A_pos = - ((omega_n2(i1,j1)*omega_n2(i2,j2)*(omega_n2(i1,j1)+omega_n2(i2,j2)) &
                    	*(1.0_rp-cos(theta_abs(i1,j1)-theta_abs(i2,j2)))) &
			            / ((omega_n2(i1,j1)+omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2)) )
		            A_neg = + ((omega_n2(i1,j1)*omega_n2(i2,j2)*(omega_n2(i2,j2)-omega_n2(i1,j1)) &
		            	*(1.0_rp+cos(theta_abs(i1,j1)-theta_abs(i2,j2)))) &
			            / ((omega_n2(i1,j1)-omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2)) )
                    !
		            B_pos = 0.5_rp/ g_star * (omega_n2(i1,j1)**2 + omega_n2(i2,j2)**2) &
		                - 0.5_rp / g_star * omega_n2(i1,j1) * omega_n2(i2,j2) * (1 - cos(theta_abs(i1,j1)-theta_abs(i2,j2))) &
		                * ( ((omega_n2(i1,j1) + omega_n2(i2,j2))**2 + g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2))&
		                /   ((omega_n2(i1,j1) + omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)+kx(i2))**2 +(ky_n2(j1)+ky_n2(j2))**2)))
                    !
		            B_neg = 0.5_rp/ g_star * (omega_n2(i1,j1)**2 + omega_n2(i2,j2)**2) &
		                + 0.5_rp / g_star * omega_n2(i1,j1) * omega_n2(i2,j2) * (1 + cos(theta_abs(i1,j1)-theta_abs(i2,j2))) &
		                * ( ((omega_n2(i1,j1) - omega_n2(i2,j2))**2 + g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2))&
		                /   ((omega_n2(i1,j1) - omega_n2(i2,j2))**2 - g_star * SQRT((kx(i1)-kx(i2))**2 +(ky_n2(j1)-ky_n2(j2))**2)))
		            IF(i1+i2 .le. n1o2p1)THEN
		                IF(j1.le.n2o2p1 .and. j1.le.n2-j2+1 )THEN
                            a_eta_pos(i1+i2-1,j2+j1-1) = a_eta_pos (i1+i2-1,j2+j1-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j2+j1-1) = a_phi_pos (i1+i2-1,j2+j1-1) - i*A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
		                IF(j1.le.n2o2p1 .and. j1.gt. n2-j2+1)THEN
                            a_eta_pos(i1+i2-1,j1-n2+j2-1) = a_eta_pos (i1+i2-1,j1-n2+j2-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j1-n2+j2-1) = a_phi_pos (i1+i2-1,j1-n2+j2-1) - i*A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
		                IF(j1.gt.n2o2p1 .and. j1+j2-n2-1.gt.n2o2p1)THEN
                            a_eta_pos(i1+i2-1,j1+j2-n2-1) = a_eta_pos (i1+i2-1,j1+j2-n2-1) + B_pos * a_eta(i1,j1) * a_eta(i2,j2)
                            a_phi_pos(i1+i2-1,j1+j2-n2-1) = a_phi_pos (i1+i2-1,j1+j2-n2-1) - i*A_pos * a_eta(i1,j1) * a_eta(i2,j2)
                        ENDIF
		            ENDIF
		            IF(i2-i1 .ge. 0)THEN
                        IF(j1 .le. n2o2p1 .and. j2-j1 .ge. n2o2p1)THEN
                            a_eta_neg (i2-i1+1,j2-j1+1) = a_eta_neg (i2-i1+1,j2-j1+1) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+1) = a_phi_neg (i2-i1+1,j2-j1+1) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                        IF(j1 .gt. n2o2p1 .and. j2-j1 .ge. 0 )THEN
                            a_eta_neg (i2-i1+1,j2-j1+1) = a_eta_neg (i2-i1+1,j2-j1+1) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+1) = a_phi_neg (i2-i1+1,j2-j1+1) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                        IF(j1 .gt. n2o2p1 .and. j2-j1 .lt. 0 )THEN
                            a_eta_neg (i2-i1+1,j2-j1+n2+1) = a_eta_neg (i2-i1+1,j2-j1+n2+1) &
                            	+ B_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                            a_phi_neg (i2-i1+1,j2-j1+n2+1) = a_phi_neg (i2-i1+1,j2-j1+n2+1) &
                            	- i*A_neg * CONJG(a_eta(i1,j1)) * a_eta(i2,j2)
                        ENDIF
                    ENDIF
	            ENDIF
	        ENDDO
	    ENDDO
    ENDIF
ENDDO
ENDDO
CALL fourier_2_space(a_eta_s,eta_s)
CALL fourier_2_space(a_eta_pos,eta_pos)
CALL fourier_2_space(a_eta_neg,eta_neg)
CALL fourier_2_space(a_phi_1,phi_1)
CALL fourier_2_space(a_phi_d,phi_d)
CALL fourier_2_space(a_phi_pos,phi_pos)
CALL fourier_2_space(a_phi_neg,phi_neg)
!
phi_d =  phi_d * eta_1
!
eta = eta_1 + eta_s + eta_pos + eta_neg
phis= phi_1 + phi_pos + phi_neg + phi_d
!
CALL space_2_fourier(eta,a_eta)
CALL space_2_fourier(phis,a_phis)
!
END SUBROUTINE initiate_NL
!
END MODULE initial_condition
