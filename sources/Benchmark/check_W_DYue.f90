program check_W_DYue
!
! This program computes the Dommermuth & Yue tests for vertical velocity accuracy
! Depending on choices of input parameters (test_x, test_y, test_xy, test_t_x ...)
! Perform the computation along x, y and/or xy (45°) direction and possibly with time integration
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
USE variables_3D
USE RF_solution
USE fourier_r2c
USE resol_HOS
!
USE runge_kutta
USE energy_calc
!
USE velocities
!
IMPLICIT NONE
!
TYPE(RF_data)       :: RF_obj
INTEGER, PARAMETER  :: RF_dealias = 1
INTEGER, PARAMETER  :: calc_case  = 1
! use 1 together with Nma=7, Mm=7 for the Dommermuth and Yue test (m1=256 and p1=15)
!     0 together with Nma=13 and Mm=9 (m1=512 and p1=21)
REAL(RP), DIMENSION(5)           :: steepness
CHARACTER(LEN=100), DIMENSION(5) :: files
INTEGER, PARAMETER               :: Nma = 7, Mm = 7 !Nma = 13, Mm = 9
! Test along x, y or 45° direction
INTEGER, PARAMETER               :: test_x = 1, test_y = 1, test_xy = 1
! Test with time stepping along x, y or 45° direction (1 for propagation, -1 for direct and inverse propagation)
INTEGER, PARAMETER               :: test_t_x = 0, test_t_y = 0, test_t_xy=0
!
INTEGER, PARAMETER               :: test_vel=0
!
INTEGER, DIMENSION(Nma,Mm,5)     :: NM
REAL(RP), DIMENSION(Nma,Mm,5,6)  :: error, log_error
INTEGER                          :: i_steep, N_loop, Npts
!
CHARACTER(LEN=100) :: filename
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: phizRF
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: phizRF_2D
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eta, phis
!
REAL(RP)                   :: pioxlen, k2, spatial_ref, thkh, angle
REAL(RP)                   :: dely, pioylen
INTEGER                    :: Mo2, part
INTEGER                    :: i1, i2, j, ideep, Nd1o2, err_x, err_y, err_xy
!
! Adaptive Time Step Runge Kutta scheme
TYPE(RK_parameters)        :: RK_param
!
REAL(RP)                   :: dt_out, dt_rk4, dt_lin, dt, dt_correc, volume, energy(4)
REAL(RP)                   :: time_cur, time_next, h_rk, h_loc
REAL(RP)                   :: error_old, error2
INTEGER                    :: idx_max(2)
!
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_eta_rk, a_phis_rk, da_eta
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_eta_ref, a_phis_ref
!
! New velocities
REAL(RP) :: error_SL
!
i_case = 0
!
g_star = 1.0_rp
depth_star = 1.0e15_rp

error = -1.0_rp
!
CALL fill_DYue_cases(steepness, files, NM, calc_case)
!
ideep = 0 ! infinite depth, deep water case
!
err_x  = 0
err_y  = 0
err_xy = 0
!
! analysis along x axis
! FIXME : do along y axis...
!
N_der(2)  = 1
n2c       = 1
N_dea(2)  = 1
i_filt(2) = 0
n2c_filt  = 1
order_max(2) = 0
i_dealias(2) = 0
!
n2       = 1
Nd2      = 1
Nd2o2p1  = 1
n2o2p1   = 1
n2m1o2m1 = -1
!
IF (test_x.EQ.1) THEN
DO i_steep = 1,5
   WRITE(*,'(A,F5.2)') 'Steepness = ', steepness(i_steep)
   !
   filename = files(i_steep)
   WRITE(*,'(A,A)') 'filename = ', files(i_steep)
   !
   CALL read_RF_data(filename, RF_obj, .TRUE.)
   !
   DO M = 2,2*Mm,2
      WRITE(*,*) M
      Mo2 = M/2
      DO N_loop = 1,Nma
         IF (NM(N_loop,M/2,i_steep) /= 0.AND.NM(N_loop,M/2,i_steep).LE.m1/2) THEN
            Npts   = 2*NM(N_loop,M/2,i_steep) ! for comparison with DYue
            WRITE(*,*) M, Npts
            n1     = Npts
            n1o2p1 = Npts/2 + 1
            n1c    = n1o2p1
            ! Memory allocation
            ALLOCATE(etapm_ext(md1, Nd2, M+1))
            !ALLOCATE(kth(md1o2p1,Nd2,M), kth_all(md1o2p1,Nd2,M))
            ALLOCATE(kth(md1o2p1,Nd2,M), kth_all(md1o2p1,Nd2o2p1,M))
            ALLOCATE(oneoj(M))
            ! Building reference solution and mesh
            !
            ! Reference solution
            CALL build_RF_reference(RF_obj, Npts)
!             CALL display_RF_data(RF_obj)
            ! Memory allocation
            !ALLOCATE(phis(Npts,1), eta(Npts,1), phizRF(Npts))
            ALLOCATE(phis(m1,m2), eta(m1,m2), phizRF(Npts))
            Nd1   = ((M+1) * n1)/2
            Nd1o2p1 = Nd1/2+1
            !
            kth = 0.0_rp
            kth_all = 0.0_rp
            kx  = 0.0_rp
            !
            N_der(1) = 2*n1o2p1
            N_der(2) = 1
            !
            ! Mesh generation
            xlen = RF_obj%lambda
            ylen = RF_obj%lambda
            dely = ylen / Npts
            x(1:Npts) = RF_obj%x(1:Npts)
            !
            ! Wave numbers
            pioxlen = TWOPI / xlen
            pioylen = TWOPI / ylen
            !
            !	wave numbers
            DO i1 = 1, Nd1o2p1
	            kx(i1)  = REAL(i1 - 1,RP) * pioxlen
            END DO
            DO i2 = 1, Nd2o2p1
	            ky(i2)  = REAL(i2 - 1,RP) * pioylen
            END DO
            !  y-wave numbers (on n2 modes)
            DO i2 = 1, n2o2p1
               ky_n2(i2) = REAL(i2 - 1,RP) * pioylen
            ENDDO
            DO i2 = 2,n2o2p1
               ky_n2(n2-i2+2) = - REAL(i2 - 1,RP) * pioylen
            END DO
            !
            IF (iseven(n2)) ky_n2(n2o2p1) = REAL(n2o2p1 - 1,RP) * pioylen
            !
            ! HOS modal coefficients of the vertical derivatives
            goomega(1,1) = 1.0_rp
            !
            i2=1
            DO i1 = 2, Nd1o2p1
	            k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	            k(i1,i2)     = SQRT(k2)
                thkh         = TANH(k(i1,i2) * depth_star)
                kth_all(i1,i2,1) = k(i1,i2)
                omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                goomega(i1,i2) = g_star / omega(i1,i2)
                IF (abs(k2) > tiny) THEN
                    c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                ELSE
                    c(i1,i2)        = SQRT(g_star * depth_star)
                END IF
                !
	            kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	            IF (M > 1) kth_all(i1,i2,2) = k2
	            DO j = 3,M
		            kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	            END DO
            END DO
            !
            DO i2 = 2, Nd2o2p1
               DO i1 = 1, Nd1o2p1
	               k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	               k(i1,i2)     = SQRT(k2)
                  thkh         = TANH(k(i1,i2) * depth_star)
                  kth_all(i1,i2,1) = k(i1,i2)
                  omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                  goomega(i1,i2) = g_star / omega(i1,i2)
                  IF (abs(k2) > tiny) THEN
                     c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                  ELSE
                     c(i1,i2)        = SQRT(g_star * depth_star)
                  END IF
                  !
	               kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	               IF (M > 1) kth_all(i1,i2,2) = k2
	               DO j = 3,M
		               kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	               END DO
               END DO
            END DO
			!
            goomega_n2(1,1) = 1.0_rp
            !
            i2=1
            DO i1 = 2, Nd1o2p1
                k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
            END DO
            !
            DO i2 = 2, n2
               DO i1 = 1, Nd1o2p1
                  k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                  omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                  goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
               END DO
            END DO
            !
            DO j=1,M
               oneoj(j) = 1.0_rp / REAL(j, RP)
            END DO
            !
            IF (RF_dealias == 1) THEN
               eta(1:Npts,1)    = RF_obj%eta_dealiased
               phis(1:Npts,1)   = RF_obj%phis_dealiased
               phizRF = RF_obj%W_dealiased  ! OK pour les erreurs en Fourier ; pas bon en spatial
            ELSE IF (RF_dealias == 2) THEN
               eta(1:Npts,1)    = RF_obj%eta_dealiased
               phis(1:Npts,1)   = RF_obj%phis_dealiased
               phizRF = RF_obj%W ! OK pour les erreurs en spatial ; pas bon en Fourier
            ELSE IF (RF_dealias == 0) THEN
               eta(1:Npts,1)    = RF_obj%eta
               phis(1:Npts,1)   = RF_obj%phis
               phizRF = RF_obj%W
                ! FIXME: test 09/2014
      			IF (iseven(n1)) THEN ! last mode of eta will be a cosine without corresponding sine for phis and W
      				CALL fourier_ini(3)
      				!
      				CALL space_2_fourier(eta,a_eta)
         			a_eta(n1o2p1,:) = 0.0_cp
         			CALL fourier_2_space(a_eta,eta)
         			!
         			CALL fourier_end(3)
     			END IF
            END IF
            !
            DEALLOCATE(RF_obj%x,RF_obj%eta,RF_obj%phis,RF_obj%W,RF_obj%eta_dealiased,RF_obj%phis_dealiased,RF_obj%W_dealiased)
            spatial_ref = 1.0_rp
            !
            ! Partial dealiasing
            part=2
            Nd1 = ((part+1) * n1) / 2
            Nd1o2p1 = Nd1/2+1
            N_dea(1)     = n1o2p1
            i_dealias(1) = 1
            order_max(1) = 2-1 ! to ensure a dealiasing at each simple product
!             N_der(1)     = ((part+1) * N) / 4
            !N_der(1)     = n1o2p1
            i_filt(1)    = 0
            n1c_filt     = n1c
            !
            CALL build_derivatives()
            CALL fourier_ini(3)
            !
            CALL space_2_fourier(eta,a_eta)
            CALL space_2_fourier(phis,a_phis)
            !
            ! calculation of the horizontal derivatives of phis and eta
            call phisxy_etaxy(a_phis, a_eta)
            ! resolution of the different order of the FS potential boundary value problem to get modes
            call HOSphis_modes_fully_dealiased(a_phis)
            !
   		    error(N_loop,M/2,i_steep,1) = MAXVAL(ABS(phizRF(1:n1) - phiz(1:n1,1)))
            !
            CALL fourier_end(3)
            !
            ! Check present results w.r.t. previous HOS results
            ! FIXME : fill in the whole matrix
            ! FIXME : create a specific calc_case for this test?
            !
            IF (calc_case == 1) THEN
				IF((error(N_loop,M/2,i_steep,1)-error(N_loop,M/2,i_steep,5))-0.05*error(N_loop,M/2,i_steep,5).GT.tiny) THEN
					print*, 'Worsening of previous HOS computation: partial dealiasing'
					print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
					err_x = 1
				ELSEIF((-error(N_loop,M/2,i_steep,1)+error(N_loop,M/2,i_steep,5))-0.05*error(N_loop,M/2,i_steep,5).GT.tiny)THEN
					IF (abs(error(N_loop,M/2,i_steep,5)-1.0_rp).GT.tiny) THEN
					print*, 'Enhancement of previous HOS computation: partial dealiasing'
					print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
					ENDIF
				ENDIF
            ENDIF
            !
            ! Total dealiasing
            part   = M
            Nd1 = ((part+1) * n1) / 2
            Nd1o2p1 = Nd1/2+1
            order_max(1) = MAX(part-1,1)
            N_dea(1)     = n1o2p1
            i_dealias(1) = 1
            !N_der(1)     = 2*n1o2p1 ! FLOOR(1.7 * n1o2p1) ! ((part+1)/2) * n1o2p1
!             N_der(1)     = ((part+1)/2) * n1o2p1
            CALL build_derivatives()
            CALL fourier_ini(3)
            CALL space_2_fourier(eta,a_eta)
            CALL space_2_fourier(phis,a_phis)
            !
            ! calculation of the horizontal derivatives of phis and eta
            call phisxy_etaxy(a_phis, a_eta)
            ! resolution of the different order of the FS potential boundary value problem to get modes
            call HOSphis_modes_fully_dealiased(a_phis)
            error(N_loop,M/2,i_steep,2) = MAXVAL(ABS(phizRF(1:n1)-phiz(1:n1,1)))
            !
            IF ((test_vel.EQ.1).AND.(M.EQ.6).AND.(Npts.EQ.64).AND.(i_steep.EQ.3)) THEN !Test the velocity calculations
                ! Initiate da_eta for energy calculation
                CALL fill_butcher_array(RK_param)
                CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
                CALL HOSvel2(25,M,a_eta,a_phis,0.0_rp) !HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current)
                CALL reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect,a_eta,da_eta,error_SL) !FIXME: last term should da_eta
                ! FIXME : compare to RF ?
                ! stop
            ENDIF
            !
			IF ((M.EQ.8).AND.(Npts.EQ.32).AND.(i_steep.EQ.3)) THEN !test the time integration...
				IF ((ABS(test_t_x).EQ.1)) THEN !test the time integration...
					!
					CALL fill_butcher_array(RK_param)
					!
					time_cur    = 0.0_rp
					T_stop_star = 1000.0_rp * RF_obj%T
					dt_out      = RF_obj%T
					!
					IF (n2 == 1) THEN
					   dt_rk4 = 2.0_rp / SQRT(PI) * SQRT(xlen / n1o2p1)
					ELSE
					   dt_rk4 = 2.0_rp / SQRT(PI) * MIN(SQRT(xlen / n1o2p1), SQRT(ylen / n2o2p1))
					END IF
					dt_lin = 10.0_rp * dt_rk4
					dt     = 0.5_rp * dt_rk4
					!
					toler = 1.e-12 !1e-7
					!
					a_eta_ref(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
					a_phis_ref(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
					!
					idx_max(1:2) = MAXLOC(ABS(a_eta_ref))
					! Initiate da_eta for energy calculation
					CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
					!
					a_eta(1:n1o2p1,1:n2)  = a_eta_ref(1:n1o2p1,1:n2)
					a_phis(1:n1o2p1,1:n2) = a_phis_ref(1:n1o2p1,1:n2)
					!
					DO WHILE (time_cur-T_stop_star <= tiny)
					   !
					   !
					   ! Output of volume and energy
					   volume = REAL(a_eta(1,1),RP)
					   energy = calc_energy(a_eta, a_phis, da_eta)
					   !
					   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
					   !
					  CALL fourier_2_space(a_eta, eta)
					  CALL fourier_2_space(a_phis,phis)
					  !
					  IF (ABS(time_cur).LT.tiny) THEN
						OPEN(66,FILE='3d_x.dat')
						WRITE(66,'(A)') 'VARIABLES="x", "y", "eta", "phis"'
						!
						OPEN(77,FILE='a_3d_x.dat')
						WRITE(77,'(A)') 'VARIABLES="kx", "ky", "a_eta", "a_phis"'
						!
						OPEN(88,FILE='error_time_x.dat')
						WRITE(88,'(A)') 'VARIABLES="time", "err amp", "err phase"'
					  ENDIF
					  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
					  DO i2 = 1, n2
						 DO i1 = 1, n1
							WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
								eta(i1,i2), phis(i1,i2)
						 END DO
					  END DO
					  !
					  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
					  DO i2 = n2o2p1+1, n2
						 DO i1 = 1, n1o2p1
							WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
								ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
						 END DO
					  END DO
					  DO i2 = 1, n2o2p1
						 DO i1 = 1, n1o2p1
							 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
								ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
						 END DO
					  END DO
					  !
					  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
						ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
						ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
						-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
					   !
					   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
					   !
					   ! Going to next time step
					   time_next = time_cur + dt_out
					   IF (time_next > T_stop_star) time_next = T_stop_star ! if last output
					   h_rk    = dt ! starting time step
					   DO WHILE(time_cur < time_next)
						  h_loc = h_rk ! local time step
						  IF (time_next < time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
						  !
						  ! Analytical integration of the linear part
						  ! starting at current time
						  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
						  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
						  !
						  ! Going to time_cur + h_loc
						  ! adaptive time step
						  !
						  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
							a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
						  error_old = error2
						  IF (error2 > toler) THEN
							dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
						  ELSE
							a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
							a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
							time_cur   = time_cur + h_loc
							dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
						  END IF
						  ! Step size ratio bounds (Mathematica)
						  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
						  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
						  !
						  ! Step size security factor
						  dt_correc = 0.98_rp * dt_correc
						  ! New step size
						  h_rk      = h_rk * dt_correc
						  ! Stability checkings
						  h_rk      = MIN(h_rk, dt_out, dt_lin)
					   END DO
					   dt = h_rk ! saving the step size for next start
					END DO
					IF (test_t_x.EQ.1) THEN
						CLOSE(66)
						CLOSE(77)
						CLOSE(88)
					ENDIF
				ENDIF
				!
				IF ((test_t_x.EQ.-1)) THEN !test the back-time integration...
					!
					CALL fill_butcher_array(RK_param)
					!
					! time_cur defined as current time
					T_stop_star = 0.0_rp
					dt_out      = -RF_obj%T
					!
					dt     = -0.5_rp * dt_rk4 !starting point identical to front propagation
					!
					! Use same time tolerance
					DO WHILE (time_cur > T_stop_star)
					   !
					   ! Output of volume and energy
					   volume = REAL(a_eta(1,1),RP)
					   energy = calc_energy(a_eta, a_phis, da_eta)
					   !
					   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
					   !
					  CALL fourier_2_space(a_eta, eta)
					  CALL fourier_2_space(a_phis,phis)
					  !
					  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
					  DO i2 = 1, n2
						 DO i1 = 1, n1
							WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
								eta(i1,i2), phis(i1,i2)
						 END DO
					  END DO
					  !
					  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
					  DO i2 = n2o2p1+1, n2
						 DO i1 = 1, n1o2p1
							WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
								ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
						 END DO
					  END DO
					  DO i2 = 1, n2o2p1
						 DO i1 = 1, n1o2p1
							 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
								ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
						 END DO
					  END DO
					  !
					  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
						ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
						ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
						-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
					   !
					   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
					   !
					   ! Going to next time step
					   time_next = time_cur + dt_out
					   IF (time_next < T_stop_star) time_next = T_stop_star ! if last output
					   h_rk    = dt ! starting time step
					   DO WHILE(time_cur > time_next)
						  h_loc = h_rk ! local time step
						  IF (time_next > time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
						  !
						  ! Analytical integration of the linear part
						  ! starting at current time
						  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
						  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
						  !
						  ! Going to time_cur + h_loc
						  ! adaptive time step
						  !
						  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
							a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
						  error_old = error2
						  IF (error2 > toler) THEN
							dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
						  ELSE
							a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
							a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
							time_cur   = time_cur + h_loc
							dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
						  END IF
						  ! Step size ratio bounds (Mathematica)
						  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
						  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
						  !
						  ! Step size security factor
						  dt_correc = 0.98_rp * dt_correc
						  ! New step size
						  h_rk      = ABS(h_rk) * dt_correc
						  ! Stability checkings
						  h_rk      = -MIN(h_rk, -dt_out, dt_lin)
					   END DO
					   dt = h_rk ! saving the step size for next start
					END DO
					!
					CLOSE(66)
					CLOSE(77)
					CLOSE(88)
				ENDIF
			ENDIF
            !
            CALL fourier_end(3)
            !
            ! Check present results w.r.t. previous HOS results
            ! FIXME : fill in the whole matrix
            ! FIXME : create a specific calc_case for this test?
            !
            IF (calc_case == 1) THEN
				IF((error(N_loop,M/2,i_steep,2)-error(N_loop,M/2,i_steep,6))-0.05*error(N_loop,M/2,i_steep,6).GT.tiny) THEN
					print*, 'Worsening of previous HOS computation: total dealiasing'
					print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
					!stop
					err_x = 1
				ELSEIF ((-error(N_loop,M/2,i_steep,2)+error(N_loop,M/2,i_steep,6)) &
					-0.05*error(N_loop,M/2,i_steep,6).GT.tiny) THEN
					IF (ABS(error(N_loop,M/2,i_steep,6)-1.0_rp).GT.tiny) THEN
						print*, 'Enhancement of previous HOS computation: total dealiasing'
						print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
					ENDIF
				ENDIF
            ENDIF
            !
            part = 4
            IF (part >= (4*Npts)/Npts-1) THEN
               ! Best dealiasing
               Nd1     = ((part+1) * n1) / 2
               Nd1o2p1 = Nd1/2+1
               order_max(1) = part-1
               i_dealias(1) = 1
               N_dea(1)     = n1o2p1
               !N_der(1)     = 2*n1o2p1
!                N_der(1)     = ((part+1)/2) * n1o2p1
               CALL build_derivatives()
               CALL fourier_ini(3)
               CALL space_2_fourier(eta,a_eta)
               CALL space_2_fourier(phis,a_phis)
               ! calculation of the horizontal derivatives of phis and eta
               call phisxy_etaxy(a_phis, a_eta)
               ! resolution of the different order of the FS potential boundary value problem to get modes
               call HOSphis_modes_fully_dealiased(a_phis)
               error(N_loop,M/2,i_steep,3) = MAXVAL(ABS(phizRF(1:n1)-phiz(1:n1,1)))
               !
               CALL fourier_end(3)
               !
            END IF

            DEALLOCATE(phis, eta, phizRF)
            DEALLOCATE(oneoj,etapm_ext)
            DEALLOCATE(kth,kth_all)
         END IF
      END DO
   END DO
END DO
!
!
WRITE(*,'(A)') 'Writing results file'
!
WRITE(*,'(A)') 'result_DYue.dat'
!
OPEN(1,FILE='result_DYue.dat')
CALL write_result(1, error(:,:,:,1))
!
WRITE(1,*)
CALL write_result(1, error(:,:,:,2))
!
WRITE(1,*)
CALL write_result(1, error(:,:,:,3))
!
CLOSE(1)
!
WRITE(*,'(A)') 'comp_DYue.dat'
!
OPEN(2,FILE='comp_DYue.dat')
!
IF (calc_case == 1) THEN
   CALL write_comp(2, error(:,:,:,1), error(:,:,:,4))
ELSE
   CALL write_comp(2, error(:,:,:,2), error(:,:,:,1))
END IF
!
DO N_loop = 1,120
   WRITE(2,'(A)', ADVANCE='NO') '_'
END DO
WRITE(2,'(A)') ''
WRITE(2,*)
!
IF (calc_case == 1) THEN
   CALL write_comp(2, error(:,:,:,2), error(:,:,:,4))
ELSE
   CALL write_comp(2, error(:,:,:,3), error(:,:,:,1))
END IF
!
IF (calc_case == 1) THEN
   DO N_loop = 1,120
      WRITE(2,'(A)', ADVANCE='NO') '_'
   END DO
   WRITE(2,'(A)') ''
   WRITE(2,*)
   !
   CALL write_comp(2, error(:,:,:,3), error(:,:,:,4))
END IF
!
CLOSE(2)
!
IF (calc_case == 1) THEN
   WRITE(*,'(A)') 'result_DYue.tex'
   !
   OPEN(1,FILE='result_DYue.tex')
   !
   WRITE(1,'(A)') '\begin{table}[!htbp]'
   WRITE(1,'(A)') '\begin{center}'
   WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
   WRITE(1,'(A)') '\hline'
   DO i_steep = 1,5
      DO N_loop = 1,Nma
         IF (NM(N_loop,1,i_steep) /= 0) THEN
            IF (N_loop == 1) THEN
               WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
            ELSE
               WRITE(1,'(X)', ADVANCE='NO')
            END IF
            WRITE(1,'(A)', ADVANCE='NO') '&'
            Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
            WRITE(1,'(I3)', ADVANCE='NO') Npts
            DO M = 2,14,2
               IF (error(N_loop,M/2,i_steep,4) > 0.0_rp) THEN
                  WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,4)
               ELSE
                  WRITE(1,'(A)', ADVANCE='NO') '&'
               END IF
            END DO
            WRITE(1,'(A)', ADVANCE='NO') '\\'
            WRITE(1,*)
         END IF
      END DO
      WRITE(1,'(A)') '\hline'
   END DO
   WRITE(1,'(A)') '\end{tabular}'
   WRITE(1,'(A)') '\end{center}'
   WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, Dommermuth et Yue}'
   WRITE(1,'(A)') '\label{tab:error-W-DYue}'
   WRITE(1,'(A)') '\end{table}'
   !
   WRITE(1,'(A)') '\begin{table}[!htbp]'
   WRITE(1,'(A)') '\begin{center}'
   WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
   WRITE(1,'(A)') '\hline'
   DO i_steep = 1,5
      DO N_loop = 1,Nma
         IF (NM(N_loop,1,i_steep) /= 0) THEN
            IF (N_loop == 1) THEN
               WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
            ELSE
               WRITE(1,'(X)', ADVANCE='NO')
            END IF
            WRITE(1,'(A)', ADVANCE='NO') '&'
            Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
            WRITE(1,'(I3)', ADVANCE='NO') Npts
            DO M = 2,14,2
               IF (error(N_loop,M/2,i_steep,1) > 0.0_rp) THEN
                  IF (error(N_loop,M/2,i_steep,4) > 1.05_rp*error(N_loop,M/2,i_steep,1)) THEN
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,1)
                  ELSE
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,1)
                  END IF
               ELSE
                  WRITE(1,'(A)', ADVANCE='NO') '&'
               END IF
            END DO
            WRITE(1,'(A)', ADVANCE='NO') '\\'
            WRITE(1,*)
         END IF
      END DO
      WRITE(1,'(A)') '\hline'
   END DO
   WRITE(1,'(A)') '\end{tabular}'
   WRITE(1,'(A)') '\end{center}'
   WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement partiel}'
   WRITE(1,'(A)') '\label{tab:error-W-simple}'
   WRITE(1,'(A)') '\end{table}'
   !
   WRITE(1,'(A)') '\begin{table}[!htbp]'
   WRITE(1,'(A)') '\begin{center}'
   WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
   WRITE(1,'(A)') '\hline'
   DO i_steep = 1,5
      DO N_loop = 1,Nma
         IF (NM(N_loop,1,i_steep) /= 0) THEN
            IF (N_loop == 1) THEN
               WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
            ELSE
               WRITE(1,'(X)', ADVANCE='NO')
            END IF
            WRITE(1,'(A)', ADVANCE='NO') '&'
            Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
            WRITE(1,'(I3)', ADVANCE='NO') Npts
            DO M = 2,14,2
               IF (error(N_loop,M/2,i_steep,2) > 0.0_rp) THEN
                  IF (error(N_loop,M/2,i_steep,4) > 1.1_rp*error(N_loop,M/2,i_steep,2)) THEN
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,2)
                  ELSE
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,2)
                  END IF
               ELSE
                  WRITE(1,'(A)', ADVANCE='NO') '&'
               END IF
            END DO
            WRITE(1,'(A)', ADVANCE='NO') '\\'
            WRITE(1,*)
         END IF
      END DO
      WRITE(1,'(A)') '\hline'
   END DO
   WRITE(1,'(A)') '\end{tabular}'
   WRITE(1,'(A)') '\end{center}'
   WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement complet}'
   WRITE(1,'(A)') '\label{tab:error-W-complet}'
   WRITE(1,'(A)') '\end{table}'
   !
   WRITE(1,'(A)') '\begin{table}[!htbp]'
   WRITE(1,'(A)') '\begin{center}'
   WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
   WRITE(1,'(A)') '\hline'
   WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
   WRITE(1,'(A)') '\hline'
   DO i_steep = 1,5
      DO N_loop = 1,Nma
         IF (NM(N_loop,1,i_steep) /= 0) THEN
            IF (N_loop == 1) THEN
               WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
            ELSE
               WRITE(1,'(X)', ADVANCE='NO')
            END IF
            WRITE(1,'(A)', ADVANCE='NO') '&'
            Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
            WRITE(1,'(I3)', ADVANCE='NO') Npts
            DO M = 2,14,2
               IF (error(N_loop,M/2,i_steep,3) > 0.0_rp) THEN
                  IF (error(N_loop,M/2,i_steep,4) > 1.05_rp*error(N_loop,M/2,i_steep,3)) THEN
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,3)
                  ELSE
                     WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,3)
                  END IF
               ELSE
                  WRITE(1,'(A)', ADVANCE='NO') '&'
               END IF
            END DO
            WRITE(1,'(A)', ADVANCE='NO') '\\'
            WRITE(1,*)
         END IF
      END DO
      WRITE(1,'(A)') '\hline'
   END DO
   WRITE(1,'(A)') '\end{tabular}'
   WRITE(1,'(A)') '\end{center}'
   WRITE(1,'(A,I2,A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement interm\''ediaire  p=',part,'}'
   WRITE(1,'(A)') '\label{tab:error-W-partiel}'
   WRITE(1,'(A)') '\end{table}'
   !
   CLOSE(1)
END IF
!
WRITE(*,'(A)') 'result_DYue.tec'
!
WHERE (error < tiny)
   log_error = 0.0_rp
ELSEWHERE
   log_error = LOG(error) / LOG(10.0_rp)
END WHERE
OPEN(3,FILE='result_DYue.tec')
WRITE(3,'(A)') 'VARIABLES="M", "N", "error", "log_1_0(error)"'
IF (calc_case == 1) THEN
   DO i_steep = 1,5
      WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="DYue',steepness(i_steep),'", I=',Nma,', J=',Mm
         DO M = 2,2*Mm,2
      DO N_loop = 1,Nma
         Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
            WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,4),log_error(N_loop,M/2,i_steep,4)
         END DO
      END DO
   END DO
END IF
DO i_steep = 1,5
   WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
      DO M = 2,2*Mm,2
   DO N_loop = 1,Nma
      Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
         WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,1),log_error(N_loop,M/2,i_steep,1)
      END DO
   END DO
END DO
DO i_steep = 1,5
   WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="comp',steepness(i_steep),'", I=',Nma,', J=',Mm
      DO M = 2,2*Mm,2
   DO N_loop = 1,Nma
      Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
         WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,2),log_error(N_loop,M/2,i_steep,2)
      END DO
   END DO
END DO
DO i_steep = 1,5
   WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
      DO M = 2,2*Mm,2
   DO N_loop = 1,Nma
      Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
         WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,3),log_error(N_loop,M/2,i_steep,3)
      END DO
   END DO
END DO
CLOSE(3)
ENDIF
!
IF (test_y.EQ.1) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Make the analysis along y direction now
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    error(:,:,:,:) = -1.0_rp
    !
    CALL fill_DYue_cases(steepness, files, NM, calc_case)
    !
    N_der(1)  = 1
    n1c       = 1
    N_dea(1)  = 1
    i_filt(1) = 0
    n1c_filt  = 1
    order_max(1) = 0
    i_dealias(1) = 0
    !
    n1       = 1
    Nd1      = 1
    Nd1o2p1  = 1
    n1o2p1   = 1
    Nd1o2    = 1
    !
    DO i_steep = 1,5
       WRITE(*,'(A,F5.2)') 'Steepness = ', steepness(i_steep)
       !
       filename = files(i_steep)
       WRITE(*,'(A,A)') 'filename = ', files(i_steep)
       !
       CALL read_RF_data(filename, RF_obj, .TRUE.)
       !
       DO M = 2,2*Mm,2
          WRITE(*,*) M
          Mo2 = M/2
          DO N_loop = 1,Nma
             IF (NM(N_loop,M/2,i_steep) /= 0.AND.NM(N_loop,M/2,i_steep).LE.m2/2) THEN
                Npts   = 2*NM(N_loop,M/2,i_steep) ! for comparison with DYue
                WRITE(*,*) M, Npts
                !
                n2     = Npts
                n2o2p1 = Npts/2 + 1
                n2c    = n2o2p1
                n2m1o2m1 = (n2-1)/2-1

                ! Memory allocation
                ALLOCATE(etapm_ext(md1, md2, M+1))
                !ALLOCATE(kth(Nd1o2p1,md2,M), kth_all(Nd1o2p1,md2,M))
                ALLOCATE(kth(Nd1o2p1,md2,M), kth_all(Nd1o2p1,md2o2p1,M))
                ALLOCATE(oneoj(M))
                ! Building reference solution and mesh
                !
                ! Reference solution
                CALL build_RF_reference(RF_obj, Npts)
    !             CALL display_RF_data(RF_obj)
                ! Memory allocation
                !ALLOCATE(phis(n1,Npts), eta(n1,Npts), phizRF(Npts))
                ALLOCATE(phis(m1,m2), eta(m1,m2), phizRF(Npts))
                Nd2   = ((M+1) * n2)/2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                !
                N_der(2) = 2*n2o2p1
                N_der(1) = 2 ! this is necessary so that it does not affect to zero in build_kth
                !
                kth = 0.0_rp
                kth_all = 0.0_rp
                kx  = 0.0_rp
                ky  = 0.0_rp
                !
                ! Mesh generation
                xlen = RF_obj%lambda
                ylen = RF_obj%lambda
                y(1:Npts) = RF_obj%x(1:Npts)
                x(1:n1) = 0.0_rp
                !
                ! Wave numbers
                pioxlen = TWOPI / xlen
                pioylen = TWOPI / ylen
                !
                !	wave numbers
                DO i1 = 1, Nd1o2p1
	                kx(i1)  = REAL(i1 - 1,RP) * pioxlen
                END DO
                DO i2 = 1, Nd2o2p1
	                ky(i2)  = REAL(i2 - 1,RP) * pioylen
                END DO
                !  y-wave numbers (on n2 modes)
                DO i2 = 1, n2o2p1
                   ky_n2(i2) = REAL(i2 - 1,RP) * pioylen
                ENDDO
                DO i2 = 2,n2o2p1
                   ky_n2(n2-i2+2) = - REAL(i2 - 1,RP) * pioylen
                END DO
                !
                IF (iseven(n2)) ky_n2(n2o2p1) = REAL(n2o2p1 - 1,RP) * pioylen
                !
                ! HOS modal coefficients of the vertical derivatives
                goomega(1,1) = 1.0_rp
                i2=1
                DO i1 = 2, Nd1o2p1
	                k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	                k(i1,i2)     = SQRT(k2)
                    thkh         = TANH(k(i1,i2) * depth_star)
                    kth_all(i1,i2,1) = k(i1,i2)
                    omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                    goomega(i1,i2) = g_star / omega(i1,i2)
                    IF (abs(k2) .GT. tiny) THEN
                       c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                    ELSE
                       c(i1,i2)        = SQRT(g_star * depth_star)
                    END IF
                    !
	                 kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	                 IF (M > 1) kth_all(i1,i2,2) = k2
	                 DO j = 3,M
		                 kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	                END DO
	            ENDDO
                DO i2 = 2, Nd2o2p1
                   DO i1 = 1, Nd1o2p1
	                   k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	                   k(i1,i2)     = SQRT(k2)
                      thkh         = TANH(k(i1,i2) * depth_star)
                      kth_all(i1,i2,1) = k(i1,i2)
                      omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                      goomega(i1,i2) = g_star / omega(i1,i2)
                      IF (abs(k2) .GT. tiny) THEN
                         c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                      ELSE
                         c(i1,i2)        = SQRT(g_star * depth_star)
                      END IF
                      !
	                   kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	                   IF (M > 1) kth_all(i1,i2,2) = k2
	                   DO j = 3,M
		                   kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	                   END DO
                   END DO
                END DO
                !
                goomega_n2(1,1) = 1.0_rp
                i2=1
                DO i1 = 2, Nd1o2p1
                    k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                    omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                    goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
                END DO
                !
                DO i2 = 2, n2
                   DO i1 = 1, Nd1o2p1
                      k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                      omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                      goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
                   END DO
                END DO
                !
                DO j=1,M
                   oneoj(j) = 1.0_rp / REAL(j, RP)
                END DO
                !
                IF (RF_dealias == 1) THEN
                   eta(1,1:Npts)    = RF_obj%eta_dealiased
                   phis(1,1:Npts)   = RF_obj%phis_dealiased
                   phizRF = RF_obj%W_dealiased  ! OK pour les erreurs en Fourier ; pas bon en spatial
                ELSE IF (RF_dealias == 2) THEN
                   eta(1,1:Npts)    = RF_obj%eta_dealiased
                   phis(1,1:Npts)   = RF_obj%phis_dealiased
                   phizRF = RF_obj%W ! OK pour les erreurs en spatial ; pas bon en Fourier
                ELSE IF (RF_dealias == 0) THEN
                   eta(1,1:Npts)    = RF_obj%eta
                   phis(1,1:Npts)   = RF_obj%phis
                   phizRF = RF_obj%W
                	! FIXME: test 09/2014
      				IF (iseven(n2)) THEN ! last mode of eta will be a cosine without corresponding sine for phis and W
      					DO i1=2,n1
                    		eta(i1,:) = eta(1,:)
                		ENDDO
      					CALL fourier_ini(3)
      					!
      					CALL space_2_fourier(eta,a_eta)
         				a_eta(:,n2o2p1) = 0.0_cp
         				CALL fourier_2_space(a_eta,eta)
         				!
                		CALL fourier_end(3)
     				END IF
                END IF
                DO i1=2,n1
                    eta(i1,:) = eta(1,:)
                    phis(i1,:) = phis(1,:)
                ENDDO
                !
                DEALLOCATE(RF_obj%x,RF_obj%eta,RF_obj%phis,RF_obj%W,RF_obj%eta_dealiased,RF_obj%phis_dealiased,RF_obj%W_dealiased)
                !
                spatial_ref = 1.0_rp
                !
                ! Partial dealiasing
                part=2
                Nd2 = ((part+1) * n2) / 2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                N_dea(2)     = n2o2p1
                i_dealias(2) = 1
                order_max(2) = 2-1 ! to ensure a dealising at each simple product
    !             N_der(1)     = ((part+1) * N) / 4
                !N_der(2)     = n2o2p1
                i_filt(2)    = 0
                n2c_filt     = n2c
                !
                CALL fourier_ini(3)
                !
                CALL build_derivatives()
    !           !
                CALL space_2_fourier(eta,a_eta)
                CALL space_2_fourier(phis,a_phis)
                !
                ! calculation of the horizontal derivatives of phis and eta
                call phisxy_etaxy(a_phis, a_eta)
                !
                ! resolution of the different order of the FS potential boundary value problem to get modes
                call HOSphis_modes_fully_dealiased(a_phis)
                !
       		    error(N_loop,M/2,i_steep,1) = MAXVAL(ABS(phizRF(1:n2) - phiz(1,1:n2)))
                !
                CALL fourier_end(3)
                !
                ! Check present results w.r.t. previous HOS results
                ! FIXME : create a specific calc_case for this test?
                !
                IF (calc_case == 1) THEN
					IF((error(N_loop,M/2,i_steep,1)-error(N_loop,M/2,i_steep,5)).GT.0.05*error(N_loop,M/2,i_steep,5)) THEN
						print*, 'Worsening of previous HOS computation: partial dealiasing'
						print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
						err_y = 1
					ELSEIF((-error(N_loop,M/2,i_steep,1)+error(N_loop,M/2,i_steep,5)) &
						-0.05*error(N_loop,M/2,i_steep,5).GT.tiny)THEN
						IF (abs(error(N_loop,M/2,i_steep,5)-1.0_rp).GT.tiny) THEN
							print*, 'Enhancement of previous HOS computation: partial dealiasing'
							print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
						ENDIF
					ENDIF
                ENDIF
                !
                ! Total dealiasing
                part   = M
                Nd2 = ((part+1) * n2) / 2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                order_max(2) = MAX(part-1,1)
                N_dea(2)     = n2o2p1
                i_dealias(2) = 1 !FIXME : changed
                !N_der(2)     = 2*n2o2p1 ! FLOOR(1.7 * n1o2p1) ! ((part+1)/2) * n1o2p1
    !             N_der(1)     = ((part+1)/2) * n1o2p1
                CALL build_derivatives()
                CALL fourier_ini(3)
                CALL space_2_fourier(eta,a_eta)
                CALL space_2_fourier(phis,a_phis)
                !
                ! calculation of the horizontal derivatives of phis and eta
                call phisxy_etaxy(a_phis, a_eta)
                ! resolution of the different order of the FS potential boundary value problem to get modes
                call HOSphis_modes_fully_dealiased(a_phis)
                error(N_loop,M/2,i_steep,2) = MAXVAL(ABS(phizRF(1:n2)-phiz(1,1:n2)))
                !
                IF ((test_vel.EQ.1).AND.(M.EQ.6).AND.(Npts.EQ.64).AND.(i_steep.EQ.3)) THEN !Test the velocity calculations
                    ! Initiate da_eta for energy calculation
                    CALL fill_butcher_array(RK_param)
                	CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
                    CALL HOSvel2(25,M,a_eta,a_phis,0.0_rp) !HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current) !CALL HOSvel2_direct(a_eta,a_phis,0.0_rp,0)
                    CALL reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect, &
                    	a_eta,da_eta,error_SL) !FIXME: last term should da_eta
                    ! FIXME : compare to RF ?
                    !stop
                ENDIF
                !
                IF ((M.EQ.8).AND.(Npts.EQ.32).AND.(i_steep.EQ.3)) THEN
					IF ((ABS(test_t_y).EQ.1)) THEN !test the time integration...
						!
						CALL fill_butcher_array(RK_param)
						!
						time_cur    = 0.0_rp
						T_stop_star = 1000.0_rp * RF_obj%T
						dt_out      = RF_obj%T
						!
						IF (n2 == 1) THEN
						   dt_rk4 = 2.0_rp / SQRT(PI) * SQRT(xlen / n1o2p1)
						ELSE
						   dt_rk4 = 2.0_rp / SQRT(PI) * MIN(SQRT(xlen / n1o2p1), SQRT(ylen / n2o2p1))
						END IF
						dt_lin = 10.0_rp * dt_rk4
						dt     = 0.5_rp * dt_rk4
						!
						toler = 1.e-12 !1e-7
						!
						a_eta_ref(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
						a_phis_ref(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
						!
						idx_max(1:2) = MAXLOC(ABS(a_eta_ref))
						! Initiate da_eta for energy calculation
						CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
						!
						a_eta(1:n1o2p1,1:n2)  = a_eta_ref(1:n1o2p1,1:n2)
						a_phis(1:n1o2p1,1:n2) = a_phis_ref(1:n1o2p1,1:n2)
						!
						DO WHILE (time_cur-T_stop_star <= tiny)
						   !
						   !
						   ! Output of volume and energy
						   volume = REAL(a_eta(1,1),RP)
						   energy = calc_energy(a_eta, a_phis, da_eta)
						   !
						   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
						   !
						  CALL fourier_2_space(a_eta, eta)
						  CALL fourier_2_space(a_phis,phis)
						  !
						  IF (ABS(time_cur).LT.tiny) THEN
							OPEN(66,FILE='3d_y.dat')
							WRITE(66,'(A)') 'VARIABLES="x", "y", "eta", "phis"'
							!
							OPEN(77,FILE='a_3d_y.dat')
							WRITE(77,'(A)') 'VARIABLES="kx", "ky", "a_eta", "a_phis"'
							!
							OPEN(88,FILE='error_time_y.dat')
							WRITE(88,'(A)') 'VARIABLES="time", "err amp", "err phase"'
						  ENDIF
						  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
						  DO i2 = 1, n2
							 DO i1 = 1, n1
								WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
									eta(i1,i2), phis(i1,i2)
							 END DO
						  END DO
						  !
						  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
						  DO i2 = n2o2p1+1, n2
							 DO i1 = 1, n1o2p1
								WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  DO i2 = 1, n2o2p1
							 DO i1 = 1, n1o2p1
								 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  !
						  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
							ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
							ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
							-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
						   !
						   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
						   !
						   ! Going to next time step
						   time_next = time_cur + dt_out
						   IF (time_next > T_stop_star) time_next = T_stop_star ! if last output
						   h_rk    = dt ! starting time step
						   DO WHILE(time_cur < time_next)
							  h_loc = h_rk ! local time step
							  IF (time_next < time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
							  !
							  ! Analytical integration of the linear part
							  ! starting at current time
							  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
							  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
							  !
							  ! Going to time_cur + h_loc
							  ! adaptive time step
							  !
							  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
								a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
							  error_old = error2
							  IF (error2 > toler) THEN
								dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
							  ELSE
								a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
								a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
								time_cur   = time_cur + h_loc
								dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
							  END IF
							  ! Step size ratio bounds (Mathematica)
							  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
							  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
							  !
							  ! Step size security factor
							  dt_correc = 0.98_rp * dt_correc
							  ! New step size
							  h_rk      = h_rk * dt_correc
							  ! Stability checkings
							  h_rk      = MIN(h_rk, dt_out, dt_lin)
						   END DO
						   dt = h_rk ! saving the step size for next start
						END DO
						IF (test_t_y.EQ.1) THEN
							CLOSE(66)
							CLOSE(77)
							CLOSE(88)
						ENDIF
					ENDIF
					!
					IF ((test_t_y.EQ.-1)) THEN !test the back-time integration...
						!
						CALL fill_butcher_array(RK_param)
						!
						! time_cur defined as current time
						T_stop_star = 0.0_rp
						dt_out      = -RF_obj%T
						!
						dt     = -0.5_rp * dt_rk4 !starting point identical to front propagation
						!
						DO WHILE (time_cur > T_stop_star)
						   !
						   ! Output of volume and energy
						   volume = REAL(a_eta(1,1),RP)
						   energy = calc_energy(a_eta, a_phis, da_eta)
						   !
						   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
						   !
						  CALL fourier_2_space(a_eta, eta)
						  CALL fourier_2_space(a_phis,phis)
						  !
						  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
						  DO i2 = 1, n2
							 DO i1 = 1, n1
								WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
									eta(i1,i2), phis(i1,i2)
							 END DO
						  END DO
						  !
						  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
						  DO i2 = n2o2p1+1, n2
							 DO i1 = 1, n1o2p1
								WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  DO i2 = 1, n2o2p1
							 DO i1 = 1, n1o2p1
								 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  !
						  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
							ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
							ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
							-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
						   !
						   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
						   !
						   ! Going to next time step
						   time_next = time_cur + dt_out
						   IF (time_next < T_stop_star) time_next = T_stop_star ! if last output
						   h_rk    = dt ! starting time step
						   DO WHILE(time_cur > time_next)
							  h_loc = h_rk ! local time step
							  IF (time_next > time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
							  !
							  ! Analytical integration of the linear part
							  ! starting at current time
							  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
							  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
							  !
							  ! Going to time_cur + h_loc
							  ! adaptive time step
							  !
							  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
								a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
							  error_old = error2
							  IF (error2 > toler) THEN
								dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
							  ELSE
								a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
								a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
								time_cur   = time_cur + h_loc
								dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
							  END IF
							  ! Step size ratio bounds (Mathematica)
							  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
							  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
							  !
							  ! Step size security factor
							  dt_correc = 0.98_rp * dt_correc
							  ! New step size
							  h_rk      = ABS(h_rk) * dt_correc
							  ! Stability checkings
							  h_rk      = -MIN(h_rk, -dt_out, dt_lin)
						   END DO
						   dt = h_rk ! saving the step size for next start
						END DO
						!
						CLOSE(66)
						CLOSE(77)
						CLOSE(88)
					ENDIF
				ENDIF
                !
                CALL fourier_end(3)
                !
                ! Check present results w.r.t. previous HOS results
                ! FIXME : fill in the whole matrix
                ! FIXME : create a specific calc_case for this test?
                !
                IF (calc_case == 1) THEN
                    !IF (i_steep == 5) THEN
                        IF((error(N_loop,M/2,i_steep,2)-error(N_loop,M/2,i_steep,6)).GT.0.05*error(N_loop,M/2,i_steep,6)) THEN
                            print*, 'Worsening previous HOS computation: total dealiasing'
                            print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
                            err_y = 1
                        ELSEIF ((-error(N_loop,M/2,i_steep,2)+error(N_loop,M/2,i_steep,6)) &
                    		-0.05*error(N_loop,M/2,i_steep,6).GT.tiny) THEN
                    		IF (ABS(error(N_loop,M/2,i_steep,6)-1.0_rp).GT.tiny) THEN
                        		print*, 'Enhancement of previous HOS computation: total dealiasing'
                        		print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
                        	ENDIF
                        ENDIF
                    !ENDIF
                ENDIF
                !
                part = 4
                IF (part >= (4*N)/N-1) THEN
                   ! Best dealiasing
                   Nd2     = ((part+1) * n2) / 2
                   Nd2o2p1 = Nd2/2+1
                   Nd2p1o2 = (Nd2+1)/2
                   order_max(2) = part-1
                   i_dealias(2) = 1
                   N_dea(2)     = n2o2p1
                   !N_der(2)     = 2*n2o2p1
    !                N_der(1)     = ((part+1)/2) * n1o2p1
                   CALL build_derivatives()
                   CALL fourier_ini(3)
                   CALL space_2_fourier(eta,a_eta)
                   CALL space_2_fourier(phis,a_phis)
                   ! calculation of the horizontal derivatives of phis and eta
                   call phisxy_etaxy(a_phis, a_eta)
                   ! resolution of the different order of the FS potential boundary value problem to get modes
                   call HOSphis_modes_fully_dealiased(a_phis)
                   error(N_loop,M/2,i_steep,3) = MAXVAL(ABS(phizRF(1:n2)-phiz(1,1:n2)))
                   !
                   CALL fourier_end(3)
                   !
                END IF

                DEALLOCATE(phis, eta, phizRF)
                DEALLOCATE(oneoj,etapm_ext)
                DEALLOCATE(kth,kth_all)
             END IF
          END DO
       END DO
    END DO
    !
    !
    WRITE(*,'(A)') 'Writing results file'
    !
    WRITE(*,'(A)') 'result_DYue_y.dat'
    !
    OPEN(1,FILE='result_DYue_y.dat')
    CALL write_result(1, error(:,:,:,1))
    !
    WRITE(1,*)
    CALL write_result(1, error(:,:,:,2))
    !
    WRITE(1,*)
    CALL write_result(1, error(:,:,:,3))
    !
    CLOSE(1)
    !
    WRITE(*,'(A)') 'comp_DYue_y.dat'
    !
    OPEN(2,FILE='comp_DYue_y.dat')
    !
    IF (calc_case == 1) THEN
       CALL write_comp(2, error(:,:,:,1), error(:,:,:,4))
    ELSE
       CALL write_comp(2, error(:,:,:,2), error(:,:,:,1))
    END IF
    !
    DO N_loop = 1,120
       WRITE(2,'(A)', ADVANCE='NO') '_'
    END DO
    WRITE(2,'(A)') ''
    WRITE(2,*)
    !
    IF (calc_case == 1) THEN
       CALL write_comp(2, error(:,:,:,2), error(:,:,:,4))
    ELSE
       CALL write_comp(2, error(:,:,:,3), error(:,:,:,1))
    END IF
    !
    IF (calc_case == 1) THEN
       DO N_loop = 1,120
          WRITE(2,'(A)', ADVANCE='NO') '_'
       END DO
       WRITE(2,'(A)') ''
       WRITE(2,*)
       !
       CALL write_comp(2, error(:,:,:,3), error(:,:,:,4))
    END IF
    !
    CLOSE(2)
    !
    IF (calc_case == 1) THEN
       WRITE(*,'(A)') 'result_DYue_y.tex'
       !
       OPEN(1,FILE='result_DYue_y.tex')
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,4) > -tiny) THEN
                      WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,4)
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, Dommermuth et Yue}'
       WRITE(1,'(A)') '\label{tab:error-W-DYue}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,1) > - tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.05_rp*error(N_loop,M/2,i_steep,1) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,1)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,1)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement partiel}'
       WRITE(1,'(A)') '\label{tab:error-W-simple}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,2) > -tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.1_rp*error(N_loop,M/2,i_steep,2) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,2)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,2)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement complet}'
       WRITE(1,'(A)') '\label{tab:error-W-complet}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,3) > -tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.05_rp*error(N_loop,M/2,i_steep,3) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,3)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,3)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A,I2,A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement interm\''ediaire  p=',part,'}'
       WRITE(1,'(A)') '\label{tab:error-W-partiel}'
       WRITE(1,'(A)') '\end{table}'
       !
       CLOSE(1)
    END IF
    !
    WRITE(*,'(A)') 'result_DYue_y.tec'
    !
    WHERE (error < tiny)
       log_error = 0.0_rp
    ELSEWHERE
       log_error = LOG(error) / LOG(10.0_rp)
    END WHERE
    OPEN(3,FILE='result_DYue_y.tec')
    WRITE(3,'(A)') 'VARIABLES="M", "N", "error", "log_1_0(error)"'
    IF (calc_case == 1) THEN
       DO i_steep = 1,5
          WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="DYue',steepness(i_steep),'", I=',Nma,', J=',Mm
             DO M = 2,2*Mm,2
          DO N_loop = 1,Nma
             Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,4),log_error(N_loop,M/2,i_steep,4)
             END DO
          END DO
       END DO
    END IF
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,1),log_error(N_loop,M/2,i_steep,1)
          END DO
       END DO
    END DO
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="comp',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,2),log_error(N_loop,M/2,i_steep,2)
          END DO
       END DO
    END DO
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,3),log_error(N_loop,M/2,i_steep,3)
          END DO
       END DO
    END DO
    CLOSE(3)
ENDIF
!
!
IF (test_xy.EQ.1) THEN
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! Make the analysis along y direction now
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    error(:,:,:,:) = -1.0_rp
    !
    CALL fill_DYue_cases(steepness, files, NM, calc_case)
    !
    DO i_steep = 1,5 !3,3 !1,5
       WRITE(*,'(A,F5.2)') 'Steepness = ', steepness(i_steep)
       !
       filename = files(i_steep)
       WRITE(*,'(A,A)') 'filename = ', files(i_steep)
       !
       CALL read_RF_data(filename, RF_obj, .TRUE.)
       !
       DO M = 2,2*Mm,2
          WRITE(*,*) M
          Mo2 = M/2
          DO N_loop = 1,Nma
             IF (NM(N_loop,M/2,i_steep) /= 0.AND.NM(N_loop,M/2,i_steep).LE.MIN(m1,m2)/2) THEN
                Npts   = 2*NM(N_loop,M/2,i_steep) ! for comparison with DYue
                WRITE(*,*) M, Npts
                !
                n1     = Npts
                n1o2p1 = Npts/2 + 1
                n1c    = n1o2p1
                !
                n2     = Npts
                n2o2p1 = Npts/2 + 1
                n2c    = n2o2p1
                n2m1o2m1 = (n2-1)/2-1

                ! Memory allocation
                ALLOCATE(etapm_ext(md1, md2, M+1))
                !ALLOCATE(kth(md1o2p1,md2,M), kth_all(md1o2p1,md2,M)) ! FIXME: trop grand?
                ALLOCATE(kth(md1o2p1,md2,M), kth_all(md1o2p1,md2o2p1,M))
                ALLOCATE(oneoj(M))
                ! Building reference solution and mesh
                !
                ! Reference solution
                CALL build_RF_reference(RF_obj, Npts)
    !             CALL display_RF_data(RF_obj)
                ! Memory allocation
                !ALLOCATE(phis(n1,Npts), eta(n1,Npts), phizRF(Npts))
                ALLOCATE(phis(m1,m2), eta(m1,m2), phizRF_2D(m1,m2))
                Nd1   = ((M+1) * n1)/2
                Nd1o2p1 = Nd1/2+1
                !
                N_der(1) = 2*n1o2p1
                !
                Nd2   = ((M+1) * n2)/2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                !
                N_der(2) = 2*n2o2p1
                !
                kth = 0.0_rp
                kth_all = 0.0_rp
                kx  = 0.0_rp
                ky  = 0.0_rp
                !
                ! Mesh generation
                angle = pi/4.0_rp !pi/2.0_rp !0.0_rp
                !
                ! ABS value necessary to have positive length and always x and y increasing
                IF((abs(angle).GT.tiny).AND.(abs(angle-pi/2.0_rp).GT.tiny)) THEN
                    xlen = ABS(RF_obj%lambda/cos(angle))
                    ylen = ABS(RF_obj%lambda/sin(angle))
                    !
                    x(1:Npts) = ABS(RF_obj%x(1:Npts)/cos(angle))
                    y(1:Npts) = ABS(RF_obj%x(1:Npts)/sin(angle))
                ELSE
                    xlen = RF_obj%lambda
                    ylen = RF_obj%lambda
                    !
                    x(1:Npts) = RF_obj%x(1:Npts)
                    y(1:Npts) = RF_obj%x(1:Npts)
                ENDIF
                !
                ! Wave numbers
                pioxlen = TWOPI / xlen
                pioylen = TWOPI / ylen
                !
                !	wave numbers
                DO i1 = 1, Nd1o2p1
	                kx(i1)  = REAL(i1 - 1,RP) * pioxlen
                END DO
                DO i2 = 1, Nd2o2p1
	                ky(i2)  = REAL(i2 - 1,RP) * pioylen
                END DO
                !  y-wave numbers (on n2 modes)
                DO i2 = 1, n2o2p1
                   ky_n2(i2) = REAL(i2 - 1,RP) * pioylen
                ENDDO
                DO i2 = 2,n2o2p1
                   ky_n2(n2-i2+2) = - REAL(i2 - 1,RP) * pioylen
                END DO
                !
                IF (iseven(n2)) ky_n2(n2o2p1) = REAL(n2o2p1 - 1,RP) * pioylen
                !
                ! HOS modal coefficients of the vertical derivatives
                goomega(1,1) = 1.0_rp
                !
                i1=1
                DO i2 = 2, Nd2o2p1
	                   k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	                   k(i1,i2)     = SQRT(k2)
                      thkh         = TANH(k(i1,i2) * depth_star)
                      kth_all(i1,i2,1) = k(i1,i2)
                      omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                      goomega(i1,i2) = g_star / omega(i1,i2)
                      IF (ABS(k2) .GT. tiny) THEN
                         c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                      ELSE
                         c(i1,i2)        = SQRT(g_star * depth_star)
                      END IF
                      !
	                   kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	                   IF (M > 1) kth_all(i1,i2,2) = k2
	                   DO j = 3,M
		                   kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	                   END DO
                END DO
                !
                DO i2 = 1, Nd2o2p1
                   DO i1 = 2, Nd1o2p1
	                   k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
	                   k(i1,i2)     = SQRT(k2)
                      thkh         = TANH(k(i1,i2) * depth_star)
                      kth_all(i1,i2,1) = k(i1,i2)
                      omega(i1,i2)       = SQRT(g_star * k(i1,i2) * thkh)
                      goomega(i1,i2) = g_star / omega(i1,i2)
                      IF (ABS(k2) .GT. tiny) THEN
                         c(i1,i2)        = omega(i1,i2) / k(i1,i2)
                      ELSE
                         c(i1,i2)        = SQRT(g_star * depth_star)
                      END IF
                      !
	                   kth_all(i1,i2,1) = kth_all(i1,i2,1) * thkh
	                   IF (M > 1) kth_all(i1,i2,2) = k2
	                   DO j = 3,M
		                   kth_all(i1,i2,j) = k2 * kth_all(i1,i2,j-2)
	                   END DO
                   END DO
                END DO
                !
                goomega_n2(1,1) = 1.0_rp
                !
                i1=1
                DO i2 = 2, n2
                      k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                      omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                      goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
                END DO
                !
                DO i2 = 1, n2
                   DO i1 = 2, Nd1o2p1
                      k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
                      omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
                      goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
                   END DO
                END DO
                !
                DO j=1,M
                   oneoj(j) = 1.0_rp / REAL(j, RP)
                END DO
                !
                !
                IF (RF_dealias == 1) THEN
					DO i2=1,n2
						eta(1:Npts,i2)       = RF_obj%eta_dealiased
						phis(1:Npts,i2)      = RF_obj%phis_dealiased
						phizRF_2D(1:Npts,i2) = RF_obj%W_dealiased ! OK pour les erreurs en Fourier ; pas bon en spatial
					ENDDO
                ELSE IF (RF_dealias == 2) THEN
					DO i2=1,n2
						eta(1:Npts,i2)       = RF_obj%eta_dealiased
						phis(1:Npts,i2)      = RF_obj%phis_dealiased
						phizRF_2D(1:Npts,i2) = RF_obj%W ! OK pour les erreurs en Fourier ; pas bon en spatial
					ENDDO
                ELSE IF (RF_dealias == 0) THEN
                   !
                	! Direct reconstruction
                	!
                	DO i1=1,n1
                    	DO i2=1,n2
                        	CALL build_RF_local(RF_obj, x(i1)*cos(angle)+y(i2)*sin(angle), eta(i1,i2), phis(i1,i2), phizRF_2D(i1,i2))
                    	ENDDO
                	ENDDO
                END IF
                !
                IF (RF_dealias == 1 .OR. RF_dealias == 2) THEN
					CALL fourier_ini(3)
					!
					CALL space_2_fourier(eta,a_eta)
					CALL space_2_fourier(phis,a_phis)
					CALL space_2_fourier(phizRF_2D,a_phiz_RF)
					!
					IF (ABS(angle-PI/4.0_rp).LT.tiny) THEN
						DO i1=1,n1o2p1
							a_eta(i1,i1)     = a_eta(i1,1)
							a_phis(i1,i1)    = a_phis(i1,1)
							a_phiz_RF(i1,i1) = a_phiz_RF(i1,1)
							DO i2=1,n2
								IF(i2.NE.i1) THEN
									a_eta(i1,i2)     = 0.0_rp
									a_phis(i1,i2)    = 0.0_rp
									a_phiz_RF(i1,i2) = 0.0_rp
								ENDIF
							ENDDO
						ENDDO
					ELSEIF (ABS(angle+PI/4.0_rp).LT.tiny) THEN
						DO i1=1,n1o2p1
							a_eta(i1,n2-i1+2)     = a_eta(i1,1)
							a_phis(i1,n2-i1+2)    = a_phis(i1,1)
							a_phiz_RF(i1,n2-i1+2) = a_phiz_RF(i1,1)
							DO i2=1,n2
								IF(i2.NE.n2-i1+2) THEN
									a_eta(i1,i2)     = 0.0_rp
									a_phis(i1,i2)    = 0.0_rp
									a_phiz_RF(i1,i2) = 0.0_rp
								ENDIF
							ENDDO
						ENDDO
					ELSEIF (ABS(angle).LT.tiny) THEN
						! Nothing to do
					ELSEIF (ABS(angle-pi/2.0_rp).LT.tiny) THEN
						DO i2=2,n2o2p1-1
							a_eta(1,i2)     = a_eta(i2,1)*0.5_rp
							a_phis(1,i2)    = a_phis(i2,1)*0.5_rp
							a_phiz_RF(1,i2) = a_phiz_RF(i2,1)*0.5_rp
							!
							a_eta(1,n2-i2+2)     = CONJG(a_eta(i2,1)*0.5_rp)
							a_phis(1,n2-i2+2)    = CONJG(a_phis(i2,1)*0.5_rp)
							a_phiz_RF(1,n2-i2+2) = CONJG(a_phiz_RF(i2,1)*0.5_rp)
						ENDDO
						i2=n2o2p1
						a_eta(1,i2)     = a_eta(i2,1)*0.5_rp
						a_phis(1,i2)    = a_phis(i2,1)*0.5_rp
						a_phiz_RF(1,i2) = a_phiz_RF(i2,1)*0.5_rp
						!
						a_eta(2:n1o2p1,:)     = 0.0_rp
						a_phis(2:n1o2p1,:)    = 0.0_rp
						a_phiz_RF(2:n1o2p1,:) = 0.0_rp
					ELSE
						print*, 'This part is not working for angles different from: pi/4, -pi/4, 0'
						STOP
					ENDIF
					!
					CALL fourier_2_space(a_eta,eta)
					CALL fourier_2_space(a_phis,phis)
					CALL fourier_2_space(a_phiz_RF,phizRF_2D)
                	!
                	CALL fourier_end(3)
                ENDIF
                !
                CALL fourier_ini(3)
                !
                CALL space_2_fourier(eta,a_eta)
                !
                ! FIXME: test 09/2014
      			IF (iseven(n1)) THEN ! last mode of eta will be a cosine without corresponding sine for phis
         			!!a_eta(n1o2p1,n2o2p1) = 0.0_cp
         			!a_eta(n1o2p1,:) = 0.0_cp
         			!a_eta(:,n2o2p1) = 0.0_cp
         			!
!         			CALL space_2_fourier(phizRF_2D,a_phiz_RF)
!         			a_phiz_RF(n1o2p1,n2o2p1) = 0.0_cp
!         			CALL fourier_2_space(a_phiz_RF,phizRF_2D)
     			END IF
                !
                CALL fourier_2_space(a_eta,eta)
                !
                CALL fourier_end(3)
				!
				DEALLOCATE(RF_obj%x,RF_obj%eta,RF_obj%phis,RF_obj%W,RF_obj%eta_dealiased,RF_obj%phis_dealiased,RF_obj%W_dealiased)
				!
                spatial_ref = 1.0_rp
                !
                ! Partial dealiasing
                part=2
                Nd1 = ((part+1) * n1) / 2
                Nd1o2p1 = Nd1/2+1
                N_dea(1)     = n1o2p1
                i_dealias(1) = 1
                order_max(1) = 2-1 ! to ensure a dealising at each simple product
                i_filt(1)    = 0
                n1c_filt     = n1c
                !
                Nd2 = ((part+1) * n2) / 2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                N_dea(2)     = n2o2p1
                i_dealias(2) = 1
                order_max(2) = 2-1 ! to ensure a dealising at each simple product
                i_filt(2)    = 0
                n2c_filt     = n2c
                !
                CALL fourier_ini(3)
                !
                CALL build_derivatives()
                !
                CALL space_2_fourier(eta,a_eta)
                CALL space_2_fourier(phis,a_phis)
                !
                ! calculation of the horizontal derivatives of phis and eta
                call phisxy_etaxy(a_phis, a_eta)
                !
                ! resolution of the different order of the FS potential boundary value problem to get modes
                call HOSphis_modes_fully_dealiased(a_phis)
                !
       		    error(N_loop,M/2,i_steep,1) = MAXVAL(ABS(phizRF_2D(1:n1,1:n2) - phiz(1:n1,1:n2)))
                !
                CALL fourier_end(3)
                !
                ! Check present results w.r.t. previous HOS results
                ! FIXME : fill in the whole matrix
                ! FIXME : create a specific calc_case for this test?
                !
                IF (calc_case == 1) THEN
					IF((error(N_loop,M/2,i_steep,1)-error(N_loop,M/2,i_steep,5)).GT.0.05*error(N_loop,M/2,i_steep,5)) THEN
						print*, 'Worsening of previous HOS computation: partial dealiasing'
						print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
						print*, 'location error', MAXLOC(ABS(phizRF_2D(1:n1,1:n2)-phiz(1:n1,1:n2)))
						err_xy = 1
					ELSEIF((-error(N_loop,M/2,i_steep,1)+error(N_loop,M/2,i_steep,5)) &
						-0.05*error(N_loop,M/2,i_steep,5).GT.tiny)THEN
						IF (abs(error(N_loop,M/2,i_steep,5)-1.0_rp).GT.tiny) THEN
							print*, 'Enhancement of previous HOS computation: partial dealiasing'
							print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,1),error(N_loop,M/2,i_steep,5)
						ENDIF
					ENDIF
                ENDIF
                !
                ! Total dealiasing
                part   = M
                Nd1 = ((part+1) * n1) / 2
                Nd1o2p1 = Nd1/2+1
                order_max(1) = MAX(part-1,1)
                N_dea(1)     = n1o2p1
                i_dealias(1) = 1
                !
                Nd2 = ((part+1) * n2) / 2
                Nd2o2p1 = Nd2/2+1
                Nd2p1o2 = (Nd2+1)/2
                order_max(2) = MAX(part-1,1)
                N_dea(2)     = n2o2p1
                i_dealias(2) = 1

                CALL build_derivatives()
                CALL fourier_ini(3)
                CALL space_2_fourier(eta,a_eta)
                CALL space_2_fourier(phis,a_phis)
                !
                ! calculation of the horizontal derivatives of phis and eta
                call phisxy_etaxy(a_phis, a_eta)
                ! resolution of the different order of the FS potential boundary value problem to get modes
                call HOSphis_modes_fully_dealiased(a_phis)
                error(N_loop,M/2,i_steep,2) = MAXVAL(ABS(phizRF_2D(1:n1,1:n2) - phiz(1:n1,1:n2)))
                !
                CALL space_2_fourier(phiz,a_phiz)
                CALL space_2_fourier(phizRF_2D,a_phiz_RF)
                !
                ! Check present results w.r.t. previous HOS results
                ! FIXME : fill in the whole matrix
                ! FIXME : create a specific calc_case for this test?
                !
                IF (calc_case == 1) THEN
                    IF((error(N_loop,M/2,i_steep,2)-error(N_loop,M/2,i_steep,6)).GT.0.05*error(N_loop,M/2,i_steep,6)) THEN
                        print*, 'Worsening of previous HOS computation: total dealiasing'
                        print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
                        print*, 'location', MAXLOC(ABS(phizRF_2D(1:n1,1:n2)-phiz(1:n1,1:n2)))
                        err_xy = 1
                    ELSEIF ((-error(N_loop,M/2,i_steep,2)+error(N_loop,M/2,i_steep,6)) &
                   		-0.05*error(N_loop,M/2,i_steep,6).GT.tiny) THEN
                   		IF (ABS(error(N_loop,M/2,i_steep,6)-1.0_rp).GT.tiny) THEN
                       		print*, 'Enhancement of previous HOS computation: total dealiasing'
                       		print*, 'N=',Npts/2,'M=',M,error(N_loop,M/2,i_steep,2),error(N_loop,M/2,i_steep,6)
                       	ENDIF
                    ENDIF
                    IF (i_steep == 5) THEN
                          IF (M.EQ.2) THEN
                            OPEN(666,FILE='3d_test.dat')
                            WRITE(666,'(A)') 'VARIABLES="x", "y", "eta", "phis", "W", "W_RF", "err"'
                            !
                            OPEN(777,FILE='a_3d_test.dat')
                            WRITE(777,'(A)') 'VARIABLES="kx", "ky", "a_eta", "a_phis", "a_W", "a_W_RF", "err"'
                          ENDIF
                          WRITE(666,'(A,I3,A,I4,A,I4)')'ZONE T = "',M,'", I=',n1,', J=',n2
                          DO i2 = 1, n2
                             DO i1 = 1, n1
                                WRITE(666,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
            	                    eta(i1,i2), phis(i1,i2), phiz(i1,i2),phizRF_2D(i1,i2), phiz(i1,i2)-phizRF_2D(i1,i2)
                             END DO
                          END DO
                          !
                          WRITE(777,'(A,I3,A,I4,A,I4)')'ZONE T = "',M,'", I=',n1o2p1,', J=',n2
                          DO i2 = n2o2p1+1, n2
                             DO i1 = 1, n1o2p1
			                    WRITE(777,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
                                    ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2)), ABS(a_phiz(i1,i2)), ABS(a_phiz_RF(i1,i2)), &
                                    	ABS(a_phiz(i1,i2)-a_phiz_RF(i1,i2))
                             END DO
                          END DO
                          DO i2 = 1, n2o2p1
                             DO i1 = 1, n1o2p1
			                    WRITE(777,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
                                    ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2)), ABS(a_phiz(i1,i2)), ABS(a_phiz_RF(i1,i2)), &
                                    	ABS(a_phiz(i1,i2)-a_phiz_RF(i1,i2))
                             END DO
                          END DO
                    ENDIF
                ENDIF
                !
                IF ((M.EQ.8).AND.(Npts.EQ.16).AND.(i_steep.EQ.3)) THEN
					IF ((ABS(test_t_xy).EQ.1)) THEN !test the time integration...
						!
						CALL fill_butcher_array(RK_param)
						!
						time_cur    = 0.0_rp
						T_stop_star = 1000.0_rp * RF_obj%T
						dt_out      = RF_obj%T
						!
						IF (n2 == 1) THEN
						   dt_rk4 = 2.0_rp / SQRT(PI) * SQRT(xlen / n1o2p1)
						ELSE
						   dt_rk4 = 2.0_rp / SQRT(PI) * MIN(SQRT(xlen / n1o2p1), SQRT(ylen / n2o2p1))
						END IF
						dt_lin = 10.0_rp * dt_rk4
						dt     = 0.5_rp * dt_rk4
						!
						toler = 1.e-12 !1e-7
						!
						a_eta_ref(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
						a_phis_ref(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
						!
						idx_max(1:2) = MAXLOC(ABS(a_eta_ref))
						! Initiate da_eta for energy calculation
						CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
						!
						a_eta(1:n1o2p1,1:n2)  = a_eta_ref(1:n1o2p1,1:n2)
						a_phis(1:n1o2p1,1:n2) = a_phis_ref(1:n1o2p1,1:n2)
						!
						DO WHILE (time_cur-T_stop_star <= tiny)
						   !
						   !
						   ! Output of volume and energy
						   volume = REAL(a_eta(1,1),RP)
						   energy = calc_energy(a_eta, a_phis, da_eta)
						   !
						   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
						   !
						  CALL fourier_2_space(a_eta, eta)
						  CALL fourier_2_space(a_phis,phis)
						  !
						  IF (ABS(time_cur).LT.tiny) THEN
							OPEN(66,FILE='3d_xy.dat')
							WRITE(66,'(A)') 'VARIABLES="x", "y", "eta", "phis"'
							!
							OPEN(77,FILE='a_3d_xy.dat')
							WRITE(77,'(A)') 'VARIABLES="kx", "ky", "a_eta", "a_phis"'
							!
							OPEN(88,FILE='error_time_xy.dat')
							WRITE(88,'(A)') 'VARIABLES="time", "err amp", "err phase"'
						  ENDIF
						  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
						  DO i2 = 1, n2
							 DO i1 = 1, n1
								WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
									eta(i1,i2), phis(i1,i2)
							 END DO
						  END DO
						  !
						  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
						  DO i2 = n2o2p1+1, n2
							 DO i1 = 1, n1o2p1
								WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky_n2(i2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  DO i2 = 1, n2o2p1
							 DO i1 = 1, n1o2p1
								 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky_n2(i2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  !
						  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
							ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
							ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
							-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
						   !
						   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
						   !
						   ! Going to next time step
						   time_next = time_cur + dt_out
						   IF (time_next > T_stop_star) time_next = T_stop_star ! if last output
						   h_rk    = dt ! starting time step
						   DO WHILE(time_cur < time_next)
							  h_loc = h_rk ! local time step
							  IF (time_next < time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
							  !
							  ! Analytical integration of the linear part
							  ! starting at current time
							  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
							  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
							  !
							  ! Going to time_cur + h_loc
							  ! adaptive time step
							  !
							  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
								a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
							  error_old = error2
							  IF (error2 > toler) THEN
								dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
							  ELSE
								a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
								a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
								time_cur   = time_cur + h_loc
								dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
							  END IF
							  ! Step size ratio bounds (Mathematica)
							  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
							  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
							  !
							  ! Step size security factor
							  dt_correc = 0.98_rp * dt_correc
							  ! New step size
							  h_rk      = h_rk * dt_correc
							  ! Stability checkings
							  h_rk      = MIN(h_rk, dt_out, dt_lin)
						   END DO
						   dt = h_rk ! saving the step size for next start
						END DO
						IF (test_t_xy.EQ.1) THEN
							CLOSE(66)
							CLOSE(77)
							CLOSE(88)
						ENDIF
					ENDIF
                	!
                	IF ((test_t_xy.EQ.-1)) THEN !test the back-time integration...
						!
						CALL fill_butcher_array(RK_param)
						!
						! time_cur defined as current time
						T_stop_star = 0.0_rp
						dt_out      = -RF_obj%T
						!
						dt     = -0.5_rp * dt_rk4 !starting point identical to front propagation
						!
						! Use same tolerance
						DO WHILE (time_cur > T_stop_star)
						   !
						   !
						   ! Output of volume and energy
						   volume = REAL(a_eta(1,1),RP)
						   energy = calc_energy(a_eta, a_phis, da_eta)
						   !
						   print*,'time cur =',time_cur,' T_stop =', T_stop_star, 'vol =', volume, 'energy = ',energy(4)
						   !
						  CALL fourier_2_space(a_eta, eta)
						  CALL fourier_2_space(a_phis,phis)
						  !
						  WRITE(66,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1,', J=',n2
						  DO i2 = 1, n2
							 DO i1 = 1, n1
								WRITE(66,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') x(i1), y(i2), &
									eta(i1,i2), phis(i1,i2)
							 END DO
						  END DO
						  !
						  WRITE(77,'(A,ES9.2,A,I4,A,I4)')'ZONE T = "',time_cur,'", I=',n1o2p1,', J=',n2
						  DO i2 = n2o2p1+1, n2
							 DO i1 = 1, n1o2p1
								WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), - ky(n2-i2+2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  DO i2 = 1, n2o2p1
							 DO i1 = 1, n1o2p1
								 WRITE(77,'(ES10.3,X,ES10.3,X,ES10.3,X,ES10.3)') kx(i1), ky(i2), &
									ABS(a_eta(i1,i2)), ABS(a_phis(i1,i2))
							 END DO
						  END DO
						  !
						  WRITE(88,'(ES10.3,X,ES10.3,X,ES10.3)') time_cur, &
							ABS(ABS(a_eta(idx_max(1),idx_max(2)))-ABS(a_eta_ref(idx_max(1),idx_max(2)))) ,&
							ATAN2(AIMAG(a_eta(idx_max(1),idx_max(2))),REAL(a_eta(idx_max(1),idx_max(2)),RP)) &
							-ATAN2(AIMAG(a_eta_ref(idx_max(1),idx_max(2))),REAL(a_eta_ref(idx_max(1),idx_max(2)),RP))
						   !
						   IF (ABS(time_cur-T_stop_star) .LT. tiny) EXIT ! output of the last zone is done
						   !
						   ! Going to next time step
						   time_next = time_cur + dt_out
						   IF (time_next < T_stop_star) time_next = T_stop_star ! if last output
						   h_rk    = dt ! starting time step
						   DO WHILE(time_cur > time_next)
							  h_loc = h_rk ! local time step
							  IF (time_next > time_cur + h_rk) h_loc = time_next - time_cur ! if last time step
							  !
							  ! Analytical integration of the linear part
							  ! starting at current time
							  a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
							  a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
							  !
							  ! Going to time_cur + h_loc
							  ! adaptive time step
							  !
							  CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, time_cur, h_loc, &
								a_phis_rk, a_eta_rk, da_eta, error2, 1.0_rp, 1.0_rp)
							  error_old = error2
							  IF (error2 > toler) THEN
								dt_correc = (error2/toler)**(-(1.0_rp/(RK_param%p-1)))
							  ELSE
								a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
								a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
								time_cur   = time_cur + h_loc
								dt_correc  = (error2/toler)**(-(1.0_rp/(RK_param%p)))
							  END IF
							  ! Step size ratio bounds (Mathematica)
							  IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
							  IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
							  !
							  ! Step size security factor
							  dt_correc = 0.98_rp * dt_correc
							  ! New step size
							  h_rk      = ABS(h_rk) * dt_correc
							  ! Stability checkings
							  h_rk      = -MIN(h_rk, -dt_out, dt_lin)
						   END DO
						   dt = h_rk ! saving the step size for next start
						END DO
						!
						CLOSE(66)
						CLOSE(77)
						CLOSE(88)
					ENDIF
                ENDIF
                !
                IF ((test_vel.EQ.1).AND.(M.EQ.6).AND.(Npts.EQ.64).AND.(i_steep.EQ.3)) THEN !Test the velocity calculations
                	! Initiate da_eta for energy calculation
                    CALL fill_butcher_array(RK_param)
                    CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
                    CALL HOSvel2(25,M,a_eta,a_phis,0.0_rp) !HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current)
                    CALL reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect,a_eta,da_eta,error_SL)
                    ! FIXME : compare to RF ?
                    !stop
                ENDIF
                !
                CALL fourier_end(3)
                !
                part = 4
                IF (part >= (4*N)/N-1) THEN
                   ! Best dealiasing
                   Nd1     = ((part+1) * n1) / 2
                   Nd1o2p1 = Nd1/2+1
                   order_max(1) = part-1
                   i_dealias(1) = 1
                   N_dea(1)     = n1o2p1
                   !
                   Nd2     = ((part+1) * n2) / 2
                   Nd2o2p1 = Nd2/2+1
                   Nd2p1o2 = (Nd2+1)/2
                   order_max(2) = part-1
                   i_dealias(2) = 1
                   N_dea(2)     = n2o2p1

                   CALL build_derivatives()
                   CALL fourier_ini(3)
                   CALL space_2_fourier(eta,a_eta)
                   CALL space_2_fourier(phis,a_phis)
                   ! calculation of the horizontal derivatives of phis and eta
                   call phisxy_etaxy(a_phis, a_eta)
                   ! resolution of the different order of the FS potential boundary value problem to get modes
                   call HOSphis_modes_fully_dealiased(a_phis)
                   error(N_loop,M/2,i_steep,3) = MAXVAL(ABS(phizRF_2D(1:n1,1:n2) - phiz(1:n1,1:n2)))
                   !
                   CALL fourier_end(3)
                   !
                END IF

                DEALLOCATE(phis, eta, phizRF_2D)
                DEALLOCATE(oneoj,etapm_ext)
                DEALLOCATE(kth,kth_all)
             END IF
          END DO
       END DO
    END DO
    !
    !
    WRITE(*,'(A)') 'Writing results file'
    !
    WRITE(*,'(A)') 'result_DYue_xy.dat'
    !
    OPEN(1,FILE='result_DYue_xy.dat')
    CALL write_result(1, error(:,:,:,1))
    !
    WRITE(1,*)
    CALL write_result(1, error(:,:,:,2))
    !
    WRITE(1,*)
    CALL write_result(1, error(:,:,:,3))
    !
    CLOSE(1)
    !
    WRITE(*,'(A)') 'comp_DYue_xy.dat'
    !
    OPEN(2,FILE='comp_DYue_xy.dat')
    !
    IF (calc_case == 1) THEN
       CALL write_comp(2, error(:,:,:,1), error(:,:,:,4))
    ELSE
       CALL write_comp(2, error(:,:,:,2), error(:,:,:,1))
    END IF
    !
    DO N_loop = 1,120
       WRITE(2,'(A)', ADVANCE='NO') '_'
    END DO
    WRITE(2,'(A)') ''
    WRITE(2,*)
    !
    IF (calc_case == 1) THEN
       CALL write_comp(2, error(:,:,:,2), error(:,:,:,4))
    ELSE
       CALL write_comp(2, error(:,:,:,3), error(:,:,:,1))
    END IF
    !
    IF (calc_case == 1) THEN
       DO N_loop = 1,120
          WRITE(2,'(A)', ADVANCE='NO') '_'
       END DO
       WRITE(2,'(A)') ''
       WRITE(2,*)
       !
       CALL write_comp(2, error(:,:,:,3), error(:,:,:,4))
    END IF
    !
    CLOSE(2)
    !
    IF (calc_case == 1) THEN
       WRITE(*,'(A)') 'result_DYue_xy.tex'
       !
       OPEN(1,FILE='result_DYue_xy.tex')
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,4) > -tiny) THEN
                      WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,4)
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, Dommermuth et Yue}'
       WRITE(1,'(A)') '\label{tab:error-W-DYue}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,1) > -tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.05_rp*error(N_loop,M/2,i_steep,1) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,1)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,1)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement partiel}'
       WRITE(1,'(A)') '\label{tab:error-W-simple}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,2) > -tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.1_rp*error(N_loop,M/2,i_steep,2) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,2)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,2)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement complet}'
       WRITE(1,'(A)') '\label{tab:error-W-complet}'
       WRITE(1,'(A)') '\end{table}'
       !
       WRITE(1,'(A)') '\begin{table}[!htbp]'
       WRITE(1,'(A)') '\begin{center}'
       WRITE(1,'(A)') '\begin{tabular}{|l|l|ccccccc|}'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '\multicolumn{9}{|c|}{$M$}\\'
       WRITE(1,'(A)') '\hline'
       WRITE(1,'(A)') '$\varepsilon$ & $N$ & 2& 4& 6& 8& 10& 12& 14 \\'
       WRITE(1,'(A)') '\hline'
       DO i_steep = 1,5
          DO N_loop = 1,Nma
             IF (NM(N_loop,1,i_steep) /= 0) THEN
                IF (N_loop == 1) THEN
                   WRITE(1,'(F5.2)', ADVANCE='NO') steepness(i_steep)
                ELSE
                   WRITE(1,'(X)', ADVANCE='NO')
                END IF
                WRITE(1,'(A)', ADVANCE='NO') '&'
                Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(1,'(I3)', ADVANCE='NO') Npts
                DO M = 2,14,2
                   IF (error(N_loop,M/2,i_steep,3) > tiny) THEN
                      IF (error(N_loop,M/2,i_steep,4) - 1.05_rp*error(N_loop,M/2,i_steep,3) > -tiny) THEN
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&\cellcolor[gray]{0.9}', error(N_loop,M/2,i_steep,3)
                      ELSE
                         WRITE(1,'(A,E10.2)', ADVANCE='NO') '&', error(N_loop,M/2,i_steep,3)
                      END IF
                   ELSE
                      WRITE(1,'(A)', ADVANCE='NO') '&'
                   END IF
                END DO
                WRITE(1,'(A)', ADVANCE='NO') '\\'
                WRITE(1,*)
             END IF
          END DO
          WRITE(1,'(A)') '\hline'
       END DO
       WRITE(1,'(A)') '\end{tabular}'
       WRITE(1,'(A)') '\end{center}'
       WRITE(1,'(A,I2,A)') '\caption{Erreur absolue maximum sur $W$, anti-repliement interm\''ediaire  p=',part,'}'
       WRITE(1,'(A)') '\label{tab:error-W-partiel}'
       WRITE(1,'(A)') '\end{table}'
       !
       CLOSE(1)
    END IF
    !
    WRITE(*,'(A)') 'result_DYue_xy.tec'
    !
    WHERE (error < tiny)
       log_error = 0.0_rp
    ELSEWHERE
       log_error = LOG(error) / LOG(10.0_rp)
    END WHERE
    OPEN(3,FILE='result_DYue_xy.tec')
    WRITE(3,'(A)') 'VARIABLES="M", "N", "error", "log_1_0(error)"'
    IF (calc_case == 1) THEN
       DO i_steep = 1,5
          WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="DYue',steepness(i_steep),'", I=',Nma,', J=',Mm
             DO M = 2,2*Mm,2
          DO N_loop = 1,Nma
             Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
                WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,4),log_error(N_loop,M/2,i_steep,4)
             END DO
          END DO
       END DO
    END IF
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,1),log_error(N_loop,M/2,i_steep,1)
          END DO
       END DO
    END DO
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="comp',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,2),log_error(N_loop,M/2,i_steep,2)
          END DO
       END DO
    END DO
    DO i_steep = 1,5
       WRITE(3,'(A,F5.2,A,I4,A,I4)') 'ZONE T="part',steepness(i_steep),'", I=',Nma,', J=',Mm
          DO M = 2,2*Mm,2
       DO N_loop = 1,Nma
          Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
             WRITE(3,'(I3,X,I3,X,ES10.3,X,ES10.3)') M,Npts,error(N_loop,M/2,i_steep,3),log_error(N_loop,M/2,i_steep,3)
          END DO
       END DO
    END DO
    CLOSE(3)
ENDIF

IF (test_x == 1) THEN
    IF (err_x == 0) THEN
        print*,'*'
        print*,'***************** Successful computation along x direction ! *****************'
        print*,'*'
    ELSE
        print*,'*'
        print*,'***************** A probable error has appeared in x direction, CHECK ! *****************'
        print*,'*'
    ENDIF
ENDIF

IF (test_y == 1) THEN
    IF (err_y == 0) THEN
        print*,'*'
        print*,'***************** Successful computation along y direction ! *****************'
        print*,'*'
    ELSE
        print*,'*'
        print*,'***************** A probable error has appeared in y direction, CHECK ! *****************'
        print*,'*'
    ENDIF
ENDIF

IF (test_xy == 1) THEN
    IF (err_xy == 0) THEN
        print*,'*'
        print*,'***************** Successful computation along xy direction ! *****************'
        print*,'*'
    ELSE
        print*,'*'
        print*,'***************** A probable error has appeared in xy direction, CHECK ! *****************'
        print*,'*'
    ENDIF
ENDIF




CONTAINS

! !
! !
! !
SUBROUTINE write_comp(unit, error, ref)
!
IMPLICIT NONE
!
INTEGER :: unit
REAL(RP), DIMENSION(Nma,Mm,5) :: error, ref
!
WRITE(2,'(A)') '  eps   N         2         4         6         8        10        12        14        16        18        20'
DO i_steep = 1,5
   IF (i_steep > 1)    WRITE(2,*)
   DO N_loop = 1,Nma
      IF (NM(N_loop,1,i_steep) /= 0) THEN
         IF (N_loop == 1) THEN
            WRITE(unit,'(F5.2,X)', ADVANCE='NO') steepness(i_steep)
         ELSE
            WRITE(unit,'(6X)', ADVANCE='NO')
         END IF
         Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
         WRITE(2,'(I3,X)', ADVANCE='NO') Npts
         DO M = 2,2*Mm,2
            IF (ref(N_loop,M/2,i_steep) > -tiny) THEN
               IF (error(N_loop,M/2,i_steep) > -tiny) THEN
                  WRITE(unit,'(F10.0)', ADVANCE='NO') &
                      (ref(N_loop,M/2,i_steep)-error(N_loop,M/2,i_steep))/error(N_loop,M/2,i_steep)*100.0_rp
               ELSE
                  WRITE(unit,'(A10)', ADVANCE='NO') '             '
               END IF
            END IF
         END DO
         WRITE(unit,*)
      END IF
   END DO
END DO

END SUBROUTINE write_comp
!
!
!
SUBROUTINE write_result(unit, error)
!
IMPLICIT NONE
!
INTEGER :: unit
REAL(RP), DIMENSION(Nma,Mm,5) :: error
!
WRITE(1,'(A)') '  eps   N      2         4         6         8        10        12        14        16        18        20'
DO i_steep = 1,5
   DO N_loop = 1,Nma
      IF (NM(N_loop,1,i_steep) /= 0) THEN
         IF (N_loop == 1) THEN
            WRITE(unit,'(F5.2,X)', ADVANCE='NO') steepness(i_steep)
         ELSE
            WRITE(unit,'(6X)', ADVANCE='NO')
         END IF
         Npts = NM(N_loop,1,i_steep) ! for comparison with DYue
         WRITE(1,'(I3,X)', ADVANCE='NO') Npts
         DO M = 2,2*Mm,2
            IF (error(N_loop,M/2,i_steep) > -tiny) THEN
               WRITE(unit,'(E10.2)', ADVANCE='NO') error(N_loop,M/2,i_steep)
            ELSE
               WRITE(unit,'(A10)', ADVANCE='NO') '          '
            END IF
         END DO
         WRITE(unit,*)
      END IF
   END DO
   WRITE(1,*)
END DO
END SUBROUTINE write_result
!
!
SUBROUTINE check_slope
!
IMPLICIT NONE
!
REAL(RP) :: slopex_max, threshold,slopey_max
!
threshold = 2.0_rp
!
slopex_max = MAXVAL(ABS(etax))
slopey_max = MAXVAL(ABS(etay))
!
IF (slopex_max > threshold) THEN
   WRITE(*,'(A)') 'Ca va trancher en x...'
   STOP
END IF
!
IF (slopey_max > threshold) THEN
  WRITE(*,'(A)') 'Ca va trancher en y...'
   STOP
!   pause
!   WRITE(*,'(A)') 'Trying to filter eta and phis...'
!   IF (n2 /= 1) THEN
!      call filtering_y(phis,0.7_rp)
!      call filtering_y(eta,0.7_rp)
!   END IF
END IF
!
END SUBROUTINE check_slope
!
!
!
SUBROUTINE display_parameters()
!
IMPLICIT NONE
!
WRITE(*,'(A)') 'General parameters'
WRITE(*,'(3X, A, F5.2)') 'Steepness       ',steepness(i_steep)
WRITE(*,'(3X, A, I4)')   'HOS order        ',M
!
WRITE(*,'(A)') 'Parameters for the x-direction'
WRITE(*,'(3X, A, I4)')  'Modes              n1=',n1
WRITE(*,'(3X, A, I4)')  'Modes used        n1c=',n1c
! IF (ctype1(1) == 1) THEN
!    WRITE(*,'(3X,A)')    'No dealiasing'
! ELSE
!    WRITE(*,'(3X,A,I2)') 'Dealiasing the products of order ', ctype1(1)
! END IF
WRITE(*,'(3X, A, I4)')  'Extended modes    Nd1=',Nd1
WRITE(*,'(3X, A, I4)')  'Dealiased modes N_dea=',N_dea(1)
WRITE(*,'(3X, A, I4)')  'Derived modes   N_der=',N_der(1)
IF (i_filt(1) == 0) THEN
   WRITE(*,'(3X,A)')    'No filtering'
ELSE
   WRITE(*,'(3X,A,I4)') 'Filtering with modes  ', n1c_filt
END IF
WRITE(*,'(3X,A,I4)')    'Maximum harmonics ', n1/2/n_lambda_x
WRITE(*,'(3X,A,I4)')    'Used harmonics    ', n1c/n_lambda_x
IF (i_filt(1) == 1) THEN
   WRITE(*,'(3X,A,I4)') 'Kept harmonics    ', n1c_filt/n_lambda_x
END IF
IF (n1c_filt/n_lambda_x < M) WRITE(*,'(A)') 'Not enough harmonics'
!
WRITE(*,'(A)') 'Parameters for the y-direction'
WRITE(*,'(3X, A, I4)')  'Modes              n2=',n2
WRITE(*,'(3X, A, I4)')  'Modes used        n2c=',n2c
! IF (ctype2(1) == 1) THEN
!    WRITE(*,'(3X,A)')    'No dealiasing'
! ELSE
!    WRITE(*,'(3X,A,I4)') 'Dealiasing the products of order ', ctype2(1)
! END IF
WRITE(*,'(3X, A, I4)')  'Extended modes   Nd2=',Nd2
WRITE(*,'(3X, A, I4)')  'Dealiased modes N_dea=',N_dea(2)
WRITE(*,'(3X, A, I4)')  'Derived modes   N_der=',N_der(2)
IF (i_filt(2) == 0) THEN
   WRITE(*,'(3X,A)')    'No filtering'
ELSE
   WRITE(*,'(3X,A,I4)') 'Filtering with modes ', n2c_filt
END IF
IF (n2 /= 1) THEN
   WRITE(*,'(3X,A,I4)')    'Maximum harmonics ', n2/2/n_lambda_y
   WRITE(*,'(3X,A,I4)')    'Used harmonics    ', n2c/n_lambda_y
   IF (i_filt(2) == 1) THEN
      WRITE(*,'(3X,A,I4)') 'Kept harmonics    ', n2c_filt/n_lambda_y
   END IF
   IF (n2c_filt/2/n_lambda_y < M) WRITE(*,'(A)') 'Not enough harmonics'
END IF
!
END SUBROUTINE display_parameters
!
!
!
SUBROUTINE fill_DYue_cases(steepness, files, NM, comp_case)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(:), INTENT(OUT)           :: steepness
CHARACTER(LEN=100), DIMENSION(:), INTENT(OUT) :: files
INTEGER, DIMENSION(:,:,:), INTENT(OUT)        :: NM
INTEGER :: comp_case

steepness = (/0.1, 0.2, 0.3, 0.35, 0.4/)
files     = (/'waverf_L628_inf_ka01_N15_30.cof ','waverf_L628_inf_ka02_N20_40.cof ',&
               'waverf_L628_inf_ka03_N25_50.cof ', 'waverf_L628_inf_ka035_N25_50.cof', &
                'waverf_L628_inf_ka04_N50_100.cof'/)
! Filling the table of N for M=2*i2 and steepness
IF (comp_case == 1) THEN
   NM(:,:,:) = 0
   NM(1,1:5,1) = 8  ! ka=0.1
   NM(2,1:5,1) = 16
   NM(1,1:5,2) = 8  ! ka=0.2
   NM(2,1:5,2) = 16
   NM(3,1:6,2) = 32
   NM(1,1:5,3) = 8  ! ka=0.3
   NM(2,1:5,3) = 16
   NM(3,1:5,3) = 32
   NM(4,1:7,3) = 64
   NM(1,1:5,4) = 8  ! ka=0.35
   NM(2,1:5,4) = 16
   NM(3,1:6,4) = 32
   NM(4,1:7,4) = 64
   NM(1,1:7,5) = 32 ! ka=0.4
   NM(2,1:7,5) = 64
   NM(3,1:7,5) = 128
   error(:,:,:,4) = -1.0_rp
   error(1,1,1,4) = 0.00075  ! ka=0.1, N=8
   error(1,2,1,4) = 0.0000068
   error(1,3,1,4) = 0.000000072
   error(1,4,1,4) = 0.0000000022
   error(1,5,1,4) = 0.000000001
   error(2,1,1,4) = 0.00075   ! ka=0.1, N=16
   error(2,2,1,4) = 0.0000068
   error(2,3,1,4) = 0.000000065
   error(2,4,1,4) = 0.00000000064
   error(2,5,1,4) = 0.000000000049
   error(1,1,2,4) = 0.0059  ! ka=0.2, N=8
   error(1,2,2,4) = 0.00022
   error(1,3,2,4) = 0.000015
   error(1,4,2,4) = 0.0000018
   error(1,5,2,4) = 0.0000013
   error(2,1,2,4) = 0.006   ! ka=0.2, N=16
   error(2,2,2,4) = 0.00022
   error(2,3,2,4) = 0.0000087
   error(2,4,2,4) = 0.00000037
   error(2,5,2,4) = 0.000000038
   error(3,1,2,4) = 0.0006  ! ka=0.2, N=32
   error(3,2,2,4) = 0.00022
   error(3,3,2,4) = 0.0000088
   error(3,4,2,4) = 0.00000035
   error(3,5,2,4) = 0.000000014
   error(3,6,2,4) = 0.00000000075
   error(1,1,3,4) = 0.019  ! ka=0.3, N=8
   error(1,2,3,4) = 0.0022
   error(1,3,3,4) = 0.00047
   error(1,4,3,4) = 0.00014
   error(1,5,3,4) = 0.00016
   error(2,1,3,4) = 0.02  ! ka=0.3, N=16
   error(2,2,3,4) = 0.0018
   error(2,3,3,4) = 0.00019
   error(2,4,3,4) = 0.000059
   error(2,5,3,4) = 0.000024
   error(3,1,3,4) = 0.02  ! ka=0.3, N=32
   error(3,2,3,4) = 0.0018
   error(3,3,3,4) = 0.00017
   error(3,4,3,4) = 0.000016
   error(3,5,3,4) = 0.0000017
   error(4,1,3,4) = 0.02  ! ka=0.3, N=64
   error(4,2,3,4) = 0.0018
   error(4,3,3,4) = 0.00017
   error(4,4,3,4) = 0.000016
   error(4,5,3,4) = 0.0000016
   error(4,6,3,4) = 0.00000021
   error(4,7,3,4) = 0.000000033
   error(1,1,4,4) = 0.031  ! ka=0.35, N=8
   error(1,2,4,4) = 0.0064
   error(1,3,4,4) = 0.0022
   error(1,4,4,4) = 0.0013
   error(1,5,4,4) = 0.0013
   error(2,1,4,4) = 0.031  ! ka=0.35, N=16
   error(2,2,4,4) = 0.0041
   error(2,3,4,4) = 0.00099
   error(2,4,4,4) = 0.00071
   error(2,5,4,4) = 0.00022
   error(3,1,4,4) = 0.031  ! ka=0.35, N=32
   error(3,2,4,4) = 0.004
   error(3,3,4,4) = 0.00053
   error(3,4,4,4) = 0.000094
   error(3,5,4,4) = 0.000095
   error(3,6,4,4) = 0.00016
   error(4,1,4,4) = 0.031  ! ka=0.35, N=64
   error(4,2,4,4) = 0.004
   error(4,3,4,4) = 0.00053
   error(4,4,4,4) = 0.000073
   error(4,5,4,4) = 0.000011
   error(4,6,4,4) = 0.0000038
   error(4,7,4,4) = 0.00068
   error(1,1,5,4) = 0.045 ! ka=0.4,N=32
   error(1,2,5,4) = 0.0079
   error(1,3,5,4) = 0.0028
   error(1,4,5,4) = 0.0081
   error(2,1,5,4) = 0.045 ! ka=0.4,N=64
   error(2,2,5,4) = 0.0079
   error(2,3,5,4) = 0.0015
   error(2,4,5,4) = 0.00035
   error(2,5,5,4) = 0.00091
   error(3,1,5,4) = 0.045 ! ka=0.4,N=128
   error(3,2,5,4) = 0.0079
   error(3,3,5,4) = 0.0015
   error(3,4,5,4) = 0.0003
   error(3,5,5,4) = 0.00089
   !
   ! Initial results from HOS to compare with...
   ! Partial dealiasing
   !
   error(:,:,:,5) = 1.0_rp
   ! ka = 0.1
   error(1,1,1,5) = 0.72E-03  ! ka=0.1, N=8
   error(1,2,1,5) = 0.64E-05
   error(1,3,1,5) = 0.62E-07
   error(1,4,1,5) = 0.36E-08
   error(1,5,1,5) = 0.32E-08
   error(1,6,1,5) = 0.32E-08
   error(1,7,1,5) = 0.32E-08
   error(2,1,1,5) = 0.74E-03  ! ka=0.1, N=16
   error(2,2,1,5) = 0.68E-05
   error(2,3,1,5) = 0.65E-07
   error(2,4,1,5) = 0.62E-09
   error(2,5,1,5) = 0.60E-11
   error(2,6,1,5) = 0.59E-13
   error(2,7,1,5) = 0.12E-14
   ! ka = 0.2
   error(1,1,2,5) = 0.59E-02  ! ka=0.2, N=8
   error(1,2,2,5) = 0.23E-03
   error(1,3,2,5) = 0.84E-05
   error(1,4,2,5) = 0.27E-05
   error(1,5,2,5) = 0.25E-05
   error(1,6,2,5) = 0.25E-05
   error(1,7,2,5) = 0.25E-05
   error(2,1,2,5) = 0.59E-02  ! ka=0.2, N=16
   error(2,2,2,5) = 0.22E-03
   error(2,3,2,5) = 0.87E-05
   error(2,4,2,5) = 0.34E-06
   error(2,5,2,5) = 0.16E-07
   error(2,6,2,5) = 0.59E-09
   error(2,7,2,5) = 0.54E-10
   error(3,1,2,5) = 0.60E-02  ! ka=0.2, N=32
   error(3,2,2,5) = 0.22E-03
   error(3,3,2,5) = 0.87E-05
   error(3,4,2,5) = 0.34E-06
   error(3,5,2,5) = 0.14E-07
   error(3,6,2,5) = 0.55E-09
   error(3,7,2,5) = 0.23E-10
   ! ka = 0.3
   error(1,1,3,5) = 0.19E-01  ! ka=0.3, N=8
   error(1,2,3,5) = 0.20E-02
   error(1,3,3,5) = 0.33E-03
   error(1,4,3,5) = 0.19E-03
   error(1,5,3,5) = 0.21E-03
   error(1,6,3,5) = 0.21E-03
   error(1,7,3,5) = 0.21E-03
   error(2,1,3,5) = 0.20E-01  ! ka=0.3, N=16
   error(2,2,3,5) = 0.18E-02
   error(2,3,3,5) = 0.16E-03
   error(2,4,3,5) = 0.24E-04
   error(2,5,3,5) = 0.28E-05
   error(2,6,3,5) = 0.30E-06
   error(2,7,3,5) = 0.67E-06
   error(3,1,3,5) = 0.20E-01  ! ka=0.3, N=32
   error(3,2,3,5) = 0.18E-02
   error(3,3,3,5) = 0.17E-03
   error(3,4,3,5) = 0.16E-04
   error(3,5,3,5) = 0.16E-05
   error(3,6,3,5) = 0.24E-06
   error(3,7,3,5) = 0.13E-06
   error(4,1,3,5) = 0.20E-01  ! ka=0.3, N=64
   error(4,2,3,5) = 0.18E-02
   error(4,3,3,5) = 0.17E-03
   error(4,4,3,5) = 0.16E-04
   error(4,5,3,5) = 0.15E-05
   error(4,6,3,5) = 0.15E-06
   error(4,7,3,5) = 0.15E-07
   ! ka = 0.35
   error(1,1,4,5) = 0.31E-01  ! ka=0.35, N=8
   error(1,2,4,5) = 0.49E-02
   error(1,3,4,5) = 0.18E-02
   error(1,4,4,5) = 0.16E-02
   error(1,5,4,5) = 0.16E-02
   error(1,6,4,5) = 0.16E-02
   error(1,7,4,5) = 0.16E-02
   error(2,1,4,5) = 0.31E-01  ! ka=0.35, N=16
   error(2,2,4,5) = 0.40E-02
   error(2,3,4,5) = 0.69E-03
   error(2,4,4,5) = 0.19E-03
   error(2,5,4,5) = 0.28E-04
   error(2,6,4,5) = 0.37E-04
   error(2,7,4,5) = 0.37E-04
   error(3,1,4,5) = 0.31E-01  ! ka=0.35, N=32
   error(3,2,4,5) = 0.40E-02
   error(3,3,4,5) = 0.53E-03
   error(3,4,4,5) = 0.79E-04
   error(3,5,4,5) = 0.28E-04
   error(3,6,4,5) = 0.26E-04
   error(3,7,4,5) = 0.16E-04
   error(4,1,4,5) = 0.31E-01  ! ka=0.35, N=64
   error(4,2,4,5) = 0.40E-02
   error(4,3,4,5) = 0.53E-03
   error(4,4,4,5) = 0.73E-04
   error(4,5,4,5) = 0.10E-04
   error(4,6,4,5) = 0.14E-05
   error(4,7,4,5) = 0.25E-06
   ! ka = 0.4
   error(1,1,5,5) = 0.45E-01 ! ka=0.4,N=32
   error(1,2,5,5) = 0.79E-02
   error(1,3,5,5) = 0.20E-02 ! M=6
   error(1,4,5,5) = 0.21E-02 ! M=8
   error(1,5,5,5) = 0.25E-02 ! M=10
   error(1,6,5,5) = 0.18E-02 ! M=12
   error(1,7,5,5) = 0.16E-02 ! M=14
   error(2,1,5,5) = 0.45E-01 ! ka=0.4,N=64
   error(2,2,5,5) = 0.79E-02
   error(2,3,5,5) = 0.15E-02 ! M=6
   error(2,4,5,5) = 0.32E-03
   error(2,5,5,5) = 0.25E-03
   error(2,6,5,5) = 0.11E-02
   error(2,7,5,5) = 0.40E-02
   error(3,1,5,5) = 0.45E-01 ! ka=0.4,N=128
   error(3,2,5,5) = 0.79E-02
   error(3,3,5,5) = 0.15E-02 ! M=6
   error(3,4,5,5) = 0.30E-03
   error(3,5,5,5) = 0.60E-04
   error(3,6,5,5) = 0.12E-04
   error(3,7,5,5) = 0.23E-04 ! 0.17E-04 (Obtained result with y and x different... take y value)
   !
   ! Initial results from HOS to compare with...
   ! Total dealiasing
   !
   error(:,:,:,6) = 1.0_rp
   ! ka = 0.1
   error(1,1,1,6) = 0.72E-03  ! ka=0.1, N=8
   error(1,2,1,6) = 0.64E-05
   error(1,3,1,6) = 0.62E-07
   error(1,4,1,6) = 0.36E-08
   error(1,5,1,6) = 0.31E-08
   error(1,6,1,6) = 0.31E-08
   error(1,7,1,6) = 0.31E-08
   error(2,1,1,6) = 0.74E-03  ! ka=0.1, N=16
   error(2,2,1,6) = 0.68E-05
   error(2,3,1,6) = 0.65E-07
   error(2,4,1,6) = 0.62E-09
   error(2,5,1,6) = 0.60E-11
   error(2,6,1,6) = 0.59E-13
   error(2,7,1,6) = 0.12E-14
   ! ka = 0.2
   error(1,1,2,6) = 0.59E-02  ! ka=0.2, N=8
   error(1,2,2,6) = 0.23E-03
   error(1,3,2,6) = 0.75E-05
   error(1,4,2,6) = 0.25E-05
   error(1,5,2,6) = 0.21E-05
   error(1,6,2,6) = 0.22E-05
   error(1,7,2,6) = 0.22E-05
   error(2,1,2,6) = 0.59E-02  ! ka=0.2, N=16
   error(2,2,2,6) = 0.22E-03
   error(2,3,2,6) = 0.87E-05
   error(2,4,2,6) = 0.34E-06
   error(2,5,2,6) = 0.14E-07
   error(2,6,2,6) = 0.54E-09
   error(2,7,2,6) = 0.11E-09
   error(3,1,2,6) = 0.60E-02  ! ka=0.2, N=32
   error(3,2,2,6) = 0.22E-03
   error(3,3,2,6) = 0.87E-05
   error(3,4,2,6) = 0.34E-06
   error(3,5,2,6) = 0.14E-07
   error(3,6,2,6) = 0.55E-09
   error(3,7,2,6) = 0.22E-10
   ! ka = 0.3
   error(1,1,3,6) = 0.19E-01  ! ka=0.3, N=8
   error(1,2,3,6) = 0.18E-02
   error(1,3,3,6) = 0.25E-03
   error(1,4,3,6) = 0.15E-03
   error(1,5,3,6) = 0.15E-03
   error(1,6,3,6) = 0.15E-03
   error(1,7,3,6) = 0.15E-03
   error(2,1,3,6) = 0.20E-01  ! ka=0.3, N=16
   error(2,2,3,6) = 0.18E-02
   error(2,3,3,6) = 0.16E-03
   error(2,4,3,6) = 0.16E-04
   error(2,5,3,6) = 0.15E-05
   error(2,6,3,6) = 0.46E-06
   error(2,7,3,6) = 0.36E-06
   error(3,1,3,6) = 0.20E-01  ! ka=0.3, N=32
   error(3,2,3,6) = 0.18E-02
   error(3,3,3,6) = 0.17E-03
   error(3,4,3,6) = 0.16E-04
   error(3,5,3,6) = 0.15E-05
   error(3,6,3,6) = 0.15E-06
   error(3,7,3,6) = 0.15E-07
   error(4,1,3,6) = 0.20E-01  ! ka=0.3, N=64
   error(4,2,3,6) = 0.18E-02
   error(4,3,3,6) = 0.17E-03
   error(4,4,3,6) = 0.16E-04
   error(4,5,3,6) = 0.15E-05
   error(4,6,3,6) = 0.15E-06
   error(4,7,3,6) = 0.15E-07
   ! ka = 0.35
   error(1,1,4,6) = 0.31E-01  ! ka=0.35, N=8
   error(1,2,4,6) = 0.38E-02
   error(1,3,4,6) = 0.13E-02
   error(1,4,4,6) = 0.83E-03
   error(1,5,4,6) = 0.84E-03
   error(1,6,4,6) = 0.84E-03
   error(1,7,4,6) = 0.84E-03
   error(2,1,4,6) = 0.31E-01  ! ka=0.35, N=16
   error(2,2,4,6) = 0.39E-02
   error(2,3,4,6) = 0.51E-03
   error(2,4,4,6) = 0.76E-04
   error(2,5,4,6) = 0.21E-04
   error(2,6,4,6) = 0.12E-04
   error(2,7,4,6) = 0.12E-04
   error(3,1,4,6) = 0.31E-01  ! ka=0.35, N=32
   error(3,2,4,6) = 0.40E-02
   error(3,3,4,6) = 0.53E-03
   error(3,4,4,6) = 0.73E-04
   error(3,5,4,6) = 0.10E-04
   error(3,6,4,6) = 0.14E-05
   error(3,7,4,6) = 0.20E-06
   error(4,1,4,6) = 0.31E-01  ! ka=0.35, N=64
   error(4,2,4,6) = 0.40E-02
   error(4,3,4,6) = 0.53E-03
   error(4,4,4,6) = 0.73E-04
   error(4,5,4,6) = 0.10E-04
   error(4,6,4,6) = 0.14E-05
   error(4,7,4,6) = 0.20E-06
   ! ka = 0.4
   error(1,1,5,6) = 0.45E-01 ! ka=0.4,N=32
   error(1,2,5,6) = 0.78E-02
   error(1,3,5,6) = 0.15E-02  ! M=6
   error(1,4,5,6) = 0.30E-03  ! M=8
   error(1,5,5,6) = 0.61E-04  ! M=10
   error(1,6,5,6) = 0.15E-04  ! M=12
   error(1,7,5,6) = 0.62E-05  ! M=14
   error(2,1,5,6) = 0.45E-01 ! ka=0.4,N=64
   error(2,2,5,6) = 0.79E-02
   error(2,3,5,6) = 0.15E-02
   error(2,4,5,6) = 0.29E-03
   error(2,5,5,6) = 0.60E-04
   error(2,6,5,6) = 0.12E-04
   error(2,7,5,6) = 0.25E-05
   error(3,1,5,6) = 0.45E-01 ! ka=0.4,N=128
   error(3,2,5,6) = 0.79E-02
   error(3,3,5,6) = 0.15E-02
   error(3,4,5,6) = 0.30E-03
   error(3,5,5,6) = 0.60E-04
   error(3,6,5,6) = 0.12E-04
   error(3,7,5,6) = 0.25E-05
ELSE
! We need to have at least N > M
! M being close to the order of perturbation
! and N the number of modes (D and Yue)
! Here we have M=2*ind
   NM(:,:,:)      = 0
   NM(1,1:Mm,1:5)  = 4  ! ka=0.1 a 0.4
   NM(2,1:Mm,1:5)  = 8
   NM(3,1:Mm,1:5)  = 16
   NM(4,1:Mm,1:5)  = 32
   NM(5,1:Mm,1:5)  = 64
   NM(6,1:Mm,1:5)  = 128
!   NM(:,:,:)      = 0
!   NM(1,1:Mm,1:5)    = 4  ! ka=0.1 a 0.4
!   NM(2,1:Mm,1:5)  = 8
!   NM(3,1:Mm,1:5)  = 12
!   NM(4,1:Mm,1:5)  = 16
!   NM(5,1:Mm,1:5)  = 20
!   NM(6,1:Mm,1:5)  = 24
!   NM(7,1:Mm,1:5)  = 32
!   NM(8,1:Mm,1:5) = 40
!   NM(9,1:Mm,1:5) = 50
!   NM(10,1:Mm,1:5) = 64
!   NM(11,1:Mm,1:5) = 80
!   NM(12,1:Mm,1:5) = 96
!   NM(13,1:Mm,1:5) = 108
!   NM(14,1:Mm,1:5) = 128
   error(:,:,:,4) = -1.0_rp
END IF
!
END SUBROUTINE fill_DYue_cases
!
SUBROUTINE build_derivatives()
!
IMPLICIT NONE
!
INTEGER :: i1, i2, j
!
! Storage for derivatives
ikx = ((0.0_rp, 0.0_rp))
!  x-derivative on n1 points (i.e. n1o2p1 modes)
DO i2 = 1, n2o2p1
   ikx(1:MIN(N_der(1),n1o2p1),i2) = i * kx(1:MIN(N_der(1),n1o2p1))
END DO
! Last mode contains cos information only and must not be part of the differentiation.
IF (iseven(n1)) ikx(n1o2p1,:) = ((0.0_rp, 0.0_rp))
! negative ky
DO i2 = 2, n2o2p1
!FIXME : check
!DO i2 = 2, n2p1o2
   ikx(1:MIN(N_der(1),n1o2p1),n2-i2+2) = i * kx(1:MIN(N_der(1),n1o2p1))
END DO
!
iky = ((0.0_rp, 0.0_rp))
! y-derivative on n1 points (i.e. n1o2p1 modes)
DO i1 = 1, n1o2p1
   iky(i1,1:MIN(N_der(2),n2o2p1)) = i * ky_n2(1:MIN(N_der(2),n2o2p1))
   ! negative ky
   DO i2 = 2, MIN(N_der(2), n2o2p1)
  ! FIXME : check
  !DO i2 = 2, MIN(N_der(2), n2p1o2)
      iky(i1,n2-i2+2) = - i * ky_n2(i2)
   END DO
   IF (iseven(n2) .AND. N_der(2)>=n2o2p1) iky(i1,n2o2p1) = ((0.0_rp, 0.0_rp)) !FIXME : check this
END DO
!
!
ikx_big = ((0.0_rp, 0.0_rp))
!  x-derivative on Nd1 points (i.e. Nd1o2p1 modes)
DO i2 = 1, Nd2
  ikx_big(1:MIN(N_der(1),Nd1o2p1),i2) = i * kx(1:MIN(N_der(1),Nd1o2p1))
END DO
! Last mode contains cos information only and must not be part of the differentiation.
IF (iseven(Nd1) .AND. N_der(1) >= Nd1o2p1) ikx_big(Nd1o2p1,:) = ((0.0_rp, 0.0_rp))
!
iky_big = ((0.0_rp, 0.0_rp))
! y-derivative on Nd1 points (i.e. Nd1o2p1 modes)
DO i1 = 1, Nd1o2p1
  ! positive ky
  iky_big(i1,1:MIN(N_der(2),Nd2o2p1)) = i * ky(1:MIN(N_der(2),Nd2o2p1))
  ! negative ky
  DO i2 = 2, MIN(N_der(2),Nd2o2p1) !(Nd2+1)/2
     iky_big(i1,Nd2-i2+2) = - i * ky(i2)
  END DO
  ! Last mode contains cos information only and must not be part of the differentiation.
  IF (iseven(Nd2) .AND. N_der(2) >= Nd2o2p1) iky_big(i1,Nd2o2p1) = ((0.0_rp, 0.0_rp)) !FIXME : check this
END DO
!
DO i2=1,Nd2o2p1
	DO i1=1,Nd1o2p1
		DO j=1,M
			kth(i1,i2,j) = kth_all(i1,i2,j)
		ENDDO
	ENDDO
ENDDO
!
! kth(N_der(1)+1:Nd1o2p1,1:Nd2o2p1,1:M) = 0.0_rp
! kth(1:Nd1o2p1,N_der(2)+1:Nd2o2p1,1:M) = 0.0_rp
!
DO i2 = 2, Nd2p1o2
   DO i1 = 1, Nd1o2p1
   	  DO j=1,M
      	kth(i1,Nd2-i2+2,1:j) = kth(i1,i2,1:j)
      ENDDO
   END DO
END DO
!
END SUBROUTINE
!
SUBROUTINE check_range(n,name)
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
END SUBROUTINE check_range
!
SUBROUTINE build_RF_local(RF_obj, x, eta, phis, W)
!
IMPLICIT NONE
!
! Input
TYPE(RF_data), INTENT(INOUT)   :: RF_obj
REAL(RP), INTENT(IN)           :: x
REAL(RP), INTENT(OUT)          :: eta, phis, W
! Local
INTEGER  :: i_loop, j_loop, N_work
REAL(RP) :: kx, keta, jk
REAL(RP) :: x_w, eta_w, phis_w, W_w
!
! Direct solution on N_space points from the modes of phi and eta
! (aliased if N_space < N_phi or N_eta)
N_work = 1
!
x_w = x
!
!DO i_loop = 1, N_work
   i_loop=1
   ! temp
   kx = RF_obj%k * x_w
   ! Free surface elevation
   eta_w = RF_obj%A(1) * 0.5_rp
   DO j_loop = 1, RF_obj%N_eta
      eta_w = eta_w + RF_obj%A(j_loop+1) * COS(j_loop*kx)
   END DO
   ! temp
   keta = RF_obj%k * eta_w
   ! Free surface potential
   phis_w = (RF_obj%C + RF_obj%B(1)) * x_w
   IF (RF_obj%inf_depth) THEN
      DO j_loop = 1, RF_obj%N_phi
         phis_w = phis_w + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
      END DO
   ELSE
      DO j_loop = 1, RF_obj%N_phi
        IF (j_loop*RF_obj%k < 500.0_rp) THEN
           phis_w = phis_w + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * COSH(j_loop*(keta+RF_obj%k))/COSH(j_loop*RF_obj%k)
        ELSE
           phis_w = phis_w + RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(j_loop*keta)
        END IF
      END DO
   END IF
   ! Vertical velocity on the free surface
   W_w = 0.0_rp
   IF (RF_obj%inf_depth) THEN
      DO j_loop = 1, RF_obj%N_phi
         W_w = W_w + REAL(j_loop, RP) * RF_obj%k * RF_obj%B(j_loop+1) * SIN(REAL(j_loop, RP)*kx) * EXP(REAL(j_loop, RP)*keta)
      END DO
   ELSE
      DO j_loop = 1, RF_obj%N_phi
         jk = REAL(j_loop, RP) * RF_obj%k
         IF (jk < 500.0_rp) THEN
            W_w = W_w + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * SINH(REAL(j_loop, RP)*keta+jk)/COSH(jk)
         ELSE
            W_w = W_w + jk * RF_obj%B(j_loop+1) * SIN(j_loop*kx) * EXP(REAL(j_loop, RP)*keta)
         END IF
      END DO
   END IF
!END DO
!
eta  = eta_w
phis = phis_w
W    = W_w

END SUBROUTINE build_RF_local
!
END PROGRAM

