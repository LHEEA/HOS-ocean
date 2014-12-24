MODULE output_post_process
!
! This module contains the input related routines
!  Subroutines :  read_input
!                 read_datum
!                 read_blank_line
!                 build_format
!                 error_message
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
!
IMPLICIT NONE
!
!
!
CONTAINS
!
!
!
SUBROUTINE init_output_post_process(i_ana,i_card)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN) :: i_ana, i_card
! Local variables
!
!INTEGER :: i1, i2
!
IF (i_ana >= 1) THEN
	!
	! Results of wave field analysis: moments of FS and wave-by-wave analysis
	OPEN(20,file='Results/Analysis.dat',status='unknown')
	WRITE(20,'(A)') 'TITLE=" Free surface analysis "'
	WRITE(20,'(A,A)') 'VARIABLES="time","ave_eta","sdev_eta","skew_eta","kurt_eta","Hs",', &
		'"H_max_up","H_max_down","crest_max","H_one3_up","H_one3_down","n_freak"'
	IF (i_ana > 1) THEN
		!
		! Free surface profile of freak waves detected
		OPEN(21,file='Results/freak_waves.dat',status='unknown')
		WRITE(21,'(A)') 'TITLE=" 3D free surface elevation of freak waves "'
		WRITE(21,'(A)') 'VARIABLES="x","y","eta"'
		!
		! Characteristics of freak events detected
		OPEN(22,file='Results/Caract_freaks.dat',status='unknown')
		WRITE(22,'(A)') 'TITLE=" Parameters of detected freak waves "'
		WRITE(22,'(A)') 'VARIABLES="time","i_freak","H_freak","x_freak","Lx_freak","y_freak","Ly_freak"'
	ENDIF
END IF
!
IF (i_card /= 0) THEN
	!
	! Velocity and pressure card
	IF (i_card == 1) THEN
    	! Tecplot output
    	OPEN(31,FILE='Results/VP_card.dat')
		WRITE(31,'(A)') 'TITLE =" Velocity and pressure field "'
    	WRITE(31,'(A)') 'VARIABLES="x","y","z","vitx","vity","vitz","Press"'
    ELSEIF (i_card == 2) THEN ! Boundary fitted coordinates
    	! Tecplot output
    	OPEN(32,FILE='Results/VP_card_fitted.dat')
    	WRITE(32,'(A)') 'TITLE =" Velocity and pressure field "'
    	WRITE(32,'(A)') 'VARIABLES="x","y","z","vitx","vity","vitz","Press"'
    ENDIF
    !
ENDIF
!
END SUBROUTINE init_output_post_process
!
!
!
SUBROUTINE output_time_step_ana(i_ana,tecplot,time,eta,ave_eta,sdev_eta,skew_eta,kurt_eta, &
		H_max_up,H_max_down,crest_max,H_one3_up,H_one3_down,&
		n_freak,H_freak,x_freak,L_freak,idx_freak,y_freak,L_freak_t,idx_freak_t)
!
! This subroutine performs the output of wave-field analysis
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)                  :: i_ana, tecplot
REAL(RP), INTENT(IN)                 :: time
REAL(RP), DIMENSION(:,:), INTENT(IN) :: eta
!
REAL(RP), INTENT(IN) :: ave_eta,sdev_eta,skew_eta,kurt_eta,H_max_up,H_max_down,crest_max,H_one3_up,H_one3_down
!
INTEGER, INTENT(IN)  :: n_freak
REAL(RP), DIMENSION(MAX(n_freak,1)), INTENT(IN) :: H_freak,x_freak,L_freak,y_freak,L_freak_t
INTEGER, DIMENSION(MAX(n_freak,1)), INTENT(IN)  :: idx_freak, idx_freak_t
! Local variables
!
INTEGER  :: i_freak, i1, i2, nptsx, nptsy, idx_tmp, idx_tmp_t
REAL(RP) :: dx, dy
!
IF (i_ana >= 1) THEN
	!
	WRITE(20,'(11(ES16.9,X),I6)') time,ave_eta,sdev_eta,skew_eta,kurt_eta, &
		4.0_rp*sdev_eta,H_max_up,H_max_down,crest_max,H_one3_up,H_one3_down,n_freak
	IF ((i_ana>1).AND.(n_freak.NE.0)) THEN
		write(22,103)'ZONE T = "t = ',time, '", I=', n_freak
		!
		DO i_freak=1,n_freak
			write(22,'(ES16.9,X,I6,X,5(ES16.9,X))') time,i_freak,H_freak(i_freak),x_freak(i_freak),L_freak(i_freak), &
				y_freak(i_freak),L_freak_t(i_freak) ! 2D-case taken into-account in input ot output routine
      	    !
    		dx    = x(2)-x(1)
    		nptsx = NINT(L_freak(i_freak)/dx)+1
    		IF (SIZE(y,1) /= 1) THEN
				dy    = y(2)-y(1)
				nptsy = NINT(L_freak_t(i_freak)/dy)+1
			ELSE
				dy    = 0.0_rp
				nptsy = 1
			ENDIF
    		IF (tecplot == 11) THEN
    			WRITE(21,103)'ZONE SOLUTIONTIME = ',time,', I=',nptsx,', J=',nptsy
    		ELSE
    			WRITE(21,103)'ZONE T = "',time,'", I=',nptsx,', J=',nptsy
    		END IF
    		IF (SIZE(y,1) /= 1) THEN
				DO i2=1,nptsy
					DO i1=1,nptsx
						idx_tmp = idx_freak(i_freak)+i1-1
						IF (idx_tmp > SIZE(x,1)) idx_tmp = idx_tmp - SIZE(x,1) ! Useful for eta
						idx_tmp_t = idx_freak_t(i_freak)+i2-1
						IF (idx_tmp_t > SIZE(y,1)) idx_tmp_t = idx_tmp_t - SIZE(y,1) ! Useful for eta
						!
						WRITE(21,'(3(ES16.9,X))') x(idx_freak(i_freak))+(i1-1)*dx, y(idx_freak_t(i_freak))+(i2-1)*dy, eta(idx_tmp,idx_tmp_t)
					ENDDO
				ENDDO
            ELSE
            	DO i2=1,nptsy
					DO i1=1,nptsx
						idx_tmp = idx_freak(i_freak)+i1-1
						IF (idx_tmp > SIZE(x,1)) idx_tmp = idx_tmp - SIZE(x,1) ! Useful for eta
						!
						WRITE(21,'(3(ES16.9,X))') x(idx_freak(i_freak))+(i1-1)*dx, 0.0_rp, eta(idx_tmp,1)
					ENDDO
				ENDDO
            ENDIF
        ENDDO
    ENDIF
ENDIF
!
103 FORMAT(A,F9.2,A,I5,A,I5)

END SUBROUTINE output_time_step_ana
!
!
!
SUBROUTINE output_time_step_card(i_card,tecplot,time,dt_out_star,zlocal,z_min,z_max,T_start,g_star,L_out,T_out,i_test,&
		imin,imax,jmin,jmax,i_zvect,vitx,vity,vitz,phit)
!
! This subroutine performs the output of velocity/pressure card
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)                  :: i_card, tecplot
REAL(RP), INTENT(IN)                 :: time, dt_out_star
!
REAL(RP), INTENT(IN)                 :: T_start,L_out,T_out,g_star,z_min,z_max
INTEGER, INTENT(IN)                  :: i_test,imin,imax,jmin,jmax,i_zvect
REAL(RP), DIMENSION(:,:), INTENT(IN) :: zlocal,vitx,vity,vitz,phit
! Local variables
!
INTEGER  :: i1,i2
REAL(RP) :: Press
!
REAL(RP)             :: tiny_sp
!
! tiny_sp is single precision: useful for inequalities check with values read from files
tiny_sp = epsilon(1.0)
!
IF (i_card /= 0) THEN
	!
	! Velocity and pressure card
	IF (i_card == 1) THEN
		!
   		! These are informations useful for eventual coupling using files VP_card
    	IF (time*T_out <= T_start+tiny_sp) THEN ! First time-step
    		OPEN(30,file='Results/data_VP_card.dat',status='unknown')
			WRITE(30,'(2(ES16.9,X))') x(imin)*L_out, x(imax)*L_out
			WRITE(30,'(2(ES16.9,X))') y(jmin)*L_out, y(jmax)*L_out
			WRITE(30,'(2(ES16.9,X))') z_min, z_max ! should be used only for (i_card == 1)
			WRITE(30,'(3(I4,X))') imax-imin+1, jmax-jmin+1, i_zvect
			WRITE(30,'(3(ES16.9,X),I6)') T_start, T_stop, dt_out_star*T_out, FLOOR((T_stop-T_start)/T_out/dt_out_star) !T_stop, T_start are dimensional quantities
			CLOSE(30)
		ENDIF
		!
		IF (i_test == 1) THEN ! First element in the z-loop
			IF (time*T_out <= T_start+tiny_sp) THEN ! First time-step
				IF (tecplot == 11) THEN
					WRITE(31,104)'ZONE SOLUTIONTIME = ',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				ELSE
					WRITE(31,104)'ZONE T = "',time*T_out,'", I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				END IF
			ELSE ! Following time-steps
				IF (tecplot == 11) THEN
					WRITE(31,104)'ZONE SOLUTIONTIME = ',time*T_out,', D=(1,2,3), I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				ELSE
					WRITE(31,104)'ZONE T = "',time*T_out,'", D=(1,2,3), I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
				END IF
			ENDIF
		ENDIF
		!
		DO i2=1,jmax-jmin+1
			DO i1=1,imax-imin+1
				Press = - g_star*zlocal(i1,i2) - 0.5_rp*(vitx(i1,i2)**2+vity(i1,i2)**2+vitz(i1,i2)**2)-phit(i1,i2)
				if(time*T_out <= T_start+tiny_sp) then
					WRITE(31,'(7(ES16.9,X))') x(i1+imin-1)*L_out, y(i2+jmin-1)*L_out, zlocal(i1,i2)*L_out, &
							vitx(i1,i2)*L_out/T_out, vity(i1,i2)*L_out/T_out, vitz(i1,i2)*L_out/T_out, &
							Press*L_out**2/T_out**2
				else
					WRITE(31,'(4(ES16.9,X))') vitx(i1,i2)*L_out/T_out, vity(i1,i2)*L_out/T_out, &
							vitz(i1,i2)*L_out/T_out, Press*L_out**2/T_out**2
				endif
			ENDDO
		ENDDO
    ELSEIF (i_card == 2) THEN
    	IF (i_test == 1) THEN ! First element in the z-loop
			IF (tecplot == 11) THEN
				WRITE(32,104)'ZONE SOLUTIONTIME = ',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
			ELSE
				WRITE(32,104)'ZONE T = "',time*T_out,', I=', imax-imin+1,' J=', jmax-jmin+1,' K=', i_zvect
			END IF
    	ENDIF
		!
		DO i2=1,jmax-jmin+1
        	DO i1=1,imax-imin+1
				Press = - g_star*zlocal(i1,i2) - 0.5_rp*(vitx(i1,i2)**2+vity(i1,i2)**2+vitz(i1,i2)**2)-phit(i1,i2)
				WRITE(32,'(7(ES16.9,X))') x(i1+imin-1)*L_out, y(i2+jmin-1)*L_out, zlocal(i1,i2)*L_out, vitx(i1,i2)*L_out/T_out, &
					vity(i1,i2)*L_out/T_out, vitz(i1,i2)*L_out/T_out, Press*L_out**2/T_out**2
			ENDDO
		ENDDO
    ENDIF
ENDIF
!
104 FORMAT(A,F9.2,A,I5,A,I5,A,I5)

END SUBROUTINE output_time_step_card
!
END MODULE output_post_process