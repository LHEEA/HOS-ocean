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
!    This program is free software: you can redistribute it and/or modify
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
USE variables_post_process
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
IF (i_ana == 1) THEN
	!
	! Results of wave field analysis: moments of FS and wave-by-wave analysis
	OPEN(20,file='Results/Analysis.dat',status='unknown')
	WRITE(20,'(A)') 'TITLE=" Free surface analysis "'
	WRITE(20,'(A,A)') 'VARIABLES="time","ave_eta","sdev_eta","skew_eta","kurt_eta","Hs",', &
		'"H_max_up","H_max_down","crest_max","trough_max","H_one3_up","H_one3_down","n_freak"'
	!
	! Free surface profile of freak waves detected
	OPEN(21,file='Results/freak_waves.dat',status='unknown')
	WRITE(21,'(A)') 'TITLE=" 3D free surface elevation of freak waves "'
	WRITE(21,'(A)') 'VARIABLES="x","y","eta"'
	!
	! Characteristics of freak events detected
	OPEN(22,file='Results/Caract_freaks.dat',status='unknown')
	WRITE(22,'(A)') 'TITLE=" Parameters of detected freak waves "'
	WRITE(22,'(A)') 'VARIABLES="time","i_freak","H_freak","x_freak","Lx_freak"'
END IF
!
IF (i_card == 1) THEN
	!
	! Velocity and pressure card
ENDIF
!
END SUBROUTINE init_output_post_process
!
!
!
SUBROUTINE output_time_step(i_ana,i_card,tecplot,time,x,eta,ave_eta,sdev_eta,skew_eta,kurt_eta, &
		H_max_up,H_max_down,crest_max,trough_max,H_one3_up,H_one3_down,n_freak,&
		H_freak,x_freak,L_freak,idx_freak)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)    :: i_ana, i_card, tecplot
REAL(RP), INTENT(IN)   :: time
REAL(RP), DIMENSION(:) :: x,eta
!
REAL(RP) :: ave_eta,sdev_eta,skew_eta,kurt_eta,H_max_up,H_max_down,crest_max,trough_max,H_one3_up,H_one3_down
!
INTEGER  :: n_freak
REAL(RP), DIMENSION(n_freak) :: H_freak,x_freak,L_freak
INTEGER, DIMENSION(n_freak)  :: idx_freak
! Local variables
!
INTEGER  :: i_freak, i1, i2, npts, idx_tmp
REAL(RP) :: dx
!
IF (i_ana == 1) THEN
	!
	WRITE(20,'(12(ES16.9,X),I6)') time,ave_eta,sdev_eta,skew_eta,kurt_eta, &
		4.0_rp*sdev_eta,H_max_up,H_max_down,crest_max,trough_max,H_one3_up,H_one3_down,n_freak
	IF(n_freak.NE.0) THEN
		write(22,103)'ZONE T = "t = ',time, '", I=', n_freak
		!
		DO i_freak=1,n_freak
      		write(22,'(ES16.9,X,I6,X,3(ES16.9,X))') time,i_freak,H_freak(i_freak),x_freak(i_freak),L_freak(i_freak)
      	    !
    		dx   = x(2)-x(1)
    		npts = NINT(L_freak(i_freak)/dx)+1
    		IF (tecplot == 11) THEN
    			WRITE(21,103)'ZONE SOLUTIONTIME = ',time,', I=',npts,', J=',1
    		ELSE
    			WRITE(21,103)'ZONE T = "',time,'", I=',npts,', J=',1
    		END IF
    		DO i2=1,1
    			DO i1=1,npts
    				idx_tmp = idx_freak(i_freak)+i1-1
    				IF (idx_tmp > SIZE(x,1)) idx_tmp = idx_tmp - SIZE(x,1) ! Useful for eta
    				WRITE(21,'(3(ES16.9,X))') x(idx_freak(i_freak))+(i1-1)*dx, 0.d0, eta(idx_tmp)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
ENDIF
!
IF (i_card == 1) THEN
	!
	! Velocity and pressure card
ENDIF
!
103 FORMAT(A,F9.2,A,I5,A,I5)
                  
END SUBROUTINE output_time_step
!
END MODULE output_post_process