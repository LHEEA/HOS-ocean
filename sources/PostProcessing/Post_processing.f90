PROGRAM Post_processing
!
! This is the main program for post-processing HOS-ocean computations
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
USE input_post_process
USE read_files
USE Analysis_wavefield
USE output_post_process
!
IMPLICIT NONE
!
REAL(RP) :: dt_out, time, time_prev
REAL(RP), ALLOCATABLE, DIMENSION(:) :: H_up, L_up, H_down, L_down, crest, trough
INTEGER, ALLOCATABLE, DIMENSION(:)  :: idx_up, idx_down, idx_crest, idx_trough
INTEGER  :: i_unit,n_waves
REAL(RP) :: H_1_3rd_up, H_1_3rd_down, H_lim
!
! For moments
REAL(RP) :: ave_eta,adev_eta,sdev_eta,var_eta,skew_eta,kurt_eta
!
! For freak waves
INTEGER :: n_freak
REAL(RP), ALLOCATABLE, DIMENSION(:) :: H_freak,L_freak,x_freak
INTEGER,ALLOCATABLE, DIMENSION(:)   :: idx_freak
!
! Read input file to define parameters of post-processing
CALL read_input('input_post_process.dat')
!
! Wave-by-wave analysis
IF (i_ana == 1) THEN
	!
	! Initialize data reading from file...
	i_unit = 101
	CALL init_read_3d(file_3d,i_unit,tecplot,n1,n2,x,y,eta,phis,dt_out)
	print*, file_3d
	print*, n1,n2,dt_out
	!
	IF (n2 > 1) THEN
		print*, 'This is limited to 2D analysis for now'
		STOP
	ENDIF
	!
	! Initialize outputs
	CALL init_output_post_process(i_ana,i_card)
	!
	! Force the time to be a multiple of dt_out... starting from closest time-step 
	time      = NINT(T_start/dt_out)*dt_out
	time_prev = 0.0_rp
	!
	DO WHILE (time <= T_stop)
		!
		! It reads the closest time in file_3d
		IF (time >= dt_out/2) CALL read_3d(i_unit,tecplot,time_prev,time,dt_out,n1,n2,eta,phis)
		!
		! Wave-by-wave analysis
		CALL wave_by_wave(eta(:,1),x,n1,n_waves,H_up,L_up,idx_up,H_down,L_down,idx_down,crest,idx_crest,trough,idx_trough)
		!
		! Evaluate H_1/3
		CALL H_onethird(H_up,n_waves,H_1_3rd_up)
		CALL H_onethird(H_down,n_waves,H_1_3rd_down)
		!
		! Evaluation of moments associated to free surface elevation
		CALL moment(n1,eta(:,1),ave_eta,adev_eta,sdev_eta,var_eta,skew_eta,kurt_eta)
		!
		! Locate the freak waves with up-crossing analysis
		H_lim = HfoHs*4.0_rp*sdev_eta
		CALL locate_freak(H_up,L_up,idx_crest,idx_up,n_waves,x,n1,H_lim,n_freak,H_freak,L_freak,x_freak,idx_freak)
		!
		! Output at the given time-step
		CALL output_time_step(i_ana,i_card,tecplot,time,x,eta(:,1),ave_eta,sdev_eta,skew_eta,kurt_eta, &
			MAXVAL(H_up),MAXVAL(H_down),MAXVAL(crest),MINVAL(trough),H_1_3rd_up,H_1_3rd_down,n_freak,H_freak,x_freak,L_freak,idx_freak)
		!
		! next time-step
		time_prev = time
		time      = time + dt_out
		!
		! Deallocate
		DEALLOCATE(H_up,L_up,H_down,L_down,crest,trough,idx_up,idx_down,idx_crest,H_freak,L_freak,x_freak,idx_freak,idx_trough)
	ENDDO
	CLOSE(i_unit)
ENDIF
! End of main program 
!
CONTAINS
!
!
!
END PROGRAM Post_processing