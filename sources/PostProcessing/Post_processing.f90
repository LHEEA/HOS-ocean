PROGRAM Post_processing
!
! This is the main program for post-processing HOS-ocean computations
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
USE input_post_process
USE read_files
USE Analysis_wavefield
USE output_post_process
USE reconstruction
USE fourier_r2c
!
IMPLICIT NONE
!
REAL(RP) :: dt_out, time, time_prev
REAL(RP), ALLOCATABLE, DIMENSION(:) :: H_up, L_up, H_down, L_down, crest, trough
INTEGER, ALLOCATABLE, DIMENSION(:)  :: idx_up, idx_down, idx_crest, idx_trough
! Transverse
REAL(RP), ALLOCATABLE, DIMENSION(:) :: L_up_t
INTEGER, ALLOCATABLE, DIMENSION(:)  :: idx_crest_t,idx_up_t
INTEGER  :: i_unit,n_waves
REAL(RP) :: H_1_3rd_up, H_1_3rd_down, H_lim
!
! For moments
REAL(RP) :: ave_eta,adev_eta,sdev_eta,var_eta,skew_eta,kurt_eta
REAL(RP), ALLOCATABLE, DIMENSION(:) :: data
!
! For freak waves
INTEGER :: n_freak
REAL(RP), ALLOCATABLE, DIMENSION(:) :: H_freak,L_freak,x_freak
INTEGER,ALLOCATABLE, DIMENSION(:)   :: idx_freak
! Transverse
REAL(RP), ALLOCATABLE, DIMENSION(:) :: L_freak_t,y_freak
INTEGER,ALLOCATABLE, DIMENSION(:)   :: idx_freak_t
!
! For SWENSE-type outputs + velocity/pressure cards
!
!INTEGER :: nz
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt
REAL(RP), ALLOCATABLE, DIMENSION(:)      :: xvect,yvect,zvect
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: vitx,vity,vitz,phit,dudt,dvdt,dwdt,zlocal
!
REAL(RP) :: dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T
INTEGER  :: imin,imax,jmin,jmax,i_xvect,i_yvect
INTEGER  :: i1,i2,i3,i_test
!
! Read input file to define parameters of post-processing
CALL read_input('input_post_process.dat')
!
! Initialize outputs
CALL init_output_post_process(i_ana,i_card)
!
! Wave-by-wave analysis
IF (i_ana /= 0) THEN
	!
	! Initialize data reading from file...
	i_unit = 101
	CALL init_read_3d(file_3d,i_unit,tecplot,n1,n2,x,y,eta,phis,dt_out)
	print*, file_3d
	print*, n1,n2,dt_out
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
		! Always compute the moments of free surface elevation
		IF (i_ana >= 1) THEN
			!
			ALLOCATE(data(n1*n2))
			DO i1=1,n1
				DO i2=1,n2
					data(i2+n2*(i1-1)) = eta(i1,i2)
				ENDDO
			ENDDO
			!
			! Evaluation of moments associated to 3D-free surface elevation
			CALL moment(n1*n2,data,ave_eta,adev_eta,sdev_eta,var_eta,skew_eta,kurt_eta)
			!
			DEALLOCATE(data)
			!
		ENDIF
		!
		IF ((n2 > 1).AND.(i_ana == 3)) THEN ! 3D analysis of 3D wave-field
			!
			! Wave-by-wave analysis for 3D wave-field
			CALL wave_by_wave_3D(eta,x,y,n1,n2,n_waves,H_up,L_up,idx_up,H_down,crest,idx_crest,L_up_t,idx_up_t,idx_crest_t)
			!
			! Evaluate H_1/3 on 3D waves
			CALL H_onethird(H_up,n_waves,H_1_3rd_up)
			CALL H_onethird(H_down,n_waves,H_1_3rd_down)
			!
			! Locate the freak waves with up-crossing analysis
			H_lim = HfoHs*4.0_rp*sdev_eta
			CALL locate_freak_3D(H_up,L_up,idx_crest,idx_up,L_up_t,idx_crest_t,idx_up_t,n_waves,x,y,n1,n2,H_lim, &
				n_freak,H_freak,L_freak,x_freak,idx_freak,L_freak_t,y_freak,idx_freak_t)
			!
			! Output at the given time-step
			CALL output_time_step_ana(i_ana,tecplot,time,eta,ave_eta,sdev_eta,skew_eta,kurt_eta,&
				MAXVAL(H_up),MAXVAL(H_down),MAXVAL(crest),H_1_3rd_up,H_1_3rd_down,&
				n_freak,H_freak,x_freak,L_freak,idx_freak,y_freak,L_freak_t,idx_freak_t)
			!
			! Deallocate
			DEALLOCATE(H_up,L_up,H_down,crest,idx_up,idx_crest,H_freak,L_freak,x_freak,idx_freak,y_freak,idx_freak_t,L_freak_t)
		ELSEIF (i_ana == 2) THEN ! 2D analysis of wave-field (if 3D treat only i2=1)
				!
				! Wave-by-wave analysis
				CALL wave_by_wave(eta(:,1),x,n1,n_waves,H_up,L_up,idx_up,H_down,L_down,idx_down,crest,idx_crest,trough,idx_trough)
				!
				! Evaluate H_1/3
				CALL H_onethird(H_up,n_waves,H_1_3rd_up)
				CALL H_onethird(H_down,n_waves,H_1_3rd_down)
				!
				! Locate the freak waves with up-crossing analysis
				H_lim = HfoHs*4.0_rp*sdev_eta
				CALL locate_freak(H_up,L_up,idx_crest,idx_up,n_waves,x,n1,H_lim,n_freak,H_freak,L_freak,x_freak,idx_freak)
				!
				! Output at the given time-step
				CALL output_time_step_ana(i_ana,tecplot,time,eta,ave_eta,sdev_eta,skew_eta,kurt_eta,&
					MAXVAL(H_up),MAXVAL(H_down),MAXVAL(crest),H_1_3rd_up,H_1_3rd_down,&
					n_freak,H_freak,x_freak,L_freak,idx_freak,0.0_rp*x_freak,0.0_rp*L_freak,0*idx_freak)
				!
				! Deallocate
				DEALLOCATE(H_up,L_up,H_down,L_down,crest,trough,idx_up,idx_down,idx_crest,H_freak,L_freak,x_freak,idx_freak,idx_trough)
		ELSEIF (i_ana == 1) THEN ! Output of results for case i_ana = 1
			n_freak = 1
			ALLOCATE(H_freak(n_freak),x_freak(n_freak),L_freak(n_freak), idx_freak(n_freak))
			ALLOCATE(y_freak(n_freak),L_freak_t(n_freak), idx_freak_t(n_freak))
			!
			! Output at the given time-step
			CALL output_time_step_ana(i_ana,tecplot,time,eta,ave_eta,sdev_eta,skew_eta,kurt_eta,&
				0.0_rp,0.0_rp,0.0_rp,0.0_rp,0.0_rp,&
				n_freak,H_freak,x_freak,L_freak,idx_freak,y_freak,L_freak_t,idx_freak_t)
			!
			DEALLOCATE(H_freak,L_freak,x_freak,idx_freak,y_freak,idx_freak_t,L_freak_t)
		ENDIF
		!
		! next time-step
		time_prev = time
		time      = time + dt_out
	ENDDO
	!
	! Close all files (including those open in output...
	CLOSE(i_unit)
	CLOSE(20)
	CLOSE(21)
	CLOSE(22)
	!
	! Deallocate variables that may be reallocated afterwards...
	DEALLOCATE(x,y,eta, phis)
ENDIF
!
! Velocities and pressure inside domain
IF (i_card /= 0) THEN
	!
	i_unit = 201
	!
	! Initialize computations of volumic informations
	! Everything is non-dimensional in file_mod
	!
	CALL recons_HOS_init(file_mod,i_unit,n1,n2,dt_out_star,T_stop_star,xlen_star,ylen_star,depth_star,g_star,L,T, &
		modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt)
	!
	n1o2p1 = n1/2+1
	n2o2p1 = n2/2+1
	!
	! Initialize Fourier (global variables have to be defined)
	m1      = n1
	m2      = n2
	Nd1     = 1
	Nd2     = 1
	Nd1o2p1 = 1
	m1o2p1  = n1o2p1
	md1o2p1 = 1
	md1     = 1
	md2     = 1
	!
	! Initialize Fourier transforms (FFTW)
	CALL fourier_ini(3)
	!
	! Check (x_min, x_max, y_min, y_max) w.r.t. domain size
	! + time window (t_min, t_max)
	CALL check_sizes(n2,x_min,x_max,y_min,y_max,T_start,T_stop,xlen_star,ylen_star,T_stop_star,L,T)
	!
	ALLOCATE(x(n1),y(n2),kx(n1o2p1),ky_n2(n2),ikx(n1o2p1,n2),iky(n1o2p1,n2),kth(n1o2p1,n2))
	!
	! Initialize mesh in physical and modal space (whole domain in HOS-ocean)
	CALL build_mesh_global(xlen_star,ylen_star,depth_star,n1,n2,x,y,kx,ky_n2,ikx,iky,kth)
	!
	! Define local meshes for zone of study
	CALL build_mesh_local(x_min,x_max,y_min,y_max,z_min,z_max,xlen_star,ylen_star,L,n1,n2,i_zvect, &
		xvect,yvect,zvect,imin,imax,jmin,jmax)
	!
	! Reconstruction of fields
	! First ALLOCATE the matrices
	i_xvect = imax-imin+1
	i_yvect = jmax-jmin+1
	!
	ALLOCATE(zlocal(i_xvect, i_yvect), vitx(i_xvect, i_yvect), vity(i_xvect, i_yvect), vitz(i_xvect, i_yvect), &
		phit(i_xvect, i_yvect), dudt(i_xvect, i_yvect), dvdt(i_xvect, i_yvect), dwdt(i_xvect, i_yvect))
	!
	! Define first time as the closest to T_start (input file)
	time      = NINT(T_start/T/dt_out_star)*dt_out_star 
	time_prev = 0.0_rp
	!
	DO WHILE (time*T <= T_stop)
		!
		write(*,'(A,ES8.1)') 'time = ',time*T
		!
		! It reads the corresponding time in file_mod (closest to time)
		IF (time >= dt_out_star/2) THEN
			CALL read_mod(file_mod,i_unit,time,dt_out_star,n1o2p1,n2, &
				modesspecx,modesspecy,modesspecz,modesspect,modesFS,modesFSt)
		ENDIF
		!
		! Make a loop over all the elements in z
		DO i3 = 1, i_zvect
			IF (i_card == 1) THEN
				! Construct the field at each zvect...
				CALL reconstruction_FFTs(modesspecx,modesspecy,modesspecz,modesspect,modesFS, &
					imin,imax,jmin,jmax,zvect(i3),depth_star,vitx,vity,vitz,phit,dudt,dvdt,dwdt)
				! Necessary for output
				DO i1=1, imax-imin+1
					DO i2=1, jmax-jmin+1
						zlocal(i1,i2) = zvect(i3)
					ENDDO
				ENDDO
			ELSEIF (i_card == 2) THEN
				! Construct the field at each zvect...
				CALL reconstruction_direct(modesspecx,modesspecy,modesspecz,modesspect,modesFS, &
					imin,imax,jmin,jmax,z_min/L,i3,i_zvect,depth_star,vitx,vity,vitz,phit,dudt,dvdt,dwdt,zlocal)
			ENDIF
			! Test to know if it is first z-element to write corresponding header
			IF (i3 == 1) THEN
				i_test = 1
			ELSE
				i_test = 0
			ENDIF
			! Output of time-step
			CALL output_time_step_card(i_card,tecplot,time,dt_out_star,zlocal,z_min,z_max,T_start,g_star,L,T,i_test,&
				imin,imax,jmin,jmax,i_zvect,vitx,vity,vitz,phit)
		ENDDO
		!
		! next time-step
		time_prev = time
		time      = time + dt_out_star
		!
	ENDDO
	!
	! Close all files (including those open in output...
	CLOSE(i_unit)
	CLOSE(30)
	CLOSE(31)
	CLOSE(32)
ENDIF
! End of main program 
!
CONTAINS
!
!
!
END PROGRAM Post_processing