MODULE output
!
! This module contains the initialization of output files and their writing at specified time-step
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
USE fourier_r2c
USE input_HOS
USE velocities
!
IMPLICIT NONE
!
!
!
CONTAINS
!
!
!
SUBROUTINE init_output(i_3d, i_a, i_vol, i_2D, i_max, i_prob, i_sw)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN) :: i_3d, i_a, i_vol, i_2D, i_max, i_prob, i_sw
! Local variables
INTEGER :: i1, i2
!
IF (i_3D == 1) THEN
    OPEN(1,file='Results/3d.dat',status='unknown')
    CALL write_input(1)
    WRITE(1,'(A)')'TITLE=" 3D free surface elevation "'
    WRITE(1,'(A)') 'VARIABLES="x","y","eta","phis"'
ENDIF
!
IF (i_a == 1) THEN
    OPEN(2,file='Results/a_3d.dat',status='unknown')
    CALL write_input(2)
    WRITE(2,'(A)')'TITLE=" 3D modes "'
    WRITE(2,'(A)') 'VARIABLES="kx","ky","a-eta","a-phis","LOG10-a-eta","LOG10-a-phis"'
ENDIF
!
IF (i_vol == 1) THEN
    OPEN(3,file='Results/vol_energy.dat',status='unknown')
    CALL write_input(3)
    WRITE(3,'(A)')'TITLE=" 3D volume and energy "'
    WRITE(3,'(A)') 'VARIABLES="t", "volume", "potential", "kinetic", "total", "dE/E_o","E_spec_dens"'
ENDIF
!
IF (i_2D == 1) THEN
    OPEN(66,file='Results/2d.dat',status='unknown')
    CALL write_input(66)
    WRITE(66,'(A)')'TITLE=" 2D free surface elevation on the center line "'
    WRITE(66,'(A)') 'VARIABLES="x","eta"'
ELSE IF (i_2D == 2) THEN ! 2D plus t thing
    OPEN(66,file='Results/2dpt.dat',status='unknown')
    CALL write_input(66)
    WRITE(66,'(A)')'TITLE=" 2D+t free surface elevation "'
    WRITE(66,'(A)') 'VARIABLES="x","t","eta"'
ENDIF
!
IF (i_max == 1) THEN
    OPEN(8,file='Results/eta_max.dat',status='unknown')
    CALL write_input(8)
    WRITE(8,'(A)')'TITLE=" max/minimum free surface elevation as a function of time"'
    WRITE(8,'(A)') 'VARIABLES="t","max-eta","min-eta","max/sqrt(t)","min/sqrt(t)"'
ENDIF
!
IF (i_case/10 == 8) THEN
    OPEN(9,file='Results/phase_shift.dat',status='unknown')
    CALL write_input(9)
    WRITE(9,'(A)')'TITLE=" phase shift (Fructus test) as a function of time"'       !
    WRITE(9,'(A)') 'VARIABLES="t","phase shift (in degrees)"'
ENDIF
!
IF (i_prob == 1) THEN
    !
    ! Wave probes location
    !
    WRITE(*,'(A)') 'Probes:'
    OPEN(55,FILE='prob.inp')
    nprobes=0
    DO
        READ(55,*,END=98)
        nprobes = nprobes + 1
        CYCLE
98    EXIT
    ENDDO
    REWIND(55)
    !
    WRITE(*,'(I3,2A)') nprobes, ' probes found in the file prob.inp'
    !
    IF (nprobes.gt.maxprobes) THEN
        PRINT*,'error: increase maxprobes in variables file'
        STOP 1
    ENDIF
    DO i1=1,nprobes
        IF (n2 == 1) THEN
                READ(55,*) xprobe(i1)
                yprobe(i1) = 0.0d0
        ELSE
                READ(55,*) xprobe(i1), yprobe(i1)
        ENDIF
    ENDDO
    CLOSE(55)
    !
    WRITE(*,'(A)') 'Probes position:'
    !
    DO i1=1,nprobes
        WRITE(*,902)'xprobe(',i1,')=',xprobe(i1),' m, yprobe(',i1,')=',yprobe(i1),' m'
    ENDDO
    !
    xprobe = xprobe / L
    yprobe = yprobe / L
    !
    !
902 FORMAT(a,2(i2,a,1es10.3,a))
    !
    OPEN(99,file='Results/probes.dat',status='unknown')
    CALL write_input(99)
    WRITE(99,'(A)') 'TITLE="Probes records versus time"'        !
    WRITE(99,'(101A)') 'VARIABLES = "time" ', ('"p'//TRIM(int2str(i1))//'" ',i1=1,nprobes)
ENDIF
!
IF (i_sw == 1) THEN
    IF ((2*n1o2p1).GE.5000) THEN
        PRINT*,'Problem in the writing of modes_HOS_swense.dat: change writing format'
        STOP 1
    ENDIF
    ! Give constants needed to computation...
    OPEN(123,file='Results/modes_HOS_SWENSE.dat',status='REPLACE', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*(2*n1o2p1))
    WRITE(123,'(5000(ES17.10,1X))',REC=1) REAL(n1,RP), REAL(n2,RP), 1.0_rp/f_out_star,T_stop_star &
          , xlen_star , ylen_star , depth_star, g_star, L, T, (0.0_rp, i1=11,2*n1o2p1)

    WRITE(123,'(5000(ES17.10,1X))',REC=2) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=3) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=4) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=5) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=6) (0.0_rp, i1=1,2*n1o2p1)
    DO i2=2,n2
        WRITE(123,'(5000(ES17.10,1X))',REC=1+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=2+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=3+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=4+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=5+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=6+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
    ENDDO
ENDIF
!
END SUBROUTINE init_output
!
!
SUBROUTINE output_time_step(i_3d, i_a, i_vol, i_2D, i_max, time, N_stop, a_eta, a_phis, da_eta, volume, energy, E_0, E_tot, &
    i_prob, dt, n_er_tot, n_rk_tot)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)               :: i_3d, i_a, i_vol, i_2D, i_max, N_stop, i_prob, n_er_tot, n_rk_tot
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_eta, a_phis, da_eta
REAL(RP), DIMENSION(m1,m2)        :: eta, phis
REAL(RP) :: time, volume, energy(3), E_0(3), E_tot, dt
! Local variables
INTEGER                        :: ii, i1, i2, it
REAL(RP)                       :: a_1, eta_mult, min_val, max_val
REAL(RP), DIMENSION(m1,m2)     :: envelope
REAL(RP), DIMENSION(m1o2p1,m2) :: abs_eta, log_eta,abs_phis,log_phis
!
eta_mult = eta_out
IF (i_3D == 1) THEN
    ! Transform modal amplitudes to physical variables
    CALL fourier_2_space(a_eta,  eta)
    CALL fourier_2_space(a_phis, phis)
    !
    IF (n2 == 1) THEN
        envelope = hilbert(a_eta)
        envelope(1:n1,1:n2) = SQRT(eta(1:n1,1:n2)**2 + envelope(1:n1,1:n2)**2)
        WHERE (ABS(envelope(1:n1,1:n2)) < 1.0E-50_rp)
                envelope(1:n1,1:n2) = 0.0_rp
        END WHERE
    ELSE
        envelope(1:n1,1:n2) = 0.0_rp
    ENDIF
    !
    ! Output of the free surface elevation versus space (size control of the output file)
    IF (ABS(time-time_restart) <= tiny) THEN
        IF (tecplot == 11) THEN
            WRITE(1,103)'ZONE SOLUTIONTIME = ',time*T_out,', I=',n1,', J=',n2
        ELSE
            WRITE(1,103)'ZONE T = "',time*T_out,'", I=',n1,', J=',n2
        ENDIF
        DO i2 = 1, n2
            DO i1 = 1, n1
                WRITE(1,102) x(i1)*L_out, y(i2)*L_out, &
                eta(i1,i2)*eta_mult*L_out, phis(i1,i2)*eta_mult*L_out**2/T_out
            ENDDO
        ENDDO
    ELSE
        IF (tecplot == 11) THEN
            WRITE(1,103)'ZONE SOLUTIONTIME = ',time*T_out,', D=(1,2), I=',n1,', J=',n2
        ELSE
            WRITE(1,103)'ZONE T = "',time*T_out,'", D=(1,2), I=',n1,', J=',n2
        ENDIF
        DO i2 = 1, n2
            DO i1 = 1, n1
            WRITE(1,104) eta(i1,i2)*eta_mult*L_out, phis(i1,i2)*eta_mult*L_out**2/T_out
            ENDDO
        ENDDO
    ENDIF
   102 FORMAT(3(ES12.5,X),ES12.5)
   103 FORMAT(A,ES12.5,A,I5,A,I5)
   104 FORMAT((ES12.5,X),ES12.5)
ENDIF
!
IF (i_2D == 1) THEN
    ! Transform modal amplitude to physical variable
    CALL fourier_2_space(a_eta,  eta)
    !
    IF (ABS(time-time_restart) <= tiny) THEN
        IF (tecplot == 11) THEN
            WRITE(66,603)'ZONE SOLUTIONTIME = ',time*T_out,', I=',n1
        ELSE
            WRITE(66,603)'ZONE T = "',time*T_out,'", I=',n1
        ENDIF
        DO i1 = 1, n1
            WRITE(66,602) x(i1)*L_out, eta(i1,MAX(1,n2/2+1))*eta_mult*L_out
        ENDDO
    ELSE
        IF (tecplot == 11) THEN
            WRITE(66,603)'ZONE SOLUTIONTIME = ',time*T_out,', D=(1), I=',n1
        ELSE
            WRITE(66,603)'ZONE T = "',time*T_out,'", D=(1), I=',n1
        ENDIF
        DO i1 = 1, n1
            WRITE(66,604) eta(i1,MAX(1,n2/2+1))*eta_mult*L_out
        ENDDO
    ENDIF
   602 FORMAT(ES12.5,X,ES12.5)
   603 FORMAT(A,ES12.5,A,I5)
   604 FORMAT(ES12.5)
ELSE IF (i_2D == 2) THEN
    ! Transform modal amplitude to physical variable
    CALL fourier_2_space(a_eta,  eta)
    !
    ! Output of the free surface elevation versus space (size control of the output file)
    IF (ABS(time-time_restart) <= tiny) THEN
        WRITE(66,'(A,I4,A,I6)')'ZONE I=',n1,', J=',N_stop
        DO i1 = 1, n1
            WRITE(66,662) x(i1)*L_out, time*T_out, eta(i1,MAX(1,n2/2+1))*eta_mult*L_out
        ENDDO
    ELSE
        DO i1 = 1, n1
            WRITE(66,662) x(i1)*L_out, time*T_out, eta(i1,MAX(1,n2/2+1))*eta_mult*L_out
        ENDDO
    ENDIF
   662 FORMAT(2(ES12.5,X),ES12.5)
ENDIF
!
IF (i_a == 1) THEN
    ! Output of the modal amplitudes
    abs_eta(1:n1o2p1,1:n2) = ABS(a_eta(1:n1o2p1,1:n2))*ABS(eta_mult)*L_out
    !
    log_eta(1:n1o2p1,1:n2) = LOG10(MAX(EPSILON(1.0_rp),abs_eta(1:n1o2p1,1:n2)))
    !
    abs_phis(1:n1o2p1,1:n2) = ABS(a_phis(1:n1o2p1,1:n2))*ABS(eta_mult)*L_out**2/T_out
    !
    log_phis(1:n1o2p1,1:n2) = LOG10(MAX(EPSILON(1.0_rp),abs_phis(1:n1o2p1,1:n2)))
    !
    IF (ABS(time-time_restart) <= tiny) THEN
        IF (tecplot == 11) THEN
            WRITE(2,103)'ZONE SOLUTIONTIME = ',time*T_out,', I=',n1o2p1,', J=',n2
        ELSE
            WRITE(2,103)'ZONE T = "',time*T_out,'", I=',n1o2p1,', J=',n2
        ENDIF
        ! Negatives ky first
        DO i2 = n2o2p1+1, n2
            DO i1 = 1, n1o2p1
                WRITE(2,202) kx(i1)/L_out, ky_n2(i2)/L_out, abs_eta(i1,i2), abs_phis(i1,i2), log_eta(i1,i2), log_phis(i1,i2)
            ENDDO
        ENDDO
        DO i2 = 1, n2o2p1
            DO i1 = 1, n1o2p1
                WRITE(2,202) kx(i1)/L_out, ky_n2(i2)/L_out, abs_eta(i1,i2), abs_phis(i1,i2), log_eta(i1,i2), log_phis(i1,i2)
            ENDDO
        ENDDO
    ELSE
        IF (tecplot == 11) THEN
            WRITE(2,103)'ZONE SOLUTIONTIME = ',time*T_out,', D=(1,2), I=',n1o2p1,', J=',n2
        ELSE
            WRITE(2,103)'ZONE T = "',time*T_out,'", D=(1,2), I=',n1o2p1,', J=',n2
        ENDIF
        DO i2 =  n2o2p1+1, n2
            DO i1 = 1, n1o2p1
                WRITE(2,204) abs_eta(i1,i2), abs_phis(i1,i2), log_eta(i1,i2), log_phis(i1,i2)
            ENDDO
        ENDDO
        DO i2 = 1, n2o2p1
            DO i1 = 1, n1o2p1
                WRITE(2,204) abs_eta(i1,i2), abs_phis(i1,i2), log_eta(i1,i2), log_phis(i1,i2)
            ENDDO
        ENDDO
   202 FORMAT(5(ES15.8,X),ES15.8)
   204 FORMAT(3(ES15.8,X),ES15.8)
    ENDIF
ENDIF
!
!************************
IF (i_vol == 1) THEN
    volume = volume * xlen_star
    IF (n2 /= 1) volume = volume * ylen_star
    IF (ABS(E_0(3)) > tiny) THEN
        WRITE(3,303) time*T_out, volume*L_out, (energy(i1)*L_out**3/T_out**2,i1=1,3), &
            ABS(energy(3)-E_0(3))/E_0(3),E_tot*L_out**3/T_out**2
    ELSE
        WRITE(3,303) time*T_out, volume*L_out, (energy(i1)*L_out**3/T_out**2,i1=1,3), &
            0.0_rp,E_tot*L_out**3/T_out**2
    ENDIF
   303 FORMAT(6(ES12.5,X),ES12.5)
ENDIF
!
IF (i_max == 1) THEN
    ! Analytical integration of the linear part
    CALL fourier_2_space(a_eta, eta)
    IF (ABS(time) > tiny) THEN
        a_1 = 1.0_rp / SQRT(time)
    ELSE
        a_1 = 1.0_rp
    ENDIF
    min_val = MINVAL(eta(1:n1,1:n2))
    max_val = MAXVAL(eta(1:n1,1:n2))
    WRITE(8,802) time*T_out, max_val*eta_mult*L_out, min_val*eta_mult*L_out, &
                max_val*eta_mult*L_out*a_1, min_val*eta_mult*L_out*a_1
   802 FORMAT(4(ES12.5,X),ES12.5)
ENDIF
!
! FIXME: make this work for xlen diff from 1
IF (i_case/10 == 8) THEN
    WRITE(9,'(ES12.5,X,ES12.5)') time * T_out, ATAN2(IMAG(a_eta(2,1)), REAL(a_eta(2,1),RP))*180.0_rp/PI
ENDIF
!
IF (i_prob == 1) THEN ! probes output
    ! Compute the probe elevation
    !
    ! Methode directe
    !
    DO ii=1,nprobes
        !
        i1 = 1
        i2 = 1
        eta_probe(ii)   =  REAL(a_eta(i1,i2),RP)
        !
        DO i2=2,n2o2p1
            eta_probe(ii)   =  eta_probe(ii)  + 2.0_rp*ABS(a_eta(i1,i2)) &
                *COS(ky_n2(i2)*yprobe(ii)+ATAN2(AIMAG(a_eta(i1,i2)),REAL(a_eta(i1,i2),RP)))
        ENDDO
        !
        DO i1=2,n1o2p1
            DO i2=1,n2
                eta_probe(ii)   =  eta_probe(ii)  + 1.0_rp*ABS(a_eta(i1,i2) * EXP(i*ky_n2(i2)*yprobe(ii))) &
                    *COS(kx(i1)*xprobe(ii)+ATAN2(AIMAG(a_eta(i1,i2) * EXP(i*ky_n2(i2)*yprobe(ii))) &
                    ,REAL(a_eta(i1,i2) * EXP(i*ky_n2(i2)*yprobe(ii)),RP)))
            ENDDO
        ENDDO
    ENDDO
    ! For the output of probes, use Hs_real and Tp_real...
    WRITE(99,'(101(ES13.5,X))') time * T_out, (eta_probe(ii)* L_out, ii=1,nprobes)
    IF (nprobes.GT.100) THEN
        PRINT*, 'Change writing format in probes.dat file'
        STOP 1
    ENDIF
ENDIF
!
IF (i_sw == 1) THEN ! time=0 has to be saved...
    it = NINT(time*T*f_out)+1
    DO i2=1,n2
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+1+6*(i2-1)) (modesspecx(i1,i2), i1=1,n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+2+6*(i2-1)) (modesspecy(i1,i2), i1=1,n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+3+6*(i2-1)) (modesspecz(i1,i2), i1=1,n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+4+6*(i2-1)) (modesspect(i1,i2), i1=1,n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+5+6*(i2-1)) (a_eta(i1,i2)   , i1=1,n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=((it)*n2*6)+6+6*(i2-1)) (da_eta(i1,i2)  , i1=1,n1o2p1)
    ENDDO
ENDIF
!
! Backup the current file (if problem occurs during writing)
CALL SYSTEM('mv -f Results/restart_FS_quantities.dat Results/restart_FS_quantities_backup.dat')
!
! Write the restart data in a file
OPEN(234,file='Results/restart_FS_quantities.dat')
! Global informations
WRITE(234,403) time, dt, n_er_tot, n_rk_tot
DO i2 = 1, n2
    DO i1 = 1, n1o2p1
        WRITE(234,404) a_eta(i1,i2), a_phis(i1,i2)
    ENDDO
ENDDO
!
CLOSE(234)
403 FORMAT(2(D25.16,X),I6,X,I6)
404 FORMAT(3(D25.16,X),D25.16)
!
END SUBROUTINE output_time_step
!
!
!
FUNCTION hilbert(x)
!
IMPLICIT NONE
! Input/output variables
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: x
REAL(RP), DIMENSION(m1,m2)      :: hilbert
!
temp_C_n(1:n1o2p1,1:n2) = x(1:n1o2p1,1:n2) * EXP(i * PIO2)
!
CALL fourier_2_space(temp_C_n, temp_R_n)
!
hilbert = temp_R_n
!
END FUNCTION hilbert
!
!
!
SUBROUTINE init_restart(i_3d, i_a, i_vol, i_2D, i_max, i_prob, i_sw, time_cur, dt, n_er_tot, n_rk_tot)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)        :: i_3d, i_a, i_vol, i_2D, i_max, i_prob, i_sw
!
! Output variables
INTEGER, INTENT(OUT)       :: n_er_tot, n_rk_tot
REAL(RP),INTENT(OUT)       :: time_cur, dt
! Local variables
INTEGER :: i1, i2
!
! Read the restart data stored in a file
OPEN(234,file='Results/restart_FS_quantities.dat')
! Global informations
READ(234,403) time_restart, dt, n_er_tot, n_rk_tot
DO i2 = 1, n2
    DO i1 = 1, n1o2p1
        READ(234,404) a_eta(i1,i2), a_phis(i1,i2)
    ENDDO
ENDDO
!
time_cur = time_restart
!
IF (i_3D == 1) THEN
    OPEN(1,file='Results/3d_restart.dat',status='unknown')
    CALL write_input(1)
    WRITE(1,'(A)')'TITLE=" 3D free surface elevation "'
    WRITE(1,'(A)') 'VARIABLES="x","y","eta","phis"'
ENDIF
!
IF (i_a == 1) THEN
    OPEN(2,file='Results/a_3d_restart.dat',status='unknown')
    CALL write_input(2)
    WRITE(2,'(A)')'TITLE=" 3D modes "'
    WRITE(2,'(A)') 'VARIABLES="kx","ky","a-eta","a-phis","LOG10-a-eta","LOG10-a-phis"'
ENDIF
!
IF (i_vol == 1) THEN
    OPEN(3,file='Results/vol_energy_restart.dat',status='unknown')
    CALL write_input(3)
    WRITE(3,'(A)')'TITLE=" 3D volume and energy "'
    WRITE(3,'(A)') 'VARIABLES="t", "volume", "potential", "kinetic", "total", "dE/E_o","E_spec_dens"'
ENDIF
!
IF (i_2D == 1) THEN
    OPEN(66,file='Results/2d_restart.dat',status='unknown')
    CALL write_input(66)
    WRITE(66,'(A)')'TITLE=" 2D free surface elevation on the center line "'
    WRITE(66,'(A)') 'VARIABLES="x","eta"'
ELSE IF (i_2D == 2) THEN ! 2D plus t thing
    OPEN(66,file='Results/2dpt_restart.dat',status='unknown')
    CALL write_input(66)
    WRITE(66,'(A)')'TITLE=" 2D+t free surface elevation "'
    WRITE(66,'(A)') 'VARIABLES="x","t","eta"'
ENDIF
!
IF (i_max == 1) THEN
    OPEN(8,file='Results/eta_max_restart.dat',status='unknown')
    CALL write_input(8)
    WRITE(8,'(A)')'TITLE=" max/minimum free surface elevation as a function of time"'
    WRITE(8,'(A)') 'VARIABLES="t","max-eta","min-eta","max/sqrt(t)","min/sqrt(t)"'
ENDIF
!
IF (i_case/10 == 8) THEN
    OPEN(9,file='Results/phase_shift_restart.dat',status='unknown')
    CALL write_input(9)
    WRITE(9,'(A)')'TITLE=" phase shift (Fructus test) as a function of time"'       !
    WRITE(9,'(A)') 'VARIABLES="t","phase shift (in degrees)"'
ENDIF
!
IF (i_prob == 1) THEN
    !
    ! Wave probes location
    !
    WRITE(*,'(A)') 'Probes:'
    OPEN(55,FILE='prob.inp')
    !
    nprobes=0
    DO
        READ(55,*,END=98)
        nprobes = nprobes + 1
        CYCLE
98      EXIT
    ENDDO
    REWIND(55)
   !
   WRITE(*,'(I3,2A)') nprobes, ' probes found in the file prob.inp'
   !
   IF (nprobes.gt.maxprobes) STOP 'error: increase maxprobes in variables file'
   DO i1=1,nprobes
      IF (n2 == 1) THEN
         READ(55,*) xprobe(i1)
         yprobe(i1) = 0.0d0
      ELSE
         READ(55,*) xprobe(i1), yprobe(i1)
      END IF
   END DO
   CLOSE(55)
   !
   WRITE(*,'(A)') 'Probes position:'
   !
   DO i1=1,nprobes
      WRITE(*,902)'xprobe(',i1,')=',xprobe(i1),' m, yprobe(',i1,')=',yprobe(i1),' m'
   END DO
   !
902 FORMAT(a,2(i2,a,1es10.3,a))
903 FORMAT("#",A,100(F10.2))
	!
    OPEN(99,file='Results/probes_restart.dat',status='unknown')
    CALL write_input(99)
    CALL write_datum(99, nprobes,             'nprobes',          'number of probes')
    WRITE(99,903)'xprobe=',(xprobe(i1),i1=1,nprobes)
    WRITE(99,903)'yprobe=',(yprobe(i1),i1=1,nprobes)
    WRITE(99,'(A)') 'TITLE="Probes records versus time"'
    WRITE(99,'(61A)') 'VARIABLES = "time" ', ('"p'//TRIM(int2str(i1))//'" ',i1=1,nprobes)
    !
    xprobe = xprobe / L
    yprobe = yprobe / L
ENDIF
!
IF (i_sw == 1) THEN
    IF ((2*n1o2p1).GE.5000) THEN
        PRINT*,'Problem in the writing of modes_HOS_swense.dat: change writing format'
        STOP
    ENDIF
    ! Give constants needed to computation...
    OPEN(123,file='Results/modes_HOS_SWENSE_restart.dat',status='REPLACE', FORM='FORMATTED', ACCESS='DIRECT',RECL=18*(2*n1o2p1))
    WRITE(123,'(5000(ES17.10,1X))',REC=1) REAL(n1,RP), REAL(n2,RP), 1.0_rp/f_out_star,T_stop_star &
          , xlen_star , ylen_star , depth_star, g_star, L, T, (0.0_rp, i1=11,2*n1o2p1)

    WRITE(123,'(5000(ES17.10,1X))',REC=2) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=3) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=4) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=5) (0.0_rp, i1=1,2*n1o2p1)
    WRITE(123,'(5000(ES17.10,1X))',REC=6) (0.0_rp, i1=1,2*n1o2p1)
    DO i2=2,n2
        WRITE(123,'(5000(ES17.10,1X))',REC=1+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=2+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=3+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=4+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=5+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
        WRITE(123,'(5000(ES17.10,1X))',REC=6+6*(i2-1)) (0.0_rp, i1=1,2*n1o2p1)
    ENDDO
ENDIF
!
CLOSE(234)
403 FORMAT(2(D25.16,X),I6,X,I6)
404 FORMAT(3(D25.16,X),D25.16)
!
END SUBROUTINE init_restart
!
!
!
END MODULE output
