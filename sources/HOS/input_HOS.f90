MODULE input_HOS
!
! This module contains the input related routines
!  Subroutines :  read_input
!                 read_datum
!                 read_blank_line
!                 build_format
!                 error_message
!                 write_input
!                 write_datum
!                 write_blank_line
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
INTERFACE read_datum
    MODULE PROCEDURE read_datum_i, read_datum_r, read_datum_c
END INTERFACE
!
INTERFACE write_datum
    MODULE PROCEDURE write_datum_i, write_datum_r, write_datum_c
END INTERFACE
!
INTEGER, PARAMETER            :: N_descr = 33
INTEGER, PARAMETER            :: N_tot   = 52
INTEGER, PARAMETER            :: len_form_read = 7
INTEGER, PARAMETER            :: len_form_write = 25
CHARACTER(LEN=len_form_read)  :: format_read(0:4)
CHARACTER(LEN=len_form_write) :: format_write(0:2)
INTEGER                       :: line_counter
CHARACTER(LEN=N_tot)          :: description
!
!
!
CONTAINS
!
!
!
SUBROUTINE read_input(filename)
!
IMPLICIT NONE
! Input variables
CHARACTER(LEN=*), INTENT(IN)  :: filename
! Local variables
INTEGER  :: unit
!
unit = 100
OPEN(unit, FILE=filename)
line_counter = 0
!
CALL read_datum(unit, i_case)          ! Choice of computed case
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, xlen)            ! Length in x-direction
CALL read_datum(unit, ylen)            ! Length in y-direction
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, T_stop)          ! Duration of the simulation
CALL read_datum(unit, f_out)           ! Sampling frequency (output)
CALL read_datum(unit, toler)           ! Tolerance of RK scheme
CALL read_datum(unit, n)               ! Dommermuth initialisation
CALL read_datum(unit, Ta)              ! Dommermuth initialisation
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, grav)            ! Gravity
CALL read_datum(unit, depth)           ! Water depth
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit,Tp_real)          ! Peak period in s
CALL read_datum(unit,Hs_real)          ! Significant wave height in m
CALL read_datum(unit, gamma)           ! JONSWAP Spectrum
CALL read_datum(unit, beta)            ! Directionality (Dysthe)
CALL read_datum(unit, random_phases)   ! Random phases generation
CALL read_blank_line(unit)
WRITE(*,*)
CALL read_datum(unit, tecplot)         ! Tecplot version
CALL read_datum(unit, i_out_dim)       ! Output: 1-dim. 0-non dim.
CALL read_datum(unit, i_3d)            ! 3d free surface quantities
CALL read_datum(unit, i_a_3d)          ! 3d modes
CALL read_datum(unit, i_2d)            ! 2d free surface, center line
CALL read_datum(unit, i_prob)          ! Wave probes in domain
CALL read_datum(unit, i_sw)            ! Swense output 1='yes',0='no'
!
CLOSE(unit)
!
END SUBROUTINE read_input
!
!
!
SUBROUTINE read_datum_i(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)  :: unit
! Output variables
INTEGER, INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('I')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",I5)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_i
!
!
!
SUBROUTINE read_datum_r(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
REAL(RP), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('R')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",ES25.16)') description(1:N_descr), input
IF (ABS(input) > tiny .AND. ABS(input) * 1.0E+16_rp < 1.0E+5_rp) THEN
    WRITE(*,'(A,A)') 'Numeric point is probably missing in current input ',description
    STOP 1
ENDIF
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_r
!
!
!
SUBROUTINE read_datum_c(unit, input)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
! Output variables
CHARACTER(LEN=*), INTENT(OUT) :: input
!
line_counter = line_counter + 1
CALL build_read_format('A')
READ(unit,format_read,ERR=101,END=101) description, input
WRITE(*,'(A," : ",A)') description(1:N_descr), input
RETURN
101 CALL error_message(description)
!
END SUBROUTINE read_datum_c
!
!
!
SUBROUTINE read_blank_line(unit)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
!
line_counter = line_counter + 1
READ(unit,*)
!
END SUBROUTINE read_blank_line
!
!
!
SUBROUTINE build_read_format(code)
!
IMPLICIT NONE
!
CHARACTER(LEN=*) :: code
!
format_read(0) = '(A'
format_read(1) = int2str(N_tot)
format_read(2) = ','
SELECT CASE (code)
    CASE('I')          ! Integer
        format_read(3) = 'I5'
    CASE('F','R')      ! Real number
        format_read(3) = 'ES25.16'
    CASE('S','C','A')  ! Character string
        format_read(3) = 'A'
END SELECT
format_read(4) = ')'
! WRITE(*,*) format_read
!
END SUBROUTINE build_read_format
!
!
!
SUBROUTINE error_message(description)
!
IMPLICIT NONE
! Input variables
CHARACTER(LEN=N_tot) :: description
!
WRITE(*,'(A,I2)') 'Error while reading the input file on line: ', line_counter
WRITE(*,'(A)') description
STOP 1

!
END SUBROUTINE error_message
!
!
!
SUBROUTINE write_input(unit)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN) :: unit
! Tecplot version
CALL write_datum(unit, i_case,            'i_case',         'Choice of computed case')
CALL write_blank_line(unit,'--- Geometry of the horizontal domain')
CALL write_datum(unit, xlen,              'xlen',           'Length in x-direction')
CALL write_datum(unit, ylen,              'ylen',           'Length in y-direction')
!
CALL write_blank_line(unit,'--- Time stuff')
CALL write_datum(unit, T_stop,            'T_stop',         'Duration of the simulation')
CALL write_datum(unit, f_out,             'f_out',          'Sampling frequency (output)')
CALL write_datum(unit, toler,             'toler',          'Tolerance of RK scheme')
CALL write_datum(unit, n,                 'n',              'Dommermuth initialisation')
CALL write_datum(unit, Ta*T,                'Ta',             'Dommermuth initialisation')
!
CALL write_blank_line(unit,'--- Physical parameters')
CALL write_datum(unit, grav,              'grav',           'Gravity')
IF (ABS(depth-1.0e15_rp) <= epsilon(1.0e15_rp)) THEN
    CALL write_datum(unit, -1.0_rp,             'depth',          'Water depth')
ELSE
    CALL write_datum(unit, depth,             'depth',          'Water depth')
ENDIF
!
CALL write_blank_line(unit,'--- Irregular waves')
CALL write_datum(unit, Tp_real,         'Tp_real',      'Peak period in s')
CALL write_datum(unit, Hs_real,         'Hs_real',      'Significant wave height in m')
CALL write_datum(unit, gamma,           'gamma',        'JONSWAP Spectrum')
CALL write_datum(unit, beta,            'beta',         'Directionality (Dysthe)')
CALL write_datum(unit, random_phases,   'random_phases','Random phases generation')
!
CALL write_blank_line(unit,'--- Output files')
CALL write_datum(unit, tecplot,           'tecplot',        'Tecplot version')
CALL write_datum(unit, i_out_dim,         'i_out_dim',      'Output: 1-dim. 0-non dim.')
CALL write_datum(unit, i_3d,              'i_3d',           '3d free surface quantities')
CALL write_datum(unit, i_a_3d,            'i_a_3d',         '3d modes')
CALL write_datum(unit, i_2d,              'i_2d',           '2d free surface, center line')
CALL write_datum(unit, i_prob,            'i_prob',         'Wave probes in domain')
CALL write_datum(unit, i_sw,              'i_sw',           'Swense output 1=yes,0=no')
!
!CALL write_blank_line(unit,'---  Swense real size (i_sw=1)')
!
CALL write_blank_line(unit,'--- Numerical parameters')
CALL write_datum(unit, n1,                'n1',             'Modes number in x-direction')
CALL write_datum(unit, n2,                'n2',             'Modes number in y-direction')
CALL write_datum(unit, M,                 'M',              'HOS nonlinearity order')
CALL write_datum(unit, p1,                'p1',             'Dealiasing in x-direction')
CALL write_datum(unit, p2,                'p2',             'Dealiasing in y-direction')
!
END SUBROUTINE write_input
!
!
!
SUBROUTINE write_datum_i(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)  :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
INTEGER, INTENT(IN) :: input
!
CALL build_write_format('I')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_i
!
!
!
SUBROUTINE write_datum_r(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
REAL(RP), INTENT(IN) :: input
!
CALL build_write_format('R')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_r
!
!
!
SUBROUTINE write_datum_c(unit, input, variable, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)   :: unit
CHARACTER(LEN=*), INTENT(IN) :: text, variable
CHARACTER(LEN=*), INTENT(IN) :: input
!
CALL build_write_format('A')
WRITE(unit,format_write) text, variable, input
!
END SUBROUTINE write_datum_c
!
!
!
SUBROUTINE write_blank_line(unit, text)
!
IMPLICIT NONE
! Input variables
INTEGER, INTENT(IN)          :: unit
CHARACTER(LEN=*), INTENT(IN) :: text
!
WRITE(unit,'("#",A)') text
!
END SUBROUTINE write_blank_line
!
!
!
SUBROUTINE build_write_format(code)
!
IMPLICIT NONE
!
CHARACTER(LEN=*) :: code
!
format_write(0) = '("#",A'//TRIM(int2str(N_descr))//',":: ",A'//TRIM(int2str(N_tot-(N_descr+3)-3))//',":: ",'
SELECT CASE (code)
    CASE('I')          ! Integer
        format_write(1) = 'I5'
    CASE('F','R')      ! Real number
        format_write(1) = 'ES25.16'
    CASE('S','C','A')  ! Character string
        format_write(1) = 'A'
END SELECT
format_write(2) = ')'
!
END SUBROUTINE build_write_format
!
!
!
END MODULE input_HOS
