MODULE input_post_process
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
INTERFACE read_datum
     MODULE PROCEDURE read_datum_i, read_datum_r, read_datum_c
END INTERFACE
!
INTEGER, PARAMETER            :: N_descr = 33
INTEGER, PARAMETER            :: N_tot   = 52
INTEGER, PARAMETER            :: len_form_read = 7
INTEGER, PARAMETER            :: len_form_write = 25
CHARACTER(LEN=len_form_read)  :: format_read(0:4)
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
CALL read_blank_line(unit)             ! --- Choice of post-processing options
CALL read_datum(unit, i_ana)           ! Wave-by-wave analysis
CALL read_datum(unit, i_card)          ! VP-card output
CALL read_datum(unit,T_start)          ! Starting time analysis
CALL read_datum(unit, T_stop)          ! Stoping time analysis
WRITE(*,*)
CALL read_blank_line(unit)             ! --- Velocities/pressure cards
CALL read_datum(unit, x_min)           ! Minimum x in VP-card
CALL read_datum(unit, x_max)           ! Maximum x in VP-card
CALL read_datum(unit, y_min)           ! Minimum y in VP-card
CALL read_datum(unit, y_max)           ! Maximum y in VP-card
CALL read_datum(unit, z_min)           ! Minimum z in VP-card
CALL read_datum(unit, z_max)           ! Maximum z in VP-card
CALL read_datum(unit, i_zvect)         ! Maximum z in VP-card
WRITE(*,*)
CALL read_blank_line(unit)             ! --- Input files
CALL read_datum(unit, tecplot)         ! Tecplot version
CALL read_datum(unit, file_3d)         ! Name of free-surface file
CALL read_datum(unit, file_mod)        ! Name of modal description
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
    STOP
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
STOP
!
END SUBROUTINE error_message
!
!
!
END MODULE input_post_process
