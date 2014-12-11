MODULE ramp
!
! defines a ramp (either in space or time or anything else)
!
!    ^
!  1 |-------------xxxxxxxxxxxxxx
!    |            x              x
!    |           x                x
!    |          x                  x
! 0.5|---------x--------------------x
!    |        x|                    |x
!    |       x |                    | x
!    |      x  |                    |  x
!  --x-----x---x--------------------x---x-------x---->
!      t_begin t_start         t_stop         t_end
! 
! t_dur = t_stop - t_start
! d_start and d_stop = lengths of the start and stop ramps
! r_start and r_stop = type of the start and stop ramps
!   0 = no ramp
!   1 = linear ramp
!   2 = polynomial ramp order 2
!   3 = polynomial ramp order 3
!   4 = polynomial ramp order 4
!   5 = polynomial ramp order 5
!
! Subroutines and functions
! * ramp_ini
! * ramp_value
! * ramp_calc
! * dramp_du_calc
! * d2ramp_du2_calc
! * draw_ramp
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
!
IMPLICIT NONE
!
TYPE :: ramp_def
	REAL(RP) :: t_begin, t_end, t_start, t_stop, t_dur ! locations
	REAL(RP) :: d_start, d_stop                        ! ramp lengths
	INTEGER  :: r_start, r_stop                        ! ramp type
END TYPE ramp_def
!
!
CONTAINS
!
!
!
FUNCTION ramp_ini(t_def, ri_def, rt_def)
!
! Initializes the ramp object
! t_def is either scalar (=t_dur), length 2 (=t_begin,t_dur) or length 3
! (t_begin, t_dur and t_end)
! ri_def is scalar (same start and stop type ramp) or length 2 (different types)
!  see the ramp definitions in ramp_calc
! rt_def is scalar (same start and stop ramp lengths) or length 2 (different lengths)
!
IMPLICIT NONE
! Input variables
REAL(RP), DIMENSION(:) :: t_def
INTEGER, DIMENSION(:)  :: ri_def
REAL(RP), DIMENSION(:) :: rt_def
! Output variables
TYPE(ramp_def) :: ramp_ini
! Local variables
INTEGER :: s
!
! Start ramp
!   type
ramp_ini%r_start = ri_def(1)
!   duration
IF (ramp_ini%r_start == 0) THEN ! no ramp_ini
   ramp_ini%d_start = 0.0_rp
ELSE ! ramp_ini
   ramp_ini%d_start = rt_def(1)
END IF
!
! Stop ramp_ini
!   different type or not ?
IF (SIZE(ri_def) > 1) THEN 
   ramp_ini%r_stop = ri_def(2)    ! type
ELSE ! same def as start ramp_ini
   ramp_ini%r_stop = ri_def(1)    ! type
END IF
!   duration
IF (ramp_ini%r_stop == 0) THEN  ! no ramp_ini
   ramp_ini%d_stop = 0.0_rp
ELSE ! ramp_ini
   IF (SIZE(rt_def) > 1) THEN
      ramp_ini%d_stop = rt_def(2) ! duration
   ELSE ! same duration as start ramp_ini
      ramp_ini%d_stop = rt_def(1) ! duration
   END IF
END IF
!
! Definition of the different ramp_ini times
s = SIZE(t_def)
!
IF (s > 3) STOP 'ramp_ini_ini: Too much times at input'
IF (s == 1) THEN ! begins at t=0
   ramp_ini%t_begin = 0.0_rp
ELSE ! begins at t=t_begin/=0
   ramp_ini%t_begin = t_def(1)
END IF
! Starting time is define when the ramp_ini is 1/2 
!   I assume it is when we are half way into the ramp_ini
ramp_ini%t_start = ramp_ini%t_begin + ramp_ini%d_start / 2
!
IF (s == 1) THEN
   ramp_ini%t_dur = t_def(1)
ELSE
   ramp_ini%t_dur = t_def(2)
END IF
!
ramp_ini%t_stop = ramp_ini%t_start + ramp_ini%t_dur
! ramp_ini%t_stop = ramp_ini%t_start + ramp_ini%t_dur + ramp_ini%d_stop / 2
!
IF (s == 3) THEN
   ramp_ini%t_end = t_def(3)
ELSE
   ramp_ini%t_end = ramp_ini%t_stop + ramp_ini%d_stop / 2
END IF
!
END FUNCTION ramp_ini
!
!
!
FUNCTION ramp_value(ramp, loc)
!
! gives the values of the ramp, first and second derivatives
!
IMPLICIT NONE
!
! Input variables
TYPE(ramp_def) :: ramp
REAL(RP)       :: loc
! Output variables
REAL(RP), DIMENSION(3) :: ramp_value
! Local variables
REAL(RP) :: u
!
IF (loc < ramp%t_begin) THEN ! before the start ramp
   ramp_value(1) = 0.0_rp
   ramp_value(2) = 0.0_rp
   ramp_value(3) = 0.0_rp
ELSE IF (ABS(ramp%d_start) .GT. tiny .AND. loc-(ramp%t_begin + ramp%d_start) <= tiny) THEN ! start ramp
   u = (loc - ramp%t_begin) / ramp%d_start
   ramp_value(1) = ramp_calc(u, ramp%r_start)
   ramp_value(2) = dramp_du_calc(u, ramp%r_start)   / ramp%d_start
   ramp_value(3) = d2ramp_du2_calc(u, ramp%r_start) / (ramp%d_start * ramp%d_start)
ELSE IF (loc-(ramp%t_stop - ramp%d_stop/2) < tiny) THEN ! 
   ramp_value = ramp_calc(1.0_rp, ramp%r_start)
   ramp_value(1) = ramp_value(1)
   ramp_value(2) = 0.0_rp
   ramp_value(3) = 0.0_rp
ELSE IF (ABS(ramp%d_stop) .GT. tiny .AND. loc-(ramp%t_stop + ramp%d_stop/2) < tiny) THEN ! start ramp
   u = - (loc - (ramp%t_stop + ramp%d_stop/2)) / ramp%d_stop
   ramp_value(1) =   ramp_calc(u, ramp%r_stop)
   ramp_value(2) = - dramp_du_calc(u, ramp%r_stop)   / ramp%d_stop
   ramp_value(3) =   d2ramp_du2_calc(u, ramp%r_stop) / (ramp%d_stop * ramp%d_stop)
ELSE
   ramp_value(1) = 0.0_rp
   ramp_value(2) = 0.0_rp
   ramp_value(3) = 0.0_rp
END IF
!
END FUNCTION ramp_value
!
!
!
FUNCTION ramp_calc(u, r_type)
!
! evaluates the ramp value
!
IMPLICIT NONE
!
! Input variables
INTEGER  :: r_type
REAL(RP) :: u
! Output variables
REAL(RP) :: ramp_calc
! Local variables
REAL(RP) :: v
!
SELECT CASE(r_type)
   CASE (0) ! no ramp
      ramp_calc = 1.0_rp
   CASE (1) ! linear ramp
      ramp_calc = u
   CASE (2) ! quadratic ramp
      IF (u < 0.5_rp) THEN
         ramp_calc = 2.0_rp * u * u
      ELSE
         v         = u - 1.0_rp
         ramp_calc = 1.0_rp - 2.0_rp * v * v
      END IF
   CASE (3) ! triic ramp
      ramp_calc = (3.0_rp - 2.0_rp * u) * u * u
   CASE (4) ! quartic ramp
      ramp_calc = (6.0_rp + u * (-8.0_rp + 3.0_rp * u)) * u * u
   CASE (5) ! quintic ramp
      ramp_calc = (10.0_rp + 3.0_rp * u * (-5.0_rp + 2.0_rp * u)) * u * u * u
   CASE (14) ! quartic ramp
      ramp_calc = (1.0_rp - (2.0_rp * u - 1.0_rp)**2)**2
   CASE (97) ! sinus ramp
      ramp_calc = SIN(PI * u)**2
   CASE (98) ! sinus ramp
      ramp_calc = SIN(PIO2 * u)
   CASE (99) ! squared sinus ramp
      ramp_calc = SIN(PIO2 * u)**2
END SELECT
!
END FUNCTION ramp_calc
!
!
!
FUNCTION dramp_du_calc(u, r_type)
!
! evaluates the value of the first derivative of the ramp
!
IMPLICIT NONE
!
! Input variables
INTEGER  :: r_type
REAL(RP) :: u
! Output variables
REAL(RP) :: dramp_du_calc
! Local variables
REAL(RP) :: v
!
SELECT CASE(r_type)
   CASE (0) ! no ramp
      dramp_du_calc = 0.0_rp
   CASE (1) ! linear ramp
      dramp_du_calc = 1.0_rp
   CASE (2) ! quadratic ramp
      IF (u < 0.5_rp) THEN
         dramp_du_calc = 4.0_rp * u
      ELSE
         v         = u - 1.0_rp
         dramp_du_calc = - 4.0_rp * v
      END IF
   CASE (3) ! triic ramp
      dramp_du_calc = (1.0_rp - u) * u * 6.0_rp
   CASE (4) ! quartic ramp
      dramp_du_calc = (1.0_rp + u * (-2.0_rp + u)) * u * 12.0_rp
   CASE (5) ! quintic ramp
      dramp_du_calc = (1.0_rp + u * (-2.0_rp + u)) * u * u * 30.0_rp
   CASE (14) ! quartic ramp
      dramp_du_calc = - 4.0_rp * (1.0_rp - (2.0_rp * u - 1.0_rp)**2) * (2.0_rp * u - 1.0_rp)
   CASE (97) ! sinus ramp
      dramp_du_calc = PI * COS(PI * u)
   CASE (98) ! sinus ramp
      dramp_du_calc = PIO2 * COS(PIO2 * u)
   CASE (99) ! squared sinus ramp
      dramp_du_calc = PIO2 * SIN(PI * u)
END SELECT
!
END FUNCTION dramp_du_calc
!
!
!
FUNCTION d2ramp_du2_calc(u, r_type)
!
! evaluates the value of the second derivative of the ramp
!
IMPLICIT NONE
!
! Input variables
INTEGER  :: r_type
REAL(RP) :: u
! Output variables
REAL(RP) :: d2ramp_du2_calc
!
SELECT CASE(r_type)
   CASE (0) ! no ramp
      d2ramp_du2_calc = 0.0_rp
   CASE (1) ! linear ramp
      d2ramp_du2_calc = 0.0_rp
   CASE (2) ! quadratic ramp
      IF (u < 0.5_rp) THEN
         d2ramp_du2_calc = 4.0_rp
      ELSE
         d2ramp_du2_calc = - 4.0_rp
      END IF
   CASE (3) ! triic ramp
      d2ramp_du2_calc = (1.0_rp - 2.0_rp * u) * 6.0_rp
   CASE (4) ! quartic ramp
      d2ramp_du2_calc = (1.0_rp + u * (-4.0_rp + 3.0_rp * u)) * 12.0_rp
   CASE (5) ! quintic ramp
      d2ramp_du2_calc = (1.0_rp + u * (-3.0_rp + 2.0_rp * u)) * u * 30.0_rp
   CASE (97) ! sinus ramp
      d2ramp_du2_calc = - PI * PI * SIN(PI * u)
   CASE (98) ! sinus ramp
      d2ramp_du2_calc = - PIO2 * PIO2 * SIN(PIO2 * u)
   CASE (99) ! squared sinus ramp
      d2ramp_du2_calc = PI * PIO2 * COS(PI * u)
END SELECT
!
END FUNCTION d2ramp_du2_calc
!
!
!
FUNCTION ramp_integ(ramp, loc)
!
! gives the values of the ramp integral
!  ramp_integ = int_0^t ramp(u) du
IMPLICIT NONE
!
! Input variables
TYPE(ramp_def) :: ramp
REAL(RP)       :: loc
! Output variables
REAL(RP) :: ramp_integ
! Local variables
REAL(RP) :: u, tmp
!
IF (loc-ramp%t_begin < tiny) THEN ! before the start ramp
   ramp_integ = 0.0_rp
ELSE IF (abs(ramp%d_start) > tiny .AND. loc-(ramp%t_begin + ramp%d_start) < tiny) THEN ! start ramp
   u = (loc - ramp%t_begin) / ramp%d_start
   ramp_integ = int_ramp_calc(u,ramp%r_start)
   ramp_integ = ramp%d_start * ramp_integ
ELSE IF (loc < ramp%t_stop - ramp%d_stop/2) THEN ! 
   ramp_integ = int_ramp_calc(1.0_rp ,ramp%r_start)
   ramp_integ = ramp%d_start * ramp_integ
   ramp_integ = ramp_integ + ramp_calc(1.0_rp ,ramp%r_start) * (loc - (ramp%t_begin + ramp%d_start))
ELSE IF (abs(ramp%d_start) > tiny .AND. loc-(ramp%t_stop + ramp%d_stop/2) < tiny) THEN ! start ramp
   ramp_integ = int_ramp_calc(1.0_rp ,ramp%r_start)
   ramp_integ = ramp%d_start * ramp_integ
   ramp_integ = ramp_integ + ramp_calc(1.0_rp ,ramp%r_start) * ((ramp%t_stop - ramp%d_stop/2) - (ramp%t_begin + ramp%d_start))
   u   = - (loc - (ramp%t_stop + ramp%d_stop/2)) / ramp%d_stop
   tmp = int_ramp_calc(1.0_rp, ramp%r_stop) - int_ramp_calc(u, ramp%r_stop)
   ramp_integ = ramp_integ + ramp%d_stop * tmp
ELSE
   ramp_integ = int_ramp_calc(1.0_rp ,ramp%r_start)
   ramp_integ = ramp%d_start * ramp_integ
   ramp_integ = ramp_integ + ramp_calc(1.0_rp ,ramp%r_start) * ((ramp%t_stop - ramp%d_stop/2) - (ramp%t_begin + ramp%d_start))
   tmp = int_ramp_calc(1.0_rp, ramp%r_stop) - int_ramp_calc(0.0_rp, ramp%r_stop)
   ramp_integ = ramp_integ + ramp%d_stop * tmp
END IF
!
END FUNCTION ramp_integ
!
!
!
FUNCTION int_ramp_calc(u, r_type)
!
! evaluates the ramp value
!
IMPLICIT NONE
!
! Input variables
INTEGER  :: r_type
REAL(RP) :: u
! Output variables
REAL(RP) :: int_ramp_calc
! Local variables
REAL(RP) :: v
!
IF (abs(u-1.0_rp) > tiny) THEN
   SELECT CASE(r_type)
      CASE (0) ! no ramp
         int_ramp_calc = u
      CASE (1) ! linear ramp
         int_ramp_calc = u * u
      CASE (2) ! quadratic ramp
         IF (u < 0.5_rp) THEN
            int_ramp_calc = 2.0_rp * u * u * u / 3.0_rp
         ELSE
            v             = u - 1.0_rp
            int_ramp_calc = u - 0.5_rp - 2.0_rp * v * v * v / 3.0_rp
         END IF
      CASE (3) ! triic ramp
         int_ramp_calc = (1.0_rp - u * 0.5_rp) * u * u * u
      CASE (4) ! quartic ramp
         int_ramp_calc = (2.0_rp + u * (-2.0_rp + 0.6_rp * u)) * u * u * u
      CASE (5) ! quintic ramp
         int_ramp_calc = (2.5_rp + u * (-3.0_rp + u)) * u * u * u * u
      CASE (14) ! quartic ramp
         int_ramp_calc = 8.0_rp * (2.0_rp /3.0_rp + u * (- 1.0_rp + 0.4_rp * u)) * u * u * u
      CASE (97) ! sinus ramp
         int_ramp_calc = (TWOPI * u - SIN(TWOPI * u)) / (4.0_rp * PI)
      CASE (98) ! sinus ramp
         int_ramp_calc = ( 1.0_rp - COS(PIO2 * u)) / PIO2
      CASE (99) ! squared sinus ramp
         int_ramp_calc = (PI * u - SIN(PI * u)) / (4.0_rp * PIO2)
   END SELECT
ELSE ! u=1
   SELECT CASE(r_type)
      CASE (0) ! no ramp
         int_ramp_calc = 1.0_rp
      CASE (1) ! linear ramp
         int_ramp_calc = 1.0_rp
      CASE (2) ! quadratic ramp
         int_ramp_calc = 0.5_rp
      CASE (3) ! triic ramp
         int_ramp_calc = 0.5_rp
      CASE (4) ! quartic ramp
         int_ramp_calc = 0.6_rp
      CASE (5) ! quintic ramp
         int_ramp_calc = 0.5_rp
      CASE (14) ! quartic ramp
         int_ramp_calc = 8.0_rp  / 15.0_rp
      CASE (97) ! sinus ramp
         int_ramp_calc = 0.5_rp
      CASE (98) ! sinus ramp
         int_ramp_calc = 1.0_rp / PIO2
      CASE (99) ! squared sinus ramp
         int_ramp_calc = 0.5_rp
   END SELECT
END IF
!
END FUNCTION int_ramp_calc
!
!
!
SUBROUTINE draw_ramp(ramp, file_name)
!
! draw the ramp between locations 0 and ramp%t_end
!
IMPLICIT NONE
!
! Input variables
TYPE(ramp_def)   :: ramp
CHARACTER(LEN=*) :: file_name
! Local variables
INTEGER, PARAMETER  :: N = 100
REAL(RP), PARAMETER :: eps = 100.0_rp * EPSILON(1.0_rp)
REAL(RP) :: loc, dloc
INTEGER  :: iloop, unit
!
unit = 101
OPEN(unit,FILE=file_name)
!
WRITE(unit, '(A)') 'VARIABLES="t", "ramp", "ramp_t", "ramp_t_t", "int ramp dt"'
! at loc=0
loc  = 0.0_rp
CALL plot_ramp(loc)
! just before the start ramp
IF (ramp%t_begin > eps) THEN
   loc  = ramp%t_begin - eps
   CALL plot_ramp(loc)
END IF
! beginning of the start ramp
loc  = ramp%t_begin
CALL plot_ramp(loc)
! start ramp
dloc = ramp%d_start / REAL(N,RP)
DO iloop = 1,N
   loc  = loc + dloc
   CALL plot_ramp(loc)
END DO
! after the start ramp
loc  = ramp%t_begin + ramp%d_start + eps
CALL plot_ramp(loc)
!
! just before the stop ramp
loc  = ramp%t_stop - ramp%d_stop/2 - eps
CALL plot_ramp(loc)
! beginning of the stop ramp
loc  = ramp%t_stop - ramp%d_stop/2
CALL plot_ramp(loc)
! stop ramp
dloc = ramp%d_stop / REAL(N,RP)
DO iloop = 1,N
   loc  = loc + dloc
   CALL plot_ramp(loc)
END DO
! just after the stop ramp
loc  = ramp%t_stop + ramp%d_stop/2 + eps
CALL plot_ramp(loc)
! end of the ramp
loc  = ramp%t_end
CALL plot_ramp(loc)
! after the end of the ramp
loc  = 1.1_rp * ramp%t_end
CALL plot_ramp(loc)
!
CONTAINS
!
SUBROUTINE plot_ramp(loc)
!
IMPLICIT NONE
! Input variables
REAL(RP) :: loc
! Local variables
REAL(RP) :: value(3)
!
value = ramp_value(ramp, loc)
WRITE(unit, '(4(ES10.3,X),ES10.3)') loc, value(1), value(2), value(3), ramp_integ(ramp, loc)
!
END SUBROUTINE plot_ramp
!
END SUBROUTINE draw_ramp
!
!
!
END MODULE ramp