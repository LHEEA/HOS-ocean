MODULE linear_wave
!
! This module contains functions allowing the use of linear dispersion relation
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
USE maths, ONLY : dichoto
!
IMPLICIT NONE
!
INTERFACE alpha_adim
   MODULE PROCEDURE alpha_adim_r, alpha_adim_v
END INTERFACE
!
CONTAINS
!
FUNCTION alpha_r(f, N, prec, h, g)
IMPLICIT NONE
REAL(RP)       :: f, prec, h, g
INTEGER        :: N
REAL(RP)       :: alpha_r
REAL(RP)       :: asympt
!
! evaluation of alpha_N
!
asympt = (N-1.0_rp)*PI
IF (xtanx(asympt-PI/4.0_rp) <= - (TWOPI*f)**2 * h / g) THEN
   alpha_r = dichoto( xtanx, - (TWOPI*f)**2 * h / g, asympt-PI/4.0_rp, asympt, prec )
   alpha_r = alpha_r / h
ELSE IF (xtanx(asympt-PI/4.0_rp) > - (TWOPI*f)**2 * h / g) THEN
   alpha_r = dichoto( xtanx, - (TWOPI*f)**2 * h / g, asympt-PI/2.0_rp+prec, asympt-PI/4.0_rp, prec )
   alpha_r = alpha_r / h
END IF
!
END FUNCTION alpha_r
!
!
!
FUNCTION alpha_adim_r(f, N, prec)
!
IMPLICIT NONE
!
REAL(RP)     :: f, prec
INTEGER      :: N
REAL(RP)     :: alpha_adim_r
REAL(RP)     :: alpha_min, alpha_max, tmp, f0, Npi
!
! evaluation of alpha_N
!
Npi       = N * PI
alpha_max = Npi
tmp       = prec
alpha_min = Npi - PI/2.0_rp + tmp
f0        = - (TWOPI*f)**2
DO WHILE(xtanx(alpha_min) > f0)
   tmp       = tmp / 2.0_rp
      WRITE(*,*) N,tmp, f0, xtanx(alpha_min)
      READ(*,*)
   alpha_min = alpha_min - tmp
END DO
alpha_adim_r = dichoto( xtanx, f0, alpha_min, alpha_max, prec )
!
END FUNCTION alpha_adim_r
!
!
!
FUNCTION alpha_adim_v(f, N, prec)
!
IMPLICIT NONE
!
REAL(RP)                     :: f, prec
INTEGER, DIMENSION(:)        :: N
REAL(RP), DIMENSION(SIZE(N)) :: alpha_adim_v
INTEGER                      :: ibou
REAL(RP)     				 :: alpha_min, alpha_max, tmp, f0, Npi
!
! evaluation of alpha_N
!
f0 = - (TWOPI*f)**2
DO ibou = 1, SIZE(N)
   !
   Npi       = N(ibou) * PI
   alpha_max = Npi
   tmp       = prec
   alpha_min = Npi - PI/2.0_rp + tmp
   DO WHILE(xtanx(alpha_min) > f0)
      tmp       = tmp / 2.0_rp
      WRITE(*,*) N(ibou),tmp, f0, xtanx(alpha_min)
      alpha_min = alpha_min - tmp
   END DO
   alpha_adim_v(ibou) = dichoto( xtanx, f0, alpha_min, alpha_max, prec )
   !
END DO
!
END FUNCTION alpha_adim_v
!
!
!
FUNCTION wave_number_r(f,h,g,prec)
IMPLICIT NONE
REAL(RP)       :: f,h,g,prec
REAL(RP)       :: wave_number_r
!
! evaluation of k
!
IF ( (TWOPI*f)**2 / g * h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   wave_number_r = dichoto(xthx,(TWOPI*f)**2 * h / g,0.0_rp,10000000.0_rp,prec)
   wave_number_r = wave_number_r / h
ELSE
   wave_number_r = (TWOPI*f)**2 / g
END IF
!
END FUNCTION wave_number_r
!
!
!
FUNCTION wave_number_adim_r(f,prec)
IMPLICIT NONE
REAL(RP)       :: f,prec
REAL(RP)       :: wave_number_adim_r
!
! evaluation of k
!
IF ( (TWOPI*f)**2 < log(HUGE(1.0_rp))-1.0_rp ) THEN
   wave_number_adim_r = dichoto(xthx,(TWOPI*f)**2,0.0_rp,1.0e6_rp,prec)
ELSE
   WRITE(*,*) 'big frequency'
   wave_number_adim_r = (TWOPI*f)**2
END IF
!
END FUNCTION wave_number_adim_r
!
!
!
FUNCTION wave_number_v(f,h,g,prec)
!
IMPLICIT NONE
REAL(RP), DIMENSION(:)       :: f
REAL(RP)                     :: h,g,prec
REAL(RP), DIMENSION(SIZE(f)) :: wave_number_v
INTEGER                		 :: ibou
!
! evaluation of k
!
DO ibou = 1, SIZE(f)
   wave_number_v(ibou) = wave_number_r(f(ibou), h, g, prec)
END DO
!
END FUNCTION wave_number_v
!
!
!
FUNCTION xthx(x)
!
IMPLICIT NONE
!
REAL(RP) :: x,xthx
!
xthx = x*tanh(x)
!
END FUNCTION xthx
!
!
!
FUNCTION xtanx(x)
!
IMPLICIT NONE
!
REAL(RP) :: x,xtanx
!
xtanx = x*tan(x)
!
END FUNCTION xtanx
!
!
!
FUNCTION phase_velocity(f,h,g)
!
IMPLICIT NONE
REAL(RP)   :: f, h, g
REAL(RP)   :: phase_velocity
REAL(RP)   :: k, prec
prec          = 1E-7_rp
k             = wave_number_r(f, h, g, prec)
phase_velocity = TWOPI * f / k
!
END FUNCTION phase_velocity
!
!
!
FUNCTION group_velocity(f,h,g)
!
IMPLICIT NONE
REAL(RP)   :: f, h, g
REAL(RP)   :: group_velocity
REAL(RP)   :: k, prec
prec = EPSILON(f)
k    = wave_number_r(f, h, g, prec)
IF ( 2.0_rp*k*h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   group_velocity = phase_velocity(f,h,g) / 2.0_rp * (1.0_rp + (2.0_rp*k*h) / SINH(2.0_rp*k*h))
ELSE
   group_velocity = phase_velocity(f,h,g) / 2.0_rp
END IF
!
END FUNCTION group_velocity
!
!
!
FUNCTION omega_seconde(f,h,g)
!
IMPLICIT NONE
REAL(RP)   :: f, h, g
REAL(RP)   :: omega_seconde
REAL(RP)   :: k, prec, v_g
prec = EPSILON(f)
k    = wave_number_r(f, h, g, prec)
v_g  = group_velocity(f,h,g)
IF ( 2.0_rp*k*h < log(HUGE(1.0_rp))-1.0_rp ) THEN
   omega_seconde = v_g/k*(1.0_rp-2.0_rp*k*h*TANH(k*h))-v_g**2/(TWOPI*f)+phase_velocity(f,h,g) / (2.0_rp*k) * (2.0_rp*k*h / &
   TANH(2.0_rp*k*h)-1.0_rp)
ELSE
   WRITE(*,*) 'not finished'
   omega_seconde = - v_g**2/(TWOPI*f)
END IF
!
END FUNCTION omega_seconde
!
!
!
END MODULE linear_wave
