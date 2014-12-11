MODULE maths
!
! This module contains useful mathematical functions
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
INTERFACE norme
   MODULE PROCEDURE norme_r,norme_d,norme_vr,norme_vd
END INTERFACE 
INTERFACE ecart_type
   MODULE PROCEDURE ecart_type_r,ecart_type_d
END INTERFACE 
INTERFACE variance
   MODULE PROCEDURE variance_r,variance_d
END INTERFACE 
INTERFACE dichoto
   MODULE PROCEDURE dichoto_r,dichoto_d
END INTERFACE 
INTERFACE phase
   MODULE PROCEDURE phase_r,phase_d !,phase_rv,phase_dv
END INTERFACE 
!
!
!
CONTAINS
!
!
!
FUNCTION norme_r(x)
!
IMPLICIT NONE
!
REAL(SP),INTENT(IN) :: x
REAL(SP)            :: norme_r
!
norme_r = SQRT( x*x )
!
END FUNCTION norme_r
!
!
!
FUNCTION norme_d(x)
!
IMPLICIT NONE
!
REAL(DP),INTENT(IN) :: x
REAL(DP)            :: norme_d
!
norme_d = SQRT( x*x )
!
END FUNCTION norme_d
!
!
!
FUNCTION norme_vr(x)
!
IMPLICIT NONE
!
REAL(SP),DIMENSION(:),INTENT(IN) :: x
REAL(SP)                         :: norme_vr
!
norme_vr = SQRT( DOT_PRODUCT( x,x ) )
!
END FUNCTION norme_vr
!
!
!
FUNCTION norme_vd(x)
!
IMPLICIT NONE
!
REAL(DP),DIMENSION(:),INTENT(IN) :: x
REAL(DP)                         :: norme_vd
!
norme_vd = SQRT( DOT_PRODUCT( x,x ) )
!
END FUNCTION norme_vd
!
!
!
FUNCTION variance_r(x)
!
IMPLICIT NONE
!
REAL(SP), DIMENSION(:),INTENT(IN) :: x
REAL(SP)                          :: variance_r
REAL(SP), DIMENSION(SIZE(x))      :: x_moy_v
REAL(SP)                          :: x_moy
INTEGER                           :: n_x
!
n_x     = SIZE( x )
IF ( n_x == 1 ) THEN
   variance_r = 1.0_sp
ELSE
   x_moy   = SUM( x ) / n_x
   x_moy_v = x_moy
   x_moy_v = x - x_moy_v
   variance_r = DOT_PRODUCT( x - x_moy_v,x - x_moy_v ) / ( n_x - 1 )
END IF
!
END FUNCTION variance_r
!
!
!
FUNCTION variance_d(x)
!
IMPLICIT NONE
!
REAL(DP), DIMENSION(:),INTENT(IN) :: x
REAL(DP)                          :: variance_d
REAL(DP), DIMENSION(SIZE(x))      :: x_moy_v
REAL(DP)                          :: x_moy
INTEGER                           :: n_x
!
n_x     = SIZE( x )
IF ( n_x == 1 ) THEN
   variance_d = 1.0_dp
ELSE
   x_moy   = SUM( x ) / n_x
   x_moy_v = x_moy
   x_moy_v = x - x_moy_v
   variance_d = DOT_PRODUCT( x_moy_v,x_moy_v ) / ( n_x - 1 )
END IF
!
END FUNCTION variance_d
!
!
!
FUNCTION ecart_type_r(x)
!
IMPLICIT NONE
!
REAL(SP), DIMENSION(:),INTENT(IN) :: x
REAL(SP)                          :: ecart_type_r
!
ecart_type_r = SQRT( variance( x ) )
!
END FUNCTION ecart_type_r
!
!
!
FUNCTION ecart_type_d(x)
!
IMPLICIT NONE
!
REAL(DP), DIMENSION(:),INTENT(IN) :: x
REAL(DP)                          :: ecart_type_d
!
ecart_type_d = SQRT( variance( x ) )
!
END FUNCTION ecart_type_d
!
!
!
FUNCTION dichoto_r(func,f0,a,b,prec)
!
!   Dichotomy to solve func(x)=f0 between a and b, with a relative accuracy prec on x
!
! inputs :
!   func => function used (declared as EXTERNAL in the calling program)
!   f0   => value of the wanted function
!   a,b  => boundaries of the search interval
!   prec => accuracy needed on x
!
IMPLICIT NONE
!
REAL(SP) :: dichoto_r,f0,a,b,prec,func
REAL(SP) :: x1,x2,x3
!
IF ( ((func(a)-f0)*(func(b)-f0)) >= 0.0_sp ) THEN
   IF ( ABS(func(a)-f0) <= REAL(tiny) ) THEN
      dichoto_r = a
   RETURN
   ELSE IF (ABS(func(b)-f0) <= REAL(tiny)) THEN
      dichoto_r = b
   RETURN
   END IF
   WRITE(*,*) 'a and b are not correctly chosen in dichoto'
END IF
!
x1 = a
x2 = b
!
DO WHILE (ABS(x1-x2) > prec*ABS(x2))
   x3=(x1+x2)/2.0_sp
   IF (ABS(func(x3)-f0) <= REAL(tiny)) THEN
      dichoto_r = x3
      RETURN
   ELSE IF ((func(a)-f0)*(func(x3)-f0) < 0.0) THEN
      x2 = x3
   ELSE
      x1 = x3
   END IF
END DO
dichoto_r = x3
!
END FUNCTION dichoto_r
!
!
!
FUNCTION dichoto_d(func,f0,a,b,prec)
!
!   Dichotomy to solve func(x)=f0 between a and b, with a relative accuracy prec on x
!
! inputs :
!   func => function used (declared as EXTERNAL in the calling program)
!   f0   => value of the wanted function
!   a,b  => boundaries of the search interval
!   prec => accuracy needed on x
!
IMPLICIT NONE
!
REAL(DP) :: dichoto_d,f0,a,b,prec,func
REAL(DP) :: x1,x2,x3
!
IF ( ((func(a)-f0)*(func(b)-f0)) >= 0.0_dp ) THEN
   IF ( ABS(func(a)-f0) <= tiny ) THEN
      dichoto_d = a
   RETURN
   ELSE IF (ABS(func(b)-f0) <= tiny) THEN
      dichoto_d = b
   RETURN
   END IF
   WRITE(*,*) 'a and b are not correctly chosen in dichoto'
END IF
!
x1 = a
x2 = b
!
DO WHILE (ABS(x1-x2) > prec*ABS(x2))
   x3=(x1+x2)/2.0_dp
   IF (ABS(func(x3)-f0) <= tiny) THEN
      dichoto_d = x3
      RETURN
   ELSE IF ((func(a)-f0)*(func(x3)-f0) < 0.0_dp) THEN
      x2 = x3
   ELSE
      x1 = x3
   END IF
END DO
!
dichoto_d = x3
!
END FUNCTION dichoto_d
!
!
!
ELEMENTAL FUNCTION phase_r(z)
!
! This function gives the phase of the complex number (between 0 and 2*PI)
!
USE type
!
IMPLICIT NONE
!
COMPLEX(SPC), INTENT(IN) :: z
COMPLEX(SPC)             :: ztemp
REAL(SP)                 :: phase_r, ptemp
!
IF ( ABS(z) <= REAL(tiny) ) THEN
   phase_r = 0.0_sp
   RETURN
ELSE
   ztemp = z / ABS(z)
   ptemp = ASIN( AIMAG( ztemp ) )
   IF ( REAL(ztemp) >= 0.0_sp ) THEN
      IF ( AIMAG(ztemp) >= 0.0_sp ) THEN
         phase_r = ptemp
         RETURN
      ELSE
         phase_r = REAL(TWOPI) + ptemp
         RETURN
      END IF
   ELSE
      phase_r = REAL(PI) - ptemp
      RETURN
   END IF
END IF
END FUNCTION phase_r
!
!
!
ELEMENTAL FUNCTION phase_d(z)
!
! This function gives the phase of the complex number (between 0 and 2*PI)
!
USE type
!
IMPLICIT NONE
!
COMPLEX(DPC), INTENT(IN) :: z
COMPLEX(DPC)             :: ztemp
REAL(DP)     :: phase_d, ptemp
!
IF ( ABS(z) <= tiny ) THEN
   phase_d = 0.0_dp
   RETURN
ELSE
   ztemp = z / ABS(z)
   ptemp = ASIN( AIMAG( ztemp ) )
   IF ( REAL(ztemp) >= 0.0_dp ) THEN
      IF ( AIMAG(ztemp) >= 0.0_dp ) THEN
         phase_d = ptemp
         RETURN
      ELSE
         phase_d = TWOPI + ptemp
         RETURN
      END IF
   ELSE
      phase_d = PI - ptemp
      RETURN
   END IF
END IF
END FUNCTION phase_d
!
!
!
FUNCTION phase_rv(z)
!
! This function gives the phase of the complex number (between 0 and 2*PI)
!
USE type
!
IMPLICIT NONE
!
COMPLEX(SPC), DIMENSION(:)   :: z
!COMPLEX(SPC)                 :: ztemp
REAL(SP), DIMENSION(SIZE(z)) :: phase_rv
!REAL(SP)                     :: ptemp
INTEGER                      :: ibou
!
DO ibou = 1, SIZE(z)
   phase_rv(ibou) = phase_r(z(ibou))
END DO
!
END FUNCTION phase_rv
!
!
!
FUNCTION phase_dv(z)
!
! This function gives the phase of the complex number (between 0 and 2*PI)
!
USE type
!
IMPLICIT NONE
!
COMPLEX(DPC), DIMENSION(:)   :: z
!COMPLEX(DPC)                 :: ztemp
REAL(DP), DIMENSION(SIZE(z)) :: phase_dv
!REAL(DP)                     :: ptemp
INTEGER                      :: ibou
!
DO ibou = 1, SIZE(z)
      phase_dv(ibou) = phase_d(z(ibou))
END DO
!
END FUNCTION phase_dv
!
!
!
END MODULE maths
