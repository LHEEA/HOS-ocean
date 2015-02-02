MODULE filters
!
! This module contains different 'filters' used in the computation for dealiasing procedure
! or filtering of highest modes
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
USE fourier_r2c
!
IMPLICIT NONE
!
!
!
CONTAINS
!
!
!
FUNCTION filter(f,ns1,ns2)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,m2) :: filter
REAL(RP), DIMENSION(m1,m2) :: f
INTEGER :: ns1, ns2,filt_rad
INTEGER :: i1,i2
REAL(RP) :: k_filt
!
! For possible radial filtering (i.e. on k instead of kx and ky) (if =1)
filt_rad = 0
!
CALL space_2_fourier(f, temp_C_n)
!
IF(filt_rad == 0)THEN
	temp_C_n(ns1+1:n1o2p1,1:n2)           = ((0.0_rp, 0.0_rp))
	temp_C_n(1:n1o2p1,ns2+1:n2-(ns2+1)+2) = ((0.0_rp, 0.0_rp))
ELSEIF(filt_rad == 1)THEN
    k_filt = max(kx(n1o2p1),ky_n2(n2o2p1))
    DO i2 = 1,n2o2p1
        DO i1 = 1,n1o2p1
    	    IF(k_abs(i1,i2) .gt. k_filt)THEN
                temp_C_n(i1,i2)           = ((0.0_rp, 0.0_rp))
                temp_C_n(i1,n2-(i2+1)+2)  = ((0.0_rp, 0.0_rp))
            ENDIF
        ENDDO
    ENDDO
ENDIF
!
CALL fourier_2_space(temp_C_n, temp_R_n)
!
filter = temp_R_n
!
END FUNCTION filter
!
!
!
FUNCTION filter_ext(f,ns1,ns2)
!
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,m2)   :: filter_ext
REAL(RP), DIMENSION(md1,md2) :: f
INTEGER :: ns1, ns2,filt_rad
INTEGER :: i1,i2
REAL(RP) :: k_filt
!
! For possible radial filtering (i.e. on k instead of kx and ky) (if =1)
filt_rad = 0
!
CALL space_2_fourier_big(f, temp_C_Nd)
!
temp_C_n = reduce_C(temp_C_Nd)
!
IF(filt_rad == 0)THEN
	temp_C_n(ns1+1:n1o2p1,1:n2)           = ((0.0_rp, 0.0_rp))
	temp_C_n(1:n1o2p1,ns2+1:n2-(ns2+1)+2) = ((0.0_rp, 0.0_rp))
ELSEIF(filt_rad == 1)THEN
    k_filt = max(kx(n1o2p1),ky_n2(n2o2p1))
    DO i2 = 1,n2o2p1
        DO i1 = 1,n1o2p1
            IF(k_abs(i1,i2) .gt. k_filt)THEN
                temp_C_n(i1,i2)           = ((0.0_rp, 0.0_rp))
                temp_C_n(i1,n2-(i2+1)+2)  = ((0.0_rp, 0.0_rp))
            ENDIF
        ENDDO
    ENDDO
ENDIF
!
CALL fourier_2_space(temp_C_n, temp_R_n)
!
filter_ext = temp_R_n
!
END FUNCTION filter_ext
!
!
!
SUBROUTINE choose_filter(ns1,ns2)
!
IMPLICIT NONE
!
INTEGER, INTENT(OUT) :: ns1, ns2
!
! Any filter ?
IF (i_filt(1) == 1) THEN
   ns1 = MAX(1,MIN(n1c_filt, n1o2p1)) ! Filtering
ELSE
   ns1 = MAX(1,MIN(n1c, n1o2p1))
END IF
!
IF (i_filt(2) == 1) THEN
   ns2 = MAX(1,MIN(n2c_filt, n2o2p1)) ! Filtering
ELSE
   ns2 = MAX(1,MIN(n2c, n2o2p1))
END IF
!
! GD change Feb 2013 : uncommented
! FIXME: check the influence, both are needed? n2 not sure...
! FIXME: influence to test on regular waves+irregular waves
! Removed 09/2014
IF ((i_case == 9).OR.(i_case == 3).OR.(i_case == 31).OR.(i_case == 32)) THEN
   !IF (iseven(n1)) ns1 = MIN(ns1, n1o2p1-1)
   !IF (iseven(n2)) ns2 = MIN(ns2, n2o2p1-1)
END IF
!
END SUBROUTINE choose_filter
!
!
!
SUBROUTINE dealias(order, var_ext)
!
IMPLICIT NONE
!
INTEGER                                         :: order
REAL(RP), DIMENSION(1:md1,1:md2), INTENT(INOUT) :: var_ext
!
! Prevent from doing 4 FFTs instead of 2 in 3D configurations
IF (((i_dealias(1) == 1).AND.(MOD(order-1,order_max(1)) == 0)).OR.((i_dealias(2) == 1).AND.(MOD(order-1,order_max(2)) == 0))) THEN
	CALL space_2_fourier_big(var_ext, temp_C_Nd)
ENDIF
!
IF (i_dealias(1) == 1) THEN
   IF (MOD(order-1,order_max(1)) == 0) THEN ! partial dealiasing along x-direction
      !CALL space_2_fourier_big(var_ext, temp_C_Nd)
      temp_C_Nd(N_dea(1)+1:Nd1o2p1,1:Nd2) = ((0.0_rp, 0.0_rp))
      !CALL fourier_2_space_big(temp_C_Nd, var_ext)
   END IF
END IF
IF (i_dealias(2) == 1) THEN
   IF (MOD(order-1,order_max(2)) == 0) THEN ! partial dealiasing along y-direction
      !CALL space_2_fourier_big(var_ext, temp_C_Nd)
      temp_C_Nd(1:Nd1o2p1,N_dea(2)+1:Nd2-(N_dea(2)+1)+2) = ((0.0_rp, 0.0_rp))
      !CALL fourier_2_space_big(temp_C_Nd, var_ext)
   END IF
END IF
!
IF (((i_dealias(1) == 1).AND.(MOD(order-1,order_max(1)) == 0)).OR.((i_dealias(2) == 1).AND.(MOD(order-1,order_max(2)) == 0))) THEN
	CALL fourier_2_space_big(temp_C_Nd, var_ext)
ENDIF
!
END SUBROUTINE dealias
!
!
!
END MODULE filters
