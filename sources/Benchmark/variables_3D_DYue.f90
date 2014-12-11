MODULE variables_3D
!
! Definition of variables for Dommermuth & Yue tests
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
USE ramp
!
IMPLICIT NONE
!
INTEGER, PARAMETER  :: n_lambda_x = 1
INTEGER, PARAMETER  :: n_lambda_y = 1
! Number of modes
INTEGER :: n1
INTEGER :: n2
! Array size
INTEGER, PARAMETER :: m1 = 256!128!512 ! max n1 in Dommermuth and Yue tests
INTEGER, PARAMETER :: m2 = 256!128!512 ! max n1 in Dommermuth and Yue tests
! HOS nonlinearity order
INTEGER :: M ! WARNING only among 1,2,3,4,5,7,8,9,11,14,15,17,19,23,29
! Dealiasing parameters
INTEGER, PARAMETER :: p1 = 14!21 must have p1 <= M
INTEGER, PARAMETER :: p2 = 14!21 must have p2 <= M
!
! INTEGER, DIMENSION(3), PARAMETER  :: ctype1 = (/p1, ((p1+1) * n1) / 4, 0/) ! p, N_der and N_filt
! INTEGER, DIMENSION(3), PARAMETER  :: ctype2 = (/p2, ((p2+1) * n2) / 4, 0/) ! p, N_der and N_filt
! Number of modes with the extended grid
INTEGER :: Nd1 ! required modes
INTEGER :: Nd2 ! required modes
! Extended array size
INTEGER, PARAMETER :: md1 = ((p1+1) * m1) / 2 ! required modes
INTEGER, PARAMETER :: md2 = ((p2+1) * m2) / 2 ! required modes
!
CHARACTER(LEN=3), PARAMETER       :: err = 'abs' ! ou 'rel'
!
!INTEGER :: n2_all
!INTEGER :: Nd2_all
!
INTEGER :: n1o2p1
INTEGER :: Nd1o2p1
INTEGER :: n2o2p1
INTEGER :: n2m1o2m1
INTEGER :: Nd2o2p1
!
INTEGER :: n2p1o2
INTEGER :: Nd2p1o2
!
! Array size
! FIXME: change implications of this change
! FIXME: check efficiency of resolution
INTEGER, PARAMETER  :: m1o2p1  = m1/2 + 1 !m1o2  = (m1+1)/2
INTEGER, PARAMETER  :: md1o2p1 = md1/2 + 1!md1o2 = (md1+1)/2
!
INTEGER, PARAMETER  :: md2o2p1 = md2/2 + 1
!
COMPLEX(CP), DIMENSION(md1o2p1,m2)  :: ikx, iky
COMPLEX(CP), DIMENSION(md1o2p1,md2) :: ikx_big, iky_big
!
! Wavenumbers
REAL(RP), DIMENSION(md1o2p1)   :: kx, kx2, kx3
REAL(RP), DIMENSION(md2o2p1) :: ky, ky2, ky3
! Spatial mesh
REAL(RP), DIMENSION(m1)    :: x
REAL(RP), DIMENSION(m2)    :: y
!
COMPLEX(CP), DIMENSION(m1o2p1,m2)      	:: a_eta, a_phis, temp_C_n, a_phiz, a_phiz_RF
COMPLEX(CP), DIMENSION(md1o2p1,md2)    	:: temp_C_Nd, temp2_C_Nd
REAL(RP), DIMENSION(m1,m2)           	:: phiz, W1
REAL(RP), DIMENSION(m1,m2)           	:: temp_R_n
REAL(RP), DIMENSION(md1,md2)         	:: temp_R_Nd, temp2_R_Nd, temp3_R_Nd
REAL(RP), DIMENSION(md1,md2)         	:: etax, phisx, etay, phisy, gradeta2
REAL(RP), DIMENSION(md1,md2)         	:: geta2phiz, phiz2, geta2phiz2
REAL(RP), DIMENSION(md1o2p1,md2o2p1)    :: omega, goomega, omega_p, omega_m, k, c
REAL(RP), DIMENSION(:,:,:), ALLOCATABLE :: etapm_ext
REAL(RP), DIMENSION(:,:,:), ALLOCATABLE :: gradeta_square_ext
REAL(RP), DIMENSION(:), ALLOCATABLE     :: oneoj
REAL(RP), DIMENSION(:,:,:), ALLOCATABLE :: kth, kth_all
! dealiasing
INTEGER :: i_dealias(2), N_dea(2), order_max(2)
! differentiation
INTEGER :: N_der(2)
! filtering
INTEGER :: i_filt(2)
!
INTEGER :: n1c, n1c_filt, n2c, n2c_filt
!
! Time marching parameters
REAL(RP) :: toler
! General
REAL(RP) :: xlen, ylen, T_stop, f_out
INTEGER  :: tecplot, i_case
! Forward speed
INTEGER  :: U_type
REAL(RP) :: U_start, U, t_min, k_min
! Pressure patch
! Pressure patch
INTEGER  :: i_patch
REAL(RP) :: patch_1, patch_2, U_pres, X_pres, P
! Physical data
REAL(RP) :: depth, grav
!
! Length and time scales
REAL(RP) :: L, T, L_out, T_out, eta_out
INTEGER  :: i_out_dim
REAL(RP) :: g_star, xlen_star, ylen_star, T_stop_star, f_out_star, depth_star
REAL(RP) :: patch_1_star, patch_2_star
!!
!! Absorbing zone
!INTEGER                    :: i_abs
!REAL(RP), DIMENSION(2)     :: abs_start, abs_strength
!REAL(RP), DIMENSION(m1,m2) :: nu
!TYPE(ramp_def)             :: ramp_abs_x, ramp_abs_y
!
! GD: add
!
REAL(RP), PARAMETER :: Ta= 0.0_rp
INTEGER,  PARAMETER :: n=2
!
REAL(RP), DIMENSION(md1o2p1,m2)    :: omega_n2,goomega_n2
REAL(RP), DIMENSION(m2)            :: ky_n2
! GD : may be useful CHECK
! used in the initialization of irreg fields
! FIXME: make it cleaner
INTEGER, PARAMETER ::  ikp = m1/4
INTEGER, PARAMETER ::  ithp = 45
REAL(RP), DIMENSION(m1o2p1,m2) :: k_abs,theta_abs
REAL(RP) :: E_tot, dth
REAL(RP),DIMENSION(ithp) :: theta_base
!
! GD: add for velocities
! FIXME: just create if output needed?
!
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: modesspec,modesspecx,modesspecy,modesspecz,modesspect,modesSL
REAL(RP), DIMENSION(md1,md2)      :: phizMm1, phiz2Mm2, phizMm2
REAL(RP), DIMENSION(m1,m2)        :: phi_SL, vitx_SL, vity_SL, vitz_SL, dphit_SL
REAL(RP), DIMENSION(m1,m2)        :: phiref_SL, vitxref_SL, vityref_SL, vitzref_SL, dphitref_SL
REAL(RP), DIMENSION(m1,m2)        :: vitx2ref_SL, vity2ref_SL, vitz2ref_SL
!
!
CONTAINS
!
!
!
FUNCTION extend_C(a)
!
IMPLICIT NONE
!
COMPLEX(CP), DIMENSION(1:m1o2p1,1:m2)   :: a
COMPLEX(CP), DIMENSION(1:md1o2p1,1:md2) :: extend_C
!
extend_C(1:Nd1o2p1,1:Nd2)           = ((0.0_rp, 0.0_rp))
extend_C(1:n1o2p1,1:n2o2p1)         = a(1:n1o2p1,1:n2o2p1)
extend_C(1:n1o2p1,Nd2-n2m1o2m1:Nd2) = a(1:n1o2p1,n2o2p1+1:n2)
!IF (iseven(n2) .AND. p2 /= 1) THEN
!   extend_C(1:n1o2p1,n2o2p1)       = extend_C(1:n1o2p1,n2o2p1) * 0.5_rp
!   !extend_C(1:n1o2p1,Nd2-n2o2p1+2) = extend_C(1:n1o2p1,n2o2p1)
!   ! FIXME : check this is equivalent : should be ok...
!   ! FIXME : does not work for n2 odd... but n2 even!
!   extend_C(1:n1o2p1,Nd2-n2m1o2m1-1) = extend_C(1:n1o2p1,n2o2p1)
!END IF
! FIXME : test may 2013... sep 2014: not useful
!IF (iseven(n1) .AND. p1 /= 1) THEN
!    extend_C(n1o2p1,1:Nd2)       = ((0.0_rp, 0.0_rp)) 
!    !extend_C(n1o2p1,1:Nd2)       = extend_C(n1o2p1,1:Nd2) * 0.5_rp
!ENDIF
! Added, sep 2014: really useful?
IF (iseven(n2) .AND. p2 /= 1) THEN
    !extend_C(1:Nd1o2p1,n2o2p1)       = ((0.0_rp, 0.0_rp))
    extend_C(1:n1o2p1,n2o2p1)       = extend_C(1:n1o2p1,n2o2p1)*0.5_rp
    extend_C(1:n1o2p1,Nd2-n2o2p1+2) = extend_C(1:n1o2p1,n2o2p1)
    IF (iseven(n1) .AND. p1 /= 1) THEN
    	extend_C(n1o2p1,n2o2p1) = extend_C(n1o2p1,n2o2p1)*2.0_rp
    	extend_C(n1o2p1,Nd2-n2o2p1+2) = extend_C(n1o2p1,n2o2p1)
    ENDIF
ENDIF
!
END FUNCTION extend_C
!
!
!
FUNCTION reduce_C(a_big)
!
IMPLICIT NONE
!
COMPLEX(CP), DIMENSION(1:md1o2p1,1:md2) :: a_big
COMPLEX(CP), DIMENSION(1:m1o2p1,1:m2)   :: reduce_C
!
reduce_C(1:n1o2p1,1:n2)        = ((0.0_rp, 0.0_rp))
reduce_C(1:n1o2p1,1:n2o2p1)    = a_big(1:n1o2p1,1:n2o2p1)
reduce_C(1:n1o2p1,n2o2p1+1:n2) = a_big(1:n1o2p1,Nd2-n2m1o2m1:Nd2)
!IF (iseven(n2) .AND. p2 /= 1) THEN
!   reduce_C(1:n1o2p1,n2o2p1) = reduce_C(1:n1o2p1,n2o2p1) * 2.0_rp
!END IF
! FIXME : check if the 2 following are necessary...
!IF (iseven(n1)) reduce_C(n1o2p1,1:n2) = 0.0_cp
!IF (iseven(n2)) reduce_C(1:n1o2p1,n2o2p1) = 0.0_cp
! FIXME : test may 2013
! Set to zero since it has to be real? (or check if it is the case?)
IF (iseven(n1) .AND. p1 /= 1) THEN
    reduce_C(n1o2p1,1:n2) = ((0.0_rp, 0.0_rp))
    !reduce_C(n1o2p1,1:n2) = reduce_C(n1o2p1,1:n2) * 2.0_rp
ENDIF
! For exact correspondance between x and y
IF (iseven(n2) .AND. p2 /= 1) THEN
    reduce_C(1:n1o2p1,n2o2p1) = ((0.0_rp, 0.0_rp))
    !reduce_C(1:n1o2p1,n2o2p1) = reduce_C(1:n1o2p1,n2o2p1) * 2.0_rp
ENDIF
!
END FUNCTION reduce_C
!
!
!
LOGICAL FUNCTION iseven(n)
!
IMPLICIT NONE
!
INTEGER :: n
!
iseven = (MOD(n,2) == 0)
!
END FUNCTION iseven
!
!
!
END MODULE variables_3D
