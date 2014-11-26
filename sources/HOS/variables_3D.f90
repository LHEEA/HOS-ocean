MODULE variables_3D
!
! This module defines the different global variables used in HOS-ocean
! Especially, this defines number of points/modes and HOS order
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
USE ramp
!
IMPLICIT NONE
!
! Number of modes
INTEGER, PARAMETER :: n1 = 256	
INTEGER, PARAMETER :: n2 = 64
! Array size
INTEGER, PARAMETER :: m1 = n1
INTEGER, PARAMETER :: m2 = n2
! HOS nonlinearity order
INTEGER, PARAMETER :: M = 3
! Dealiasing parameters
INTEGER, PARAMETER :: p1 = M ! must have p1 <= M
INTEGER, PARAMETER :: p2 = M ! must have p2 <= M
!
INTEGER, DIMENSION(3), PARAMETER  :: ctype1 = (/p1, 2*(n1/2+1), 0/) ! p, N_der and N_filt
INTEGER, DIMENSION(3), PARAMETER  :: ctype2 = (/p2, 2*(n2/2+1), 0/) ! p, N_der and N_filt
! Number of modes with the extended grid
INTEGER, PARAMETER :: Nd1 = ((p1+1) * n1) / 2 ! required modes
INTEGER, PARAMETER :: Nd2 = ((p2+1) * n2) / 2 ! required modes
! Extended array size
INTEGER, PARAMETER :: md1 = ((p1+1) * m1) / 2 ! required modes
INTEGER, PARAMETER :: md2 = ((p2+1) * m2) / 2 ! required modes
!
CHARACTER(LEN=3), PARAMETER       :: err = 'abs' ! ou 'rel'
INTEGER, PARAMETER  			  :: RF_dealias = 1
!
!
INTEGER, PARAMETER  :: n2_all  = n2
INTEGER, PARAMETER  :: Nd2_all = Nd2
INTEGER, PARAMETER  :: n1o2p1   = n1/2+1
INTEGER, PARAMETER  :: n2o2p1   = n2/2+1
INTEGER, PARAMETER  :: n2p1o2   = (n2+1)/2
INTEGER, PARAMETER  :: n2m1o2m1 = (n2-1)/2-1
INTEGER, PARAMETER  :: Nd1o2p1  = Nd1/2+1
INTEGER, PARAMETER  :: Nd2o2p1  = Nd2/2+1
INTEGER, PARAMETER  :: Nd2p1o2   = (Nd2+1)/2
! Array size
! FIXME: change implications of this change
! FIXME: check efficiency of resolution
INTEGER, PARAMETER  :: m1o2p1   = m1/2+1
INTEGER, PARAMETER  :: md1o2p1  = md1/2+1
! FIXME: maybe not necessary (only used in kth_vel)
INTEGER, PARAMETER  :: md2o2p1 = md2/2 + 1
!
COMPLEX(CP), DIMENSION(md1o2p1,n2)  :: ikx, iky
COMPLEX(CP), DIMENSION(md1o2p1,Nd2)	:: ikx_big, iky_big
!
! Wavenumbers
REAL(RP), DIMENSION(md1o2p1)   	:: kx, kx2, kx3
REAL(RP), DIMENSION(Nd2o2p1) 	:: ky, ky2, ky3
! Spatial mesh
REAL(RP), DIMENSION(m1)    :: x
REAL(RP), DIMENSION(n2)    :: y
!
COMPLEX(CP), DIMENSION(m1o2p1,n2)      		:: a_eta, a_phis, temp_C_n
COMPLEX(CP), DIMENSION(md1o2p1,Nd2)    		:: temp_C_Nd, temp2_C_Nd
REAL(RP), DIMENSION(m1,n2)           		:: phiz, W1, phis, eta
REAL(RP), DIMENSION(m1,n2)           		:: temp_R_n
REAL(RP), DIMENSION(md1,Nd2)         		:: temp_R_Nd, temp2_R_Nd, temp3_R_Nd
REAL(RP), DIMENSION(md1,Nd2,M+1)     		:: etapm_ext
REAL(RP), DIMENSION(md1,Nd2,MAX((M+1)/2,5)) :: gradeta_square_ext
REAL(RP), DIMENSION(md1,Nd2)         		:: etax, phisx, etay, phisy, gradeta2
REAL(RP), DIMENSION(md1,Nd2)         		:: geta2phiz, phiz2, geta2phiz2
REAL(RP), DIMENSION(M)               		:: oneoj
REAL(RP), DIMENSION(md1o2p1,Nd2,MAX(M,2))	:: kth
REAL(RP), DIMENSION(md1o2p1,Nd2)       		:: k4
REAL(RP), DIMENSION(md1o2p1,Nd2o2p1)   		:: omega, omega_p, omega_m, k, c, goomega
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
! GD : may be useful CHECK
! used in the initialization of irreg fields
! FIXME: make it cleaner
INTEGER, PARAMETER ::  ikp = n1/4
INTEGER, PARAMETER ::  ithp = 45
REAL(RP), DIMENSION(m1o2p1,n2) :: k_abs,theta_abs
REAL(RP) :: E_tot, dth
REAL(RP),DIMENSION(ithp) :: theta_base
!
! GD : some parameters are necessary... FIXME: true?
REAL(RP) :: depth, grav
!
! Length and time scales
REAL(RP) :: L, T, L_out, T_out, eta_out
INTEGER  :: i_out_dim
REAL(RP) :: g_star, xlen_star, ylen_star, T_stop_star, f_out_star, depth_star
!
! Volume and energy
REAL(RP) :: volume, E_o(4)
!
! Output numbers
INTEGER :: i_3d, i_a_3d, i_2d, i_prob
!
! Irregular waves
INTEGER  :: i_sw ! Outputs...
INTEGER  :: n                    		! Spreading parameter
REAL(RP) :: E_cible,gamma,beta,Ta,Tp_real,Hs_real,Tp
!
REAL(RP), DIMENSION(md1o2p1,n2)  :: omega_n2,goomega_n2
REAL(RP), DIMENSION(n2)          :: ky_n2
!
! GD: add for velocities
! FIXME: just create if output needed?
!
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: modesspec,modesspecx,modesspecy,modesspecz,modesspect
REAL(RP), DIMENSION(md1,md2)      :: phizMm1, phiz2Mm2, phizMm2
REAL(RP), DIMENSION(m1,m2)        :: phi_SL, vitx_SL, vity_SL, vitz_SL, dphit_SL
REAL(RP), DIMENSION(m1,m2)        :: phiref_SL, vitxref_SL, vityref_SL, vitzref_SL, dphitref_SL
REAL(RP), DIMENSION(m1,m2)        :: vitx2ref_SL, vity2ref_SL, vitz2ref_SL
!
! GD: add for probes
!
INTEGER, PARAMETER :: maxprobes=5
INTEGER :: nprobes
REAL(RP), DIMENSION(maxprobes) :: xprobe, yprobe, eta_probe
!
CONTAINS
!
FUNCTION extend_C(a)
!
IMPLICIT NONE
!
COMPLEX(CP), DIMENSION(1:m1o2p1,1:n2)   :: a
COMPLEX(CP), DIMENSION(1:md1o2p1,1:Nd2) :: extend_C
!
extend_C(1:Nd1o2p1,1:Nd2)           = 0.0_cp
extend_C(1:n1o2p1,1:n2o2p1)         = a(1:n1o2p1,1:n2o2p1)
extend_C(1:n1o2p1,Nd2-n2m1o2m1:Nd2) = a(1:n1o2p1,n2o2p1+1:n2)
!IF (iseven(n2) .AND. p2 /= 1) THEN
!   extend_C(1:n1o2p1,n2o2p1)       = extend_C(1:n1o2p1,n2o2p1) * 0.5_rp
!   !extend_C(1:n1o2p1,Nd2-n2o2p1+2) = extend_C(1:n1o2p1,n2o2p1)
!   ! FIXME : check this is equivalent : should be ok...
!   ! FIXME : does not work for n2 odd... but n2 even!
!   extend_C(1:n1o2p1,Nd2-n2m1o2m1-1) = extend_C(1:n1o2p1,n2o2p1)
!END IF
!IF (iseven(n2)) extend_C(1:Nd1o2p1,n2o2p1) = 0.0_cp
!! FIXME : test may 2013
!IF (iseven(n1) .AND. p1 /= 1) THEN
!!    extend_C(n1o2p1,2:Nd2)       = extend_C(n1o2p1,2:Nd2) * 0.5_rp
!    extend_C(n1o2p1,1:Nd2)       = 0.0_rp
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
COMPLEX(CP), DIMENSION(1:md1o2p1,1:Nd2) :: a_big
COMPLEX(CP), DIMENSION(1:m1o2p1,1:n2)   :: reduce_C
!
reduce_C(1:n1o2p1,1:n2)        = 0.0_cp
reduce_C(1:n1o2p1,1:n2o2p1)    = a_big(1:n1o2p1,1:n2o2p1)
reduce_C(1:n1o2p1,n2o2p1+1:n2) = a_big(1:n1o2p1,Nd2-n2m1o2m1:Nd2)
!IF (iseven(n2) .AND. p2 /= 1) THEN
!   reduce_C(1:n1o2p1,n2o2p1) = reduce_C(1:n1o2p1,n2o2p1) * 2.0_rp
!END IF
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
