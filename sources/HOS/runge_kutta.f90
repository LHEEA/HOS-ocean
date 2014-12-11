MODULE runge_kutta
!
! This module performs the time integration with Runge-Kutta scheme: 5(4) from Cash & Karp
! The linear part of equations is integrated analytically
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
USE resol_HOS
USE variables_3D
USE fourier_r2c
USE velocities
!
IMPLICIT NONE
!
TYPE :: RK_parameters
   INTEGER :: s=6
   REAL(RP), DIMENSION(6,6) :: A ! Butcher array
   REAL(RP), DIMENSION(6)   :: b ! Butcher array: increment factor for slope
   REAL(RP), DIMENSION(6)   :: c ! Butcher array: increment factor for abscissa
   REAL(RP), DIMENSION(6)   :: e ! Butcher array: weighting factors for the error estimate
   INTEGER                  :: p
END TYPE
!
CONTAINS
!
SUBROUTINE fill_butcher_array(RK_param)
!
IMPLICIT NONE
!
TYPE(RK_parameters), INTENT(INOUT) :: RK_param
!
WRITE(*,*) '5(4) Runge Kutta Fehlberg scheme (Cash & Karp)'

! RK_param%s      = 6
! A
RK_param%A(2,1) = 1.0_rp     / 5.0_rp
RK_param%A(3,1) = 3.0_rp     / 40.0_rp
RK_param%A(4,1) = 3.0_rp     / 10.0_rp
RK_param%A(5,1) =-11.0_rp    / 54.0_rp
RK_param%A(6,1) = 1631.0_rp  / 55296.0_rp
RK_param%A(3,2) = 9.0_rp     / 40.0_rp
RK_param%A(4,2) =-9.0_rp     / 10.0_rp
RK_param%A(5,2) = 5.0_rp     / 2.0_rp
RK_param%A(6,2) = 175.0_rp   / 512.0_rp
RK_param%A(4,3) = 6.0_rp     / 5.0_rp
RK_param%A(5,3) =-70.0_rp    / 27.0_rp     !***(5,4)
RK_param%A(6,3) = 575.0_rp   / 13824.0_rp  !***(6,5)
RK_param%A(5,4) = 35.0_rp    / 27.0_rp       
RK_param%A(6,4) = 44275.0_rp / 110592.0_rp  !***(6,5)
RK_param%A(6,5) = 253.0_rp   / 4096_rp
! e
RK_param%e(1)   = 2825.0_rp / 27648.0_rp
RK_param%e(3)   = 18575.0_rp / 48384.0_rp
RK_param%e(4)   = 13525.0_rp / 55296.0_rp
RK_param%e(5)   = 277.0_rp / 14336.0_rp
RK_param%e(6)   = 1.0_rp / 4.0_rp
! c
RK_param%c(2)   = 1.0_rp / 5.0_rp
RK_param%c(3)   = 3.0_rp / 10.0_rp
RK_param%c(4)   = 3.0_rp / 5.0_rp
RK_param%c(5)   = 1.0_rp
RK_param%c(6)   = 7.0_rp / 8.0_rp
! b
RK_param%b(1)   = 37.0_rp / 378.0_rp
RK_param%b(3)   = 250.0_rp / 621.0_rp
RK_param%b(4)   = 125.0_rp / 594.0_rp
RK_param%b(6)   = 512.0_rp / 1771.0_rp
! order p
RK_param%p      = 5
!

! WRITE(*,*) '4(3) Runge Kutta Fehlberg scheme'
! ! A
! RK_param%A      = 0.0_rp
! RK_param%A(2,1) = 1.0_rp / 2.0_rp
! RK_param%A(3,2) = 1.0_rp / 2.0_rp
! RK_param%A(4,3) = 1.0_rp
! RK_param%A(5,1) =-1.0_rp
! RK_param%A(5,2) = 2.0_rp
! ! b
! RK_param%b(1)   = 1.0_rp / 6.0_rp
! RK_param%b(2)   = 1.0_rp / 3.0_rp
! RK_param%b(3)   = RK_param%b(2)
! RK_param%b(4)   = RK_param%b(1)
! ! c
! RK_param%c      = 0.0_rp
! RK_param%c(2)   = 1.0_rp / 2.0_rp
! RK_param%c(3)   = 1.0_rp / 2.0_rp
! RK_param%c(4)   = 1.0_rp
! RK_param%c(5)   = 1.0_rp
! ! e
! RK_param%e      = 0.0_rp
! RK_param%e(1)   = 1.0_rp / 6.0_rp
! RK_param%e(2)   = 2.0_rp / 3.0_rp
! RK_param%e(5)   = 1.0_rp / 6.0_rp
! ! order p
! RK_param%p      = 3
!
!WRITE(*,*) '3(2) Runge Kutta Fehlberg scheme'
!RK_A      = 0.0
!RK_A(2,1) = 1.0 / 2.0
!RK_A(3,1) = - 1.0
!RK_A(3,2) = 2.0
!RK_b(1)   = 1.0 / 6.0
!RK_b(2)   = 2.0 / 3.0
!RK_b(3)   = 1.0 / 6.0
!RK_c      = 0.0
!RK_c(2)   = 1.0 / 2.0
!RK_c(3)   = 1.0
!RK_e      = 0.0
!RK_e(1)   = 1.0 / 2.0
!RK_e(2)   = 1.0 / 2.0
!RK_p      = 3
!
END SUBROUTINE fill_butcher_array
!
!
!
SUBROUTINE RK_adapt_2var_3D_in_mo_lin(iCPUtime, RK_param, t_o, h, var_1, var_2, d_var_2, err, scale_1, scale_2)
!
IMPLICIT NONE
!
! Input variables
INTEGER, INTENT(IN)                       :: iCPUtime
TYPE(RK_parameters), INTENT(IN)           :: RK_param
REAL(RP), INTENT(IN)                      :: h, t_o
! REAL(RP), INTENT(INOUT)                   :: vol
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(INOUT) :: var_2, var_1
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(OUT)   :: d_var_2 ! d_var_1
REAL(RP), OPTIONAL                        :: err, scale_1, scale_2
REAL(RP), DIMENSION(m1,m2)              :: d_var2_r
! Local variables
COMPLEX(CP), DIMENSION(m1o2p1,m2)            :: var_1m, var_2m
COMPLEX(CP), DIMENSION(m1o2p1,m2,2)          :: erm
COMPLEX(CP), DIMENSION(m1o2p1,m2,RK_param%s) :: k_var_1, k_var_2
INTEGER                   :: jloop, ns, i1, i2
REAL(RP)                  :: time, dt
! REAL(RP)                  :: vol_m
COMPLEX(CP)               :: tmp1, tmp2
COMPLEX(CP)               :: exp_p_exp, exp_m_exp
COMPLEX(CP)               :: da_eta_dt, da_phis_dt
REAL(RP)                  :: ti, tf
!
!	CPU times inlet
IF (iCPUtime.EQ.1) THEN
	PRINT*,'entering subroutine runge4'
	CALL CPU_TIME(ti)
ENDIF
!
k_var_2(:,:,2:RK_param%s-1) = ((0.0_rp, 0.0_rp))
k_var_1(:,:,2:RK_param%s-1) = ((0.0_rp, 0.0_rp))
time = t_o
!
 CALL solveHOS_lin(var_1, k_var_1(:,:,1), var_2, k_var_2(:,:,1), time)
!
DO jloop = 2, RK_param%s
   ! free surface modes (elevation and potential)
   var_1m = ((0.0_rp, 0.0_rp))
   var_2m = ((0.0_rp, 0.0_rp))
   ! Current (jloop) Runge Kutta time step increment
   dt = h * RK_param%c(jloop)
   !
   ! remaining modes
   DO i2=1,n2
      ! x and y varying modes
      DO i1=1,n1o2p1
         tmp1 = ((0.0_rp, 0.0_rp))
         tmp2 = ((0.0_rp, 0.0_rp))
         DO ns = 1, jloop-1 ! explicit RK
            CALL calc_exps(t_o + h * RK_param%c(ns), t_o, i1, i2)
            da_eta_dt  = k_var_2(i1,i2,ns) * goomega_n2(i1,i2)
            da_phis_dt = k_var_1(i1,i2,ns)
            tmp2 = tmp2 + (    exp_p_exp * da_eta_dt - i * exp_m_exp * da_phis_dt) * RK_param%A(jloop, ns)
            tmp1 = tmp1 + (i * exp_m_exp * da_eta_dt +     exp_p_exp * da_phis_dt) * RK_param%A(jloop, ns)
         END DO
         CALL calc_exps(t_o, t_o + dt, i1, i2)
         da_eta_dt  = var_2(i1,i2) * goomega_n2(i1,i2) + h * tmp2
         da_phis_dt = var_1(i1,i2)                     + h * tmp1
         var_2m(i1,i2) = (   exp_p_exp * da_eta_dt - i * exp_m_exp * da_phis_dt) / goomega_n2(i1,i2)
         var_1m(i1,i2) = i * exp_m_exp * da_eta_dt +     exp_p_exp * da_phis_dt
      END DO
   END DO
   !
   time = t_o + dt
   CALL solveHOS_lin(var_1m, k_var_1(:,:,jloop), var_2m, k_var_2(:,:,jloop), time)

   IF (abs(dt-h) < tiny) THEN ! evaluate deta_dt
     CALL fourier_2_space(k_var_2(:,:,jloop),d_var2_r)
     d_var2_r = d_var2_r + W1
     CALL space_2_fourier(d_var2_r,d_var_2)
   ENDIF
    !
END DO
!
!
!
IF (PRESENT(err)) THEN
   !
   ! estimate of lower order
   erm(:,:,1:2) = ((0.0_rp, 0.0_rp))

   DO i2=1,n2
      ! x and y varying modes
      DO i1=1,n1o2p1
         tmp1 = ((0.0_rp, 0.0_rp))
         tmp2 = ((0.0_rp, 0.0_rp))
         DO jloop = 1, RK_param%s
            CALL calc_exps(0.0_rp, t_o + h * RK_param%c(jloop), i1, i2)
            da_eta_dt  = k_var_2(i1,i2,jloop) * goomega_n2(i1,i2)
            da_phis_dt = k_var_1(i1,i2,jloop)
            tmp2 = tmp2 + (    exp_p_exp * da_eta_dt - i * exp_m_exp * da_phis_dt) * (RK_param%b(jloop)-RK_param%e(jloop))
            tmp1 = tmp1 + (i * exp_m_exp * da_eta_dt +     exp_p_exp * da_phis_dt) * (RK_param%b(jloop)-RK_param%e(jloop))
         END DO
         erm(i1,i2,1) = h * tmp2
         erm(i1,i2,2) = h * tmp1
      END DO
   END DO
   !
   erm(:,:,1) = ABS(erm(:,:,1)) / scale_1
   erm(:,:,2) = ABS(erm(:,:,2)) / scale_2
   !
   err = MAXVAL(ABS(erm(1:n1o2p1,1:n2,1:2)))
END IF
!
! free surface modes (elevation and potential)
var_1m = ((0.0_rp, 0.0_rp))
var_2m = ((0.0_rp, 0.0_rp))
dt     = h

DO i2=1,n2
   ! x and y varying modes
   DO i1=1,n1o2p1
      tmp1 = ((0.0_rp, 0.0_rp))
      tmp2 = ((0.0_rp, 0.0_rp))
      DO jloop = 1, RK_param%s
         CALL calc_exps(t_o + h * RK_param%c(jloop), t_o, i1, i2)
         da_eta_dt  = k_var_2(i1,i2,jloop) * goomega_n2(i1,i2)
         da_phis_dt = k_var_1(i1,i2,jloop)
         tmp2 = tmp2 + (    exp_p_exp * da_eta_dt - i * exp_m_exp * da_phis_dt) * RK_param%b(jloop)
         tmp1 = tmp1 + (i * exp_m_exp * da_eta_dt +     exp_p_exp * da_phis_dt) * RK_param%b(jloop)
      END DO
      CALL calc_exps(t_o, t_o + dt, i1, i2)
      da_eta_dt  = var_2(i1,i2) * goomega_n2(i1,i2) + h * tmp2
      da_phis_dt = var_1(i1,i2)                     + h * tmp1
      var_2m(i1,i2) = (   exp_p_exp * da_eta_dt - i * exp_m_exp * da_phis_dt) / goomega_n2(i1,i2)
      var_1m(i1,i2) = i * exp_m_exp * da_eta_dt +     exp_p_exp * da_phis_dt
   END DO
END DO
!
var_2(:,:)   = var_2m(:,:)
var_1(:,:)   = var_1m(:,:)
!
!	CPU times outlet
IF (iCPUtime.EQ.1) THEN
   CALL CPU_TIME(tf)
	WRITE(*,910)'quitting subroutine runge4, total CPU time: ',tf-ti,'s'
ENDIF
910  FORMAT(a,1ES11.4,a)
!
!
!
CONTAINS
!
!
!
SUBROUTINE calc_exps(t_1, t_2, i1, i2)
!
IMPLICIT NONE
!
REAL(RP) :: t_1, t_2
INTEGER  :: i1, i2
!
COMPLEX(CP) :: exp_pip, exp_mim
REAL(RP)    :: psi
!
! FIXME : check that this is omega that should be used
psi   = omega_n2(i1,i2) * (t_2 - t_1)
exp_pip = EXP(+ i * psi)
exp_mim = EXP(- i * psi)
!
exp_p_exp = 0.5_rp * (exp_pip + exp_mim)
exp_m_exp = 0.5_rp * (exp_pip - exp_mim)
!
END SUBROUTINE calc_exps
!
END SUBROUTINE RK_adapt_2var_3D_in_mo_lin
!
!
!
END MODULE runge_kutta
