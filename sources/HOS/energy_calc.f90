MODULE energy_calc
!
! This module is related to the evaluation of free-surface wave fields energy.
! 
! contains  : one function named calc_energy
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
!
IMPLICIT NONE
!
CONTAINS
!
!
!
FUNCTION calc_energy(a_eta, a_phis, da_eta)
!
IMPLICIT NONE
!
! returns the wavefield energy (potential and kinetic energy)
!
! inputs    :  a_eta   = free surface elevation (fourier amplitudes)
!              a_phis  = free surface potential (fourier amplitudes)
!              da_eta  = time derivative of the free surface elevation (fourier amplitudes)
!
! needs     : non-dimensional gravity, g_star
!
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: a_eta, a_phis, da_eta
REAL(RP), DIMENSION(4)                        :: calc_energy
!
 calc_energy(:)=0.0_rp
! Potential energy
!
!  int eta * eta dx (gravity related)
!   eta on the extended grid
temp_C_Nd = extend_C(a_eta)
CALL fourier_2_space_big(temp_C_Nd, temp_R_Nd)
!   eta * eta
temp_R_Nd(1:Nd1,1:Nd2) = temp_R_Nd(1:Nd1,1:Nd2) * temp_R_Nd(1:Nd1,1:Nd2)
!   integral
CALL space_2_fourier_big(temp_R_Nd, temp_C_Nd)
!
calc_energy(1) = g_star * REAL(temp_C_Nd(1,1))
!
! Kinetic energy
!
!  int deta_dt * phis dx
!
!   deta_dt on big mesh
temp_C_Nd = extend_C(da_eta)
CALL fourier_2_space_big(temp_C_Nd, temp_R_Nd)
!
!   phis on big mesh
temp_C_Nd = extend_C(a_phis)
CALL fourier_2_space_big(temp_C_Nd, temp2_R_Nd)
!
!   product deta_dt times phis
temp_R_Nd(1:Nd1,1:Nd2) = temp_R_Nd(1:Nd1,1:Nd2) * temp2_R_Nd(1:Nd1,1:Nd2)
!   integral
CALL space_2_fourier_big(temp_R_Nd, temp_C_Nd)
!
calc_energy(3) = REAL(temp_C_Nd(1,1))
calc_energy(4) = calc_energy(1) + calc_energy(2) + calc_energy(3)
calc_energy(1:4) = 0.5_rp * calc_energy(1:4)
!
RETURN
!
END FUNCTION calc_energy
!
END MODULE energy_calc
