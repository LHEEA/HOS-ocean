MODULE variables_3d
!
! This module defines the different common variables for post-processing
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
INTEGER :: n1,n2
REAL(RP), ALLOCATABLE, DIMENSION(:)   :: x,y
REAL(RP), ALLOCATABLE, DIMENSION(:,:) :: eta, phis
!
REAL(RP), ALLOCATABLE, DIMENSION(:)      :: kx,ky_n2
REAL(RP), ALLOCATABLE, DIMENSION(:,:)    :: kth
COMPLEX(CP), ALLOCATABLE, DIMENSION(:,:) :: ikx,iky
!
! Input file
INTEGER               :: i_ana, i_card, tecplot, i_zvect
REAL(RP)              :: T_start, T_stop, x_min, x_max, y_min, y_max, z_min, z_max
CHARACTER(LEN=100)    :: file_3d, file_mod
!
INTEGER, PARAMETER    :: n_hdr = 34 ! Headerlines in '3d.dat' including line variables
REAL(RP), PARAMETER   :: HfoHs = 2.0_rp ! Freak wave height on Hs threshold for detection
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
INTEGER :: n1o2p1,n2o2p1
! test for fourier
INTEGER :: m1,m2,Nd1,Nd2,Nd1o2p1,m1o2p1,md1o2p1,md1,md2
!
CONTAINS
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
END MODULE variables_3d