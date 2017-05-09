MODULE random_numbers
!
! This module defines a function for effective random number generation
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
USE variables_3D, only: n1o2p1,n2
USE type
USE iso_fortran_env, only: int64
!USE IFPORT
!
CONTAINS
!
SUBROUTINE init_random_seed()
!
IMPLICIT NONE
INTEGER, ALLOCATABLE :: seed(:)
INTEGER :: i, n, un, istat, dt(8), pid
INTEGER(int64) :: t
!
CALL random_seed(size = n)
ALLOCATE(seed(n))
! First try if the OS provides a random number generator
OPEN(newunit=un, file="/dev/urandom", access="stream", &
    form="unformatted", action="read", status="old", iostat=istat)
IF (istat == 0) THEN
    READ(un) seed
    CLOSE(un)
ELSE
    ! Fallback to XOR:ing the current time and pid. The PID is
    ! useful in case one launches multiple instances of the same
    ! program in parallel.
    CALL SYSTEM_CLOCK(t)
    IF (t == 0) THEN
        CALL DATE_AND_TIME(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24_int64 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
    ENDIF
    pid = GETPID()
    t = IEOR(t, INT(pid, KIND(t)))
    DO i = 1, n
        seed(i) = lcg(t)
    ENDDO
ENDIF
!
CALL random_seed(put=seed)
!
END SUBROUTINE init_random_seed
!

FUNCTION lcg(s)
! This simple PRNG might not be good enough for real work, but is
! sufficient for seeding a better PRNG.
INTEGER :: lcg
INTEGER(int64) :: s
!
IF (s == 0) THEN
    s = 104729
ELSE
    s = mod(s, 4294967296_int64)
ENDIF
!
s = MOD(s * 279470273_int64, 4294967291_int64)
lcg = INT(MOD(s, INT(HUGE(0), int64)), KIND(0))

END FUNCTION lcg
!
SUBROUTINE init_not_random_seed()
!
! This functions create a seed which is the same for every run at a given n1 and n2
!
IMPLICIT NONE
INTEGER, ALLOCATABLE :: seed(:)
INTEGER :: i, n, un, istat, dt(8), pid
INTEGER(int64) :: t
!
CALL random_seed(size = n)
ALLOCATE(seed(n))
DO i=1,n
    seed(i)=n1o2p1*n2*i
ENDDO
print*,'n',n
!
CALL random_seed(put=seed)
!
END SUBROUTINE init_not_random_seed
!
!
!
END MODULE random_numbers