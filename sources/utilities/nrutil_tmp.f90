MODULE nrutil_tmp
!
! This module contains utility functions to swap two elements
! Prevents from the creation of new element in main routines
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
!
IMPLICIT NONE
!
INTERFACE swap
   MODULE PROCEDURE swap_i,swap_r,swap_rv,swap_c, &
        swap_cv,swap_cm,swap_z,swap_zv,swap_zm, &
        masked_swap_rs,masked_swap_rv,masked_swap_rm
END INTERFACE
!
CONTAINS
!
SUBROUTINE swap_i(a,b)
!
! Swap the contents of a and b.
!
INTEGER, INTENT(INOUT) :: a,b
INTEGER :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_i
!
SUBROUTINE swap_r(a,b)
REAL(RP), INTENT(INOUT) :: a,b
REAL(RP) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_r
!
SUBROUTINE swap_rv(a,b)
REAL(RP), DIMENSION(:), INTENT(INOUT) :: a,b
REAL(RP), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_rv
!
SUBROUTINE swap_c(a,b)
COMPLEX(CP), INTENT(INOUT) :: a,b
COMPLEX(CP) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_c
!
SUBROUTINE swap_cv(a,b)
COMPLEX(CP), DIMENSION(:), INTENT(INOUT) :: a,b
COMPLEX(CP), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_cv
!
SUBROUTINE swap_cm(a,b)
COMPLEX(SPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
COMPLEX(SPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_cm
!
SUBROUTINE swap_z(a,b)
COMPLEX(SPC), INTENT(INOUT) :: a,b
COMPLEX(SPC) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_z
!
SUBROUTINE swap_zv(a,b)
COMPLEX(SPC), DIMENSION(:), INTENT(INOUT) :: a,b
COMPLEX(SPC), DIMENSION(SIZE(a)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zv
!
SUBROUTINE swap_zm(a,b)
COMPLEX(DPC), DIMENSION(:,:), INTENT(INOUT) :: a,b
COMPLEX(DPC), DIMENSION(size(a,1),size(a,2)) :: dum
dum=a
a=b
b=dum
END SUBROUTINE swap_zm
!
SUBROUTINE masked_swap_rs(a,b,mask)
REAL(RP), INTENT(INOUT) :: a,b
LOGICAL, INTENT(IN) :: mask
REAL(RP) :: swp
if (mask) then
   swp=a
   a=b
   b=swp
end if
END SUBROUTINE masked_swap_rs
!
SUBROUTINE masked_swap_rv(a,b,mask)
REAL(RP), DIMENSION(:), INTENT(INOUT) :: a,b
LOGICAL, DIMENSION(:), INTENT(IN) :: mask
REAL(RP), DIMENSION(size(a)) :: swp
where (mask)
   swp=a
   a=b
   b=swp
end where
END SUBROUTINE masked_swap_rv
!
SUBROUTINE masked_swap_rm(a,b,mask)
REAL(RP), DIMENSION(:,:), INTENT(INOUT) :: a,b
LOGICAL, DIMENSION(:,:), INTENT(IN) :: mask
REAL(RP), DIMENSION(size(a,1),size(a,2)) :: swp
where (mask)
   swp=a
   a=b
   b=swp
end where
END SUBROUTINE masked_swap_rm
!
!
!
FUNCTION arth(first,increment,n)
INTEGER, PARAMETER :: NPAR_ARTH=16, NPAR2_ARTH=8
REAL(RP), INTENT(IN) :: first, increment
INTEGER, INTENT(IN) :: n
REAL(RP), DIMENSION(n) :: arth
INTEGER      :: k,k2
REAL(RP) :: temp
IF (n>0) arth(1)=first
IF(n<=NPAR_ARTH) THEN
   DO k=2,n
      arth(k)=arth(k-1)+increment
   END DO
ELSE
   DO k=2,NPAR2_ARTH
      arth(k)=arth(k-1)+increment
   END DO
   temp=increment*NPAR2_ARTH
   k=NPAR2_ARTH
   DO
      IF (k>=n) EXIT
      k2=k+k
      arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
      temp=temp+temp
      k=k2
   END DO
END IF
END FUNCTION arth
!
END MODULE nrutil_tmp
