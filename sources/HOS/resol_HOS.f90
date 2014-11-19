MODULE resol_HOS
!
! This module enables the computation of free-surface boundary conditions
! This includes particularly the computation of vertical velocity with HOS method
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
USE variables_3D
USE fourier_r2c
USE filters
!
!
!
CONTAINS
!
!
!
SUBROUTINE phisxy_etaxy(a_phisrk,a_etark)
!
IMPLICIT NONE
!% INPUT VARIABLES
!	free surface potential
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_phisrk
!	free surface elevation
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_etark
!% LOCAL VARIABLES
INTEGER  :: j, iCPUtime
REAL(RP) :: ti, tf
!
iCPUtime = 0
!	CPU times inlet
if (iCPUtime.eq.1) then
	print*,'entering subroutine phisxy_etaxy'
	call CPU_TIME(ti)
endif
!
! the powers of eta (eta^m/m!) (on a 2N grid)
etapm_ext(:,:,1) = 1.0_rp
!
! m=1
temp_C_Nd = extend_C(a_etark)
!
CALL fourier_2_space_big(temp_C_Nd, etapm_ext(:,:,1+1))
!
! m=2 to M
DO j=2,M
   etapm_ext(1:Nd1,1:Nd2,j+1) = etapm_ext(1:Nd1,1:Nd2,j-1+1) * oneoj(j) * etapm_ext(1:Nd1,1:Nd2,1+1)
   CALL dealias(j, etapm_ext(:,:,j+1))
END DO
!
!   calculation of etax in the Fourier domain (on a (p+1)/2 N grid)
temp_C_Nd = extend_C(a_etark)
temp_C_Nd(1:Nd1o2p1,1:Nd2) = ikx_big(1:Nd1o2p1,1:Nd2) * temp_C_Nd(1:Nd1o2p1,1:Nd2)
CALL fourier_2_space_big(temp_C_Nd, etax)
!
!   calculation of etay in the Fourier domain (on a (p+1)/2 N grid)
temp_C_Nd = extend_C(a_etark)
temp_C_Nd(1:Nd1o2p1,1:Nd2) = iky_big(1:Nd1o2p1,1:Nd2) * temp_C_Nd(1:Nd1o2p1,1:Nd2)
CALL fourier_2_space_big(temp_C_Nd, etay)
! 
! grad(eta)^2 (on a (p+1)/2 N grid)
gradeta2(1:Nd1,1:Nd2) = etax(1:Nd1,1:Nd2) * etax(1:Nd1,1:Nd2) + etay(1:Nd1,1:Nd2) * etay(1:Nd1,1:Nd2)
CALL dealias(2,gradeta2)
!
!   calculation of phisx in the Fourier domain (on a (p+1)/2 N grid)
temp_C_Nd = extend_C(a_phisrk)
temp_C_Nd(1:Nd1o2p1,1:Nd2) = ikx_big(1:Nd1o2p1,1:Nd2) * temp_C_Nd(1:Nd1o2p1,1:Nd2)
CALL fourier_2_space_big(temp_C_Nd, phisx)
!
!   calculation of phisy in the Fourier domain (on a (p+1)/2 N grid)
temp_C_Nd = extend_C(a_phisrk)
temp_C_Nd(1:Nd1o2p1,1:Nd2) = iky_big(1:Nd1o2p1,1:Nd2) * temp_C_Nd(1:Nd1o2p1,1:Nd2)
CALL fourier_2_space_big(temp_C_Nd, phisy)
!
!	CPU times outlet
if (iCPUtime.eq.1) then
	call CPU_TIME(tf)
	write(*,910)'quitting subroutine phisxy_etaxy, total CPU time: ',tf-ti,'s'
endif
910  format(a,1ES11.4,a)
!
END SUBROUTINE phisxy_etaxy
!
!
! 
SUBROUTINE HOSphis_modes_fully_dealiased(a_phisrk)
!
IMPLICIT NONE

! input variables
!	free surface potential
COMPLEX(CP) , DIMENSION(m1o2p1,m2) :: a_phisrk
!
! local variables
!
!	order m of phi and i^th derivatives at order m
!
COMPLEX(CP), DIMENSION(md1o2p1,md2,M) :: aHOS
REAL(RP), DIMENSION(md1,md2,M)      :: phizm
REAL(RP), DIMENSION(md1,md2)        :: phizi
REAL(RP), DIMENSION(md1,md2)        :: phim_ext
!
INTEGER :: j, q, i_m, j2, ns1, ns2, Mo2
REAL(RP), DIMENSION(md1o2p1,md2) :: kth_local
!
Mo2 = M / 2
! Initializations
phizm(1:Nd1,1:Nd2,1:M) = 0.0_rp
!
aHOS(:,:,1) = extend_C(a_phisrk)
!
! Calculation at each order i_m, problem at next order: phim at order i_m+1 and phizm at order i_m in powers of eta
DO i_m = 1, M
   !
   phim_ext(1:Nd1,1:Nd2) = 0.0_rp
   !
   DO j = i_m,1,-1
      j2  = i_m - j + 1
      !
      kth_local = build_kth(j2,j)
      ! FIXME : test bug may 2013
      !kth_local(1:Nd1o2p1,1:Nd2) = kth(1:Nd1o2p1,1:Nd2,j)
      !
      ! Direct FFT: j^th z-derivatives of phim at order i_m-j+1 for all j
      temp_C_Nd(1:Nd1o2p1,1:Nd2) = aHOS(1:Nd1o2p1,1:Nd2,j2)
      temp_C_Nd(1:Nd1o2p1,1:Nd2) = kth_local(1:Nd1o2p1,1:Nd2) * temp_C_Nd(1:Nd1o2p1,1:Nd2)
      CALL fourier_2_space_big(temp_C_Nd, phizi)
      !
      ! Construction of phim at order i_m+1 (see the product at order i_m-j+1+j=i_m+1)
	   phim_ext(1:Nd1,1:Nd2)  = phim_ext(1:Nd1,1:Nd2)  - phizi(1:Nd1,1:Nd2) * etapm_ext(1:Nd1,1:Nd2,j+1)
      ! Construction of phizm components at order i_m in eta (cf above)
      phizm(1:Nd1,1:Nd2,i_m) = phizm(1:Nd1,1:Nd2,i_m) + phizi(1:Nd1,1:Nd2) * etapm_ext(1:Nd1,1:Nd2,j)
   END DO
   !
   ! De-aliasing of phizm at order i_m
   IF (i_m > 1) THEN
      CALL dealias(i_m, phizm(:,:,i_m))
   END IF
   ! De-aliasing of phim at order i_m+1
   IF (i_m < M) THEN
      CALL dealias(i_m+1, phim_ext)
      ! Storing it for next i_m
      CALL space_2_fourier_big(phim_ext, temp_C_Nd)
      aHOS(1:Nd1o2p1,1:Nd2,i_m+1) = temp_C_Nd(1:Nd1o2p1,1:Nd2)
   END IF
   ! Storing W1 = phizm at order 1
   IF (i_m == 1) temp2_R_Nd(1:Nd1,1:Nd2) = phizi(1:Nd1,1:Nd2) 
   !
END DO
!
!Assembling of phiz at order M-2, M-1 and M
temp_R_Nd(1:Nd1,1:Nd2) = 0.0_rp
! velocity calculation
phizMm2(1:Nd1,1:Nd2)    = 0.0_rp
phizMm1(1:Nd1,1:Nd2)    = 0.0_rp
DO j=1,M-2
   phizMm2(1:Nd1,1:Nd2) = phizMm2(1:Nd1,1:Nd2) + phizm(1:Nd1,1:Nd2,j)
END DO
DO j=1,M-1
   phizMm1(1:Nd1,1:Nd2) = phizMm1(1:Nd1,1:Nd2) + phizm(1:Nd1,1:Nd2,j)
END DO
temp_R_Nd(1:Nd1,1:Nd2) = phizMm1(1:Nd1,1:Nd2) + phizm(1:Nd1,1:Nd2,M)
! Dealiasing
!  everything has been dealiased step by step previously so it's OK
!
! Filtering if required
CALL choose_filter(ns1,ns2)
phiz = filter_ext(temp_R_Nd, ns1,ns2)
W1   = filter_ext(temp2_R_Nd,ns1,ns2)
!
! Assembling of gradeta2*phiz at order M (i.e. with phiz at order M-2)
!  with correct p dealiasing
geta2phiz(1:Nd1,1:Nd2) = 0.0_rp
DO j=1,M-3
   temp_R_Nd(1:Nd1,1:Nd2) = gradeta2(1:Nd1,1:Nd2) * phizm(1:Nd1,1:Nd2,j)
   CALL dealias(j+2,temp_R_Nd)
   geta2phiz(1:Nd1,1:Nd2) = geta2phiz(1:Nd1,1:Nd2) + temp_R_Nd(1:Nd1,1:Nd2)
END DO
IF (M >= 3) THEN
   temp_R_Nd(1:Nd1,1:Nd2) = gradeta2(1:Nd1,1:Nd2) * phizm(1:Nd1,1:Nd2,MAX(M-2,1))
   CALL dealias(M,temp_R_Nd)
   geta2phiz(1:Nd1,1:Nd2) = geta2phiz(1:Nd1,1:Nd2) + temp_R_Nd(1:Nd1,1:Nd2)
END IF
!
! Assembling of phiz2 at order M and 
!  gradeta2*phiz at order M (i.e. with phiz at order M-2) 
!  with correct p dealiasing
!  in two steps
phiz2(1:Nd1,1:Nd2)      = 0.0_rp
geta2phiz2(1:Nd1,1:Nd2) = 0.0_rp
temp_R_Nd(1:Nd1,1:Nd2)  = 0.0_rp
! velocity calculation
phiz2Mm2(1:Nd1,1:Nd2)   = 0.0_rp
! Step 1: phiz2 up to order M-2 and gradeta2*phiz2 at order M
! squares
DO j = 1, (M-2)/2
   temp2_R_Nd(1:Nd1,1:Nd2) = phizm(1:Nd1,1:Nd2,j) * phizm(1:Nd1,1:Nd2,j)
   CALL dealias(2*j,temp2_R_Nd)
   temp_R_Nd(1:Nd1,1:Nd2)  = temp_R_Nd(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
   temp2_R_Nd(1:Nd1,1:Nd2) =  gradeta2(1:Nd1,1:Nd2) * temp2_R_Nd(1:Nd1,1:Nd2)
   CALL dealias(2*j+2,temp2_R_Nd)
   geta2phiz2(1:Nd1,1:Nd2) = geta2phiz2(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
END DO
! double products
DO j = 2, M-3
   DO q = 1, MIN(M-2-j,j-1)
      temp2_R_Nd(1:Nd1,1:Nd2) = 2.0_rp * phizm(1:Nd1,1:Nd2,j) * phizm(1:Nd1,1:Nd2,q)
      CALL dealias(j+q,temp2_R_Nd)
      temp_R_Nd(1:Nd1,1:Nd2)  = temp_R_Nd(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
      temp2_R_Nd(1:Nd1,1:Nd2) =  gradeta2(1:Nd1,1:Nd2) * temp2_R_Nd(1:Nd1,1:Nd2)
      CALL dealias(j+q+2,temp2_R_Nd)
      geta2phiz2(1:Nd1,1:Nd2) = geta2phiz2(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
   END DO
END DO
! velocity calculation
! FIXME: check this is OK
phiz2Mm2(1:Nd1,1:Nd2) = temp_R_Nd(1:Nd1,1:Nd2)
! Step 2: phiz2 up to order M
! squares
IF (M >= 2) THEN
   temp2_R_Nd(1:Nd1,1:Nd2) = phizm(1:Nd1,1:Nd2,Mo2) * phizm(1:Nd1,1:Nd2,Mo2)
   CALL dealias(2*Mo2,temp2_R_Nd)
   temp_R_Nd(1:Nd1,1:Nd2)  = temp_R_Nd(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
END IF
! double products
DO j = 2, M-1
   DO q = MAX(0,MIN(M-2-j,j-1))+1, MIN(M-j,j-1)
      temp2_R_Nd(1:Nd1,1:Nd2) = 2.0_rp * phizm(1:Nd1,1:Nd2,j) * phizm(1:Nd1,1:Nd2,q)
      CALL dealias(j+q,temp2_R_Nd)
      temp_R_Nd(1:Nd1,1:Nd2)  = temp_R_Nd(1:Nd1,1:Nd2) + temp2_R_Nd(1:Nd1,1:Nd2)
   END DO
END DO
!
phiz2(1:Nd1,1:Nd2) = temp_R_Nd(1:Nd1,1:Nd2)
!
! this looks surprising but the cases M=2, 3 and 4 are correctly treated
! (cf work done in december 2007)
!
CONTAINS
!
FUNCTION build_kth(q,j)
!
IMPLICIT NONE
!
INTEGER :: q, j
REAL(RP), DIMENSION(md1o2p1,md2) :: build_kth
!
INTEGER :: N_der_local, i1, i2 
!
!build_kth = kth(:,:,j)
!FIXME; change feb.2013 OK?
build_kth = 0.d0
build_kth(1:Nd1o2p1,1:Nd2) = kth(1:Nd1o2p1,1:Nd2,j)
!
! in x-direction
! FIXME : check dependence
!Nd_local = 2*n1
!! Nd_local = MIN(q, M-q+1) * n1
!! Nd_local = MIN(Nd_local, FLOOR(4.0*n1))
!! Nd_local = MIN(q, (M+1)/2) * n1
!!
!N_der_local = Nd_local/2 + 1
N_der_local = N_der(1)
!
IF (iseven(N_der_local)) N_der_local = N_der_local - 1
!
build_kth(N_der_local+1:Nd1o2p1,1:Nd2o2p1) = 0.0_rp
!
! in y-direction
IF (n2/=1) THEN
!   Nd_local = 2*n2
!   ! Nd_local = MIN(q, M-q+1) * n2
!   ! Nd_local = MIN(q, (M+1)/2) * n2
!   !
!   N_der_local = Nd_local/2+1
    N_der_local = N_der(2)
   !
   IF (iseven(N_der_local)) N_der_local = N_der_local - 1
   !
   build_kth(1:Nd1o2p1,N_der_local+1:Nd2o2p1) = 0.0_rp
   !
   ! second half of the matrix
   DO i2 = 2, Nd2p1o2
      DO i1 = 1, Nd1o2p1
         build_kth(i1,Nd2-i2+2) = build_kth(i1,i2)
      END DO
   END DO
END IF
!
END FUNCTION build_kth
!
END SUBROUTINE
!
!
!
SUBROUTINE solveHOS_lin(a_phisrk, da_phisrk, a_etark, da_etark, time)
!
IMPLICIT NONE
!% INPUT VARIABLES
REAL(RP) :: time
!       free surface elevation
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN)  :: a_etark
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(OUT) :: da_etark
!       free surface potential
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN)  :: a_phisrk
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(OUT) :: da_phisrk
!% LOCAL VARIABLES
REAL(RP), DIMENSION(m1,m2)                   :: deta, dphis
!
INTEGER :: ns1, ns2, i1, i2
!       de-aliased intermediate products
REAL(RP)                          :: dphi_dn
!       extended FSBCs terms
REAL(RP), DIMENSION(m1,m2) :: eta_temp
!
! Evaluating the elevation
 CALL fourier_2_space(a_etark,eta_temp)
!
! calculation of the powers of eta and
!  the horizontal derivatives of phis and eta
 CALL phisxy_etaxy(a_phisrk, a_etark)
!
! resolution of the different order of the FS potential boundary value problem to get modes
 CALL HOSphis_modes_fully_dealiased(a_phisrk)
!
 CALL choose_filter(ns1,ns2)
! Linear case
IF (M == 1) THEN
      deta(1:n1,1:n2)  = 0.0_rp
      dphis(1:n1,1:n2) = 0.0_rp !- nu(1:n1,1:n2) * phiz(1:n1,1:n2)
ELSE ! Nonlinear case
   ! Evaluation of the RHS of the FSBCs, in physical space
   ! and dealiasing
   !
   !  temp_R_Nd stands for deta_dt-phiz  (on a (p+1)N/2 grid)
   temp_R_Nd(1:Nd1,1:Nd2) = - etax(1:Nd1,1:Nd2) * phisx(1:Nd1,1:Nd2) - etay(1:Nd1,1:Nd2) * phisy(1:Nd1,1:Nd2)
   CALL dealias(2,temp_R_Nd)
   !
   temp_R_Nd(1:Nd1,1:Nd2) = temp_R_Nd(1:Nd1,1:Nd2) + geta2phiz(1:Nd1,1:Nd2)
   !
   !  temp2_R_Nd stands for dphis_dt (on a (p+1)N/2 grid)
   temp2_R_Nd(1:Nd1,1:Nd2) = - phisx(1:Nd1,1:Nd2) * phisx(1:Nd1,1:Nd2) - phisy(1:Nd1,1:Nd2) * phisy(1:Nd1,1:Nd2)
   CALL dealias(2,temp2_R_Nd)
   !
   temp2_R_Nd(1:Nd1,1:Nd2) = temp2_R_Nd(1:Nd1,1:Nd2) + geta2phiz2(1:Nd1,1:Nd2) + phiz2(1:Nd1,1:Nd2)
   !
   ! Filter ?
   ! phiz already filtered in HOS_phi...
   !
   deta     = filter_ext(temp_R_Nd,ns1,ns2)
   !
   dphis    = filter_ext(temp2_R_Nd,ns1,ns2)
   !
   temp_R_n = filter(eta_temp,ns1,ns2)
   !
  IF(abs(Ta).GT.tiny) THEN
      DO i2=1,n2
         DO i1=1,n1
            dphi_dn      = phiz(i1,i2) + deta(i1,i2)
            deta(i1,i2)  = (dphi_dn - W1(i1,i2))* (1.0_rp - exp(-(time/Ta)**n)) 
            ! only the nonlinear part as A is then build as a o(|grad eta|^1)
            dphis(i1,i2) = (0.5_rp * dphis(i1,i2))* (1.0_rp - exp(-(time/Ta)**n)) ! - nu(i1,i2) * dphi_dn)* (1.0_rp - exp(-(time/Ta)**n))
         END DO
      END DO
   ELSE
    DO i2=1,n2
      DO i1=1,n1
         dphi_dn      = phiz(i1,i2) + deta(i1,i2)
         deta(i1,i2)  = dphi_dn - W1(i1,i2)
! only the nonlinear part as A is then build as a o(|grad eta|^1)
         dphis(i1,i2) = 0.5_rp * dphis(i1,i2) ! - nu(i1,i2) * dphi_dn
      END DO
    END DO
   ENDIF
END IF
!
 CALL space_2_fourier(deta,  da_etark)
 CALL space_2_fourier(dphis, da_phisrk)
!
da_phisrk(1,1) = da_phisrk(1,1) - a_etark(1,1)
!
END SUBROUTINE solveHOS_lin
!
!
!
END MODULE resol_HOS

