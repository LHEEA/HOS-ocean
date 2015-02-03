MODULE velocities
!
! This module computes the modal amplitudes of velocity components + dphi/dt + ...
! This uses adaptation of H2-operator used in DNO to HOS formalism
! Or a direct inversion of the full matrices may be done
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
USE resol_HOS
USE fourier_r2c
USE RF_solution

CONTAINS
!
! start HOSvel2_bis *************************************************
!
!     ======================================
!      SUBROUTINE HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current,Ta,n,iCPUtime)
        SUBROUTINE HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current)
!     ======================================

IMPLICIT NONE
!
!% INPUT VARIABLES
!   Number of levels for eta2
INTEGER, INTENT(IN)  :: meta2, M_HOSvel
REAL(RP), INTENT(IN) :: time_current
!
COMPLEX(CP), DIMENSION(m1o2p1,m2) , INTENT(IN) :: a_eta_l, a_phis_l
!% LOCAL VARIABLES
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_deta_l, a_dphis_l, as2
REAL(RP), DIMENSION(m1,m2)        :: eta3, deta3, dphis3, eta2

REAL(RP), DIMENSION(md1,md2,M_HOSvel+1) :: etapmext1, etapmext2
REAL(RP), DIMENSION(md1o2p1,md2o2p1,M_HOSvel) :: kth_vel, kth2_vel
REAL(RP) :: k2
REAL(RP), DIMENSION(M_HOSvel) :: oneoj_vel
REAL(RP), DIMENSION(md1,md2) :: etax2, etay2
REAL(RP), DIMENSION(md1,md2) :: phimx_ext2, phimy_ext2, phimz_ext2

REAL(RP), DIMENSION(m1,m2) :: detaphiz_d, phimxetax, phiz_etax, phiz_etay, phisx_d, phisy_d

REAL(RP), DIMENSION(m1,m2)        :: phim, phimx, phimy, phimz, phit, phimx2, phimy2, phimz2
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_phim, a_phimx, a_phimy, a_phimz, a_phit
!
REAL(RP), DIMENSION(m1,m2) :: Press_SL, CCSL
!
REAL(RP), DIMENSION(md1,md2)        :: phi_int, vitx_int, vity_int, vitz_int, phit_int
COMPLEX(CP), DIMENSION(md1o2p1,md2) :: a_phi_int, a_vitx_int, a_vity_int, a_vitz_int, a_phit_int
!
REAL(RP), DIMENSION(md1,md2) :: phim_ext, phimx_ext, phimy_ext, phimz_ext, phit_ext
REAL(RP), DIMENSION(md1,md2) :: phim2_ext, phimx2_ext, phimy2_ext, phimz2_ext, phit2_ext
!
COMPLEX(CP), DIMENSION(md1o2p1,md2,M_HOSvel) :: aHOS, aHOSx, aHOSy, aHOSz, aHOSt

REAL(RP), DIMENSION(md1,md2) :: phizi, phizix, phiziy, phiziz, phizit
COMPLEX(CP), DIMENSION(md1o2p1,md2) :: a_phizi, a_phizix, a_phiziy, a_phiziz, a_phizit
!
INTEGER :: i1, j, iHOS, jm1, j2, i2, ii1
!
REAL(RP), DIMENSION(m1,m2) :: Press_bis, CCSL_bis, etax_tmp, etay_tmp
!
! One uses mode as input
!
CALL fourier_2_space(a_eta_l,  eta3)
!
! Compute da_eta
!
CALL solveHOS_lin(a_phis_l, a_dphis_l, a_eta_l, a_deta_l, time_current)
!
! FIXME: check where it comes from in resol_HOS
a_dphis_l(1,1) = a_dphis_l(1,1)+a_eta_l(1,1)
!
CALL fourier_2_space(a_deta_l,  deta3)
CALL fourier_2_space(a_dphis_l,  dphis3)
!
DO i2=1,n2
    DO i1=1,n1
        deta3(i1,i2)  = (deta3(i1,i2) + W1(i1,i2))
        dphis3(i1,i2) = (dphis3(i1,i2) - eta3(i1,i2))
    ENDDO
ENDDO
!
! Definition of new variables needed on M_HOSvel
!
kth_vel(1,1,1:M_HOSvel) = 0.0_rp
kth2_vel(1,1,1:M_HOSvel) = 0.0_rp
!
DO i1=2, Nd1o2p1
    k2           = kx(i1) * kx(i1)
    kth_vel(i1,1,1)  = SQRT(k2)*TANH(SQRT(k2)*depth_star)
    kth2_vel(i1,1,1)  = SQRT(k2)/TANH(SQRT(k2)*depth_star)
    IF (M_HOSvel > 1) THEN
        kth_vel(i1,1,2)  = k2
        kth2_vel(i1,1,2) = k2
    ENDIF
    DO j = 3,M_HOSvel
        kth_vel(i1,1,j) = k2 * kth_vel(i1,1,j-2)
        kth2_vel(i1,1,j) = k2 * kth2_vel(i1,1,j-2)
    ENDDO
ENDDO
!
DO i2 = 2, Nd2o2p1
    k2          = ky(i2) * ky(i2)
    kth_vel(1,i2,1) = SQRT(k2)*TANH(SQRT(k2)*depth_star)
    kth2_vel(1,i2,1) = SQRT(k2)/TANH(SQRT(k2)*depth_star)
    IF (M_HOSvel > 1) THEN
        kth_vel(1,i2,2) = k2*TANH(SQRT(k2)*depth_star)
        kth2_vel(1,i2,2) = k2/TANH(SQRT(k2)*depth_star)
    ENDIF
    DO j = 3,M_HOSvel
        kth_vel(1,i2,j)  = k2 * kth_vel(1,i2,j-2)
        kth2_vel(1,i2,j) = k2 * kth2_vel(1,i2,j-2)
    ENDDO
    DO i1 = 2, Nd1o2p1
        k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
        kth_vel(i1,i2,1) = SQRT(k2)*TANH(SQRT(k2)*depth_star)
        kth2_vel(i1,i2,1)= SQRT(k2)/TANH(SQRT(k2)*depth_star)
        IF (M_HOSvel > 1) THEN
            kth_vel(i1,i2,2)  = k2
            kth2_vel(i1,i2,2) = k2
        ENDIF
        DO j = 3,M_HOSvel
            kth_vel(i1,i2,j)  = k2 * kth_vel(i1,i2,j-2)
            kth2_vel(i1,i2,j) = k2 * kth2_vel(i1,i2,j-2)
        ENDDO
    ENDDO
ENDDO
!
! the powers of eta (eta^m/m!) (on a 2N grid)
! m=1
DO j=1,M_HOSvel
    oneoj_vel(j) = 1.0_rp / REAL(j, RP)
ENDDO
!
! the powers of eta (eta^m/m!) (on a 2N grid)
etapmext1(:,:,1) = 1.0_rp
!
! m=1
temp_C_Nd = extend_C(a_eta_l) !FIXME: temp_C_Nd global variable
!
CALL fourier_2_space_big(temp_C_Nd, etapmext1(:,:,1+1))
!
! m=2 to M
DO j=2,M_HOSvel
    etapmext1(1:Nd1,1:Nd2,j+1) = etapmext1(1:Nd1,1:Nd2,j-1+1) * oneoj_vel(j) * etapmext1(1:Nd1,1:Nd2,1+1)
    CALL dealias(j, etapmext1(:,:,j+1))
ENDDO
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        etax2(i1,i2) = etax(i1,i2) * etax(i1,i2)
        etay2(i1,i2) = etay(i1,i2) * etay(i1,i2)
    ENDDO
ENDDO

CALL dealias(2,etax2)
CALL dealias(2,etay2)
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        etapmext2(i1,i2,1:M_HOSvel+1)=etapmext1(i1,i2,1:M_HOSvel+1)
    ENDDO
ENDDO
!
! initializations
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        phimx_ext(i1,i2)= (phisx(i1,i2) - etax(i1,i2) * phizMm1(i1,i2))
        phimy_ext(i1,i2)= (phisy(i1,i2) - etay(i1,i2) * phizMm1(i1,i2))
        phimx_ext2(i1,i2)= (phisx(i1,i2)*phisx(i1,i2) + etax2(i1,i2)*phiz2Mm2(i1,i2) &
           - 2.0_rp*phisx(i1,i2)*etax(i1,i2) * phizMm1(i1,i2))
        phimy_ext2(i1,i2)= (phisy(i1,i2)*phisy(i1,i2) + etay2(i1,i2)*phiz2Mm2(i1,i2) &
           - 2.0_rp*phisy(i1,i2)*etay(i1,i2) * phizMm1(i1,i2))
        phimz_ext2(i1,i2) = phiz2(i1,i2)
        temp_R_Nd(i1,i2) =   (etax2(i1,i2)+ etay2(i1,i2)) * phiz2Mm2(i1,i2) + phiz2(i1,i2) &
           - phisx(i1,i2) * etax(i1,i2) * phizMm1(i1,i2) - phisy(i1,i2) * etay(i1,i2) * phizMm1(i1,i2) !detaphiz(i1,i2)
        temp2_R_Nd(i1,i2) = ((etax(i1,i2)*etax(i1,i2) + etay(i1,i2)*etay(i1,i2)) * phizMm2(i1,i2) &
           - phisx(i1,i2) * etax(i1,i2) - phisy(i1,i2) * etay(i1,i2)) !phimxetax_ext(i1,i2)
    ENDDO
ENDDO
!
detaphiz_d = filter_ext(temp_R_Nd,n1o2p1,n2o2p1)
phimxetax  = filter_ext(temp2_R_Nd,n1o2p1,n2o2p1)
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        temp_R_Nd(i1,i2) = etax(i1,i2) * phizMm1(i1,i2)  !phiz_etax_ext(i1,i2)
        temp2_R_Nd(i1,i2) = etay(i1,i2) * phizMm1(i1,i2) !phiz_etay_ext(i1,i2)
    ENDDO
ENDDO
!
phiz_etax = filter_ext(temp_R_Nd,n1o2p1,n2o2p1)
phiz_etay = filter_ext(temp2_R_Nd,n1o2p1,n2o2p1)
phisx_d   = filter_ext(phisx,n1o2p1,n2o2p1)
phisy_d   = filter_ext(phisy,n1o2p1,n2o2p1)
!
DO i1=1,n1
    DO i2=1,n2
        phimx(i1,i2) = (phisx_d(i1,i2) - phiz_etax(i1,i2))
        phimy(i1,i2) = (phisy_d(i1,i2) - phiz_etay(i1,i2))
        phit(i1,i2)  = dphis3(i1,i2) - (detaphiz_d(i1,i2))
        phimz(i1,i2) = phiz(i1,i2)
    ENDDO
ENDDO
!
phimx  = filter_ext(phimx_ext,n1o2p1,n2o2p1)
phimy  = filter_ext(phimy_ext,n1o2p1,n2o2p1)
phimx2 = filter_ext(phimx_ext2,n1o2p1,n2o2p1)
phimy2 = filter_ext(phimy_ext2,n1o2p1,n2o2p1)
phimz2 = filter_ext(phimz_ext2,n1o2p1,n2o2p1)
!
etax_tmp = filter_ext(etax,n1o2p1,n2o2p1)
etay_tmp = filter_ext(etay,n1o2p1,n2o2p1)
!
DO i1=1,n1
    DO i2=1,n2
        phiref_SL(i1,i2)  = phim(i1,i2)
        vitxref_SL(i1,i2) = phimx(i1,i2)
        vityref_SL(i1,i2) = phimy(i1,i2)
        vitzref_SL(i1,i2) = phimz(i1,i2)
        dphitref_SL(i1,i2) = phit(i1,i2)
        ! Comment to have exact CDSL... FIXME ?
        !phimx2(i1,i2) = phimx(i1,i2)*phimx(i1,i2)
        !phimy2(i1,i2) = phimy(i1,i2)*phimy(i1,i2)
        !phimz2(i1,i2) = phimz(i1,i2)*phimz(i1,i2)
        !
        vitx2ref_SL(i1,i2) = phimx2(i1,i2)
        vity2ref_SL(i1,i2) = phimy2(i1,i2)
        vitz2ref_SL(i1,i2) = phimz2(i1,i2)

        Press_SL(i1,i2) = - eta3(i1,i2) - phit(i1,i2) +(- 0.5_rp*(phimx2(i1,i2) + phimy2(i1,i2) + phimz2(i1,i2)))    ! DFSBC OK if phimx2.NE.phimx*phimx
        CCSL(i1,i2) = deta3(i1,i2) - phimxetax(i1,i2) - phimz(i1,i2)       ! CCSL OK
        ! Approximation (no dealiasing)
        Press_bis(i1,i2) = - eta3(i1,i2) - phit(i1,i2) - 0.5_rp*(phimx(i1,i2)**2 + phimy(i1,i2)**2 + phimz(i1,i2)**2)
        CCSL_bis(i1,i2)  = deta3(i1,i2) + (phimx(i1,i2)*etax_tmp(i1,i2)+phimy(i1,i2)*etay_tmp(i1,i2)) - phimz(i1,i2) ! KFSBC OK
    ENDDO
ENDDO
!
WRITE(*,*) 'KFSBC/DFSBC =', MAX(MAXVAL(ABS(CCSL(1:n1,1:n2))),MAXVAL(ABS(Press_SL(1:n1,1:n2)))), 'approx',&
    MAX(MAXVAL(ABS(CCSL_bis(1:n1,1:n2))),MAXVAL(ABS(Press_bis(1:n1,1:n2))))
!
CALL space_2_fourier(phim,a_phim)
CALL space_2_fourier(phimx,a_phimx)
CALL space_2_fourier(phimy,a_phimy)
CALL space_2_fourier(phimz,a_phimz)
CALL space_2_fourier(phit,a_phit)

a_phi_int = extend_C(a_phim)
a_vitx_int = extend_C(a_phimx)
a_vity_int = extend_C(a_phimy)
a_vitz_int = extend_C(a_phimz)
a_phit_int = extend_C(a_phit)

CALL fourier_2_space_big(a_phi_int,phi_int)
CALL fourier_2_space_big(a_vitx_int,vitx_int)
CALL fourier_2_space_big(a_vity_int,vity_int)
CALL fourier_2_space_big(a_vitz_int,vitz_int)
CALL fourier_2_space_big(a_phit_int,phit_int)
!
phim_ext(:,:)  = phi_int(:,:)
phimx_ext(:,:) = vitx_int(:,:)
phimy_ext(:,:) = vity_int(:,:)
phimz_ext(:,:) = vitz_int(:,:)
phit_ext(:,:)  = phit_int(:,:)

DO ii1=1,meta2-1

    DO i1 = 1,Nd1
        DO i2 = 1,Nd2
            etapmext1(i1,i2,:)=etapmext2(i1,i2,:)
        ENDDO
    ENDDO

    DO i1=1,n1
        DO i2=1,n2
            eta2(i1,i2)=(eta3(i1,i2))*(meta2-1-(ii1-1))/(meta2)
        ENDDO
    ENDDO
    !
    ! inverse FFT: determination of the modes of eta2
    !
    CALL space_2_fourier(eta2,as2)

    etapmext2(:,:,1) = 1.0_rp
    !
    ! m=1
    temp_C_Nd = extend_C(as2) !FIXME: temp_C_Nd global variable
    !
    CALL fourier_2_space_big(temp_C_Nd, etapmext2(:,:,1+1))
    !
    ! m=2 to M
    DO j=2,M_HOSvel
        etapmext2(1:Nd1,1:Nd2,j+1) = etapmext2(1:Nd1,1:Nd2,j-1+1) * oneoj_vel(j) * etapmext2(1:Nd1,1:Nd2,1+1)
        CALL dealias(j, etapmext2(:,:,j+1))
    ENDDO
    !
    ! calculation at each order iHOS
    DO iHOS=1,M_HOSvel
        ! de-aliasing of aHOS(:,:,iHOS)
        CALL dealias(1,phim_ext)
        CALL dealias(1,phimx_ext)
        CALL dealias(1,phimy_ext)
        CALL dealias(1,phimz_ext)
        CALL dealias(1,phit_ext)
        !
        CALL space_2_fourier_big(phim_ext,aHOS(:,:,iHOS))
        CALL space_2_fourier_big(phimx_ext,aHOSx(:,:,iHOS))
        CALL space_2_fourier_big(phimy_ext,aHOSy(:,:,iHOS))
        CALL space_2_fourier_big(phimz_ext,aHOSz(:,:,iHOS))
        CALL space_2_fourier_big(phit_ext,aHOSt(:,:,iHOS))

        DO i1 = 1, Nd1
            DO i2=1, Nd2
                phim_ext(i1,i2)  = 0.0_rp
                phimx_ext(i1,i2)  = 0.0_rp
                phimy_ext(i1,i2)  = 0.0_rp
                phimz_ext(i1,i2)  = 0.0_rp
                phit_ext(i1,i2)   = 0.0_rp
                !
                phim2_ext(i1,i2) = 0.0_rp
                phimx2_ext(i1,i2) = 0.0_rp
                phimy2_ext(i1,i2) = 0.0_rp
                phimz2_ext(i1,i2) = 0.0_rp
                phit2_ext(i1,i2)  = 0.0_rp
            ENDDO
        ENDDO

        DO j = iHOS,2,-1
            jm1 = j - 1
            j2  = iHOS - jm1
            !
            ! FIXME : optimize
            a_phizi  = ((0.0_rp, 0.0_rp))
            a_phizix = ((0.0_rp, 0.0_rp))
            a_phiziy = ((0.0_rp, 0.0_rp))
            a_phiziz = ((0.0_rp, 0.0_rp))
            a_phizit = ((0.0_rp, 0.0_rp))
            DO i2 = 1, MIN(N_der(2),Nd2o2p1)
                DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                    a_phizi(i1,i2)  = aHOS(i1,i2,j2)  * kth_vel(i1,i2,j)
                    a_phizix(i1,i2) = aHOSx(i1,i2,j2) * kth_vel(i1,i2,j)
                    a_phiziy(i1,i2) = aHOSy(i1,i2,j2) * kth_vel(i1,i2,j)
                    a_phiziz(i1,i2) = aHOSz(i1,i2,j2) * kth2_vel(i1,i2,j)
                    a_phizit(i1,i2) = aHOSt(i1,i2,j2) * kth_vel(i1,i2,j)
                ENDDO
            ENDDO
            DO i2 = 2, MIN(N_der(2),Nd2o2p1)
                DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                    a_phizi(i1,Nd2-i2+2)  = aHOS(i1,Nd2-i2+2,j2)  * kth_vel(i1,i2,j)
                    a_phizix(i1,Nd2-i2+2) = aHOSx(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
                    a_phiziy(i1,Nd2-i2+2) = aHOSy(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
                    a_phiziz(i1,Nd2-i2+2) = aHOSz(i1,Nd2-i2+2,j2) * kth2_vel(i1,i2,j)
                    a_phizit(i1,Nd2-i2+2) = aHOSt(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
                ENDDO
            ENDDO
            !
            CALL fourier_2_space_big(a_phizi,phizi)
            CALL fourier_2_space_big(a_phizix,phizix)
            CALL fourier_2_space_big(a_phiziy,phiziy)
            CALL fourier_2_space_big(a_phiziz,phiziz)
            CALL fourier_2_space_big(a_phizit,phizit)

            ! change : new etapmext goes from 1 to mHOS+1
            DO i1=1,Nd1
                DO i2=1,Nd2
                    phim_ext(i1,i2)  = phim_ext(i1,i2) - phizi(i1,i2) * etapmext1(i1,i2,j+1)
                    phimx_ext(i1,i2) = phimx_ext(i1,i2) - phizix(i1,i2) * etapmext1(i1,i2,j+1)
                    phimy_ext(i1,i2) = phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext1(i1,i2,j+1)
                    phimz_ext(i1,i2) = phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext1(i1,i2,j+1)
                    phit_ext(i1,i2)  = phit_ext(i1,i2) - phizit(i1,i2) * etapmext1(i1,i2,j+1)
                    !
                    phim2_ext(i1,i2) = phim2_ext(i1,i2) + phizi(i1,i2) * etapmext2(i1,i2,j+1)
                    phimx2_ext(i1,i2) = phimx2_ext(i1,i2) + phizix(i1,i2) * etapmext2(i1,i2,j+1)
                    phimy2_ext(i1,i2) = phimy2_ext(i1,i2) + phiziy(i1,i2) * etapmext2(i1,i2,j+1)
                    phimz2_ext(i1,i2) = phimz2_ext(i1,i2) + phiziz(i1,i2) * etapmext2(i1,i2,j+1)
                    phit2_ext(i1,i2) = phit2_ext(i1,i2) + phizit(i1,i2) * etapmext2(i1,i2,j+1)
                ENDDO
            ENDDO
            ! dealias(j+1) FIXME: test 09/2014
        ENDDO

        !
        ! FIXME: optimize
        a_phizi  = ((0.0_rp, 0.0_rp))
        a_phizix = ((0.0_rp, 0.0_rp))
        a_phiziy = ((0.0_rp, 0.0_rp))
        a_phiziz = ((0.0_rp, 0.0_rp))
        a_phizit = ((0.0_rp, 0.0_rp))
        DO i2 = 1, MIN(N_der(2),Nd2o2p1)
            DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                a_phizi(i1,i2)  = aHOS(i1,i2,iHOS)  * kth_vel(i1,i2,1)
                a_phizix(i1,i2) = aHOSx(i1,i2,iHOS) * kth_vel(i1,i2,1)
                a_phiziy(i1,i2) = aHOSy(i1,i2,iHOS) * kth_vel(i1,i2,1)
                a_phiziz(i1,i2) = aHOSz(i1,i2,iHOS) * kth2_vel(i1,i2,1)
                a_phizit(i1,i2) = aHOSt(i1,i2,iHOS) * kth_vel(i1,i2,1)
            ENDDO
        ENDDO
        DO i2 = 2, MIN(N_der(2),Nd2o2p1)
            DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                a_phizi(i1,Nd2-i2+2)  = aHOS(i1,Nd2-i2+2,iHOS)  * kth_vel(i1,i2,1)
                a_phizix(i1,Nd2-i2+2) = aHOSx(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
                a_phiziy(i1,Nd2-i2+2) = aHOSy(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
                a_phiziz(i1,Nd2-i2+2) = aHOSz(i1,Nd2-i2+2,iHOS) * kth2_vel(i1,i2,1)
                a_phizit(i1,Nd2-i2+2) = aHOSt(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
            ENDDO
        ENDDO

        CALL fourier_2_space_big(a_phizi, phizi)
        CALL fourier_2_space_big(a_phizix, phizix)
        CALL fourier_2_space_big(a_phiziy, phiziy)
        CALL fourier_2_space_big(a_phiziz, phiziz)
        CALL fourier_2_space_big(a_phizit, phizit)

        ! change : new etapmext goes from 1 to mHOS+1
        DO i1=1,Nd1
            DO i2=1,Nd2
                phim_ext(i1,i2)   = phim_ext(i1,i2) - phizi(i1,i2) * etapmext1(i1,i2,1+1)
                phimx_ext(i1,i2)  = phimx_ext(i1,i2) - phizix(i1,i2) * etapmext1(i1,i2,1+1)
                phimy_ext(i1,i2)  = phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext1(i1,i2,1+1)
                phimz_ext(i1,i2)  = phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext1(i1,i2,1+1)
                phit_ext(i1,i2)   = phit_ext(i1,i2) - phizit(i1,i2) * etapmext1(i1,i2,1+1)
                !
                phim2_ext(i1,i2)  = phim2_ext(i1,i2) + phizi(i1,i2) * etapmext2(i1,i2,1+1)
                phimx2_ext(i1,i2) = phimx2_ext(i1,i2) + phizix(i1,i2) * etapmext2(i1,i2,1+1)
                phimy2_ext(i1,i2) = phimy2_ext(i1,i2) + phiziy(i1,i2) * etapmext2(i1,i2,1+1)
                phimz2_ext(i1,i2) = phimz2_ext(i1,i2) + phiziz(i1,i2) * etapmext2(i1,i2,1+1)
                phit2_ext(i1,i2)  = phit2_ext(i1,i2) + phizit(i1,i2) * etapmext2(i1,i2,1+1)
                !   construction of velocities
                phi_int(i1,i2)  = phi_int(i1,i2)  + phim_ext(i1,i2)  + phim2_ext(i1,i2)
                vitx_int(i1,i2) = vitx_int(i1,i2) + phimx_ext(i1,i2) + phimx2_ext(i1,i2)
                vity_int(i1,i2) = vity_int(i1,i2) + phimy_ext(i1,i2) + phimy2_ext(i1,i2)
                vitz_int(i1,i2) = vitz_int(i1,i2) + phimz_ext(i1,i2) + phimz2_ext(i1,i2)
                phit_int(i1,i2) = phit_int(i1,i2) + phit_ext(i1,i2) + phit2_ext(i1,i2)
            ENDDO
        ENDDO
        ! dealias(2) !FIXME: added 09/2014
    ENDDO

    DO i1=1,Nd1
       DO i2=1,Nd2
            phim_ext(i1,i2)  = phi_int(i1,i2)
            phimx_ext(i1,i2) = vitx_int(i1,i2)
            phimy_ext(i1,i2) = vity_int(i1,i2)
            phimz_ext(i1,i2) = vitz_int(i1,i2)
            phit_ext(i1,i2)  = phit_int(i1,i2)
        ENDDO
    ENDDO
ENDDO

DO i1=1,Nd1
    DO i2=1,Nd2
        etapmext1(i1,i2,1:M_HOSvel+1)=etapmext2(i1,i2,1:M_HOSvel+1)
    ENDDO
ENDDO

CALL HOSvel(M_HOSvel,phi_int,vitx_int,vity_int,vitz_int,phit_int,etapmext1, kth_vel, kth2_vel)


END SUBROUTINE
!
! end HOSvel2_bis ***************************************************
!
! start HOSvel *************************************************
!
!     ======================================
        SUBROUTINE HOSvel(M_HOSvel,phi_int,vitx_int,vity_int,vitz_int,phit_int,etapmext1, kth_vel, kth2_vel)
!     ======================================

IMPLICIT NONE

!% LOCAL VARIABLES
!   order m of phi and i^th derivatives at order m
REAL(RP), DIMENSION(md1,md2), INTENT(IN) :: phi_int,vitx_int,vity_int,vitz_int,phit_int
REAL(RP), DIMENSION(md1,md2,M_HOSvel+1), INTENT(IN) :: etapmext1
REAL(RP), DIMENSION(md1o2p1,md2o2p1,M_HOSvel), INTENT(IN) :: kth_vel, kth2_vel

REAL(RP), DIMENSION(md1,md2) :: phi_ext, vitx_ext, vity_ext, vitz_ext, dphit_ext
COMPLEX(CP), DIMENSION(md1o2p1,md2,M_HOSvel) :: aHOS, aHOSx, aHOSy, aHOSz, aHOSt
REAL(RP), DIMENSION(md1,md2) :: phizi, phizix, phiziy, phiziz, phizit

COMPLEX(CP), DIMENSION(md1o2p1,md2) :: a_phizi, a_phizix, a_phiziy, a_phiziz, a_phizit

REAL(RP), DIMENSION(md1,md2) :: phim_ext, phimx_ext, phimy_ext, phimz_ext, phit_ext
REAL(RP), DIMENSION(m1,m2) :: B, Bx, By, Bz, Bt

INTEGER :: i1, j, iHOS, jm1, j2, i2, M_HOSvel
!
! initializations
!
DO i1=1,Nd1
    DO i2=1,Nd2
        phi_ext(i1,i2) = phi_int(i1,i2)
        vitx_ext(i1,i2)= vitx_int(i1,i2)
        vity_ext(i1,i2)= vity_int(i1,i2)
        vitz_ext(i1,i2)= vitz_int(i1,i2)
        dphit_ext(i1,i2)= phit_int(i1,i2)
        !
        phim_ext(i1,i2) = phi_int(i1,i2)
        phimx_ext(i1,i2)= vitx_int(i1,i2)
        phimy_ext(i1,i2)= vity_int(i1,i2)
        phimz_ext(i1,i2)= vitz_int(i1,i2)
        phit_ext(i1,i2)= phit_int(i1,i2)
    ENDDO
ENDDO
!
! calculation at each order iHOS
!
DO iHOS=1,M_HOSvel

!   de-aliasing of aHOS(:,:,iHOS)
! FIXME: possible to optimize
  CALL dealias(1,phim_ext)
  CALL dealias(1,phimx_ext)
  CALL dealias(1,phimy_ext)
  CALL dealias(1,phimz_ext)
  CALL dealias(1,phit_ext)
  !
  CALL space_2_fourier_big(phim_ext,aHOS(:,:,iHOS))
  CALL space_2_fourier_big(phimx_ext,aHOSx(:,:,iHOS))
  CALL space_2_fourier_big(phimy_ext,aHOSy(:,:,iHOS))
  CALL space_2_fourier_big(phimz_ext,aHOSz(:,:,iHOS))
  CALL space_2_fourier_big(phit_ext,aHOSt(:,:,iHOS))

    DO i1 = 1, Nd1
        DO i2=1,Nd2
            phim_ext(i1,i2)  = 0.0_rp
            phimx_ext(i1,i2) = 0.0_rp
            phimy_ext(i1,i2) = 0.0_rp
            phimz_ext(i1,i2) = 0.0_rp
            phit_ext(i1,i2)  = 0.0_rp
        ENDDO
    ENDDO

    DO j = iHOS,2,-1
        jm1 = j - 1
        j2  = iHOS - jm1
        ! FIXME: optimize
        a_phizi  = ((0.0_rp, 0.0_rp))
        a_phizix = ((0.0_rp, 0.0_rp))
        a_phiziy = ((0.0_rp, 0.0_rp))
        a_phiziz = ((0.0_rp, 0.0_rp))
        a_phizit = ((0.0_rp, 0.0_rp))
        DO i2 = 1, MIN(N_der(2),Nd2o2p1)
            DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                a_phizi(i1,i2)  = aHOS(i1,i2,j2)  * kth_vel(i1,i2,j)
                a_phizix(i1,i2) = aHOSx(i1,i2,j2) * kth_vel(i1,i2,j)
                a_phiziy(i1,i2) = aHOSy(i1,i2,j2) * kth_vel(i1,i2,j)
                a_phiziz(i1,i2) = aHOSz(i1,i2,j2) * kth2_vel(i1,i2,j)
                a_phizit(i1,i2) = aHOSt(i1,i2,j2) * kth_vel(i1,i2,j)
            ENDDO
        ENDDO
        DO i2 = 2, MIN(N_der(2),Nd2o2p1)
            DO i1 = 1, MIN(N_der(1),Nd1o2p1)
                a_phizi(i1,Nd2-i2+2)  = aHOS(i1,Nd2-i2+2,j2)  * kth_vel(i1,i2,j)
                a_phizix(i1,Nd2-i2+2) = aHOSx(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
                a_phiziy(i1,Nd2-i2+2) = aHOSy(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
                a_phiziz(i1,Nd2-i2+2) = aHOSz(i1,Nd2-i2+2,j2) * kth2_vel(i1,i2,j)
                a_phizit(i1,Nd2-i2+2) = aHOSt(i1,Nd2-i2+2,j2) * kth_vel(i1,i2,j)
            ENDDO
        ENDDO
        !
        CALL fourier_2_space_big(a_phizi,phizi)
        CALL fourier_2_space_big(a_phizix,phizix)
        CALL fourier_2_space_big(a_phiziy,phiziy)
        CALL fourier_2_space_big(a_phiziz,phiziz)
        CALL fourier_2_space_big(a_phizit,phizit)
        ! change : new etapmext goes from 1 to mHOS+1
        DO i1=1,Nd1
            DO i2=1,Nd2
                phim_ext(i1,i2)  = phim_ext(i1,i2)  - phizi(i1,i2) * etapmext1(i1,i2,j+1)
                phimx_ext(i1,i2) = phimx_ext(i1,i2) - phizix(i1,i2) * etapmext1(i1,i2,j+1)
                phimy_ext(i1,i2) = phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext1(i1,i2,j+1)
                phimz_ext(i1,i2) = phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext1(i1,i2,j+1)
                phit_ext(i1,i2)  = phit_ext(i1,i2)  - phizit(i1,i2) * etapmext1(i1,i2,j+1)
            ENDDO
        ENDDO
         ! dealias(j+1) FIXME: test 09/2014
    ENDDO
    ! FIXME: optimize
    a_phizi(:,:)  = ((0.0_rp, 0.0_rp))
    a_phizix(:,:) = ((0.0_rp, 0.0_rp))
    a_phiziy(:,:) = ((0.0_rp, 0.0_rp))
    a_phiziz(:,:) = ((0.0_rp, 0.0_rp))
    a_phizit(:,:) = ((0.0_rp, 0.0_rp))
    !
    DO i2 = 1, MIN(N_der(2),Nd2o2p1)
        DO i1 = 1, MIN(N_der(1),Nd1o2p1)
            a_phizi(i1,i2)  = aHOS(i1,i2,iHOS)  * kth_vel(i1,i2,1)
            a_phizix(i1,i2) = aHOSx(i1,i2,iHOS) * kth_vel(i1,i2,1)
            a_phiziy(i1,i2) = aHOSy(i1,i2,iHOS) * kth_vel(i1,i2,1)
            a_phiziz(i1,i2) = aHOSz(i1,i2,iHOS) * kth2_vel(i1,i2,1)
            a_phizit(i1,i2) = aHOSt(i1,i2,iHOS) * kth_vel(i1,i2,1)
        ENDDO
    ENDDO
    DO i2 = 2, MIN(N_der(2),Nd2o2p1)
        DO i1 = 1, MIN(N_der(1),Nd1o2p1)
            a_phizi(i1,Nd2-i2+2)  = aHOS(i1,Nd2-i2+2,iHOS)  * kth_vel(i1,i2,1)
            a_phizix(i1,Nd2-i2+2) = aHOSx(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
            a_phiziy(i1,Nd2-i2+2) = aHOSy(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
            a_phiziz(i1,Nd2-i2+2) = aHOSz(i1,Nd2-i2+2,iHOS) * kth2_vel(i1,i2,1)
            a_phizit(i1,Nd2-i2+2) = aHOSt(i1,Nd2-i2+2,iHOS) * kth_vel(i1,i2,1)
        ENDDO
    ENDDO

    CALL fourier_2_space_big(a_phizi, phizi)
    CALL fourier_2_space_big(a_phizix, phizix)
    CALL fourier_2_space_big(a_phiziy, phiziy)
    CALL fourier_2_space_big(a_phiziz, phiziz)
    CALL fourier_2_space_big(a_phizit, phizit)
    ! change : new etapmext goes from 1 to mHOS+1
    DO i1=1,Nd1
        DO i2=1,Nd2
            phim_ext(i1,i2) = phim_ext(i1,i2)  - phizi(i1,i2)  * etapmext1(i1,i2,1+1)
            phimx_ext(i1,i2)= phimx_ext(i1,i2) - phizix(i1,i2) * etapmext1(i1,i2,1+1)
            phimy_ext(i1,i2)= phimy_ext(i1,i2) - phiziy(i1,i2) * etapmext1(i1,i2,1+1)
            phimz_ext(i1,i2)= phimz_ext(i1,i2) - phiziz(i1,i2) * etapmext1(i1,i2,1+1)
            phit_ext(i1,i2) = phit_ext(i1,i2)  - phizit(i1,i2) * etapmext1(i1,i2,1+1)
            !   construction of velocities / potential
            phi_ext(i1,i2)   = phi_ext(i1,i2)   + phim_ext(i1,i2)
            vitx_ext(i1,i2)  = vitx_ext(i1,i2)  + phimx_ext(i1,i2)
            vity_ext(i1,i2)  = vity_ext(i1,i2)  + phimy_ext(i1,i2)
            vitz_ext(i1,i2)  = vitz_ext(i1,i2)  + phimz_ext(i1,i2)
            dphit_ext(i1,i2) = dphit_ext(i1,i2) + phit_ext(i1,i2)
        ENDDO
    ENDDO
ENDDO

B  = filter_ext(phi_ext,n1o2p1,n2o2p1)
Bx = filter_ext(vitx_ext,n1o2p1,n2o2p1)
By = filter_ext(vity_ext,n1o2p1,n2o2p1)
Bz = filter_ext(vitz_ext,n1o2p1,n2o2p1)
Bt = filter_ext(dphit_ext,n1o2p1,n2o2p1)

CALL space_2_fourier(B,modesspec)
CALL space_2_fourier(Bx,modesspecx)
CALL space_2_fourier(By,modesspecy)
CALL space_2_fourier(Bz,modesspecz)
CALL space_2_fourier(Bt,modesspect)

END SUBROUTINE
!
! end HOSvel ***************************************************
!
!     ============================================================
SUBROUTINE reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect,modesSL, modesSLt, error)
!     ============================================================

IMPLICIT NONE
!% INPUT VARIABLES
COMPLEX(CP), DIMENSION(m1o2p1,m2), INTENT(IN) :: modesspec,modesspecx,modesspecy,modesspecz,modesspect,modesSL,modesSLt
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: phi_l, vitx_l, vity_l, vitz_l, dphit_l
REAL(RP), DIMENSION(m1,m2) :: eta_SL, eta_l
COMPLEX(CP) :: coeff, coeff2
REAL(RP) :: k_n2, error
INTEGER :: i1,i2,ii,ii2
!
! Test for Laplacian of phi
!
REAL(RP), DIMENSION(m1,m2) :: accx_SL, accz_SL, accy_SL
! For CCSL
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_etax, a_etay
REAL(RP), DIMENSION(m1,m2) :: etax_l, etay_l, CCSL, CDSL, deta_l
!
! FIXME: check that modesSL are the modes, should be OK
CALL fourier_2_space(modesSL,eta_l)
CALL fourier_2_space(modesSLt,deta_l)
!
DO i1=1,n1o2p1
    DO i2=1,n2
        a_etax(i1,i2) = modesSL(i1,i2)*ikx(i1,i2)
        a_etay(i1,i2) = modesSL(i1,i2)*iky(i1,i2)
    ENDDO
ENDDO
CALL fourier_2_space(a_etax,etax_l)
CALL fourier_2_space(a_etay,etay_l)
!
phi_SL(1:n1,1:n2)   = 0.0_rp
vitx_SL(1:n1,1:n2)  = 0.0_rp
vity_SL(1:n1,1:n2)  = 0.0_rp
vitz_SL(1:n1,1:n2)  = 0.0_rp
dphit_SL(1:n1,1:n2) = 0.0_rp
eta_SL(1:n1,1:n2)   = 0.0_rp
!
accx_SL(1:n1,1:n2)  = 0.0_rp
accy_SL(1:n1,1:n2)  = 0.0_rp
accz_SL(1:n1,1:n2)  = 0.0_rp
!
! Direct method
!
DO ii=1,n1
    DO ii2=1,n2
        !
        ! Take into account finite depth...
        !
        i1 = 1
        i2 = 1
        k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
        ! constant mode
        !
        phi_SL(ii,ii2)   = 0.0_rp
        vitx_SL(ii,ii2)  = REAL(modesspecx(i1,i2),RP)
        vity_SL(ii,ii2)  = REAL(modesspecy(i1,i2),RP)
        vitz_SL(ii,ii2)  = REAL(modesspecz(i1,i2),RP)
        dphit_SL(ii,ii2) = REAL(modesspect(i1,i2),RP)
        eta_SL(ii,ii2)   = REAL(modesSL(i1,i2),RP)
        !
        accx_SL(ii,ii2) = 0.0_rp ! it's a sine serie... but constant mode?
        accy_SL(ii,ii2) = 0.0_rp ! i*ky_n2(i2)*modesspecy(i1,i2)
        accz_SL(ii,ii2) = 0.0_rp ! modesspecz(i1,i2)
        !
        ! i1=1
        DO i2=2,n2o2p1
            k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
            IF ((k_n2*(eta_l(ii,ii2)+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
                coeff = COSH(k_n2*(eta_l(ii,ii2)+depth_star))/COSH(k_n2*depth_star)
                coeff2= SINH(k_n2*(eta_l(ii,ii2)+depth_star))/SINH(k_n2*depth_star)
            ELSE
                coeff = EXP(k_n2*eta_l(ii,ii2))
                coeff2= coeff
            ENDIF
            !
            phi_SL(ii,ii2)   =  phi_SL(ii,ii2)  + 0.0_rp
            vitx_SL(ii,ii2)  =  vitx_SL(ii,ii2) + 2.0_rp*ABS(modesspecx(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecx(i1,i2)),REAL(modesspecx(i1,i2),RP)))
            vity_SL(ii,ii2)  =  vity_SL(ii,ii2) + 2.0_rp*ABS(modesspecy(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecy(i1,i2)),REAL(modesspecy(i1,i2),RP)))
            vitz_SL(ii,ii2)  =  vitz_SL(ii,ii2) + 2.0_rp*ABS(modesspecz(i1,i2) * coeff2) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecz(i1,i2)),REAL(modesspecz(i1,i2),RP)))
            dphit_SL(ii,ii2) = dphit_SL(ii,ii2) + 2.0_rp*ABS(modesspect(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspect(i1,i2)),REAL(modesspect(i1,i2),RP)))
            eta_SL(ii,ii2)   =  eta_SL(ii,ii2)  + 2.0_rp*ABS(modesSL(i1,i2)) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesSL(i1,i2)),REAL(modesSL(i1,i2),RP)))
            !
            accx_SL(ii,ii2)  =  accx_SL(ii,ii2) + 0.0_rp
            accy_SL(ii,ii2)  =  accy_SL(ii,ii2) - 2.0_rp*ABS(modesspecy(i1,i2) * coeff * i * ky_n2(i2)) &
                *SIN(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecy(i1,i2)),REAL(modesspecy(i1,i2),RP)))
            accz_SL(ii,ii2)  =  accz_SL(ii,ii2) + 2.0_rp*ABS(modesspecz(i1,i2) * coeff2 * k_n2 &
                / TANH(k_n2*(eta_l(ii,ii2)+depth_star)))*COS(ky_n2(i2)*y(ii2) &
                + ATAN2(AIMAG(modesspecz(i1,i2)),REAL(modesspecz(i1,i2),RP)))
        ENDDO
        ! With the previous computation, mode n2o2p1 is computed twice if n2 even but should be computed only once...
        IF (iseven(n2)) THEN
            i2=n2o2p1
            phi_SL(ii,ii2)   =  phi_SL(ii,ii2)  + 0.0_rp
            vitx_SL(ii,ii2)  =  vitx_SL(ii,ii2) - 1.0_rp*ABS(modesspecx(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecx(i1,i2)),REAL(modesspecx(i1,i2),RP)))
            vity_SL(ii,ii2)  =  vity_SL(ii,ii2) - 1.0_rp*ABS(modesspecy(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecy(i1,i2)),REAL(modesspecy(i1,i2),RP)))
            vitz_SL(ii,ii2)  =  vitz_SL(ii,ii2) - 1.0_rp*ABS(modesspecz(i1,i2) * coeff2) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecz(i1,i2)),REAL(modesspecz(i1,i2),RP)))
            dphit_SL(ii,ii2) = dphit_SL(ii,ii2) - 1.0_rp*ABS(modesspect(i1,i2) * coeff) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspect(i1,i2)),REAL(modesspect(i1,i2),RP)))
            eta_SL(ii,ii2)   =  eta_SL(ii,ii2)  - 1.0_rp*ABS(modesSL(i1,i2)) &
                *COS(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesSL(i1,i2)),REAL(modesSL(i1,i2),RP)))
            !
            accx_SL(ii,ii2)  =  accx_SL(ii,ii2) + 0.0_rp
            accy_SL(ii,ii2)  =  accy_SL(ii,ii2) + 1.0_rp*ABS(modesspecy(i1,i2) * coeff * i * ky_n2(i2)) &
                *SIN(ky_n2(i2)*y(ii2)+ATAN2(AIMAG(modesspecy(i1,i2)),REAL(modesspecy(i1,i2),RP)))
            accz_SL(ii,ii2)  =  accz_SL(ii,ii2) - 1.0_rp*ABS(modesspecz(i1,i2) * coeff2 * k_n2 &
                / TANH(k_n2*(eta_l(ii,ii2)+depth_star)))*COS(ky_n2(i2)*y(ii2) &
                + ATAN2(AIMAG(modesspecz(i1,i2)),REAL(modesspecz(i1,i2),RP)))
        ENDIF
        !
        DO i1=2,n1o2p1
            DO i2=1,n2
                k_n2 = SQRT(kx(i1)**2+ky_n2(i2)**2)
                IF ((k_n2*(eta_l(ii,ii2)+depth_star).LT.50.).AND.(k_n2*depth_star.LT.50.)) THEN
                    coeff = COSH(k_n2*(eta_l(ii,ii2)+depth_star))/COSH(k_n2*depth_star) * EXP(i*ky_n2(i2)*y(ii2))
                    coeff2= SINH(k_n2*(eta_l(ii,ii2)+depth_star))/SINH(k_n2*depth_star) * EXP(i*ky_n2(i2)*y(ii2))
                ELSE
                    coeff = EXP(k_n2*eta_l(ii,ii2)) * EXP(i*ky_n2(i2)*y(ii2))
                    coeff2= coeff
                ENDIF
                ! FIXME; why this *kx(i1) ???
                phi_l(i1,i2)   = modesspec(i1,i2)  * coeff
                vitx_l(i1,i2)  = modesspecx(i1,i2) * coeff
                vity_l(i1,i2)  = modesspecy(i1,i2) * coeff
                vitz_l(i1,i2)  = modesspecz(i1,i2) * coeff2
                dphit_l(i1,i2) = modesspect(i1,i2) * coeff
                !
                phi_SL(ii,ii2)   =  phi_SL(ii,ii2)  &
                + 1.0_rp*ABS(phi_l(i1,i2))*COS(kx(i1)*x(ii)+ATAN2(AIMAG(phi_l(i1,i2)),REAL(phi_l(i1,i2),RP)))
                vitx_SL(ii,ii2)  =  vitx_SL(ii,ii2) &
                + 1.0_rp*ABS(vitx_l(i1,i2))*COS(kx(i1)*x(ii)+ATAN2(AIMAG(vitx_l(i1,i2)),REAL(vitx_l(i1,i2),RP)))
                vity_SL(ii,ii2)  =  vity_SL(ii,ii2) &
                + 1.0_rp*ABS(vity_l(i1,i2))*COS(kx(i1)*x(ii)+ATAN2(AIMAG(vity_l(i1,i2)),REAL(vity_l(i1,i2))))
                vitz_SL(ii,ii2)  =  vitz_SL(ii,ii2) &
                + 1.0_rp*ABS(vitz_l(i1,i2))*COS(kx(i1)*x(ii)+ATAN2(AIMAG(vitz_l(i1,i2)),REAL(vitz_l(i1,i2),RP)))
                dphit_SL(ii,ii2) = dphit_SL(ii,ii2) &
                + 1.0_rp*ABS(dphit_l(i1,i2))*COS(kx(i1)*x(ii)+ATAN2(AIMAG(dphit_l(i1,i2)),REAL(dphit_l(i1,i2),RP)))
                eta_SL(ii,ii2)   =  eta_SL(ii,ii2)  + 1.0_rp*ABS(modesSL(i1,i2) * EXP(i*ky_n2(i2)*y(ii2))) &
                *COS(kx(i1)*x(ii)+ATAN2(AIMAG(modesSL(i1,i2) * EXP(i*ky_n2(i2)*y(ii2))) &
                ,REAL(modesSL(i1,i2) * EXP(i*ky_n2(i2)*y(ii2)),RP)))
                !
                accx_SL(ii,ii2)  =  accx_SL(ii,ii2) - 1.0_rp*ABS(kx(i1)*vitx_l(i1,i2)) &
                *SIN(kx(i1)*x(ii)+ATAN2(AIMAG(vitx_l(i1,i2)),REAL(vitx_l(i1,i2),RP)))
                accy_SL(ii,ii2)  =  accy_SL(ii,ii2) + 1.0_rp*ABS(i * ky_n2(i2) * vity_l(i1,i2)) &
                *COS(kx(i1)*x(ii)+ATAN2(AIMAG(vity_l(i1,i2)*i*ky_n2(i2)),REAL(vity_l(i1,i2)*i*ky_n2(i2),RP)))
                accz_SL(ii,ii2)  =  accz_SL(ii,ii2) + 1.0_rp*ABS(k_n2/TANH(k_n2*(eta_l(ii,ii2)+depth_star))*vitz_l(i1,i2)) &
                *COS(kx(i1)*x(ii)+ATAN2(AIMAG(vitz_l(i1,i2)),REAL(vitz_l(i1,i2),RP)))
            ENDDO
        ENDDO
        !
        ! Evaluate free surface boundary conditions
        ! FIXME: change variable names
        CDSL(ii,ii2) = - g_star*eta_SL(ii,ii2) - 0.5_rp*(vitx_SL(ii,ii2)**2+vity_SL(ii,ii2)**2+vitz_SL(ii,ii2)**2) - dphit_SL(ii,ii2)
        CCSL(ii,ii2) =  deta_l(ii,ii2) - vitz_SL(ii,ii2) + vitx_SL(ii,ii2)*etax_l(ii,ii2) + vity_SL(ii,ii2)*etay_l(ii,ii2)
    ENDDO
ENDDO
!
IF (n2==1) THEN
    accy_SL = 0.0_rp
ENDIF
!
WRITE(*,*) 'KFSBC/DFSBC reconstructed', MAX(MAXVAL(ABS(CCSL(1:n1,1:n2))),MAXVAL(ABS(CDSL(1:n1,1:n2))))
!
! Compute error
!
IF (n1.NE.1) THEN
    error = MAXVAL(ABS((vitxref_SL(1:n1,1:n2) - vitx_SL(1:n1,1:n2))))/MAXVAL(ABS(vitxref_SL(1:n1,1:n2)))
ELSE
    error = 0.0_rp
ENDIF
error = MAX(error,MAXVAL(ABS((vitzref_SL(1:n1,1:n2) - vitz_SL(1:n1,1:n2))))/MAXVAL(ABS(vitzref_SL(1:n1,1:n2))))
error = MAX(error,MAXVAL(ABS((dphitref_SL(1:n1,1:n2) - dphit_SL(1:n1,1:n2))))/MAXVAL(ABS(dphitref_SL(1:n1,1:n2))))
IF (n2.NE.1) THEN
    error = MAX(error,MAXVAL(ABS(vityref_SL(1:n1,1:n2) - vity_SL(1:n1,1:n2)))/MAXVAL(ABS((vityref_SL(1:n1,1:n2)))))
ENDIF
!
END SUBROUTINE
!
! end reconstruction

!
! start HOSvel2_direct *************************************************
!
!     ======================================
!      SUBROUTINE HOSvel2(meta2,M_HOSvel,a_eta_l,a_phis_l,time_current,Ta,n,iCPUtime)
        SUBROUTINE HOSvel_direct(a_eta_l,a_phis_l,time_current)
!     ======================================

IMPLICIT NONE !double precision (a-h,o-z)

!
!% INPUT VARIABLES
!   Number of levels for eta2
REAL(RP), INTENT(IN) :: time_current
!
COMPLEX(CP), DIMENSION(m1o2p1,m2) , INTENT(IN) :: a_eta_l, a_phis_l
!% LOCAL VARIABLES
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_deta_l, a_dphis_l
REAL(RP), DIMENSION(m1,m2)        :: eta3, deta3, dphis3, phis3
!
REAL(RP), DIMENSION(md1,md2) :: etax2,etay2
REAL(RP), DIMENSION(md1,md2) :: phimx_ext,phimy_ext,phimx_ext2,phimy_ext2,phimz_ext2
REAL(RP), DIMENSION(m1,m2)   :: detaphiz_d,phimxetax,phiz_etax,phiz_etay,phisx_d,phisy_d,phimx,phimy,phimz,phit
REAL(RP), DIMENSION(m1,m2)   :: phimx2,phimy2,phimz2,etax_tmp,etay_tmp,Press_bis,CCSL_bis,CCSL,Press_SL
!
INTEGER :: i1,i2
!
! One uses mode as input
!
CALL fourier_2_space(a_eta_l,   eta3)
CALL fourier_2_space(a_phis_l,  phis3)
!FIXME: as2 not defined
!CALL fourier_2_space(as2, phis3)
!
! Compute da_eta
!
CALL solveHOS_lin(a_phis_l, a_dphis_l, a_eta_l, a_deta_l, time_current)
!
! FIXME: check where it comes from in resol_HOS
a_dphis_l(1,1) = a_dphis_l(1,1)+a_eta_l(1,1)
!
CALL fourier_2_space(a_deta_l,  deta3)
CALL fourier_2_space(a_dphis_l,  dphis3)
!
DO i2=1,n2
    DO i1=1,n1
        deta3(i1,i2)  = (deta3(i1,i2) + W1(i1,i2))
        dphis3(i1,i2) = (dphis3(i1,i2) - eta3(i1,i2))
    ENDDO
ENDDO
!
! de-aliasing of terms involved in multiple products
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        etax2(i1,i2) = etax(i1,i2) * etax(i1,i2)
        etay2(i1,i2) = etay(i1,i2) * etay(i1,i2)
    ENDDO
ENDDO

CALL dealias(2,etax2)
CALL dealias(2,etay2)
!
! initializations
!
DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        phimx_ext(i1,i2)= (phisx(i1,i2) - etax(i1,i2) * phizMm1(i1,i2))!*(1.0_rp - exp(-(time_current/Ta)**n))
        phimy_ext(i1,i2)= (phisy(i1,i2) - etay(i1,i2) * phizMm1(i1,i2))!*(1.0_rp - exp(-(time_current/Ta)**n))
        phimx_ext2(i1,i2)= (phisx(i1,i2)*phisx(i1,i2) + etax2(i1,i2)*phiz2Mm2(i1,i2) &
           - 2.0_rp*phisx(i1,i2)*etax(i1,i2) * phizMm1(i1,i2))!*(1.0_rp - exp(-(time_current/Ta)**n))
        phimy_ext2(i1,i2)= (phisy(i1,i2)*phisy(i1,i2) + etay2(i1,i2)*phiz2Mm2(i1,i2) &
           - 2.0_rp*phisy(i1,i2)*etay(i1,i2) * phizMm1(i1,i2))!*(1.0_rp - exp(-(time_current/Ta)**n))
        phimz_ext2(i1,i2) = phiz2(i1,i2)!*(1.0_rp - exp(-(time_current/Ta)**n))
        temp_R_Nd(i1,i2) =   (etax2(i1,i2)+ etay2(i1,i2)) * phiz2Mm2(i1,i2) + phiz2(i1,i2) &
           - phisx(i1,i2) * etax(i1,i2) * phizMm1(i1,i2) - phisy(i1,i2) * etay(i1,i2) * phizMm1(i1,i2) !detaphiz(i1,i2)
        temp2_R_Nd(i1,i2) = ((etax(i1,i2)*etax(i1,i2) + etay(i1,i2)*etay(i1,i2)) * phizMm2(i1,i2) &
           - phisx(i1,i2) * etax(i1,i2) - phisy(i1,i2) * etay(i1,i2)) !phimxetax_ext(i1,i2)
    ENDDO
ENDDO

detaphiz_d = filter_ext(temp_R_Nd,n1o2p1,n2o2p1)
phimxetax  = filter_ext(temp2_R_Nd,n1o2p1,n2o2p1)

DO i1 = 1,Nd1
    DO i2 = 1,Nd2
        temp_R_Nd(i1,i2) = etax(i1,i2) * phizMm1(i1,i2)  !phiz_etax_ext(i1,i2)
        temp2_R_Nd(i1,i2) = etay(i1,i2) * phizMm1(i1,i2) !phiz_etay_ext(i1,i2)
    ENDDO
ENDDO

phiz_etax = filter_ext(temp_R_Nd,n1o2p1,n2o2p1)
phiz_etay = filter_ext(temp2_R_Nd,n1o2p1,n2o2p1)
phisx_d   = filter_ext(phisx,n1o2p1,n2o2p1)
phisy_d   = filter_ext(phisy,n1o2p1,n2o2p1)

DO i1=1,n1
    DO i2=1,n2
        ! FIXME: phis3 not defined
        !phim(i1,i2)  = phis3(i1,i2)
        phimx(i1,i2) = (phisx_d(i1,i2) - phiz_etax(i1,i2))
        phimy(i1,i2) = (phisy_d(i1,i2) - phiz_etay(i1,i2))
        phit(i1,i2)  = dphis3(i1,i2) - (detaphiz_d(i1,i2))
        phimz(i1,i2) = phiz(i1,i2)
    ENDDO
ENDDO

phimx  = filter_ext(phimx_ext,n1o2p1,n2o2p1)
phimy  = filter_ext(phimy_ext,n1o2p1,n2o2p1)
phimx2 = filter_ext(phimx_ext2,n1o2p1,n2o2p1)
phimy2 = filter_ext(phimy_ext2,n1o2p1,n2o2p1)
phimz2 = filter_ext(phimz_ext2,n1o2p1,n2o2p1)

etax_tmp = filter_ext(etax,n1o2p1,n2o2p1)
etay_tmp = filter_ext(etay,n1o2p1,n2o2p1)
! FIXME : vity2ref_SL ???? why was it commented?
DO i1=1,n1
    DO i2=1,n2
        phiref_SL(i1,i2)  = phis3(i1,i2) !phim(i1,i2)
        vitxref_SL(i1,i2) = phimx(i1,i2)
        vityref_SL(i1,i2) = phimy(i1,i2)
        vitzref_SL(i1,i2) = phimz(i1,i2)
        dphitref_SL(i1,i2) = phit(i1,i2)
        ! Comment to have exact CDSL... FIXME ?
        !phimx2(i1,i2) = phimx(i1,i2)*phimx(i1,i2)
        !phimy2(i1,i2) = phimy(i1,i2)*phimy(i1,i2)
        !phimz2(i1,i2) = phimz(i1,i2)*phimz(i1,i2)
        !
        vitx2ref_SL(i1,i2) = phimx2(i1,i2)
        vity2ref_SL(i1,i2) = phimy2(i1,i2)
        vitz2ref_SL(i1,i2) = phimz2(i1,i2)

        Press_SL(i1,i2) = - eta3(i1,i2) - phit(i1,i2) +(- 0.5_rp*(phimx2(i1,i2) + phimy2(i1,i2) + phimz2(i1,i2))) ! CDSL OK si phimx2.NE.phimx*phimx
        CCSL(i1,i2) = deta3(i1,i2) - phimxetax(i1,i2) - phimz(i1,i2)       ! CCSL OK
        ! Approximation (no dealiasing)
        Press_bis(i1,i2) = - eta3(i1,i2) - phit(i1,i2) - 0.5_rp*(phimx(i1,i2)**2 + phimy(i1,i2)**2 + phimz(i1,i2)**2)
        CCSL_bis(i1,i2)  = deta3(i1,i2) + (phimx(i1,i2)*etax_tmp(i1,i2)+phimy(i1,i2)*etay_tmp(i1,i2)) - phimz(i1,i2)       ! CCSL OK
    ENDDO
ENDDO

WRITE(*,*) 'KFSBC/DFSBC =', MAX(MAXVAL(ABS(CCSL(1:n1,1:n2))),MAXVAL(ABS(Press_SL(1:n1,1:n2)))), 'approx',&
    MAX(MAXVAL(ABS(CCSL_bis(1:n1,1:n2))),MAXVAL(ABS(Press_bis(1:n1,1:n2))))
!
CALL Velocity_modes_direct(phis3,dphitref_SL,eta3,kx,ky_n2,x,y,depth_star, &
        modesspec,modesspecx,modesspecy,modesspecz,modesspect)
!
END SUBROUTINE HOSvel_direct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Velocity_modes_direct(phis,dphisdt,eta,kx,ky,x,y,depth, &
    modesphi,modesphix,modesphiy,modesphiz,modesphit)
!
!USE Variables_3D, ONLY : m1,m2,m1o2p1
IMPLICIT NONE
!
REAL(RP), DIMENSION(m1,1)   :: phis,dphisdt,eta
REAL(RP), DIMENSION(1)      :: ky,y
REAL(RP), DIMENSION(m1)     :: x
REAL(RP), DIMENSION(m1o2p1) :: kx
REAL(RP) :: depth
!
COMPLEX(CP), DIMENSION(m1*1,m1*1) :: M
COMPLEX(CP), DIMENSION(m1o2p1,1) :: modesphi, modesphix, modesphiy, modesphiz, modesphit
!
! Check m2=1
!
IF (m2.NE.1) THEN
    PRINT*,'direct inversion only working (RAM memory use) in 2D'
    stop
ENDIF
!
! Create the matrix
!
CALL CreateMatrix(kx,ky,x,y,eta,depth,M)
!!
!! Invert the matrix
!!
!CALL InvertMatrix(M,n1*n2)
!!
!! Compute the modal amplitudes
!!
!CALL modal_amplitude(M,phis,modesspec)
!
! Solve AX=B for phis
!
CALL solve_system(M,phis,modesphi)
!!
!! Compute derivatives
!!
!DO i1=1,n1o2p1
!   DO i2=1,n2
!       k = SQRT(kx(i1)**2+ky(i2)**2)
!       kthk = k*TANH(k*depth)
!       modesphiz(i1,i2) = modesphi(i1,i2)*kthk
!       modesphix(i1,i2) = modesphi(i1,i2)*i*kx(i1)
!   ENDDO
!ENDDO
!IF (n2.EQ.1) THEN
!   modesphiy = 0.0_rp
!ELSE
!   DO i1=1,n1o2p1
!       modesphiy(i1,1:n2) = modesphi(i1,1:n2)*i*ky(1:n2)
!   ENDDO
!ENDIF
!
! Use velocities on FS
!
CALL CreateMatrix(kx,ky,x,y,eta,depth,M)
CALL solve_system(M,vitxref_SL,modesphix)
!
CALL CreateMatrix2(kx,ky,x,y,eta,depth,M)
CALL solve_system(M,vitzref_SL,modesphiz)
!
!modesphiz(1,1) = 0.0_rp
!
IF(n2.EQ.1)THEN
    modesphiy=((0.0_rp, 0.0_rp))
ELSE
    !
    CALL CreateMatrix(kx,ky,x,y,eta,depth,M)
    CALL solve_system(M,vityref_SL,modesphiy)
    !
ENDIF
!
! Create the matrix
!
CALL CreateMatrix(kx,ky,x,y,eta,depth,M)
!
!
! Solve AX=B for dphisdt
!
CALL Solve_system(M,dphisdt,modesphit)
!
END SUBROUTINE Velocity_modes_direct
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
SUBROUTINE Solve_system(A,B,X)
!
IMPLICIT NONE
!
COMPLEX(CP), DIMENSION(m1*1,m1*1), INTENT(IN) :: A
COMPLEX(CP), DIMENSION(m1*1,m1*1) :: A_save
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: B
COMPLEX(CP), DIMENSION(m1*1) :: work,work2
COMPLEX(CP), DIMENSION(m1o2p1,1), INTENT(OUT) :: X
!
INTEGER, DIMENSION(m1*1)  :: ipiv
!
INTEGER :: i1,i2,index1,info
!
! Solve AX=B for dphisdt
!
DO i2=1,1
    DO i1=1,n1
        index1 = i1+n1*(i2-1)
        work(index1) = B(i1,i2)
        work2(index1)= work(index1)
    ENDDO
ENDDO

A_save(1:n1*1,1:n1*1)=A(1:n1*1,1:n1*1)

CALL zgesv(n1*1,1,A_save(1:n1*1,1:n1*1),n1*1,ipiv, work(1:n1*1), n1*1, info)
!
IF(info/=0)THEN
  PRINT *, 'Problems with inversion of the matrix',info
  STOP
ENDIF
!!
!! Evaluate error of inversion
!!
!DO index1=1,n1*n2
!   error(index1) = SUM(A(index1,1:n1*n2)*work(1:n1*n2))-work2(index1)
!ENDDO
!!
!PRINT*,'error max=', MAXVAL(ABS(error)),MAXVAL(ABS(work2))
!
! One just save half of the matrix (symetry)
!
i1=1
i2=1
index1=i1+n1*(i2-1)
X(i1,i2) = work(index1)
!
DO i1=2,n1o2p1
    index1=i1+n1*(i2-1)
    X(i1,i2) = 2.0_rp*work(index1)
ENDDO
!
i1=1
DO i2=2,1
    index1=i1+n1*(i2-1)
    X(i1,i2) = work(index1)
ENDDO
!
DO i2=2,1
    DO i1=2,n1o2p1
        index1=i1+n1*(i2-1)
        X(i1,i2) = 2.0_rp*work(index1)
    ENDDO
ENDDO
IF(iseven(n1)) THEN
    i1=n1o2p1
    DO i2=1,1
        index1=i1+n1*(i2-1)
        X(i1,i2) = work(index1)
    ENDDO
ENDIF
!
END SUBROUTINE solve_system
!
SUBROUTINE CreateMatrix(kx,ky,x,y,eta,h,M)
! By G. Ducrozet.
USE Variables_3D, ONLY : n1,n2,n1o2p1
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1), INTENT(IN) :: x
REAL(RP), DIMENSION(m1o2p1), INTENT(IN) :: kx
REAL(RP), DIMENSION(1), INTENT(IN) :: y, ky
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: eta
COMPLEX(CP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
REAL(RP), INTENT(IN) :: h
REAL(RP) :: k
!
M = ((0.0_rp, 0.0_rp))
!
IF(h.LE.0.0_rp) THEN ! Infinite depth
    DO i2=1,n2
        DO i1=1,n1
            index1=i1+n1*(i2-1)
            DO ii2=1,n2
                DO ii1=1,n1o2p1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT(kx(ii1)**2+ky(ii2)**2)
                    M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                ENDDO
                DO ii1=n1o2p1+1,n1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT((-kx(n1-ii1+2))**2+ky(ii2)**2)
                    M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ELSE
    DO i2=1,n2
        DO i1=1,n1
            index1=i1+n1*(i2-1)
            DO ii2=1,n2
                DO ii1=1,n1o2p1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT(kx(ii1)**2+ky(ii2)**2)
                    IF(k*h.LT.50.0_rp) THEN
                        M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*COSH(k*(eta(i1,i2)+h))/COSH(k*h)
                    ELSE
                        M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                    ENDIF
                ENDDO
                DO ii1=n1o2p1+1,n1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT((-kx(n1-ii1+2))**2+ky(ii2)**2)
                    IF(k*h.LT.50.0_rp) THEN
                        M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*COSH(k*(eta(i1,i2)+h))/COSH(k*h)
                    ELSE
                        M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDIF
!
END SUBROUTINE CreateMatrix
!
SUBROUTINE CreateMatrix2(kx,ky,x,y,eta,h,M)
! By G. Ducrozet.
USE Variables_3D, ONLY : n1,n2,n1o2p1
IMPLICIT NONE
!
INTEGER :: i1,i2,ii1,ii2,index1,index2
REAL(RP), DIMENSION(m1), INTENT(IN) :: x
REAL(RP), DIMENSION(m1o2p1), INTENT(IN) :: kx
REAL(RP), DIMENSION(1), INTENT(IN) :: y, ky
REAL(RP), DIMENSION(m1,1), INTENT(IN) :: eta
COMPLEX(CP), DIMENSION(m1*1,m1*1), INTENT(OUT) :: M
REAL(RP), INTENT(IN) :: h
REAL(RP) :: k
!
M = ((0.0_rp, 0.0_rp))
!
IF(h.LE.0.0_rp) THEN ! Infinite depth
    DO i2=1,n2
        DO i1=1,n1
            index1=i1+n1*(i2-1)
            DO ii2=1,n2
                DO ii1=1,n1o2p1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT(kx(ii1)**2+ky(ii2)**2)
                    M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                ENDDO
                DO ii1=n1o2p1+1,n1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT((-kx(n1-ii1+2))**2+ky(ii2)**2)
                    M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ELSE
    DO i2=1,n2
        DO i1=1,n1
            index1=i1+n1*(i2-1)
            DO ii2=1,n2
                DO ii1=1,n1o2p1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT(kx(ii1)**2+ky(ii2)**2)
                    IF(k*h.LT.50.0_rp) THEN
                        M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*SINH(k*(eta(i1,i2)+h))/SINH(k*h)
                    ELSE
                        M(index1,index2) = exp(i*(kx(ii1)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                    ENDIF
                ENDDO
                DO ii1=n1o2p1+1,n1
                    index2=ii1+n1*(ii2-1)
                    k=SQRT((-kx(n1-ii1+2))**2+ky(ii2)**2)
                    IF(k*h.LT.50.0_rp) THEN
                        M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*SINH(k*(eta(i1,i2)+h))/SINH(k*h)
                    ELSE
                        M(index1,index2) = exp(i*(-kx(n1-ii1+2)*x(i1)+ky(ii2)*y(i2)))*EXP(k*eta(i1,i2))
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    ENDDO
ENDIF
!
! Specific treatment of k=0
!
M(:,1) = ((1.0_rp, 0.0_rp))
!
END SUBROUTINE CreateMatrix2
!
SUBROUTINE InvertMatrix(M,rank)
!
! Inverse real matrix with LAPACK
!
IMPLICIT NONE
INTEGER :: rank, lwork, info
REAL(RP), DIMENSION(rank,rank) :: M
REAL(RP), DIMENSION(5*rank**2) :: work
INTEGER, DIMENSION(rank) :: ipiv
rank = SIZE(M,1)
CALL dgetrf(rank,rank,M,rank,ipiv,info)
IF(info/=0)THEN
    PRINT *, 'Problems with L-U of the matrix',info
    STOP
ENDIF
lwork=rank
CALL dgetri(rank,M,rank,ipiv,work,lwork,info)
IF(info/=0)THEN
  PRINT *, 'Problems with inversion of the matrix',info
  STOP
ENDIF
END SUBROUTINE InvertMatrix
!
END MODULE velocities
!
