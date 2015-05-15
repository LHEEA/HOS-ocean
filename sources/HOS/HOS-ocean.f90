PROGRAM HOS_ocean
!
! This is the main program for HOS-ocean computations
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
USE output
USE Runge_Kutta
USE energy_calc
USE fourier_r2c
USE ramp
USE input_HOS
USE initial_condition
USE filters
USE linear_wave
!
IMPLICIT NONE
! Time marching parameters
INTEGER, PARAMETER  :: i_adapt    = 1
!
REAL(RP), DIMENSION(m1,m2)        :: eta_ref, phis_ref
COMPLEX(CP), DIMENSION(m1o2p1,m2) :: a_phis_rk, a_eta_rk, da_eta
!
REAL(RP)                          :: delx, dely, pioxlen, pioylen, k2, thkh
INTEGER                           :: Mo2
INTEGER                           :: i1, i2, j, iCPUtime, jj
REAL(RP)                          :: energy(3)
!
REAL(RP) :: t_i, t_f, time_cur, t_tot, time_next, t_i_indiv, t_f_indiv
REAL(RP) :: dt_out, dt_rk4, dt, h_rk, h_loc, dt_lin
INTEGER  :: n_rk, n_rk_tot, n_error, n_er_tot
INTEGER  :: n_hour, n_min, n_sec
REAL(RP) :: error, dt_correc, error_old, eta_scale, phis_scale
REAL(RP) :: kp_real
!
INTEGER, DIMENSION(2) :: idx
!
! Adaptive Time Step Runge Kutta scheme
TYPE(RK_parameters)        :: RK_param
! H2 operator
REAL(RP) :: error_H2
!
! Input file
CALL read_input('input_HOS.dat')
!
CALL initiate_parameters()
!
! Scales
IF(i_case == 3 .or. i_case==31 .or. i_case == 32) THEN
    !
    IF(ABS(depth-1.0e15_rp) <= epsilon(1.0e15_rp)) THEN
        kp_real = (TWOPI/Tp_real)**2/grav
    ELSE
        kp_real = wave_number_r(1/Tp_real,depth,grav,1.0E-15_rp)
    ENDIF
    !
    L = 1.0_rp/kp_real
    T = Tp_real/TWOPI
    ! Outputs
    IF (i_out_dim == 1 .OR. i_out_dim == 3) THEN ! dimensional output
        L_out = L
        T_out = T
    ELSE
        L_out = 1.d0
        T_out = 1.d0
    ENDIF
ELSE
    IF (ABS(depth-1.0e15_rp) <= epsilon(1.0e15_rp)) THEN
       ! Infinite depth
       L = 1.0_rp
       T = 1.0_rp/SQRT(grav)
    ELSE
       ! Finite depth
       L = depth
       T = SQRT(depth / grav)
    ENDIF
    IF (i_out_dim == 1 .OR. i_out_dim == 3) THEN ! dimensional output
        L_out = L
        T_out = T
    ELSE
        L_out = 1.0_rp
        T_out = 1.0_rp
    ENDIF
ENDIF
!
! Going to non dimensional form (1/2)
g_star        = grav / (L/T**2)
xlen_star     = xlen / L
ylen_star     = ylen / L
T_stop_star   = T_stop / T
f_out_star    = f_out * T
depth_star    = depth / L
Ta            = Ta / T
!
! Specific scaling for irreg. cases
IF(i_case.EQ.3 .or. i_case .EQ. 31 .or. i_case .EQ. 32) THEN
    ! FIXME: clarify non-dimensionalization
    E_cible = E_cible / L **2
ENDIF
!
! Normalisation factor for elevation eta
!
eta_out = 1.0_rp
!
! x-direction parameters
order_max(1) = ctype1(1)-1
n1c          = n1o2p1
! Derivative mode
IF (ctype1(2) /= 0) THEN
    N_der(1) = ctype1(2)
ELSE
    N_der(1) = n1c
ENDIF
!
! Dealiasing
IF (ctype1(1) == 1) THEN
    i_dealias(1) = 0
    N_dea(1)     = Nd1o2p1
ELSE
    i_dealias(1) = 1
    N_dea(1)     = n1o2p1 ! dealiased modes = proper order p dealiasing
ENDIF
!
! Filtering mode
IF (ctype1(3) /= 0) THEN
    i_filt(1) = 1
    n1c_filt  = MIN(ctype1(3),n1o2p1)
    n1c       = MAX(1,MIN(n1c_filt,n1o2p1))
ELSE
    i_filt(1) = 0
    n1c_filt  = n1c
ENDIF
!
! y-direction parameters
order_max(2) = ctype2(1)-1
n2c          = n2o2p1
! Derivative mode
IF (n2 == 1) THEN
    N_der(2) = 1
ELSE
    ! dealiased modes
    IF (ctype2(2) /= 0) THEN
        N_der(2) = ctype2(2)
    ELSE
        N_der(2) = n2c
    ENDIF
ENDIF
!
! Dealiasing
IF (ctype2(1) == 1) THEN
    i_dealias(2) = 0
    N_dea(2)     = Nd2o2p1
ELSE
    i_dealias(2) = 1
    N_dea(2)     = n2o2p1 ! dealiased modes
ENDIF
!
! Filtering mode
IF (ctype2(3) /= 0) THEN
    i_filt(2) = 1
    IF (n2 == 1) THEN
        n2c_filt  = 1
    ELSE
        n2c_filt  = MIN(ctype2(3),n2o2p1)
        n2c       = MAX(1,MIN(n2c_filt,n2o2p1))
    ENDIF
ELSE
    i_filt(2) = 0
    n2c_filt  = n2c
ENDIF
! profiling
iCPUtime = 0
!
CALL display_parameters()
!
! Specify length of domain
!
pioxlen = TWOPI / xlen_star
!
IF (n2 == 1) THEN
    pioylen = 0.0_rp
ELSE
    pioylen = TWOPI / ylen_star
ENDIF
!
!   mesh generation
!
delx = xlen_star / n1
DO i1 = 1,n1
    x(i1) = (i1 - 1) * delx
ENDDO
!
IF (n2 == 1) THEN
    dely = 2.0_rp * delx ! 2D case but see the choice of dt
ELSE
    dely = ylen_star / n2
ENDIF
DO i2 = 1,n2
    y(i2) = (i2 - 1) * dely
ENDDO
!
!   wave numbers
DO i1 = 1, Nd1o2p1
    kx(i1)  = REAL(i1 - 1,RP) * pioxlen
ENDDO
DO i2 = 1, Nd2o2p1
    ky(i2)  = REAL(i2 - 1,RP) * pioylen
ENDDO
!  y-wave numbers (on n2 modes)
DO i2 = 1, n2o2p1
    ky_n2(i2) = REAL(i2 - 1,RP) * pioylen
ENDDO
DO i2 = 2,n2o2p1
    ky_n2(n2-i2+2) = - REAL(i2 - 1,RP) * pioylen
ENDDO
!
IF (iseven(n2)) ky_n2(n2o2p1) = REAL(n2o2p1 - 1,RP) * pioylen
! Storage for derivatives
ikx = 0.0_cp
!  x-derivative on n1 points (i.e. n1o2p1 modes)
DO i2 = 1, n2o2p1
    ikx(1:MIN(N_der(1),n1o2p1),i2) = i * kx(1:MIN(N_der(1),n1o2p1))
ENDDO
! Last mode contains cos information only and must not be part of the differentiation.
IF (iseven(n1)) ikx(n1o2p1,:) = 0.0_cp
! negative ky
DO i2 = 2, n2o2p1
    ikx(1:MIN(N_der(1),n1o2p1),n2-i2+2) = i * kx(1:MIN(N_der(1),n1o2p1))
ENDDO
!
iky = 0.0_cp
! y-derivative on n1 points (i.e. n1o2p1 modes)
DO i1 = 1, n1o2p1
    iky(i1,1:MIN(N_der(2),n2o2p1)) = i * ky_n2(1:MIN(N_der(2),n2o2p1))
    ! negative ky
    DO i2 = 2, MIN(N_der(2), n2o2p1)
        iky(i1,n2-i2+2) = - i * ky_n2(i2)
    ENDDO
    IF (iseven(n2) .AND. N_der(2)>=n2o2p1) iky(i1,n2o2p1) = 0.0_cp
ENDDO
!
!
ikx_big = 0.0_cp
!  x-derivative on Nd1 points (i.e. Nd1o2p1 modes)
DO i2 = 1, Nd2
    ikx_big(1:MIN(N_der(1),Nd1o2p1),i2) = i * kx(1:MIN(N_der(1),Nd1o2p1))
ENDDO
! Last mode contains cos information only and must not be part of the differentiation.
IF (iseven(Nd1) .AND. N_der(1) >= Nd1o2p1) ikx_big(Nd1o2p1,:) = 0.0_cp
!
iky_big = 0.0_cp
! y-derivative on Nd1 points (i.e. Nd1o2p1 modes)
DO i1 = 1, Nd1o2p1
    ! positive ky
    iky_big(i1,1:MIN(N_der(2),Nd2o2p1)) = i * ky(1:MIN(N_der(2),Nd2o2p1))
    ! negative ky
    DO i2 = 2, MIN(N_der(2),Nd2o2p1) !(Nd2+1)/2
        iky_big(i1,Nd2-i2+2) = - i * ky(i2)
    ENDDO
    ! Last mode contains cos information only and must not be part of the differentiation.!FIXME: check
    IF (iseven(Nd2) .AND. N_der(2) >= Nd2o2p1) iky_big(i1,Nd2o2p1) = 0.0_cp
ENDDO
!
! k_abs may be useful... CHECK
! some variables useful for initialization
! FIXME : make it cleaner
DO i2=1,n2o2p1
    theta_abs(1,i2) = 0.0_rp
    DO i1=2,n1o2p1
        k_abs(i1,i2)=SQRT(kx(i1)*kx(i1)+(ky_n2(i2)*ky_n2(i2)))
        theta_abs(i1,i2)=ATAN2(ky_n2(i2),kx(i1))
        IF (theta_abs(i1,i2) .LT. 0.0_rp) THEN
                theta_abs(i1,i2)=theta_abs(i1,i2) + 2.0_rp*PI
        ENDIF
    ENDDO
ENDDO
DO i2=2,n2o2p1
    theta_abs(1,n2-i2+1) = 0.0_rp
    DO i1=2,n1o2p1
        k_abs(i1,n2-i2+2)=SQRT(kx(i1)*kx(i1)+(ky_n2(i2)*ky_n2(i2)))
        theta_abs(i1,n2-i2+2)=ATAN2(-ky_n2(i2),kx(i1))
        IF (theta_abs(i1,n2-i2+2) .LT. 0.0_rp) THEN
                theta_abs(i1,n2-i2+2)=theta_abs(i1,n2-i2+2) + 2.0_rp*PI
        ENDIF
    ENDDO
ENDDO
!------ Polar angular basis
IF(n2 == 1) THEN
    dth = 1.0
    theta_base=0.0_rp
ELSE
    dth = pi / REAL(ithp-1,RP)
    DO jj=1,ithp/2+1
        theta_base(jj)=dth * (jj-1)
    ENDDO
    DO jj=floor(real(ithp,RP)/2.0)+2,ithp
        theta_base(jj)=3.0_rp * pi / 2.0_rp + dth * (jj-floor(real(ithp,RP)/2.0)-2)
    ENDDO
ENDIF
! HOS modal coefficients of the vertical derivatives
DO i2 = 1, Nd2o2p1
    DO i1 = 1, Nd1o2p1
      k2           = kx(i1) * kx(i1) + ky(i2) * ky(i2)
      k(i1,i2)     = SQRT(k2)
        thkh         = TANH(k(i1,i2) * depth_star)
        kth(i1,i2,1) = k(i1,i2)
      k4(i1,i2)    = k2*k2
        omega(i1,i2) = SQRT(g_star * k(i1,i2) * thkh)
        IF (ABS(k2) > tiny) THEN
         goomega(i1,i2) = g_star / omega(i1,i2)
                c(i1,i2)       = omega(i1,i2) / k(i1,i2)
        ELSE
         goomega(i1,i2) = 1.0_rp
                c(i1,i2)       = SQRT(g_star * depth_star)
        ENDIF
        !
        kth(i1,i2,1) = kth(i1,i2,1) * thkh
        IF (M > 1) kth(i1,i2,2) = k2
        DO j = 3,M
                kth(i1,i2,j) = k2 * kth(i1,i2,j-2)
        ENDDO
    ENDDO
ENDDO
!
goomega(1,1)    = 1.0_rp
!
DO i2 = 1, n2
    DO i1 = 1, Nd1o2p1
        k2                 = kx(i1) * kx(i1) + ky_n2(i2) * ky_n2(i2)
        omega_n2(i1,i2)    = SQRT(g_star * SQRT(k2) * TANH(SQRT(k2)*depth_star))
        IF (ABS(k2) > tiny) THEN
            goomega_n2(i1,i2)  = g_star / omega_n2(i1,i2)
        ELSE
            goomega_n2(i1,i2) = 1.0_rp
        ENDIF
    ENDDO
ENDDO
!
goomega_n2(1,1) = 1.0_rp
!
k4(N_der(1)+1:Nd1o2p1,1:Nd2o2p1) = 0.0_rp
k4(1:Nd1o2p1,N_der(2)+1:Nd2o2p1) = 0.0_rp
!
DO i2 = 2, Nd2p1o2
    DO i1 = 1, Nd1o2p1
        k4(i1, Nd2-i2+2)   = k4(i1,i2)
        kth(i1,Nd2-i2+2,:) = kth(i1,i2,:)
    ENDDO
ENDDO
!
CALL display_velocities()
!
! Initial conditions
!
! Evaluates the modal amplitudes (FT)
CALL fourier_ini(3)
CALL fill_butcher_array(RK_param)
!
Mo2 = M / 2
!
DO j=1,M
    oneoj(j) = 1.0_rp / REAL(j, RP)
ENDDO
!
! Time step control
!
dt_out = 1.0_rp / f_out_star
IF (n2 == 1) THEN
    dt_rk4 = 2.0_rp / SQRT(PI) * SQRT(xlen / n1o2p1)
ELSE
    dt_rk4 = 2.0_rp / SQRT(PI) * MIN(SQRT(xlen / n1o2p1), SQRT(ylen / n2o2p1))
ENDIF
dt_lin = 10.0_rp * dt_rk4
dt     = 0.5_rp * dt_rk4
!
WRITE(*,*) dt_out, dt_rk4, dt_lin
!
! Initial solution
CALL initiate(time_cur, eta, phis, RK_param)
!
! keeping only n1c and n2c modes
CALL space_2_fourier(eta,  a_eta)
CALL space_2_fourier(phis, a_phis)
!
! zero forcing above n1c
a_eta( n1c+1:n1o2p1,:) = 0.0_cp
a_phis(n1c+1:n1o2p1,:) = 0.0_cp
! zero forcing above n2c
a_eta(:, n2c+1:n2-(n2c+1)+2) = 0.0_cp
a_phis(:,n2c+1:n2-(n2c+1)+2) = 0.0_cp
!
! making sure the volume is zero, FIXME: necessary?
IF (ABS(a_eta(1,1)) > tiny) THEN
    WRITE(*,'(A)') 'Nonzero initial volume'
    a_eta(1,1) = 0.0_cp
    a_phis(1,1) = 0.0_cp
ENDIF
!
! Back to physical space
CALL fourier_2_space(a_eta,  eta)
CALL fourier_2_space(a_phis, phis)
! Storing t=0 as reference
IF (i_case.EQ.1 .OR. i_case.EQ.9) THEN ! dimensional initialization
    eta_ref(1:n1,1:n2)  = eta(1:n1,1:n2) / L
    phis_ref(1:n1,1:n2) = phis(1:n1,1:n2) / (L**2 / T)
ELSE ! nondimensional initialization
    eta_ref(1:n1,1:n2)  = eta(1:n1,1:n2)
    phis_ref(1:n1,1:n2) = phis(1:n1,1:n2)
ENDIF
!
! Analytical integration of the linear part
! Evaluates the modal amplitudes (FT)
CALL space_2_fourier(eta_ref,  a_eta)
CALL space_2_fourier(phis_ref, a_phis)
!
! guessing deta_dt at t=0 for energy output at t=0
CALL RK_adapt_2var_3D_in_mo_lin(0, RK_param, 0.0_rp, 0.0_rp, a_phis, a_eta, da_eta)
!****************
energy = calc_energy(a_eta, a_phis, da_eta)
E_o = energy
!
CALL CPU_TIME(t_i)
n_rk      = 0
n_rk_tot  = 0
n_error   = 0
n_er_tot  = 0
time_cur  = 0.0_rp
error_old = toler
!
! Output files
!
CALL init_output(i_3d=i_3d, i_a=i_a_3d, i_vol=1, i_2D=i_2d, i_max=0, i_prob=i_prob, i_sw=i_sw)
!
!
IF (err == 'rel') THEN
    eta_scale     = MAXVAL(ABS(a_eta))  ! relative error
    phis_scale    = MAXVAL(ABS(a_phis)) ! relative error
ELSE IF (err == 'abs') THEN
    eta_scale     = 1.0_rp ! absolute error
    phis_scale    = 1.0_rp ! absolute error
ENDIF
!
DO WHILE (time_cur <= T_stop_star)
    !
    n_er_tot = n_error + n_er_tot ! number or time steps with too large error (wrong time steps)
    n_rk_tot = n_rk_tot + n_rk    ! total number of time steps
    CALL CPU_TIME(t_f)
    ! Dsplays the CPU time for one time steps
    IF (time_cur <= 3000 * dt_out) THEN ! only for the first 100 steps
        WRITE(*,'(A,F6.2,2(X,I4),3(X,F9.3))') 'CPU time for Output time step ',t_f-t_i, n_rk, n_error, time_cur*T_out, dt_rk4, dt
    ENDIF
    ! Guess of the simulation's duration
    if (time_cur <= 100*dt_out) then ! only for the first 10 steps
        t_tot  = (t_f - t_i) * T_stop_star * f_out_star
        n_hour = FLOOR(t_tot / 3600)
        n_min  = FLOOR((t_tot - 3600 * n_hour) / 60)
        n_sec  = FLOOR(t_tot - 60 * (n_min + 60 * n_hour))
        WRITE(*,'(A,I3,A,I2,A,I2,A)') 'Expected time ',n_hour,' h ',n_min, ' min ',n_sec,' sec.'
    end if
    t_i = t_f
    !
    ! Output of volume and energy
    volume = REAL(a_eta(1,1),RP)
    energy = calc_energy(a_eta, a_phis, da_eta)
    !
    ! Test addition velocity comp.
    !
    IF (i_sw == 1) THEN !Test the velocity calculations
        ! Output of errors induced by H2 operator on free surface quantities
        PRINT*, '.......... H2 operator ..........'
        CALL HOSvel2(25,M,a_eta,a_phis,time_cur)
        CALL reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect,a_eta,da_eta,error_H2)
        IF(error_H2.GT.0.1_rp.AND.m2.EQ.1) THEN
            !
            ! Direct method
            !
            PRINT*, '.......... Direct method ..........'
            CALL HOSvel_direct(a_eta,a_phis,time_cur)
            CALL reconstruction_SL(modesspec,modesspecx,modesspecy,modesspecz,modesspect,a_eta,da_eta,error_H2)
        ELSEIF (error_H2.GT.0.1_rp.AND.m2.GT.1) THEN
            PRINT*, 'Max relative error greater than 10% at free surface ...', error_H2
        ENDIF
       ! FIXME : compare to RF ?
    ENDIF
    !
    ! Rough estimation of peak period
    idx = MAXLOC(abs(a_eta))
    Tp = TWOPI/omega_n2(idx(1),idx(2))
    PRINT*,'Hs_out=',4.0_rp*SQRT(energy(3))* L_out,', T_peak=',Tp * T_out
    PRINT*,'***************************'
    CALL output_time_step(i_3d=i_3d, i_a=i_a_3d, i_vol=1, i_2D=i_2d, i_max=0, &
                         time=time_cur, N_stop=FLOOR(T_stop_star / dt_out), &
                         a_eta=a_eta, a_phis=a_phis, da_eta= da_eta, volume=volume, energy=energy, E_0=E_o, E_tot=E_tot, &
                         i_prob=i_prob)
    !
    IF (ABS(time_cur-T_stop_star) <= tiny) EXIT ! output of the last zone is done
    !
    ! Going to next time step
    time_next = time_cur + dt_out
    IF (time_next > T_stop_star) time_next = T_stop_star ! if last output
    h_rk    = dt ! starting time step
    n_rk    = 0  ! local number of time steps
    n_error = 0  ! local number of wrong time steps
    DO WHILE(time_cur < time_next)
        CALL CPU_TIME(t_i_indiv)
        !
        CALL check_slope
        h_loc = h_rk ! local time step
        IF (time_next - time_cur < h_rk) h_loc = time_next - time_cur ! if last time step
        !
        ! Analytical integration of the linear part
        ! starting at current time
        a_eta_rk(1:n1o2p1,1:n2)  = a_eta(1:n1o2p1,1:n2)
        a_phis_rk(1:n1o2p1,1:n2) = a_phis(1:n1o2p1,1:n2)
        !
        ! Going to time_cur + h_loc
        IF (i_adapt == 1) THEN ! adaptive time step    !
            CALL RK_adapt_2var_3D_in_mo_lin(iCPUtime, RK_param, time_cur, h_loc, &
            a_phis_rk, a_eta_rk, da_eta, error, phis_scale, eta_scale)
            error_old = error
            IF (error > toler) THEN
                dt_correc = (error/toler)**(-(1.0_rp/(RK_param%p-1)))
                n_error   = n_error + 1
            ELSE
                a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
                a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
                time_cur   = time_cur + h_loc
                n_rk       = n_rk + 1
                dt_correc  = (error/toler)**(-(1.0_rp/(RK_param%p)))
            ENDIF
            ! Step size ratio bounds (Mathematica)
            IF (dt_correc > 4.0_rp)   dt_correc = 4.0_rp
            IF (dt_correc < 0.125_rp) dt_correc = 0.125_rp
            !
            ! Step size security factor
            dt_correc = 0.98_rp * dt_correc
            ! New step size
            h_rk      = h_rk * dt_correc
            ! Stability checkings
            h_rk      = MIN(h_rk, dt_out, dt_lin)
        ELSE IF (i_adapt == 0) THEN   ! fixed time step
            CALL RK_adapt_2var_3D_in_mo_lin(iCPUtime, RK_param, time_cur, h_loc, a_phis_rk, a_eta_rk, da_eta)
            a_eta(1:n1o2p1,1:n2)  = a_eta_rk(1:n1o2p1,1:n2)
            a_phis(1:n1o2p1,1:n2) = a_phis_rk(1:n1o2p1,1:n2)
            time_cur = time_cur + h_loc
            n_error  = 0
            n_rk     = n_rk + 1
        ENDIF
!
        CALL CPU_TIME(t_f_indiv)
        IF (time_cur <= 10 * dt_out) THEN
            IF ((t_f_indiv-t_i_indiv).gt.2.0d-2) THEN
                WRITE(*,'(A,F6.2,3(X,ES11.4))') 'CPU time for individual time step ',t_f_indiv-t_i_indiv,h_rk,error,time_cur*T_out
            ENDIF
        ENDIF
    ENDDO
    dt = h_rk ! saving the step size for next start
ENDDO
!
CLOSE(1)
CLOSE(2)
CLOSE(3)
CLOSE(4)
CLOSE(5)
CLOSE(6)
CLOSE(7)
CLOSE(8)
CLOSE(9)
CLOSE(10)
CLOSE(11)
CLOSE(12)
CLOSE(130)
!
CONTAINS
!
!
SUBROUTINE check_slope
!
IMPLICIT NONE
!
REAL(RP) :: slopex_max, threshold,slopey_max
!
threshold = 10.0_rp
!
slopex_max = MAXVAL(ABS(etax))
slopey_max = MAXVAL(ABS(etay))
! WRITE(*,*) slopex_max
!
IF (slopex_max > threshold) THEN
    WRITE(*,'(A,2(ES11.3,X))') 'Ca va trancher en x...', time_cur, slopex_max
    STOP
ENDIF
!
IF (slopey_max > threshold) THEN
    WRITE(*,'(A)') 'Ca va trancher en y...'
    STOP
ENDIF
!
END SUBROUTINE check_slope
!
!
!
SUBROUTINE display_velocities()
!
IMPLICIT NONE
!
WRITE(*,'(A)') 'Parameters for the x-direction'
IF (depth < 1.0E14_rp) THEN
    WRITE(*,'(3X, A, ES12.5)')  'Shallow water velocity             c_0      =',SQRT(g_star*depth_star)*L/T
ENDIF
!
END SUBROUTINE display_velocities
!
!
!
SUBROUTINE display_parameters()
!
IMPLICIT NONE
!
WRITE(*,'(A)') 'General parameters'
WRITE(*,'(3X, A, I2)')   'HOS order        ',M
!
WRITE(*,'(A)') 'Parameters for the x-direction'
WRITE(*,'(3X, A, I5)')  'Modes              n1=',n1
WRITE(*,'(3X, A, I5)')  'Modes used        n1c=',n1c
IF (ctype1(1) == 1) THEN
    WRITE(*,'(3X,A)')    'No dealiasing'
ELSE
    WRITE(*,'(3X,A,I2)') 'Dealiasing the products of order ', ctype1(1)
ENDIF
WRITE(*,'(3X, A, I5)')  'Extended modes    Nd1=',Nd1
WRITE(*,'(3X, A, I5)')  'Dealiased modes N_dea=',N_dea(1)
WRITE(*,'(3X, A, I5)')  'Derived modes   N_der=',N_der(1)
IF (i_filt(1) == 0) THEN
    WRITE(*,'(3X,A)')    'No filtering'
ELSE
    WRITE(*,'(3X,A,I5)') 'Filtering with modes  ', n1c_filt
ENDIF
!
WRITE(*,'(A)') 'Parameters for the y-direction'
WRITE(*,'(3X, A, I5)')  'Modes              n2=',n2
WRITE(*,'(3X, A, I5)')  'Modes used        n2c=',n2c
IF (ctype2(1) == 1) THEN
    WRITE(*,'(3X,A)')    'No dealiasing'
ELSE
    WRITE(*,'(3X,A,I2)') 'Dealiasing the products of order ', ctype2(1)
ENDIF
WRITE(*,'(3X, A, I5)')  'Extended modes   Nd2=',Nd2
WRITE(*,'(3X, A, I5)')  'Dealiased modes N_dea=',N_dea(2)
WRITE(*,'(3X, A, I5)')  'Derived modes   N_der=',N_der(2)
IF (i_filt(2) == 0) THEN
    WRITE(*,'(3X,A)')    'No filtering'
ELSE
    WRITE(*,'(3X,A,I5)') 'Filtering with modes ', n2c_filt
ENDIF
!
END SUBROUTINE display_parameters
!
END PROGRAM HOS_ocean
