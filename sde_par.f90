MODULE sde_par
    
	implicit none

    integer, parameter :: dp1 = selected_real_kind(15, 307)
    integer, parameter :: dp = kind(0.d0)           ! double precision
    
    integer, parameter :: N1 = 10000.0_dp, N = 2*(10**8), R1 = 100.0_dp, L = N/R1
    integer, parameter :: H_0 = 500.0_dp, A_1 = 700.0_dp, A_2 = 800.0_dp, K_1 = 35*(H_0 + A_1), K_2 = 150*A_2
	
    real(kind=dp1) :: T = 200.0_dp, beta1 = 0.75_dp, beta2 = 0.365_dp, beta3 = 0.4015_dp

    real(kind=dp1) :: mu_h = 0.01428571_dp, z_1 = 14.6_dp, tz_1 = 36.5_dp, tz_2 = 40.15_dp,	tpi_h = 0.0009_dp
    real(kind=dp1) :: pi_h = 0.0009_dp,	pi_v = 0.03_dp,	tpi_v = 0.49_dp, mu_1 = 0.16666667_dp, mu_2 = 0.4_dp
    real(kind=dp1) :: at_1 = 281.561_dp, at_2 = 255.5_dp, e_1 = 2.4455_dp, e_2 = 2.555_dp, w_1 = 0.9_dp
    real(kind=dp1) :: w_2 = 0.9_dp,	eta_1 = 0.1_dp	
	
	real(kind=dp1) :: Dt_reso, Dt, r_1, r_2    
    real(kind=dp1) :: alpha_1, tsigma_2, alpha_2, tsigma_1, sigma_11, sigma_12, talpha_11, talpha_12, talpha_2

    real(dp) :: X_temp(5,1), dW(3,R1), F(5,1), X(5,1), Wdisc(3,1), G(5,3), X_exact_temp(5,1), X_exact(5,1), F1(5,1)		!double precision
    real(dp) :: X1((L/100)+1,N1), X2((L/100)+1,N1), X3((L/100)+1,N1), X4((L/100)+1,N1), X5((L/100)+1,N1)
    real(dp) :: solu_exact(L+1,5), esp_x(L+1,5), var_x(L+1,5), t_h(L+1,1), t_h_apr((L/100)+1,1)

END MODULE sde_par