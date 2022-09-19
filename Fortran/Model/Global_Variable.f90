!
!   Global_Variable.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
module Global_Variable
    implicit none
! Parameters
integer, parameter :: iter_max = 1d4, iter_IR = 101
real(8), parameter :: TOL = 1d-6, TOL1 = 0.003, TOL2 = 1d-2

! Control
integer, parameter :: calibration = 0, Task = 3
integer :: Np, iunit
character(100) :: Filename, FMT
integer :: Model, Step, SPP, flag
real(8) :: sup, t1, t2

! ======================================================================
! Model parameters
! ======================================================================

real(8), parameter :: theta = 0.9d0, sigma0 = 2d0               ! Preference and ES
real(8), parameter :: delta = 0.1d0                             ! Depreciation
real(8), parameter :: h_HL_star = 0d0, lambda = 0.5d0           ! RD and human capital
real(8), parameter :: epsilon_H = 0.3d0, epsilon_L = 0.7d0      ! Population
real(8) :: a_HL = 5d0
real(8) :: rho
real(8) :: p, sigma                                             ! Elasticity
real(8) :: A, B_N, b_h, B_NH, B_NL                              ! Productivity
real(8) :: mu_N, mu_I, mu_h, mu_hH, mu_hL                       ! RD and human capital
real(8) :: alpha, alpha_H, alpha_L, lambda_h                    ! RD and human capital
real(8) :: shock

! ======================================================================
! Variables
! ======================================================================

! State variables
real(8) :: N, I_tilde, S_tilde, h_HL, gamma_HL 
real(8) :: h, h_H, h_L, k, k_H, k_L

! Statistic 
real(8) :: d_Eta_I, d_Gamma_I, d_Gamma_S, d_I
real(8) :: g_wN, g_wI, g_wS, g_wH
real(8) :: gamma_IL, gamma_SH, gamma_SL, beta_S, t_tilde
real(8) :: a_I, a_S

! Non-state variables
real(8) :: y, pi, c, c_H, c_L, V, V_H, V_L,                     & ! Output
           L, ll_H, ll_L, L_H, L_L,                             & ! Output
           Eta, Gamma, Gamma_H, Gamma_L,                        & ! Productivity
           s_K, s_L, s_LH, s_LL,                                & ! Share
           R, rr, w, w_H, w_L,                                  & ! Price        
           V_NH, V_NL, V_II,                                    & ! Profit
           p_N, p_I,                                            & ! Patent
           epsilon_N, epsilon_I,                                & ! R&D
           g_N, g_I, g_h, g_hH, g_hL, g,                        & ! Growth rate
           g_k, g_L, g_LH, g_LL                                   ! Transition

! ======================================================================
! Policy
! ======================================================================

real(8) :: tau_k, tau_N, tau_I, tau_h
real(8) :: tau_hH, tau_hL 
real(8) :: Trans


! ======================================================================
! Calibration
! ======================================================================

! Input
integer, parameter :: Nparam0 = 14
real(8) :: param0(Nparam0,iter_max)
integer, parameter :: Nparam1 = 18
real(8) :: param1(Nparam1,iter_max)
integer, parameter :: Nparam = 10, Nmoment = 10
real(8) :: moment_data(Nmoment), moment_model(Nmoment), moment_weight(Nmoment)
real(8) :: L_N0

real(8) :: rr0, rr1, w0, w1, w_H0, w_H1, y0, y1
real(8) :: s_K0, s_K1, s_LL0, s_LL1, s_LH0, s_LH1
real(8) :: V_NH0, V_NH1, V_NL0, V_NL1, V_II0, V_II1

! ======================================================================
! Transition
! ======================================================================

! Transition path
real(8), parameter :: T = 100d0, lambda_T = 1d0, rho_ts = 0.05d0
real(8) :: dt
real(8) :: I0, k0, h_HL0, S0
real(8) :: I1, k1, h_HL1, S1
real(8) :: L0, L_H0, L_L0, g_N0, g_I0, epsilon_N0, epsilon_I0
real(8) :: L1, L_H1, L_L1, g_N1, g_I1, epsilon_N1, epsilon_I1
real(8) :: g_h0, g_hH0, g_hL0, ll_H0, ll_L0
real(8) :: g_h1, g_hH1, g_hL1, ll_H1, ll_L1
real(8), dimension(iter_IR) :: ts_t, ts_dt
real(8), dimension(iter_IR) :: ts_I0, ts_k0, ts_h_HL0
real(8), dimension(iter_IR) :: ts_I1, ts_k1, ts_h_HL1
real(8), dimension(iter_IR) :: ts_L0, ts_L_L0, ts_L_H0 
real(8), dimension(iter_IR) :: ts_L1, ts_L_L1, ts_L_H1 
real(8), dimension(iter_IR) :: ts_g_N0, ts_g_I0, ts_g_hH0, ts_g_hL0
real(8), dimension(iter_IR) :: ts_g_N1, ts_g_I1, ts_g_hH1, ts_g_hL1 
real(8), dimension(iter_IR) :: ts_a_dk, ts_aH_dk, ts_aL_dk
real(8), dimension(iter_IR) :: ts_S, ts_p_N, ts_p_I
real(8), dimension(iter_IR) :: ts_g_k, ts_g_L, ts_g_LH, ts_g_LL
real(8), dimension(iter_IR) :: ts_V_NH0, ts_V_NL0, ts_V_II0
real(8), dimension(iter_IR) :: ts_V_NH1, ts_V_NL1, ts_V_II1

! ======================================================================
! Simulation
! ======================================================================

real(8), dimension(iter_IR) :: N_IR,                     & ! Innovation
                               I_tilde_IR,               & ! Automation level
                               S_tilde_IR,               & ! Allocation
                               h_HL_IR,                  & ! Human capital gap
                               h_IR,                     & ! Human capital
                               h_H_IR,                   & ! HS capital
                               h_L_IR,                   & ! LS human capital
                               k_IR,                     & ! Capital
                               k_H_IR,                   & ! HS labor capital
                               k_L_IR,                   & ! LS labor capital
                               L_IR,                     & ! Labor
                               L_H_IR,                   & ! HS labor
                               L_L_IR,                   & ! LS labor
                               y_IR,                     & ! Output
                               c_IR,                     & ! Consumption
                               c_H_IR,                   & ! HS labor consumption
                               c_L_IR,                   & ! LS labor consumption
                               V_IR,                     & ! Welfare
                               V_H_IR,                   & ! HS labor walfare
                               V_L_IR,                   & ! LS labor walfare
                               Eta_IR,                   & ! Capital productivity
                               Gamma_IR,                 & ! Labor productivity
                               Gamma_H_IR,               & ! HS labor productivity
                               Gamma_L_IR,               & ! LS labor productivity
                               s_K_IR,                   & ! Capital share
                               s_L_IR,                   & ! Labor share
                               s_LH_IR,                  & ! HS labor share
                               s_LL_IR,                  & ! LS labor share
                               rr_IR,                    & ! Real interest rate
                               w_IR,                     & ! Wage
                               w_H_IR,                   & ! HS wage
                               w_L_IR,                   & ! LS wage
                               p_N_IR,                   & ! Innovation value
                               p_I_IR,                   & ! Automation value
                               epsilon_N_IR,             & ! Innovation scientist
                               epsilon_I_IR,             & ! Automation scientist
                               g_N_IR,                   & ! Innovation rate
                               g_I_IR,                   & ! Automation rate
                               g_h_IR,                   & ! Human capital growth rate
                               g_hH_IR,                  & ! HS human capital growth rate
                               g_hL_IR,                  & ! LS human capital growth rate
                               g_IR                        ! Growth rate


contains


subroutine File_Create0
    implicit none
        
    open(2, file='Simulate.txt', status='replace')
        FMT = '(A10,10A10)'
        write(2,FMT) 'number', 'p', 'sigma', 'A', 'B_N', 'b_h',             &
                               'mu_N', 'mu_I', 'mu_h', 'alpha', 'shock'
    close(2)

    open(2, file='Policy.txt', status='replace')
        FMT = '(A10,5A10)'
        write(2,FMT) 'number', 'tau_k', 'tau_N', 'tau_I', 'tau_h'
    close(2)
    
    open(2, file='BGP1.txt', status='replace')
        FMT = '(8A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'y', 'k', 'L', 'c', 'V'
    close(2)

    open(2, file='BGP2.txt', status='replace')
        FMT = '(9A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'I_tilde', 'Eta', 'Gamma', 'r', 'w', 's_L'
    close(2)
    
    open(2, file='BGP3.txt', status='replace')
        FMT = '(12A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'RD', 'g', 'g_N', 'g_I', 'g_h', 'p_N', 'p_I', 'eps_N', 'eps_I'
    close(2)

    open(2, file='BGP_SPP1.txt', status='replace')
        FMT = '(8A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'y', 'k', 'L', 'c', 'V'
    close(2)

    open(2, file='BGP_SPP2.txt', status='replace')
        FMT = '(9A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'I_tilde', 'Eta', 'Gamma', 'r', 'w', 's_L'
    close(2)

    open(2, file='BGP_SPP3.txt', status='replace')
        FMT = '(12A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'RD', 'g', 'g_N', 'g_I', 'g_h', 'p_N', 'p_I', 'eps_N', 'eps_I'
    close(2)
    
end subroutine File_Create0


subroutine Read_Param0
    implicit none

    rho = moment_data(1)-theta*(moment_data(2)+moment_data(3))
    p = param0(1,iunit); sigma =param0(2,iunit)*(1d0-p)+p 
    A = param0(3,iunit)**(1d0/(1d0-p))
    B_N = param0(4,iunit); b_h = param0(5,iunit)
    mu_N = param0(6,iunit); mu_I = param0(7,iunit); mu_h = param0(8,iunit)
    alpha = param0(9,iunit); shock = param0(10,iunit)
    tau_k = param0(11,iunit); tau_N = param0(12,iunit); tau_I = param0(13,iunit)
    tau_h = param0(14,iunit)

    
    write(*,*) '=============================================================================='
    FMT = '(6A10)'
    write(*,FMT) 'N', 'p', 'sigma', 'A', 'B_N', 'b_h'
    FMT = '(I10,5F10.2)'
    write(*,FMT) iunit, p, sigma, A, B_N, b_h
    FMT = '(5A10)'
    write(*,FMT) '', 'mu_N', 'mu_I', 'mu_h', 'alpha'
    FMT = '(A10,4F10.2)'
    write(*,FMT) '', mu_N, mu_I, mu_h, alpha
    write(*,*) ''

    open(2, file='Simulate.txt', status='old', position='append')
        FMT = '(I10,10F10.4)'
        write(2,FMT) iunit, param0(1:10,iunit)
    close(2)

    open(2, file='Policy.txt', status='old', position='append')
        FMT = '(I10,4F10.4)'
        write(2,FMT) iunit, param0(11:14,iunit)
    close(2)

endsubroutine Read_Param0


subroutine Save_BGP0(t)
    implicit none
integer, intent(in) :: t
    
    FMT = '(6A10)'
    write(*,FMT) '', 'y', 'k', 'L', 'c', 'V'
    FMT = '(A10,5F10.4)'
    write(*,FMT) '', y, k, L, c, V
    FMT = '(6A10)'
    write(*,FMT) 'I_tilde', 'Eta', 'Gamma', 'r', 'w', 's_L'
    FMT = '(6F10.4)'
    write(*,FMT) I_tilde, Eta, Gamma, rr, w, s_L
    FMT = '(6A10)'
    write(*,FMT) '', 'RD', 'g', 'g_N', 'g_I', 'g_h'
    FMT = '(A10,5F10.4)'
    write(*,FMT) '', w*(epsilon_N+epsilon_I)/y, g, B_N*g_N, B_N*g_I, b_h*g_h

    call Save_Result0(t)
    if (SPP == 0) then
        call File_Write0
    elseif (SPP == 1) then
        call File_Write_SPP0
    endif

end subroutine Save_BGP0


subroutine Read_Result0(t)
    implicit none
integer, intent(in) :: t

    N = N_IR(t); I_tilde = I_tilde_IR(t); h = h_IR(t)
    k = k_IR(t); L = L_IR(t); y = y_IR(t)
    c = c_IR(t); V = V_IR(t)
    Eta = Eta_IR(t); Gamma = Gamma_IR(t)
    s_K = s_K_IR(t); s_L = s_L_IR(t); rr = rr_IR(t); w = w_IR(t)
    p_N = p_N_IR(t); p_I = p_I_IR(t)
    epsilon_N = epsilon_N_IR(t); epsilon_I = epsilon_I_IR(t)
    g = g_IR(t); g_N = g_N_IR(t); g_I = g_I_IR(t); g_h = g_h_IR(t)
    
end subroutine Read_Result0


subroutine Save_Result0(t)
    implicit none
integer, intent(in) :: t

    N_IR(t) = N; I_tilde_IR(t) = I_tilde; h_IR(t) = h
    k_IR(t) = k; L_IR(t) = L; y_IR(t) = y
    c_IR(t) = c; V_IR(t) = V
    Eta_IR(t) = Eta; Gamma_IR(t) = Gamma
    s_K_IR(t) = s_K; s_L_IR(t) = s_L; rr_IR(t) = rr; w_IR(t) = w 
    p_N_IR(t) = p_N; p_I_IR(t) = p_I
    epsilon_N_IR(t) = epsilon_N; epsilon_I_IR(t) = epsilon_I
    g_IR(t) = g; g_N_IR(t) = g_N; g_I_IR(t) = g_I; g_h_IR(t) = g_h
    
end subroutine Save_Result0


subroutine File_Write0
    implicit none
    
    open(2, file='BGP1.txt', status='old', position='append')
        FMT = '(3I10,5F10.4)'
        write(2,FMT) iunit, Model, flag, y, k, L, c, V
    close(2)

    open(2, file='BGP2.txt', status='old', position='append')
        FMT = '(3I10,6F10.4)'
        write(2,FMT) iunit, Model, flag, I_tilde, Eta, Gamma, rr, w, s_L
    close(2)
        
    open(2, file='BGP3.txt', status='old', position='append')
        FMT = '(3I10,9F10.4)'
        write(2,FMT) iunit, Model, flag, w*(epsilon_N+epsilon_I)/y, g, B_N*g_N, B_N*g_I, b_h*g_h, p_N, p_I, epsilon_N, epsilon_I
    close(2)
    
end subroutine File_Write0


subroutine File_Write_SPP0
    implicit none

    open(2, file='BGP_SPP1.txt', status='old', position='append')
        FMT = '(3I10,5F10.4)'
        write(2,FMT) iunit, Model, flag, y, k, L, c, V
    close(2)

    open(2, file='BGP_SPP2.txt', status='old', position='append')
        FMT = '(3I10,6F10.4)'
        write(2,FMT) iunit, Model, flag, I_tilde, Eta, Gamma, rr, w, s_L
    close(2)
    
    open(2, file='BGP_SPP3.txt', status='old', position='append')
        FMT = '(3I10,9F10.4)'
        write(2,FMT) iunit, Model, flag, w*(epsilon_N+epsilon_I)/y, g, B_N*g_N, B_N*g_I, b_h*g_h, p_N, p_I,epsilon_N, epsilon_I
    close(2)
    
end subroutine File_Write_SPP0


subroutine Export_Result0
    implicit none
integer :: j
    
    write(Filename,'("IR",I3,".txt")') iunit
    open(unit=2, file=Filename, status='replace')
        write(2,'(6A10)') 't', 'N', 'I_tilde', 'h','Eta', 'Gamma'
        do j = 1,iter_IR
            write(2,'(6F10.4)') ts_t(j), N_IR(j), I_tilde_IR(j), h_IR(j),      &
                                         Eta_IR(j), Gamma_IR(j)
        end do
    close(2)

    write(Filename,'("IR_outout",I3,".txt")') iunit
    open(unit=2, file=Filename, status='replace')
        write(2,'(6A10)') 't', 'y', 'k', 'L', 'c', 'V'
                                   
        do j = 1,iter_IR
            write(2,'(6F10.4)') ts_t(j), y_IR(j), k_IR(j), L_IR(j), c_IR(j), V_IR(j)
        end do
    close(2)

    write(Filename,'("IR_price",I3,".txt")') iunit
    open(unit=2, file=Filename, status='replace')
        write(2,'(4A10)') 't', 'r', 'w', 's_L'
        do j = 1,iter_IR
            write(2,'(4F10.4)') ts_t(j), rr_IR(j), w_IR(j), s_L_IR(j)
        end do
    close(2)

    write(Filename,'("IR_RD",I3,".txt")') iunit
    open(unit=2, file=Filename, status='replace')
        write(2,'(5A10)') 't', 'p_N', 'p_I', 'eps_N', 'eps_I'
        do j = 1,iter_IR
            write(2,'(5F10.4)') ts_t(j), p_N_IR(j), p_I_IR(j), epsilon_N_IR(j), epsilon_I_IR(j)
        end do
    close(2)
   
    write(Filename,'("IR_growth",I3,".txt")') iunit
    open(unit=2, file=Filename, status='replace')
        write(2,'(5A10)') 't', 'g', 'g_N', 'g_I', 'g_h'
        do j = 1,iter_IR
            write(2,'(5F10.4)') ts_t(j), g_IR(j), B_N*g_N_IR(j), B_N*g_I_IR(j), b_h*g_h_IR(j)
        end do
    close(2)
    
end subroutine Export_Result0


subroutine File_Create1
    implicit none
        
    open(2, file='Simulate.txt', status='replace')
        FMT = '(2A10,13A10)'
        write(2,FMT) 'number', 'Model', 'p', 'sigma', 'A', 'B_NH', 'B_NL', 'b_h',           &       
                                        'mu_N', 'mu_I', 'mu_h', 'alpha_H', 'alpha_L',       &
                                        'lambda_h', 'shock'
    close(2)

    open(2, file='Policy.txt', status='replace')
        FMT = '(A10,5A10)'
        write(2,FMT) 'number', 'tau_k', 'tau_N', 'tau_I', 'tau_hH', 'tau_hL'
    close(2)
    
    open(2, file='BGP1.txt', status='replace')
        FMT = '(12A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'y', 'k', 'L_H', 'L_L', 'c_H', 'c_L', 'V_H', 'V_L', 'V_HL'
    close(2)

    open(2, file='BGP2.txt', status='replace')
        FMT = '(16A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'I_tilde', 'S_tilde', 'h_HL', 'Eta', 'Gamma_H', 'Gamma_L',     &
                     'r', 'w_H', 'w_L', 'w_HL', 's_LH', 's_LL', 's_L'
    close(2)
    
    open(2, file='BGP3.txt', status='replace')
        FMT = '(13A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'RD', 'g', 'g_N', 'g_I', 'g_hH', 'g_hL', 'p_N', 'p_I', 'eps_N', 'eps_I'
    close(2)
    
    open(2, file='BGP_SPP1.txt', status='replace')
        FMT = '(12A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'y', 'k', 'L_H', 'L_L', 'c_H', 'c_L', 'V_H', 'V_L', 'V_HL'
    close(2)

    open(2, file='BGP_SPP2.txt', status='replace')
        FMT = '(16A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'I_tilde', 'S_tilde', 'h_HL', 'Eta', 'Gamma_H', 'Gamma_L',     &
                     'r', 'w_H', 'w_L', 'w_HL', 's_LH', 's_LL', 's_L'
    close(2)
    
    open(2, file='BGP_SPP3.txt', status='replace')
        FMT = '(13A10)'
        write(2,FMT) 'number', 'Model', 'flag',                                     &
                     'RD', 'g', 'g_N', 'g_I', 'g_hH', 'g_hL', 'p_N', 'p_I', 'eps_N', 'eps_I'
    close(2)

    open(2, file='Decompose.txt', status='replace')
        FMT = '(10A10)'
        write(2,FMT) 'number', 'Model', 'flag', 'share', '', 'deepen', '', 'gap', '', 'log(w)'
    close(2)
    
end subroutine File_Create1


subroutine Read_Param1
    implicit none

    rho = moment_data(1)-theta*(moment_data(2)+moment_data(3))
    p = param1(1,iunit); sigma = param1(2,iunit)*(1d0-p)+p 
    A = param1(3,iunit)**(1d0/(1d0-p))
    B_NH = param1(4,iunit); B_NL = param1(5,iunit); b_h = param1(6,iunit)
    mu_N = param1(7,iunit); mu_I = param1(8,iunit); mu_h = param1(9,iunit)
    alpha_H = param1(10,iunit); alpha_L = param1(11,iunit); lambda_h = param1(12,iunit)
    shock = param1(13,iunit)  
    tau_k = param1(14,iunit); tau_N = param1(15,iunit); tau_I = param1(16,iunit)
    tau_hH = param1(17,iunit); tau_hL = param1(18,iunit)
    
    write(*,*) '=============================================================================='
    FMT = '(7A10)'
    write(*,FMT) 'N', 'p', 'sigma', 'A', 'B_NH', 'B_NL', 'b_h'
    FMT = '(I10,6F10.2)'
    write(*,FMT) iunit, p, sigma, A, B_NH, B_NL, b_h
    FMT = '(7A10)'
    write(*,FMT) 'mu_N', 'mu_I', 'mu_h', 'h_HL*', 'alpha_H', 'alpha_L', 'lambda_h'
    FMT = '(7F10.2)'
    write(*,FMT)  mu_N, mu_I, mu_h, h_HL_star, alpha_H, alpha_L, lambda_h
    write(*,*) ''

    open(2, file='Simulate.txt', status='old', position='append')
        FMT = '(2I10,13F10.4)'
        write(2,FMT) iunit, Model, param1(1:13,iunit)
    close(2)

    open(2, file='Policy.txt', status='old', position='append')
        FMT = '(I10,5F10.4)'
        write(2,FMT) iunit, param1(14:18,iunit)
    close(2)

end subroutine Read_Param1


subroutine Save_BGP1(t)
    implicit none
integer, intent(in) :: t

    FMT = '(7A10)'
    write(*,FMT) '', 'y', 'k', 'L_H', 'L_L', 'eps_N', 'eps_I'
    FMT = '(A10,6F10.4)'
    write(*,FMT) '', y, k, L_H, L_L, epsilon_N, epsilon_I
    FMT = '(7A10)'
    write(*,FMT) '', '', 'c_H', 'c_L', 'V_H', 'V_L', 'V_HL'
    FMT = '(2A10,5F10.4)'
    write(*,FMT) '', '', c_H, c_L, V_H, V_L, V
    FMT = '(7A10)'
    write(*,FMT) '', 'I_tilde', 'S_tilde', 'h_HL', 'Eta', 'Gamma_H', 'Gamma_L'
    FMT = '(A10,6F10.4)'
    write(*,FMT) '', I_tilde, S_tilde, h_HL, Eta, Gamma_H, Gamma_L
    FMT = '(7A10)'
    write(*,FMT) '', 'r', 'w_H', 'w_L', 'w_HL', 's_K', 's_L'
    FMT = '(A10,6F10.4)'
    write(*,FMT) '', rr, w_H, w_L, w, s_K, 1d0-s_K
    FMT = '(7A10)'
    write(*,FMT) '', 'RD', 'g', 'g_N', 'g_I', 'g_hH', 'g_hL'
    FMT = '(A10,6F10.4)'
    write(*,FMT) '', w_H*(epsilon_N+epsilon_I)/y, 100d0*g,100d0*B_NH*g_N, 100d0*B_NH*g_I, 100d0*b_h*g_hH, 100d0*b_h*g_hL

    call Save_Result1(t)
    if (SPP == 0) then
        call File_Write1
    elseif (SPP == 1) then
        call File_Write_SPP1
    endif

end subroutine Save_BGP1


subroutine Read_Result1(t)
    implicit none
integer, intent(in) :: t

    N = N_IR(t); I_tilde = I_tilde_IR(t); S_tilde = S_tilde_IR(t)
    h_H = h_H_IR(t); h_L = h_L_IR(t); h_HL =  h_HL_IR(t)   
    k = k_IR(t); k_H = k_H_IR(t); k_L = k_L_IR(t)
    L = L_IR(t); L_H = L_H_IR(t); L_L = L_L_IR(t)
    y = y_IR(t); c = c_IR(t); c_H = c_H_IR(t); c_L = c_L_IR(t)
    V = V_IR(t); V_H = V_H_IR(t); V_L = V_L_IR(t)
    Eta = Eta_IR(t); Gamma = Gamma_IR(t); Gamma_H = Gamma_H_IR(t); Gamma_L = Gamma_L_IR(t)
    s_K = s_K_IR(t); s_L = s_L_IR(t); s_LH = s_LH_IR(t); s_LL = s_LL_IR(t)
    rr = rr_IR(t); w = w_IR(t); w_H = w_H_IR(t); w_L = w_L_IR(t)
    p_N = p_N_IR(t); p_I = p_I_IR(t)
    epsilon_N = epsilon_N_IR(t); epsilon_I = epsilon_I_IR(t)
    g = g_IR(t); g_N = g_N_IR(t); g_I = g_I_IR(t)
    g_h = g_h_IR(t); g_hL = g_hL_IR(t); g_hH = g_hH_IR(t) 
    
end subroutine Read_Result1


subroutine Save_Result1(t)
    implicit none
integer, intent(in) :: t

    N_IR(t) = N; I_tilde_IR(t) = I_tilde; S_tilde_IR(t) = S_tilde
    h_H_IR(t) = h_H; h_L_IR(t) = h_L; h_HL_IR(t) = h_HL
    k_IR(t) = k; k_H_IR(t) = k_H; k_L_IR(t) = k_L
    L_IR(t) = L; L_H_IR(t) = L_H; L_L_IR(t) = L_L 
    y_IR(t) = y; c_IR(t) = c; c_H_IR(t) = c_H; c_L_IR(t) = c_L
    V_IR(t) = V; V_H_IR(t) = V_H; V_L_IR(t) = V_L
    Eta_IR(t) = Eta; Gamma_IR(t) = Gamma; Gamma_H_IR(t) = Gamma_H; Gamma_L_IR(t) = Gamma_L
    s_K_IR(t) = s_K; s_L_IR(t) = s_L; s_LH_IR(t) = s_LH; s_LL_IR(t) = s_LL
    rr_IR(t) = rr; w_IR(t) = w; w_H_IR(t) = w_H; w_L_IR(t) = w_L; 
    p_N_IR(t) = p_N; p_I_IR(t) = p_I
    epsilon_N_IR(t) = epsilon_N; epsilon_I_IR(t) = epsilon_I
    g_IR(t) = g; g_N_IR(t) = g_N; g_I_IR(t) = g_I
    g_h_IR(t) = g_h; g_hL_IR(t) = g_hL; g_hH_IR(t) = g_hH
    
end subroutine Save_Result1


subroutine File_Write1
    implicit none
    
    open(2, file='BGP1.txt', status='old', position='append')
        FMT = '(3I10,9F10.4)'
        write(2,FMT) iunit, Model, flag, y, k, L_H, L_L, c_H, c_L, V_H, V_L, V
    close(2)

    open(2, file='BGP2.txt', status='old', position='append')
        FMT = '(3I10,13F10.4)'
        write(2,FMT) iunit, Model, flag, I_tilde, S_tilde, h_HL, Eta, Gamma_H, Gamma_L,      &
                                         rr, w_H, w_L, w, 1d0-s_K-s_LL, s_LL, 1d0-s_K
    close(2)
        
    open(2, file='BGP3.txt', status='old', position='append')
        FMT = '(3I10,10F10.4)'
        write(2,FMT) iunit, Model, flag, w_H*(epsilon_N+epsilon_I)/y,                        &
                                         100d0*g, 100d0*B_NH*g_N, 100d0*B_NH*g_I, 100d0*b_h*g_hH, 100d0*b_h*g_hL, p_N, p_I, epsilon_N, epsilon_I
    close(2)
    
end subroutine File_Write1


subroutine File_Write_SPP1
    implicit none
    
    open(2, file='BGP_SPP1.txt', status='old', position='append')
        FMT = '(3I10,12F10.4)'
        write(2,FMT) iunit, Model, flag, y, k, L_H, L_L, c_H, c_L, V_H, V_L, V
    close(2)

    open(2, file='BGP_SPP2.txt', status='old', position='append')
        FMT = '(3I10,13F10.4)'
        write(2,FMT) iunit, Model, flag, I_tilde, S_tilde, h_HL, Eta, Gamma_H, Gamma_L,  &
                                         rr, w_H, w_L, w, 1d0-s_K-s_LL, s_LL, 1d0-s_K
    close(2)
        
    open(2, file='BGP_SPP3.txt', status='old', position='append')
        FMT = '(3I10,10F10.4)'
        write(2,FMT) iunit, Model, flag, w_H*(epsilon_N+epsilon_I)/y,                    &
                                         100d0*g, 100d0*B_NH*g_N, 100d0*B_NH*g_I, 100d0*b_h*g_hH, 100d0*b_h*g_hL, p_N, p_I, epsilon_N, epsilon_I
    close(2)
    
end subroutine File_Write_SPP1


subroutine Export_Result1
    implicit none
integer :: j
    
    write(Filename,'("IR",I3,"_",I1,".txt")') iunit, Model
    open(unit=2, file=Filename, status='replace')
        write(2,'(10A10)') 't', 'N', 'I', 'S', 'h_HL', 'h_H', 'h_L',                &
                               'Eta', 'Gamma_H', 'Gamma_L'
        do j = 1,iter_IR
            write(2,'(10F10.4)') ts_t(j), N_IR(j), I_tilde_IR(j), S_tilde_IR(j),    &
                                          h_HL_IR(j), h_H_IR(j), h_L_IR(j),         &
                                          Eta_IR(j), Gamma_H_IR(j), Gamma_L_IR(j)
        end do
    close(2)

    write(Filename,'("IR_outout",I3,"_",I1,".txt")') iunit, Model
    open(unit=2, file=Filename, status='replace')
        write(2,'(10A10)') 't', 'y', 'k', 'L_H', 'L_L', 'c_H', 'c_L', 'V_H', 'V_L', 'V_HL'
                                   
        do j = 1,iter_IR
            write(2,'(10F10.4)') ts_t(j), y_IR(j), k_IR(j), L_H_IR(j), L_L_IR(j),   &
                                          c_H_IR(j), c_L_IR(j),                     &
                                          V_H_IR(j), V_L_IR(j), V_IR(j)
                                         
        end do
    close(2)

    write(Filename,'("IR_price",I3,"_",I1,".txt")') iunit, Model
    open(unit=2, file=Filename, status='replace')
        write(2,'(7A10)') 't', 'r', 'w_H',' w_L', 'w_HL', 's_LH', 's_LL'
        do j = 1,iter_IR
            write(2,'(7F10.4)') ts_t(j), rr_IR(j), w_H_IR(j), w_L_IR(j), w_IR(j),     &
                                         s_LH_IR(j), s_LL_IR(j)
        end do
    close(2)

    write(Filename,'("IR_RD",I3,"_",I1,".txt")') iunit, Model
    open(unit=2, file=Filename, status='replace')
        write(2,'(5A10)') 't', 'p_N', 'p_I', 'eps_N', 'eps_I'
        do j = 1,iter_IR
            write(2,'(5F10.4)') ts_t(j), p_N_IR(j), p_I_IR(j), epsilon_N_IR(j), epsilon_I_IR(j)
        end do
    close(2)
    
    write(Filename,'("IR_growth",I3,"_",I1,".txt")') iunit, Model
    open(unit=2, file=Filename, status='replace')
        write(2,'(6A10)') 't', 'g', 'g_N', 'g_I', 'g_hH', 'g_hL'
        do j = 1,iter_IR
            write(2,'(6F10.4)') ts_t(j), 100d0*g_IR(j), 100d0*B_NH*g_N_IR(j), 100d0*B_NH*g_I_IR(j),  &
                                         100d0*b_h*g_hH_IR(j), 100d0*b_h*g_hL_IR(j)
        end do
    close(2)
        
end subroutine Export_Result1

end module Global_Variable

