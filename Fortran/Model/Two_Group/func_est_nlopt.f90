!
!   func_est_nlopt.f90
!   Growth_Model
!
!   Created by Rong Fan on 4/22/20.
!   Copyright 2020 FanRong. All rights reserved.
!
subroutine func_est_nlopt(L_N,Nx,x,grad,need_gradient,f_data)
use Global_Variable
use Simulation
    implicit none
integer, intent(in) :: need_gradient
integer, intent(in) :: Nx
real(8), intent(in) :: x(Nx), grad(Nx)
real(8), intent(out) :: L_N
real(8), intent(inout) :: f_data(*)


	if (need_gradient .ne. 0) then
        STOP 'ERROR: nloptfcn: gradient requested, but have not programmed this'
    endif

if (calibration == 1) then 
    
    call Read_Param
    call BGP(0)
    rr0 = rr; g_N0 = g_N; g_h0 = g_h
    s_K0 = s_K; w0 = w; w_H0 = w_H; y0 = y
    epsilon_N0 = epsilon_N; epsilon_I0 = epsilon_I
    ll_H0 = ll_H; ll_L0 = ll_L
    
    mu_I = shock*mu_I
    call BGP(0)
    rr1 = rr; g_N1 = g_N; g_h1 = g_h
    s_K1 = s_K; w1 = w; w_H1 = w_H; y1 = y
    epsilon_N1 = epsilon_N; epsilon_I1 = epsilon_I
    ll_H1 = ll_H; ll_L1 = ll_L

    moment_model(1) = rr0                                       ! interest rate
    moment_model(2) = B_NH*g_N0                                 ! growth rate (RD 1980)
    moment_model(3) = b_h*g_h0                                  ! growth rate (human capital 1980)
    moment_model(4) = 1d0-s_K0                                  ! labor share (1980)
    moment_model(5) = w0                                        ! wage premium (1980)
    moment_model(6) = w_H0/y0*(epsilon_N0+epsilon_I0)           ! RD/GDP ratio (1980)
    moment_model(7) = 1d0-s_K1                                  ! labor share (2005)
    moment_model(8) = w1                                        ! wage premium (1980)
    moment_model(9) = log(1d0-ll_H1)-log(1d0-ll_H0)             ! Change of training time (skilled)
    moment_model(10) = log(1d0-ll_L1)-log(1d0-ll_L0)            ! Change of training time (unskilled)

    L_N = 1d6 
    if (flag == 1) L_N = sum(moment_weight*(moment_model-moment_data)**2d0)
    if (L_N < L_N0) then 

        call Read_Param
        write(*,*) '=============================================================================='
        FMT = '(7A10)'
        write(*,FMT) 'N', 'p', 'sigma', 'A', 'B_NH', 'B_NL', 'b_h'
        FMT = '(I10,6F10.2)'
        write(*,FMT) iunit, p, sigma, A, B_NH, B_NL, b_h
        FMT = '(7A10)'
        write(*,FMT) 'mu_N', 'mu_I', 'mu_h', 'alpha_H', 'alpha_L','lambda_h', 'shock'
        FMT = '(7F10.2)'
        write(*,FMT)  mu_N, mu_I, mu_h, alpha_H, alpha_L, lambda_h, shock
        write(*,*) ''

        call BGP(0)
        call Save_BGP1(1)
        write(*,*) ''  

        mu_I = shock*mu_I
        call BGP(0)
        call Save_BGP1(iter_IR)

        call File_Write
        open(2, file='Simulate.txt', status='old', position='append')
            FMT = '(2I10,13F10.4)'
            write(2,FMT) iunit, Model, param1(1:13,iunit)
        close(2)

        open(2, file='Policy.txt', status='old', position='append')
            FMT = '(I10,5F10.4)'
            write(2,FMT) iunit, param1(14:18,iunit)
        close(2)

        L_N0 = L_N
        iunit = iunit+1

    endif 

elseif (calibration == 0) then
    select case (Model)
        case (1)          
            g_N = x(1)
            g_I = x(2)
            g_hH = g_hH0
            g_hL = g_hL0
            g = B_NH*g_N+b_h*g_hH 
        case (2)
            if (Step == 1) then 
                g_N = x(1)
                g_I = x(2)
                g = B_NH*g_N+b_h*g_hH 
            elseif (Step == 2) then 
                g_hH = x(1)
                g_hL = g_hL0
            endif           
        case (3) 
            if (Step == 1) then 
                g_N = x(1)
                g_I = x(2)
                g = B_NH*g_N+b_h*g_hH 
            elseif (Step == 2) then 
                g_hH = x(1)
                g_hL = x(2)
            endif    
    end select 

        epsilon_N = (g_N*mu_N)**(1d0/lambda)
        epsilon_I = (g_I*mu_I)**(1d0/lambda)
        mu_hH = mu_hH_h(h_HL)
        ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
        if (calibration == 0 .and. Task == 2) then 
            L_H = epsilon_H*ll_H-epsilon_N-epsilon_I
        else
            L_H = epsilon_H*ll_H
        endif
        mu_HL = mu_h
        ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
        L_L = epsilon_L*ll_L

        d_Eta_I = 1d0/I_tilde/(sigma-1d0)
        d_Gamma_I = B_NL/(1d0-exp(B_NL*(S_tilde-I_tilde)*(sigma-1d0)))
        d_Gamma_S = B_NL/(1d0-exp(B_NL*(I_tilde-S_tilde)*(sigma-1d0)))
        d_I = s_K/(1-p)/sigma
        a_S = 1d0/((sigma-1d0)/sigma*(s_LL+s_LH)/s_LH*d_Gamma_S+B_NH-B_NL)

        ! With interection between human capital 
        ! g_wN = B_NH*g_N+s_LL/(s_LL+s_LH)*d_I*b_h*(g_hH-g_hL)+d_I*(g_k-g_L)
        g_wN = B_NH*g_N+s_LL/(s_LL+s_LH)*d_I*b_h*(g_hH-g_hL)
        g_wI = (d_I*d_Eta_I+s_LL/(s_LL+s_LH)*(1d0-d_I)*d_Gamma_I)*(g_I-g_N)
        g_wS = (B_NH-B_NL)/a_S*((g_LL-g_LH+b_h*(g_hL-g_hH))/sigma-(sigma-1d0)/sigma*d_Gamma_I*(g_I-g_N))

        if (SPP == 0) then  
            ! Solve patent value sequentially      
            ! V_NH0 = (p*(w_H/A)**(1d0-sigma)*dt+V_NH1)/((rr-g+(sigma-1d0)*B_NH*g_N)*dt+1d0)
            ! V_NL0 = (p*(w_L/A*gamma_HL)**(1d0-sigma)*dt+V_NL1)/((rr-g+(sigma-1d0)*B_NL*g_N)*dt+1d0)
            ! V_II0 = (p*(R/A)**(1d0-sigma)*dt+V_II1)/((rr-g)*dt+1d0)

            ! Solve current patent value
            V_NH0 = p*(w_H/A)**(1d0-sigma)/(rr-g+(sigma-1d0)*B_NH*g_N)
            V_NL0 = p*(w_L/A*gamma_HL)**(1d0-sigma)/(rr-g+(sigma-1d0)*B_NL*g_N)
            V_II0 = p*(R/A)**(1d0-sigma)/(rr-g)
            gamma_IL = exp(B_NL*(1d0-I_tilde)*(1d0-sigma))
            gamma_SH = exp(B_NH*(1d0-S_tilde)*(1d0-sigma))
            gamma_SL = exp(B_NL*(1d0-S_tilde)*(1d0-sigma))
            t_tilde = i-(1d0-S_tilde)/g_N
            beta_S = exp(-(rr-g)*(1d0-S_tilde)/g_N)
            if (t_tilde .ge. 1.5d0) then 
                V_NH = ts_V_NH0(t_tilde)
                V_NL = ts_V_NL0(t_tilde)
            else 
                V_NH = ts_V_NH0(1)
                V_NL = ts_V_NL0(1)
            endif 
            p_N = (1d0-tau_N)*(V_NH0+beta_S*(gamma_SL*V_NL-gamma_SH*V_NH))-(1d0-tau_I)*V_II0
            p_I = (1d0-tau_I)*V_II0-(1d0-tau_N)*gamma_IL*V_NL0
        elseif (SPP == 1) then 
            g_wH = g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS+b_h*g_hH
            V_NH0 = (1-p)/(sigma-1d0)*(w_H/A)**(1d0-sigma)*y/(R/(1d0-p)-delta-g_wH)
            V_NL0 = (1-p)/(sigma-1d0)*(w_L/A*gamma_HL)**(1d0-sigma)*y/(R/(1d0-p)-delta-g_wH)
            V_II0 = (1-p)/(sigma-1d0)*(R/A)**(1d0-sigma)*y/(R/(1d0-p)-delta-g_wH)
            p_N = s_LL*(B_NH-B_NL)/(R/(1d0-p)-delta-g_wH)+V_NH0-V_II0
            p_I = V_II0-gamma_IL*V_NL0
        endif 

    select case (Model)
        case (1)          
            L_N = (lambda*p_N*epsilon_N**(lambda-1d0)/mu_N-s_LH/L_H)**2d0                   &
                + (lambda*p_I*epsilon_I**(lambda-1d0)/mu_I-s_LH/L_H)**2d0
        case (2)
            if (SPP == 0 .and. Step == 1) then
                L_N = (lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H)**2d0               &
                    + (lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H)**2d0 
            elseif (SPP == 0 .and. Step == 2) then
                L_N = (b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H+b_h*g_hH            &
                    - (1d0+tau_hH)*rr+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS)**2d0  
            elseif (SPP == 1 .and. Step == 1) then 
                L_N = (lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H)**2d0               &
                    + (lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H)**2d0          
            elseif (SPP == 1 .and. Step == 2) then 
                L_N = (b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_L/epsilon_H            &
                    - R/(1d0-p)-delta+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS+b_h*g_hH)**2d0 
            endif 
        case (3) 
            if (SPP == 0 .and. Step == 1) then
                L_N = (lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H)**2d0               &
                    + (lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H)**2d0              
            elseif (SPP == 0 .and. Step == 2) then
                L_N = (b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H+b_h*g_hH            &
                    - (1d0+tau_hH)*rr+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS)**2d0                 &
                    + (b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*ll_L+b_h*g_hL            &
                    - (1d0+tau_hL)*rr+g_wN+g_wI+s_LL/(s_LH+s_LL)*g_wS)**2d0  
            elseif (SPP == 1 .and. Step == 1) then 
                L_N = (lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H)**2d0               &
                    + (lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H)**2d0    
            elseif (SPP == 1 .and. Step == 2) then 
                L_N = (b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_L/epsilon_H            &
                    - R/(1d0-p)-delta+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS+b_h*g_hH)**2d0        &
                    + (b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*ll_L                     &
                    - R/(1d0-p)-delta+g_wN+g_wI+s_LL/(s_LH+s_LL)*g_wS+b_h*g_hL)**2d0    
            endif 
    end select 
endif 

contains

subroutine Read_Param
    implicit none

    rho = moment_data(1)-theta*(moment_data(2)+moment_data(3))
    p = x(1); sigma = sigma0*(1d0-p)+p 
    A = x(2)**(1d0/(1d0-p))
    B_NH = x(3); B_NL = 1d0; b_h = 1d0
    mu_N = x(4); mu_I = x(5); mu_h = x(6)
    alpha_H = x(7); alpha_L = x(8)
    lambda_h = x(9); shock = x(10)
    param1([1,3:4,7:13],iunit) = x
    param1(2,iunit) = sigma0
    param1(5,iunit) = B_NL 
    param1(6,iunit) = b_h

end subroutine Read_Param


subroutine File_Write
    implicit none

    open(2, file='L_N.txt', status='old', position='append')
        FMT = '(I10,11F10.4)'
        write(2,FMT) iunit, L_N, moment_model
    close(2)
    
end subroutine File_Write

end subroutine
