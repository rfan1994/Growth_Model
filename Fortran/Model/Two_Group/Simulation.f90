!   Simulation.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
module Simulation
use Global_Variable
use My_Function
use Burkardt_fsolve

    implicit none
integer :: i 
    
contains


! ======================================================================
! Static Equilibirum
! ======================================================================

subroutine Static(SPP0)
    implicit none
integer, intent(in) :: SPP0
real(8) :: x(2), fvec(2)

    SPP = SPP0
    do i = 1,99
        I_tilde = 0.01*i  
        call fsolve(Static_Equation,1,x,fvec,TOL,flag)
        call Save_BGP1(1)
    enddo 
    
contains

subroutine Static_Equation(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)

    S_tilde = x(1)
    if (S_tilde < I_tilde .or. S_tilde > 1d0) then
        fvec = -1d6 
        return 
    endif

    call Firm_Problem_R
    call RD_Problem

    fvec(1) = w-exp((B_NH-B_NL)*S_tilde+b_h*h_HL)

end subroutine Static_Equation

endsubroutine Static


! ======================================================================
! BGP
! ======================================================================

subroutine BGP(SPP0)
    implicit none
integer, intent(in) :: SPP0
real(8) :: x(5), fvec(5)

    SPP = SPP0
    if (rho+delta > A) then 
        select case (Model)
            case (1)
                do i = 1,9
                    x(1) = 0.1d0*i; x(2) = 0.1d0*(i+10d0)/2d0          
                    x(3) = moment_data(2)/B_NH
                    call fsolve(BGP_Equation1,3,x,fvec,TOL,flag)
                    if (flag == 1) exit
                enddo 
            case (2)
                do i = 1,9
                    x(1) = 0.1d0*i; x(2) = 0.1d0*(i+10d0)/2d0          
                    x(3) = moment_data(2)/B_NH
                    x(4) = h_HL_star
                    call fsolve(BGP_Equation2,4,x,fvec,TOL,flag)
                    if (flag == 1) exit
                enddo 
            case (3)   
                do i = 1,9
                    x(1) = 0.1d0*i; x(2) = 0.1d0*(i+10d0)/2d0           
                    x(3) = moment_data(2)/B_NH; x(4) = moment_data(3)/b_h
                    x(5) = h_HL_star
                    call fsolve(BGP_Equation3,5,x,fvec,TOL,flag)
                    if (flag == 1) exit
                enddo 
        end select

        N = 1d0
        k = w_H*L_H/R*s_K/S_LH 
        call Firm_Problem     
        if (SPP == 0) then
            k_H = a_HL*k/(a_HL*epsilon_H+epsilon_L)
            k_L = k/(a_HL*epsilon_H+epsilon_L)            
            Trans = tau_N*(s_LH+s_LL)*y+tau_I*s_K*y                      &
                  + tau_hH*w_H*(1d0-ll_H)*epsilon_H+tau_hL*w_L*(1d0-ll_L)*epsilon_L
            c_H = (rr-g)*k_H+(w_H*L_H+pi)/epsilon_H-tau_hH*w_H*(1d0-ll_H)+Trans
            c_L = (rr-g)*k_L+w_L*L_L/epsilon_L-tau_hL*w_L*(1d0-ll_L)+Trans 
            c = epsilon_H*c_H+epsilon_L*c_L
        elseif (SPP == 1) then 
            c = y-(delta+g)*k
            c_H = c 
            c_L = c
        endif 
        V_H = u(c_H)/(rho-(1d0-theta)*g)
        V_L = u(c_L)/(rho-(1d0-theta)*g)
        V = (V_H/V_L)**(1d0/(1d0-theta))-1d0

    elseif (rho+delta .le. A) then 
        I_tilde = 1d0
        N = 1d0 
        k = 1d0 
        call Firm_Problem
        if (SPP == 0) then
            k_H = a_HL*k/(a_HL*epsilon_H+epsilon_L)
            k_L = k/(a_HL*epsilon_H+epsilon_L)
            Trans = tau_N*(s_LH+s_LL)*y+tau_I*s_K*y                      &
                  + tau_hH*w_H*(1d0-ll_H)*epsilon_H+tau_hL*w_L*(1d0-ll_L)*epsilon_L
            c_H = (rr-g)*k_H+(w_H*L_H+pi)/epsilon_H-tau_hH*w_H*(1d0-ll_H)+Trans
            c_L = (rr-g)*k_L+w_L*L_L/epsilon_L-tau_hL*w_L*(1d0-ll_L)+Trans 
            c = epsilon_H*c_H+epsilon_L*c_L        
        elseif (SPP == 1) then
            g = (R/(1-p)-delta-rho)/theta
            c = y-(delta+g)*k
            c_H = c 
            c_L = c
        endif 
        V_H = u(c_H)/(rho-(1d0-theta)*g)
        V_L = u(c_L)/(rho-(1d0-theta)*g)
        V = (V_H/V_L)**(1d0/(1d0-theta))-1d0
    endif 
    
contains

subroutine BGP_Equation1(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)

    I_tilde = x(1)
    S_tilde = x(2)
    g_N = x(3)
    g_I = g_N
    g_hH = g_hH0
    g_hL = g_hL0
    h_HL = h_HL0
    g = B_NH*g_N+b_h*g_hH

    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    L_H = epsilon_H*ll_H
    mu_hL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    if (I_tilde < 0d0 .or. I_tilde > 1d0 .or.           &
        S_tilde < I_tilde .or. S_tilde > 1d0 .or.       &
        epsilon_N < 0d0 .or. L_H < 0d0) then
        fvec = -1d6 
        return 
    endif

    if (SPP == 0) then
        rr = rho+theta*g
        R = rr+delta
    elseif (SPP == 1) then
        rr = rho+theta*g
        R = (rr+delta)*(1-p)
    endif 
    call Firm_Problem_R
    call RD_Problem

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H
    if (I_tilde < TOL) fvec(2) = max(fvec(2),0d0)
    fvec(3) = w-exp((B_NH-B_NL)*S_tilde+b_h*h_HL)
    if (I_tilde > 1d0-TOL) fvec(3) = max(fvec(3),0d0)

end subroutine BGP_Equation1


subroutine BGP_Equation2(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)

    I_tilde = x(1)
    S_tilde = x(2)
    g_N = x(3)
    g_I = g_N
    g_hH = g_hH0
    g_hL = g_hL0
    h_HL = x(4)
    g = B_NH*g_N+b_h*g_hH

    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    L_H = epsilon_H*ll_H
    mu_hL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    if (I_tilde < 0d0 .or. I_tilde > 1d0 .or.           &
        S_tilde < I_tilde .or. S_tilde > 1d0 .or.       &
        ll_H < 0d0 .or. ll_H > 1d0 .or.                 &
        epsilon_N < 0d0 .or. L_H < 0d0) then
        fvec = -1d6 
        return 
    endif

    if (SPP == 0) then
        rr = rho+theta*g
        R = rr+delta
    elseif (SPP == 1) then
        rr = rho+theta*g
        R = (rr+delta)*(1-p)
    endif 
    call Firm_Problem_R
    call RD_Problem

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H
    if (I_tilde < TOL) fvec(2) = max(fvec(2),0d0)
    fvec(3) = w-exp((B_NH-B_NL)*S_tilde+b_h*h_HL)
    if (I_tilde > 1d0-TOL) fvec(3) = max(fvec(3),0d0)
    if (SPP == 0) then
        fvec(4) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H      &
                + b_h/mu_hH*(1d0-ll_H)**alpha_H                         &
                - (1d0+tau_hH)*rr+B_NH*g_N   
    elseif (SPP == 1) then 
        fvec(4) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_H/epsilon_H-rr+g     
    endif 
        
end subroutine BGP_Equation2


subroutine BGP_Equation3(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)

    I_tilde = x(1)
    S_tilde = x(2)
    g_N = x(3)
    g_I = g_N
    g_hH = x(4)
    g_hL = g_hH
    h_HL = x(5)
    g = B_NH*g_N+b_h*g_hH

    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    if (calibration == 0 .and. Task == 3) then 
        L_H = epsilon_H*ll_H
    else 
        L_H = epsilon_H*ll_H-epsilon_N-epsilon_I
    endif 
    mu_hL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    if (I_tilde < 0d0 .or. I_tilde > 1d0 .or.           &
        S_tilde < I_tilde .or. S_tilde > 1d0 .or.       &
        ll_H < 0d0 .or. ll_H > 1d0 .or.                 &
        ll_L < 0d0 .or. ll_L > 1d0 .or.                 &
        epsilon_N < 0d0 .or. L_H < 0d0) then
        fvec = -1d6 
        return 
    endif

    if (SPP == 0) then
        rr = rho+theta*g
        R = rr+delta
    elseif (SPP == 1) then
        rr = rho+theta*g
        R = (rr+delta)*(1-p)
    endif 
    call Firm_Problem_R
    call RD_Problem

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H
    if (I_tilde < TOL) fvec(2) = max(fvec(2),0d0)
    fvec(3) = w-exp((B_NH-B_NL)*S_tilde+b_h*h_HL)
    if (I_tilde > 1d0-TOL) fvec(3) = max(fvec(3),0d0)
    if (SPP == 0) then
        fvec(4) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H      &
                + b_h/mu_hH*(1d0-ll_H)**alpha_H                         &
                - (1d0+tau_hH)*rr+B_NH*g_N                                                    
        fvec(5) = b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*ll_L      &
                + b_h/mu_hL*(1d0-ll_L)**alpha_L                         &
                - (1d0+tau_hL)*rr+B_NH*g_N  
    elseif (SPP == 1) then
        fvec(4) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_H/epsilon_H-rr+g     
        fvec(5) = b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*L_L/epsilon_L-rr+g     
    endif
                                             
end subroutine BGP_Equation3


end subroutine BGP


! ======================================================================
! Transition
! ======================================================================


subroutine Transition
    implicit none
integer :: iter
real(8) :: kmin, kmax, k_dot, g_dot
integer :: p0, p1

    ! Read
    I0 = I_tilde_IR(1); I1 = I_tilde_IR(iter_IR)
    k0 = k_IR(1); k1 = k_IR(iter_IR) 
    h_HL0 = h_HL_IR(1); h_HL1 = h_HL_IR(iter_IR) 
    g_N0 = g_N_IR(1); g_N1 = g_N_IR(iter_IR)
    g_I0 = g_I_IR(1); g_I1 = g_I_IR(iter_IR)
	g_hH0 = g_hH_IR(1); g_hH1 = g_hH_IR(iter_IR)
    g_hL0 = g_hL_IR(1); g_hL1 = g_hL_IR(iter_IR)
    L_H0 = L_H_IR(1); L_H1 = L_H_IR(iter_IR)
    L_L0 = L_L_IR(1); L_L1 = L_L_IR(iter_IR)

    ! Initialize
    dt = T**lambda_T/real(iter_IR-1,8)
    do i = 1,iter_IR
        ts_t(i) = (dt*real(i-1,8))**(1d0/lambda_T)
    enddo
    do i = 1,iter_IR-1
        ts_dt(i) = ts_t(i+1)-ts_t(i)
    enddo
    ts_dt(iter_IR) = ts_dt(iter_IR-1)

    ! State 
    ts_I0 = I1; ts_I1 = I1
    ts_k0 = k1; ts_k1 = k1
    ts_h_HL0 = h_HL1; ts_h_HL1 = h_HL1
    call Growth_Rate_T
    
    do iter = 1,iter_max
		do i = 1,(iter_IR-1)
            ts_g_k(i) = (ts_k0(i+1)-ts_k0(i))/ts_k0(i)/ts_dt(i)
            ts_g_LH(i) = (ts_L_H0(i+1)-ts_L_H0(i))/ts_L_H0(i)/ts_dt(i)
            ts_g_LL(i) = (ts_L_L0(i+1)-ts_L_L0(i))/ts_L_L0(i)/ts_dt(i)
        enddo
        ts_g_k(iter_IR) = 0d0
        ts_g_LH(iter_IR) = 0d0
        ts_g_LL(iter_IR) = 0d0
 
        ! HH problem
        call Policy_Function
        ! RD problem
        call Growth_Rate  

        ! Update I
        ts_I1(1) = I0
        do i = 2,iter_IR-1
            g_dot = max(ts_g_I1(i-1)-ts_g_N1(i-1),0d0)
            I_tilde = ts_I0(i-1)+g_dot*ts_dt(i-1)
            ts_I1(i) = I_tilde
            if (I_tilde .le. min(I0,I1)) then
                ts_I1(i) = min(I0,I1)
            elseif (I_tilde .ge. max(I0,I1)) then
                ts_I1(i) = max(I0,I1)
            else
                ts_I1(i) = I_tilde
            endif
        enddo
        
        ! Update k
        ts_k1(1) = k0
        do i = 2,iter_IR-1
            g_dot = min(ts_a_dk(i-1),0d0)
            k = ts_k0(i-1)+g_dot*ts_dt(i-1)
            if (k .le. min(k0,k1)) then
                ts_k1(i) = min(k0,k1)
            elseif (k .ge. max(k0,k1)) then
                ts_k1(i) = max(k0,k1)
            else
                ts_k1(i) = k
            endif
        enddo

        ! Update h_HL
        ts_h_HL1(1) = h_HL0
        do i = 2,iter_IR-1
            g_dot = ts_g_hH1(i-1)-ts_g_hL1(i-1)
            h_HL = ts_h_HL0(i-1)+g_dot*ts_dt(i-1)  
            if (h_HL .le. min(h_HL0,h_HL1)) then
                ts_h_HL1(i) = min(h_HL0,h_HL1)
            elseif (k .ge. max(h_HL0,h_HL1)) then
                ts_h_HL1(i) = max(h_HL0,h_HL1)
            else
                ts_h_HL1(i) = h_HL
            endif    
        enddo

        ! Check convergence
        sup = max(maxval(abs(ts_I1-ts_I0)),         &
                  maxval(abs(ts_k1-ts_k0)),         &    
                  maxval(abs(ts_h_HL1-ts_h_HL0)),   &                    
                  maxval(abs(ts_g_N1-ts_g_N0)),     &
                  maxval(abs(ts_g_I1-ts_g_I0)),     &
                  maxval(abs(ts_g_hH1-ts_g_hH0)),   &
                  maxval(abs(ts_g_hL1-ts_g_hL0)))
        print*, iter, sup
        if (sup<TOL1) exit

        ts_I0 = ts_I0+rho_ts*(ts_I1-ts_I0)
        ts_k0 = ts_k0+rho_ts*(ts_k1-ts_k0)
        ts_h_HL0 = ts_h_HL0+rho_ts*(ts_h_HL1-ts_h_HL0)
        ts_g_N0 = ts_g_N0+rho_ts*(ts_g_N1-ts_g_N0)
        ts_g_I0 = ts_g_I0+rho_ts*(ts_g_I1-ts_g_I0)
        ts_g_hH0 = ts_g_hH0+rho_ts*(ts_g_hH1-ts_g_hH0)
        ts_g_hL0 = ts_g_hL0+rho_ts*(ts_g_hL1-ts_g_hL0)
        ts_L_H0 = ts_L_H0+rho_ts*(ts_L_H1-ts_L_H0)
        ts_L_L0 = ts_L_L0+rho_ts*(ts_L_L1-ts_L_L0)
        ts_V_NH0 = ts_V_NH0+rho_ts*(ts_V_NH1-ts_V_NH0)
        ts_V_NL0 = ts_V_NL0+rho_ts*(ts_V_NL1-ts_V_NL0)
        ts_V_II0 = ts_V_II0+rho_ts*(ts_V_II1-ts_V_II0)
    enddo

    N_IR(1) = 1d0; h_H_IR = h_HL0; h_L_IR = 0d0
    call Value_Function
    do i = 2,iter_IR-1
        I_tilde = ts_I0(i)
        k = ts_k0(i)
        h_HL = ts_h_HL0(i) 
        S_tilde = ts_S(i) 
        L_H = ts_L_H0(i)
        L_L = ts_L_L0(i)
        epsilon_N = (g_N*mu_N/(1d0-tau_N))**(1d0/lambda) 
        epsilon_I = (g_I*mu_I/(1d0-tau_I))**(1d0/lambda)
        call Firm_Problem
        call RD_Problem          
        g_N = ts_g_N0(i)
        g_I = ts_g_I0(i)
        g_hH = ts_g_hH0(i)
        g_hL = ts_g_hL0(i)
        g = B_NH*g_N+b_h*g_hH      
        c_H = c_H_IR(i)
        c_L = c_L_IR(i)
        V_H = V_H_IR(i)
        V_L = V_L_IR(i)
        V = V_IR(i)     
        p_N = ts_p_N(i)
        p_I = ts_p_I(i)
        N = N_IR(i-1)+g_N*ts_dt(i-1)
        h_H = h_H_IR(i-1)+g_hH*ts_dt(i-1)
        h_L = h_L_IR(i-1)+g_hL*ts_dt(i-1)
        call Save_Result1(i)
    enddo
    
end subroutine Transition


subroutine Growth_Rate_T
    implicit none
real(8) :: x(4), fvec(4)

    I_tilde = I1
    k = k1
    h_HL = h_HL1
    g_k = 0d0
    g_LH = 0d0
    g_LL = 0d0

    select case (Model)
        case (1)
            x(1) = g_N1
            x(2) = g_I1
            call fsolve(RD_Equation_T1,2,x,fvec,TOL,flag)
        case (2)
            x(1) = g_N1
            x(2) = g_I1
            x(3) = g_hH1
            call fsolve(RD_Equation_T2,3,x,fvec,TOL,flag)
        case (3)           
            x(1) = g_N1
            x(2) = g_I1
            x(3) = g_hH1
            x(4) = g_hL1
            call fsolve(RD_Equation_T3,4,x,fvec,TOL,flag)
    end select

    ts_g_N0 = g_N; ts_g_N1 = g_N
    ts_g_I0 = g_I; ts_g_I1 = g_I
    ts_g_hH0 = g_hH; ts_g_hH1 = g_hH
    ts_g_hL0 = g_hL; ts_g_hL1 = g_hL
    ts_L_H0 = L_H; ts_L_H1 = L_H
    ts_L_L0 = L_L; ts_L_L1 = L_L
    ts_S = S_tilde

    ts_V_NH0 = V_NH; ts_V_NH1 = V_NH
    ts_V_NL0 = V_NL; ts_V_NL1 = V_NL
    ts_V_II0 = V_II; ts_V_II1 = V_II
    ts_p_N = p_N; ts_p_I = p_I
        
contains

subroutine RD_Equation_T1(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)
real(8) :: d_Eta_I, d_Gamma_I, d_Gamma_S, d_I
real(8) :: beta_NS

    g_N = x(1)
    g_I = x(2)
    g_hH = g_hH0
    g_hL = g_hL0
    g = B_NH*g_N+b_h*g_hH
    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    L_H = epsilon_H*ll_H
    mu_HL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    S_tilde = Brent_Root(S_tilde_Equation,I_tilde+TOL,1d0-TOL,TOL,flag)
    call Firm_Problem
    call RD_Problem

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H

end subroutine RD_Equation_T1


subroutine RD_Equation_T2(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)
real(8) :: d_Eta_I, d_Gamma_I, d_Gamma_S, d_I
real(8) :: g_wN, g_wI, g_wS

    g_N = x(1)
    g_I = x(2)
    g_hH = x(3)
    g_hL = g_hL0
    g = B_NH*g_N+b_h*g_hH
    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    L_H = epsilon_H*ll_H
    mu_HL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    S_tilde = Brent_Root(S_tilde_Equation,I_tilde+TOL,1d0-TOL,TOL,flag)  
    call Firm_Problem
    call RD_Problem

    d_Eta_I = 1d0/I_tilde/(sigma-1d0)
    d_Gamma_I = B_NL/(1d0-exp(B_NL*(S_tilde-I_tilde)*(sigma-1d0)))
    d_Gamma_S = B_NL/(1d0-exp(B_NL*(I_tilde-S_tilde)*(sigma-1d0)))
    d_I = s_K/(1-p)/sigma
    a_S = (sigma-1d0)/sigma*(s_LL+s_LH)/s_LH*d_Gamma_S+B_NH-B_NL

    g_L = s_LL/(s_LH+s_LL)*g_LL+s_LH/(s_LH+s_LL)*g_LH
    g_wN = B_NH*g_N+s_LL/(s_LL+s_LH)*d_I*b_h*(g_hH-g_hL)+d_I*(g_k-g_L)
    g_wI = (d_I*d_Eta_I+s_LL/(s_LL+s_LH)*(1d0-d_I)*d_Gamma_I)*(g_I-g_N)
    g_wS = (B_NH-B_NL)/a_S*((g_LL-g_LH+b_h*(g_hL-g_hH))/sigma-(sigma-1d0)/sigma*d_Gamma_I*(g_I-g_N))

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H
    if (SPP == 0) then
        fvec(3) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H                      &
                + b_h/mu_hH*(1d0-ll_H)**alpha_H                                         &
                - (1d0+tau_hH)*rr+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS  
    elseif (SPP == 1) then 
        fvec(3) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_L/epsilon_H             &
                - R/(1d0-p)-delta+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS+b_h*g_hH    
    endif     

end subroutine RD_Equation_T2


subroutine RD_Equation_T3(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)
real(8) :: d_Eta_I, d_Gamma_I, d_Gamma_S, d_I
real(8) :: g_wN, g_wI, g_wS

    g_N = x(1)
    g_I = x(2)
    g_hH = x(3)
    g_hL = x(4)
    g = B_NH*g_N+b_h*g_hH
    epsilon_N = (g_N*mu_N)**(1d0/lambda)
    epsilon_I = (g_I*mu_I)**(1d0/lambda)
    mu_hH = mu_hH_h(h_HL)
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    L_H = epsilon_H*ll_H
    mu_HL = mu_h
    ll_L = 1d0-(g_hL*mu_hL)**(1d0/alpha_L)
    L_L = epsilon_L*ll_L

    S_tilde = Brent_Root(S_tilde_Equation,I_tilde+TOL,1d0-TOL,TOL,flag)  
    call Firm_Problem
    call RD_Problem

    d_Eta_I = 1d0/I_tilde/(sigma-1d0)
    d_Gamma_I = B_NL/(1d0-exp(B_NL*(S_tilde-I_tilde)*(sigma-1d0)))
    d_Gamma_S = B_NL/(1d0-exp(B_NL*(I_tilde-S_tilde)*(sigma-1d0)))
    d_I = s_K/(1-p)/sigma
    a_S = (sigma-1d0)/sigma*(s_LL+s_LH)/s_LH*d_Gamma_S+B_NH-B_NL

    g_L = s_LL/(s_LH+s_LL)*g_LL+s_LH/(s_LH+s_LL)*g_LH
    g_wN = B_NH*g_N+s_LL/(s_LL+s_LH)*d_I*b_h*(g_hH-g_hL)+d_I*(g_k-g_L)
    g_wI = (d_I*d_Eta_I+s_LL/(s_LL+s_LH)*(1d0-d_I)*d_Gamma_I)*(g_I-g_N)
    g_wS = (B_NH-B_NL)/a_S*((g_LL-g_LH+b_h*(g_hL-g_hH))/sigma-(sigma-1d0)/sigma*d_Gamma_I*(g_I-g_N))

    fvec(1) = lambda/mu_N*p_N*epsilon_N**(lambda-1d0)-s_LH/L_H
    fvec(2) = lambda/mu_I*p_I*epsilon_I**(lambda-1d0)-s_LH/L_H
    if (SPP == 0) then
        fvec(3) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*ll_H                      &
                + b_h/mu_hH*(1d0-ll_H)**alpha_H                                         &
                - (1d0+tau_hH)*rr+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS    
        fvec(4) = b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*ll_L                      &
                + b_h/mu_hL*(1d0-ll_L)**alpha_L                                         &
                - (1d0+tau_hL)*rr+g_wN+g_wI+s_LL/(s_LH+s_LL)*g_wS
    elseif (SPP == 1) then 
        fvec(3) = b_h/mu_hH*alpha_H*(1d0-ll_H)**(alpha_H-1d0)*L_H/epsilon_H             &
                - R/(1d0-p)-delta+g_wN+g_wI-s_LH/(s_LH+s_LL)*g_wS+b_h*g_hH    
        fvec(4) = b_h/mu_hL*alpha_L*(1d0-ll_L)**(alpha_L-1d0)*ll_L                      &
                - R/(1d0-p)-delta+g_wN+g_wI+s_LL/(s_LH+s_LL)*g_wS+b_h*g_hL
    endif

end subroutine RD_Equation_T3


real(8) function S_tilde_Equation(S_tilde0)
    implicit none
real(8), intent(in) :: S_tilde0
real(8) :: gamma_HL_S

    S_tilde = S_tilde0
    call Firm_Problem
    gamma_HL_S = exp((B_NH-B_NL)*S_tilde+b_h*h_HL)  
    S_tilde_Equation = w-gamma_HL_S

end function S_tilde_Equation

end subroutine Growth_Rate_T


subroutine Policy_Function
    implicit none 
real(8) :: dt, g_c, c0, c1
real(8) :: c_H0, c_H1, c_L0, c_L1

    c1 = c_IR(iter_IR)
    c_H1 = c_H_IR(iter_IR)
    c_L1 = c_L_IR(iter_IR)
    ts_a_dk(iter_IR) = 0d0

    do i = iter_IR-1,1,-1
        dt = ts_dt(i)
        I_tilde = ts_I0(i+1)
        k = ts_k0(i+1)
        h_HL = ts_h_HL0(i+1)
        L_H = ts_L_H0(i+1)
        L_L = ts_L_L0(i+1) 
        S_tilde = Brent_Root(S_tilde_Equation,I_tilde+TOL,1d0-TOL,TOL,flag)   
        call Firm_Problem

        g_N = ts_g_N0(i)
        g_hH = ts_g_hH0(i)
        g = B_NH*g_N+b_h*g_hH  
        if (SPP == 0) then  
            g_c = (rr-rho)/theta-g 
        elseif (SPP == 1) then 
            g_c = (R/(1-p)-delta-rho)/theta-g 
        endif 
        c0 = c1/(1d0+g_c*dt)
        c_H0 = c_H1/(1d0+g_c*dt)
        c_L0 = c_L1/(1d0+g_c*dt)
        
        I_tilde = ts_I0(i)
        k = ts_k0(i)
        h_HL = ts_h_HL0(i)
        L_H = ts_L_H0(i)
        L_L = ts_L_L0(i)
        S_tilde = Brent_Root(S_tilde_Equation,I_tilde+TOL,1d0-TOL,TOL,flag)
        call Firm_Problem
        ts_a_dk(i) = y-(delta+g)*k-c0

        ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
        ll_L = L_L/epsilon_L
        Trans = rr*k+tau_N*(s_LH+s_LL)*y+tau_I*s_K*y                          &
              + tau_hH*w_H*(1d0-ll_H)*epsilon_H+tau_hL*w_L*(1d0-ll_L)*epsilon_L
        ts_aH_dk(i) = (rr-g)*k_H+(w_H*L_H+pi)/epsilon_H-tau_hH*w_H*(1d0-ll_H)+Trans-c_H
        ts_aL_dk(i) = (rr-g)*k_L+w_L*L_L/epsilon_L-tau_hL*w_L*(1d0-ll_L)+Trans-c_L
        c_IR(i) = c0 
        c_H_IR(i) = c_H0 
        c_L_IR(i) = c_L0 
        ts_S(i) = S_tilde

        ! Update 
        c1 = c0
        c_H1 = c_H0
        c_L1 = c_L0
    enddo

contains

real(8) function S_tilde_Equation(S_tilde0)
    implicit none
real(8), intent(in) :: S_tilde0
real(8) :: gamma_HL_S

    S_tilde = S_tilde0
    call Firm_Problem
    gamma_HL_S = exp((B_NH-B_NL)*S_tilde+b_h*h_HL)  
    S_tilde_Equation = w-gamma_HL_S

end function S_tilde_Equation

end subroutine Policy_Function


subroutine Growth_Rate
    implicit none
real(8), allocatable :: x(:), x_range(:,:)
real(8), parameter :: rmin = 0.8d0, rmax = 1.5d0

    do i = iter_IR-1,1,-1
        I_tilde = ts_I0(i)
        S_tilde = ts_S(i)
        k = ts_k0(i)
        h_HL = ts_h_HL0(i)
        L_H = ts_L_H0(i)
        L_L = ts_L_L0(i) 
        call Firm_Problem

        g_k = ts_g_k(i)
        g_LH = ts_g_LH(i)
        g_LL = ts_g_LL(i)
        g_L = s_LL/(s_LH+s_LL)*g_LL+s_LH/(s_LH+s_LL)*g_LH
        V_NH1 = ts_V_NH1(i+1)
        V_NL1 = ts_V_NL1(i+1)
        V_II1 = ts_V_II1(i+1)  
        
        select case (Model)
            case (1)       
                allocate(x(2),x_range(2,2))   
                x(1) = g_N0; x_range(1,1) = rmin*g_N0; x_range(2,1) = rmax*g_N1 
                x(2) = g_I0; x_range(1,2) = rmin*g_I0; x_range(2,2) = rmax*g_I1 
                Np = 2; flag = 1; call nlopt(6,x,x_range)
                deallocate(x,x_range)
            case (2)
                allocate(x(3),x_range(2,3))   
                x(1) = g_N0; x_range(1,1) = rmin*g_N0; x_range(2,1) = rmax*g_N1 
                x(2) = g_I0; x_range(1,2) = rmin*g_I0; x_range(2,2) = rmax*g_I1
                x(3) = g_hH0; ; x_range(1,3) = rmin*g_hH0; x_range(2,3) = rmax*g_hH1
                Np = 3; flag = 1; call nlopt(6,x,x_range)
                deallocate(x,x_range)
            case (3)  
                allocate(x(4),x_range(2,4))      
                x(1) = g_N0; x_range(1,1) = rmin*g_N0; x_range(2,1) = rmax*g_N1 
                x(2) = g_I0; x_range(1,2) = rmin*g_I0; x_range(2,2) = rmax*g_I1
                x(3) = g_hH0; ; x_range(1,3) = rmin*g_hH0; x_range(2,3) = rmax*g_hH1
                x(4) = g_hL0; ; x_range(1,4) = g_hL0; x_range(2,4) = g_hL1
                Np = 4; flag = 1; call nlopt(6,x,x_range)
                deallocate(x,x_range)
        end select

        ! Updata
        ts_g_N1(i) = g_N
        ts_g_I1(i) = g_I
        ts_g_hH1(i) = g_hH
        ts_g_hL1(i) = g_hL
        ts_L_H1(i) = L_H
        ts_L_L1(i) = L_L

        ts_V_NH1(i) = V_NH0
        ts_V_NL1(i) = V_NL0
        ts_V_II1(i) = V_II0
        ts_p_N(i) = p_N
        ts_p_I(i) = p_I
    enddo

end subroutine Growth_Rate


subroutine Value_Function 
    implicit none 
real(8) :: V_H0, V_H1, V_L0, V_L1

    do i = 2,iter_IR
        k_H_IR(i) = k_H_IR(i-1)+ts_aH_dk(i)*ts_dt(i-1)
        k_L_IR(i) = k_L_IR(i-1)+ts_aL_dk(i)*ts_dt(i-1)
    enddo

    I_tilde = ts_I0(iter_IR)
    S_tilde = ts_S(iter_IR)
    k = ts_k0(iter_IR)
    h_HL = ts_h_HL0(iter_IR)
    L_H = ts_L_H0(iter_IR)
    L_L = ts_L_L0(iter_IR)
    call Firm_Problem

    g_N = ts_g_N0(iter_IR)
    g_hH = ts_g_hH0(iter_IR)
    g = B_NH*g_N+b_h*g_hH 
    ll_H = 1d0-(g_hH*mu_hH)**(1d0/alpha_H)
    ll_L = L_L/epsilon_L
    Trans = tau_N*(s_LH+s_LL)*y+tau_I*s_K*y                      &
          + tau_hH*w_H*(1d0-ll_H)*epsilon_H+tau_hL*w_L*(1d0-ll_L)*epsilon_L
    c_H = (rr-g)*k_H+(w_H*L_H+pi)/epsilon_H-tau_hH*w_H*(1d0-ll_H)+Trans
    c_L = (rr-g)*k_L+w_L*L_L/epsilon_L-tau_hL*w_L*(1d0-ll_L)+Trans 
    V_H_IR(iter_IR) = u(c_H)/(rho-(1d0-theta)*g)
    V_L_IR(iter_IR) = u(c_L)/(rho-(1d0-theta)*g)
    V_IR(iter_IR) = (u(c_H)/u(c_L))**(1d0/(1d0-theta))-1d0

    V_H1 = V_H_IR(iter_IR)
    V_L1 = V_L_IR(iter_IR)
    do i = iter_IR-1,1,-1  
        dt = ts_dt(i) 
        c_H = c_H_IR(i)  
        c_L = c_L_IR(i) 
        g_N = ts_g_N0(i)
        g_hH = ts_g_hH0(i)
        g_hL = ts_g_hL0(i)
        g = B_NH*g_N+b_h*g_hH     
        V_H0 = (u(c_H)*dt+V_H1)/((rho-(1d0-theta)*g)*dt+1d0)
        V_L0 = (u(c_L)*dt+V_L1)/((rho-(1d0-theta)*g)*dt+1d0)
        V_H_IR(i) = V_H0
        V_L_IR(i) = V_L0
        V_IR(i) = (V_H0/V_L0)**(1d0/(1d0-theta))-1d0

        ! Update
        V_H1 = V_H0
        V_L1 = V_L0
    enddo 

end subroutine Value_Function


subroutine Firm_Problem_R
    implicit none

    ! Productivity
    Eta = I_tilde**(1d0/(sigma-1d0))
    gamma_HL = exp(B_NH-B_NL+b_h*h_HL)
    Gamma_H = ((1d0-exp(-B_NH*(1d0-S_tilde)*(sigma-1d0)))/(B_NH*(sigma-1d0)))**(1d0/(sigma-1d0))   
    Gamma_L = ((exp(-B_NL*(1d0-S_tilde)*(sigma-1d0))-exp(-B_NL*(1d0-I_tilde)*(sigma-1d0)))/(B_NL*(sigma-1d0)))**(1d0/(sigma-1d0))/gamma_HL

    ! Price
    w = (Gamma_H/Gamma_L)**((sigma-1d0)/sigma)*(L_L/L_H)**(1d0/sigma)
    s_K = (1d0-p)*(R/A/Eta)**(1d0-sigma)
    s_LH = ((1d0-p)-s_K)*w*L_H/(L_L+w*L_H)
    s_LL = ((1d0-p)-s_K)*L_L/(L_L+w*L_H)
    w_H = A*Gamma_H*(s_LH/(1d0-p))**(1d0/(1d0-sigma))
    w_L = A*Gamma_L*(s_LL/(1d0-p))**(1d0/(1d0-sigma))

end subroutine Firm_Problem_R


subroutine Firm_Problem
    implicit none
real(8) :: d_Eta_I, d_Gamma_I, d_Gamma_S

    ! Productivity
    Eta = I_tilde**(1d0/(sigma-1d0))
    gamma_HL = exp(B_NH-B_NL+b_h*h_HL)
    Gamma_H = ((1d0-exp(-B_NH*(1d0-S_tilde)*(sigma-1d0)))/(B_NH*(sigma-1d0)))**(1d0/(sigma-1d0))   
    Gamma_L = ((exp(-B_NL*(1d0-S_tilde)*(sigma-1d0))-exp(-B_NL*(1d0-I_tilde)*(sigma-1d0)))/(B_NL*(sigma-1d0)))**(1d0/(sigma-1d0))/gamma_HL

    if (I_tilde .le. TOL) then
        ! Output
        k = 0d0
        y = A/(1d0-p)*((Gamma_L*L_L)**((sigma-1d0)/sigma)+(Gamma_H*L_H)**((sigma-1d0)/sigma))**(sigma/(sigma-1d0))
        R = 0d0
        rr = rho+theta*g
        w_H = A*Gamma_H*(y*(1d0-p)/(A*Gamma_H*L_H))**(1d0/sigma)
        w_L = A*Gamma_L*(y*(1d0-p)/(A*Gamma_L*L_L))**(1d0/sigma)
        w = w_H/w_L
        s_K = R*k/y
        s_LH = w_H*L_H/y
        s_LL = w_L*L_L/y 
        pi = (1d0-tau_I)*p/(1d0-p)*s_K*y+(1d0-tau_N)*p/(1d0-p)*(s_LH+s_LL)*y
    elseif (I_tilde .ge. 1d0-TOL) then
        ! Output
        k = 1d0
        y = A/(1d0-p)*Eta*k
        R = A*Eta
        if (SPP == 0) then 
            rr = R-delta
        elseif (SPP == 1) then 
            rr = R/(1d0-p)-delta
        endif 
        w_H = 0d0
        w_L = 0d0
        w = 0d0
        s_K = R*k/y
        s_LH = w_H*L_H/y
        s_LL = w_L*L_L/y 
        pi = (1d0-tau_I)*p/(1d0-p)*s_K*y+(1d0-tau_N)*p/(1d0-p)*(s_LH+s_LL)*y
    elseif (TOL < I_tilde .and. I_tilde < 1d0-TOL) then
        ! Output
        y = A/(1d0-p)*((Eta*k)**((sigma-1d0)/sigma)+(Gamma_L*L_L)**((sigma-1d0)/sigma)+(Gamma_H*L_H)**((sigma-1d0)/sigma))**(sigma/(sigma-1d0))
        R = A*Eta*(y*(1d0-p)/(A*Eta*k))**(1d0/sigma)
        if (SPP == 0) then 
            rr = R-delta
        elseif (SPP == 1) then 
            rr = R/(1d0-p)-delta
        endif 
        w_H = A*Gamma_H*(y*(1d0-p)/(A*Gamma_H*L_H))**(1d0/sigma)
        w_L = A*Gamma_L*(y*(1d0-p)/(A*Gamma_L*L_L))**(1d0/sigma)
        w = w_H/w_L
        s_K = R*k/y
        s_LH = w_H*L_H/y
        s_LL = w_L*L_L/y
        pi = (1d0-tau_I)*p/(1d0-p)*s_K*y+(1d0-tau_N)*p/(1d0-p)*(s_LH+s_LL)*y
    endif 

    d_Eta_I = 1d0/I_tilde/(sigma-1d0)
    d_Gamma_I = B_NL/(1d0-exp(B_NL*(S_tilde-I_tilde)*(sigma-1d0)))
    d_Gamma_S = B_NL/(1d0-exp(B_NL*(I_tilde-S_tilde)*(sigma-1d0)))
    a_I = d_Gamma_I/(s_K/(1d0-eta)*(d_Eta_I-S_LL/(s_LH+s_LL)*d_Gamma_I))
    a_S = (sigma-1d0)/sigma*(s_LL+s_LH)/s_LH*d_Gamma_S+B_NH-B_NL   

end subroutine Firm_Problem


subroutine RD_Problem
    implicit none
 real(8) :: gamma_IL, gamma_SH, gamma_SL, beta_S

    gamma_IL = exp(B_NL*(1d0-I_tilde)*(1d0-sigma))
    gamma_SH = exp(B_NH*(1d0-S_tilde)*(1d0-sigma))
    gamma_SL = exp(B_NL*(1d0-S_tilde)*(1d0-sigma))
    beta_S = exp(-(rr-g)*(1d0-S_tilde)/g_N)
    if (SPP == 0) then
        V_NH = p*(w_H/A)**(1d0-sigma)/(rr-g+(sigma-1d0)*B_NH*g_N)
        V_NL = p*(w_L/A*gamma_HL)**(1d0-sigma)/(rr-g+(sigma-1d0)*B_NL*g_N)
        V_II = p*(R/A)**(1d0-sigma)/(rr-g)
        p_N = (1d0-tau_N)*(V_NH+beta_S*(gamma_SL*V_NL-gamma_SH*V_NH))-(1d0-tau_I)*V_II
        p_I = (1d0-tau_I)*V_II-(1d0-tau_N)*gamma_IL*V_NL     
    elseif (SPP == 1) then
        V_NH = (1-p)/(sigma-1d0)*(w_H/A)**(1d0-sigma)/(R/(1d0-p)-delta-g)
        V_NL = (1-p)/(sigma-1d0)*(w_L/A*gamma_HL)**(1d0-sigma)/(R/(1d0-p)-delta-g)
        V_II = (1-p)/(sigma-1d0)*(R/A)**(1d0-sigma)/(R/(1d0-p)-delta-g)
        p_N = s_LL*(B_NH-B_NL)/(R/(1d0-p)-delta-g)+V_NH-V_II
        p_I = V_II-gamma_IL*V_NL
    endif
    
end subroutine RD_Problem


! ======================================================================
! Functions
! ======================================================================


real(8) function u(c)
    implicit none
real(8), intent(in) :: c
     u = c**(1d0-theta)/(1d0-theta)
end function u


real(8) function c1(dV)
    implicit none
real(8), intent(in) :: dV
    c1 = dV**(-1d0/theta)
end function c1


real(8) function mu_hH_h(h_HL)
    implicit none
real(8), intent(in) :: h_HL
    mu_hH_h = mu_h*exp(lambda_h*(h_HL-h_HL_star))
    mu_hL = mu_h
end function mu_hH_h


real(8) function h_HL_mu(mu_hH)
    implicit none 
real(8), intent(in) :: mu_hH 
    h_HL_mu = h_HL_star+log(mu_hH/mu_h)/lambda_h
end function  h_HL_mu


end module Simulation
