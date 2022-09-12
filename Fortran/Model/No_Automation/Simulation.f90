!
!   Simulation.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
module Simulation
use Global_Variable
use My_Function
    implicit none
    
contains

! ======================================================================
! Transition
! ======================================================================


subroutine Transition
    implicit none
integer :: i, iter

    ! Read
    g_N0 = g_N_IR(1); g_N1 = g_N_IR(iter_IR)
    g_h0 = g_h_IR(1); g_h1 = g_h_IR(iter_IR)
    L0 = L_IR(1); L1 = L_IR(iter_IR)
    epsilon_N0 = epsilon_N_IR(1); epsilon_N1 = epsilon_N_IR(iter_IR)
    
    ! Initialize
    dt = T0/real(iter_IR-1,8)
    do i = 1,iter_IR
        ts_t(i) = dt*real(i-1,8)
        ts_dt(i) = dt
    enddo

    ts_g_N0(i) = g_N1
    ts_g_h(i) = g_h1
    
    do iter = 1,iter_max
        ! RD problem
        call Patent_Value
        
        ! Check convergence
        sup = maxval(abs(ts_g_N1-ts_g_N0))
        if (real(int(real(iter,8)/50d0),8)==real(iter,8)/50d0) print*, iter, sup  
        if (sup < TOL1) exit

        ts_g_N0 = ts_g_N0+rho_ts*(ts_g_N1-ts_g_N0)
        
    enddo
    
    N_IR(1) = 1d0; h_IR(1) = 0d0
    do i = 2,iter_IR
        g_N = ts_g_N0(i)
        g_h = ts_g_h(i)
        V_N = ts_V_N(i)
        epsilon_N = mu_N*g_N
        L = 1d0-mu_h*g_h-epsilon_N
        call Firm_Problem
        if (SPP == 0) then
            S_N = (1d0-p)*B_N
        elseif (SPP == 1) then
            S_N = B_N
        endif
        N = N_IR(i-1)+g_N*ts_dt(i-1)
        h = h_IR(i-1)+g_h*ts_dt(i-1)
        g = B_N*g_N+b_h*g_h+(y-y_IR(i-1))/ts_dt(i-1)
        rr = rho+theta*g
        c = w*(1d0-mu_h*g_h)
        call Save_Result0(i)
    enddo
    call Export_Result0

end subroutine Transition


subroutine Patent_Value
    implicit none
integer :: i
   
    do i = iter_IR,1,-1
        g_N = ts_g_N0(i)
        g_h = ts_g_h(i)
        g = B_N*g_N+b_h*g_h
        rr = rho+theta*g
        if (i == iter_IR) then
            ts_V_N(i) = (1d0-p)*B_N/(rr-g+(sigma-1)*B_N*g_N)
        elseif (i < iter_IR) then
            ts_V_N(i) = ((1d0-p)*B_N*ts_dt(i)+ts_V_N(i+1))/        &
                        ((rr-g+(sigma-1)*B_N*g_N)*ts_dt(i)+1d0)
        endif
        V_N = ts_V_N(i)
        L = (1d0-p)*mu_N/V_N
        g_N = (rho/theta*(mu_h/b_h)-(1d0-theta)/theta-L)/               &
              ((1d0-theta)/theta*(mu_h/b_h)+(mu_N/B_N))/B_N
        g_h = (1d0/theta*((b_h/mu_h)-rho)+(1d0-theta)/theta*B_N*g_N)/b_h
        
        ! Updata
        ts_g_N1(i) = g_N
        ts_g_h(i) = g_h
    enddo
    
end subroutine Patent_Value


! ======================================================================
! Steady State
! ======================================================================

subroutine BGP(SPP0)
    implicit none
integer, intent(in) :: SPP0

    SPP = SPP0
    if (SPP == 1) p = 0d0

    I_tilde = 0d0
    if ((B_N/mu_N) < (b_h/mu_h)**2d0/(rho/theta-(1d0-theta)/theta*(b_h/mu_h))) then
        rr = b_h/mu_h
        g_N = 0d0
        g_h = (1d0/theta*((b_h/mu_h)-rho))/b_h
    elseif ((b_h/mu_h) < sigma/(theta+sigma-1d0)*rho-(1d0-theta)/(theta+sigma-1d0)*(B_N/mu_N)) then
        rr = (sigma-1d0)/(theta+sigma-1d0)*rho+theta/(theta+sigma-1d0)*(B_N/mu_N)
        g_N = (1d0/(theta+sigma-1d0)*((B_N/mu_N)-rho))/B_N
        g_h = 0d0
    else
        rr = (b_h/mu_h)+((rho/theta-(1d0-theta)/theta*(b_h/mu_h))*((B_N/mu_N)-(b_h/mu_h)))/ &
                        (sigma*(b_h/mu_h)+(1d0-theta)/theta*((B_N/mu_N)-(b_h/mu_h)))
        g_N = (rr-(b_h/mu_h))/B_N
        g_h = ((rr-rho)/theta-B_N*g_N)/b_h
    endif
    g = B_N*g_N+b_h*g_h
    epsilon_N = g_N*mu_N
    epsilon_I = g_I*mu_I
    L = 1d0-epsilon_N-epsilon_I
    call Firm_Problem
    if (SPP == 0) then
        S_N = (1d0-p)*B_N
    elseif (SPP == 1) then
        S_N = B_N
    endif
    
    V_N = S_N/(rr-g+(sigma-1d0)*B_N*g_N)
    c = w*L+pi
    Welfare = u(c)/(rho-(1d0-theta)*g)

    if (SPP == 0) then
        call File_Write
    elseif (SPP == 1) then
        call File_Write_SPP
    endif

end subroutine BGP



subroutine Firm_Problem
    implicit none

    ! Productivity
    Gamma = ((1d0-exp(-B_N*(sigma-1d0)))/(B_N*(sigma-1d0)))**(1d0/(sigma-1d0))
    y = A/(1d0-p)*Gamma*L
    pi = p*y
    w = A*Gamma
    s_L = (1d0-p)

end subroutine Firm_Problem



real(8) function u(c)
    implicit none
real(8), intent(in) :: c
     u = c**(1d0-theta)/(1d0-theta)
end function u

end module Simulation
