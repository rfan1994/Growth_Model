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
use Burkardt_fsolve

    implicit none
    
contains

! ======================================================================
! Transition
! ======================================================================


subroutine Transition
    implicit none
integer :: iter, i
real(8) :: kmin, kmax, g_dot
integer :: p0, p1

    ! Read
    I0 = I_tilde_IR(1); I1 = I_tilde_IR(iter_IR)
    k0 = k_IR(1); k1 = k_IR(iter_IR)   
    g_N0 = g_N_IR(1); g_N1 = g_N_IR(iter_IR)
    L0 = L_IR(1); L1 = L_IR(iter_IR)
    epsilon_N0 = epsilon_N_IR(1); epsilon_N1 = epsilon_N_IR(iter_IR)
    epsilon_I0 = epsilon_I_IR(1); epsilon_I1 = epsilon_I_IR(iter_IR)
    rr0 = rr_IR(1); rr1 = rr_IR(iter_IR)

    ! Initialize
    dt = T/real(iter_IR-1,8)
    do i = 1,iter_IR
        ts_t(i) = dt*real(i-1,8)
        ts_dt(i) = dt
    enddo
    
    ! State
    ts_I0 = I1
    ts_k0 = k1
    call Patent_Value_T

    do iter = 1,iter_max
        ! RD problem
        call Patent_Value
        ! HH problem
        call Policy_Function
           
        ! Update I
        ts_I1(1) = I0
        do i = 2,iter_IR
            g_dot = ts_g_I1(i)-ts_g_N1(i)
            I_tilde = ts_I1(i-1)+g_dot*ts_dt(i-1)
            if (I_tilde .le. min(I0,I1)) then
                ts_I1(i) = min(I0,I1)
            elseif (I_tilde .ge. max(I0,I1)+TOL2) then
                ts_I1(i) = max(I0,I1)+TOL2
            else
                ts_I1(i) = I_tilde
            endif
        enddo

        ! Update k
        ts_k1(1) = k0
        do i = 2,iter_IR
            k = ts_k1(i-1)+ts_a_dk(i)*ts_dt(i-1)
            if (k .le. min(k0,k1)) then
                ts_k1(i) = min(k0,k1)
            elseif (k .ge. max(k0,k1)) then
                ts_k1(i) = max(k0,k1)
            else
                ts_k1(i) = k
            endif
        enddo

        ! Check convergence
        sup = max(maxval(abs(ts_I1-ts_I0)),         &
                  maxval(abs(ts_k1-ts_k0)),         & 
                  maxval(abs(ts_g_N1-ts_g_N0)))
        if (real(int(real(iter,8)/50d0),8)==real(iter,8)/50d0) print*, iter, sup
        if (sup < TOL1) exit

        ts_I0 = ts_I0+rho_ts*(ts_I1-ts_I0)
        ts_k0 = ts_k0+rho_ts*(ts_k1-ts_k0)
        ts_g_N0 = ts_g_N0+rho_ts*(ts_g_N1-ts_g_N0)
        ts_g_I0 = ts_g_I0+rho_ts*(ts_g_I1-ts_g_I0)       
        ts_L0 = ts_L0+rho_ts*(ts_L1-ts_L0)
        ts_V_NN0 = ts_V_NN0+rho_ts*(ts_V_NN1-ts_V_NN0)
        ts_V_NI0 = ts_V_NI0+rho_ts*(ts_V_NI1-ts_V_NI0)
        ts_V_II0 = ts_V_II0+rho_ts*(ts_V_II1-ts_V_II0)
        
    enddo

    N_IR(1) = 1d0; h_IR = 0d0
    do i = 2,iter_IR
        I_tilde = ts_I0(i)
        k = ts_k0(i)
        L = ts_L0(i)
        g_N = ts_g_N0(i)
        g_I = ts_g_I0(i)
        c = ts_a_c(i)
        Welfare = ts_V(i)
        call Firm_Problem
        call Decompose('rr')
        call RD_Problem
        V_N = ts_V_N(i)
        V_I = ts_V_I(i)
        N = N_IR(i-1)+g_N*ts_dt(i-1)
        g = B_N*g_N+(y-y_IR(i-1))/ts_dt(i-1)
        epsilon_N = mu_N*g_N
        epsilon_I = mu_I*g_I
        call Save_Result0(i)
    enddo
    call Export_Result0

end subroutine Transition


subroutine Patent_Value_T
    implicit none
real(8) :: x(2), fvec(2)

    I_tilde = I1
    k = k1

    x(1) = epsilon_N1; x(2) = epsilon_I1
    call fsolve(RD_Equation,2,x,fvec,TOL,flag)

    ts_g_N0 = g_N; ts_g_N1 = g_N
    ts_g_I0 = g_I; ts_g_I1 = g_I
    ts_L0 = L; ts_L1 = L

    ts_V_NN0 = V_NN_T; ts_V_NN1 = V_NN_T
    ts_V_NI0 = V_NI_T; ts_V_NI1 = V_NI_T
    ts_V_II0 = V_II_T; ts_V_II1 = V_II_T
    ts_V_N = V_NN_T-V_II_T
    ts_V_I = V_II_T-V_NI_T

contains
 
subroutine RD_Equation(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)
real(8) :: gamma_I

    if (x(1) < 0d0 .or. x(2) < 0d0 .or.    &
        x(1)+x(2) > 1d0) then 
        fvec = -1d6 
        return 
    endif

    epsilon_N = x(1)
    epsilon_I = x(2)
    L = 1d0-epsilon_N-epsilon_I

    call Firm_Problem
    g_N = 1d0/mu_N*epsilon_N**lambda
    g_I = 1d0/mu_I*epsilon_I**lambda
    g = B_N*g_N

    gamma_I = exp(B_N*(1d0-I_tilde))
    if (SPP == 0) then
        S_NN = p*(w/A)**(1d0-sigma)
        S_NI = p*(w/A*gamma_I)**(1d0-sigma)
        S_II = p*(R/A)**(1d0-sigma)
    elseif (SPP == 1) then
        S_NN = (w/A)**(1d0-sigma)
        S_NI = (w/A*gamma_I)**(1d0-sigma)
        S_II = (R/A)**(1d0-sigma)
    endif

    V_NN_T = S_NN/(rr-g+(sigma-1d0)*B_N*g_N)
    V_NI_T = S_NI/(rr-g+(sigma-1d0)*B_N*g_N)
    V_II_T = S_II/(rr-g)
    V_N = V_NN_T-V_II_T
    V_I = V_II_T-V_NI_T

    fvec(1) = lambda/mu_N*V_N*epsilon_N**(lambda-1d0)-s_L/L
    fvec(2) = lambda/mu_I*V_I*epsilon_I**(lambda-1d0)-s_L/L

end subroutine RD_Equation

end subroutine Patent_Value_T



subroutine Patent_Value
    implicit none
integer :: i
real(8) :: V_NN0, V_NI0, V_II0
real(8) :: V_NN1, V_NI1, V_II1
real(8) :: dt, x(2), fvec(2)

    do i = iter_IR-1,1,-1
        dt = ts_dt(i)
        I_tilde = ts_I0(i)
        k = ts_k0(i)

        V_NN1 = ts_V_NN0(i+1)
        V_NI1 = ts_V_NI0(i+1)
        V_II1 = ts_V_II0(i+1)

        x(1) = epsilon_N1; x(2) = epsilon_I1
        call fsolve(RD_Equation,2,x,fvec,TOL,flag)
        
        ! Updata
        ts_g_N1(i) = g_N
        ts_g_I1(i) = g_I
        ts_L1(i) = L
        ts_V_NN1(i) = V_NN0
        ts_V_NI1(i) = V_NI0
        ts_V_II1(i) = V_II0
        ts_V_N(i) = V_NN0-V_II0
        ts_V_I(i) = V_II0-V_NI0
    enddo

contains
 
subroutine RD_Equation(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)
real(8) :: gamma_I

    if (x(1) < 0d0 .or. x(2) < 0d0 .or.    &
        x(1)+x(2) > 1d0) then 
        fvec = -1d6 
        return 
    endif

    epsilon_N = x(1)
    epsilon_I = x(2)
    L = 1d0-epsilon_N-epsilon_I

    call Firm_Problem
    call Decompose('rr')
    g_N = 1d0/mu_N*epsilon_N**lambda
    g_I = 1d0/mu_I*epsilon_I**lambda
    g = B_N*g_N

    gamma_I = exp(B_N*(1d0-I_tilde))

    if (SPP == 0) then
        S_NN = p*(w/A)**(1d0-sigma)
        S_NI = p*(w/A*gamma_I)**(1d0-sigma)
        S_II = p*(R/A)**(1d0-sigma)
    elseif (SPP == 1) then
        S_NN = (w/A)**(1d0-sigma)
        S_NI = (w/A*gamma_I)**(1d0-sigma)
        S_II = (R/A)**(1d0-sigma)
    endif

    V_NN0 = (S_NN*dt+V_NN1)/((rr-g+(sigma-1d0)*B_N*g_N)*dt+1d0)
    V_NI0 = (S_NI*dt+V_NI1)/((rr-g+(sigma-1d0)*B_N*g_N)*dt+1d0)
    V_II0 = (S_II*dt+V_II1)/((rr-g)*dt+1d0)
    V_N = V_NN0-V_II0
    V_I = V_II0-V_NI0

    V_NN_T = S_NN/(rr-g+(sigma-1d0)*B_N*g_N)
    V_NI_T = S_NI/(rr-g+(sigma-1d0)*B_N*g_N)
    V_II_T = S_II/(rr-g)
    V_N = V_NN_T-V_II_T
    V_I = V_II_T-V_NI_T

    fvec(1) = lambda/mu_N*V_N*epsilon_N**(lambda-1d0)-s_L/L
    fvec(2) = lambda/mu_I*V_I*epsilon_I**(lambda-1d0)-s_L/L

end subroutine RD_Equation

end subroutine Patent_Value



subroutine Value_Function_T
    implicit none
integer :: ki
real(8) :: c0

    ! Grid point
    kmin = r_k0*min(k0,k1)
    kmax = r_k1*max(k0,k1)
    do ki = 1,Nk
        k_grid(ki) = kmin+((kmax-kmin)**(1d0/lambda_k)*real(ki-1,8)/real(Nk-1,8))**lambda_k
    enddo
    do ki = 1,Nk-1
        dk(ki) = k_grid(ki+1)-k_grid(ki)
    enddo 
    
    I_tilde = I1
    k = k1
    L = L1
    call Firm_Problem

    g_N = g_N1
    g = B_N*g_N
    
    do ki = 1,Nk
        k = k_grid(ki)
        c0 = (rr-g)*k+w*L+pi
        V_T(ki) = u(c0)/(rho-(1d0-theta)*g)
    enddo

    k = k1
    ts_a_dk(iter_IR) = 0d0
    ts_a_c(iter_IR) = (rr-g)*k+w*L+pi

end subroutine Value_Function_T



subroutine Value_Function
    implicit none
integer :: iter, ki, i
integer :: p0, p1
real(8), dimension(Nk) :: V0, V1, a_dk, a_c, Vchange
real(8) :: dVf, dVb, c0, cf, cb, mu, muf, mub
real(8) :: dVt, dt

    V0 = V_T
    V1 = V_T

    do i = iter_IR-1,1,-1
        dt = ts_dt(i)
        I_tilde = ts_I0(i)
        k = ts_k0(i)
        L = ts_L0(i)
        call Firm_Problem

        g_N = ts_g_N0(i)
        g = B_N*g_N
               
        do iter = 1,iter_max
            do ki = 1,Nk
                k = k_grid(ki)
                if (ki == 1) then
                    dVf = (V0(ki+1)-V0(ki))/dk(ki)
                    dVb = 0d0
                    cf = c1(dVf)
                    cb = c1(dVb)
                    muf = (rr-g)*k+w*L+pi-cf
                    mub = (rr-g)*k+w*L+pi-cb
                elseif (ki > 1 .and. ki < NK) then
                    dVf = (V0(ki+1)-V0(ki))/dk(ki)
                    dVb = (V0(ki)-V0(ki-1))/dk(ki-1)
                    cf = c1(dVf)
                    cb = c1(dVb)
                    muf = (rr-g)*k+w*L+pi-cf
                    mub = (rr-g)*k+w*L+pi-cb
                elseif (ki == Nk) then
                    dVf = 0d0
                    dVb = (V0(ki)-V0(ki-1))/dk(ki-1)
                    cf = c1(dVf)
                    cb = c1(dVb)
                    muf = (rr-g)*k+w*L+pi-cf
                    mub = (rr-g)*k+w*L+pi-cb
                endif
            
                dVt = (V1(ki)-V0(ki))/dt
                if (muf > 0d0) then
                    Vchange(ki) = u(cf)+dVf*muf+dVt-(rho-(1d0-theta)*g)*V0(ki)
                    a_dk(ki) = muf
                    a_c(ki) = cf
                elseif (mub < 0d0) then
                    Vchange(ki) = u(cb)+dVb*mub+dVt-(rho-(1d0-theta)*g)*V0(ki)
                    a_dk(ki) = mub
                    a_c(ki) = cb
                else
                    c0 = (rr-g)*k+w*L+pi
                    Vchange(ki) = u(c0)+dVt-(rho-(1d0-theta)*g)*V0(ki)
                    a_dk(ki) = 0d0
                    a_c(ki) = c0
                endif
            enddo
            
            sup = maxval(abs(Vchange))
            if (sup < TOL) exit    
            V0 = V0+Vchange*rho_V
        enddo

        V1 = V0
        k = ts_k0(i)
        p0 = 1; p1 = Nk
        call Interpolate(Nk,k_grid,k,p0,p1)
        ts_a_dk(i) = a_dk(p0)
        ts_a_c(i) = a_c(p0)

    enddo

end subroutine Value_Function


subroutine Policy_Function
    implicit none 
integer :: i
real(8) :: dt, g_c, c0, c1, V0, V1

    c1 = c_IR(iter_IR)
    V1 = V_IR(iter_IR)
    ts_a_c(iter_IR) = c1
    ts_a_dk(iter_IR) = 0d0
    ts_V(iter_IR) = V1

    do i = iter_IR-1,1,-1
        dt = ts_dt(i)
        I_tilde = ts_I0(i+1)
        k = ts_k0(i+1)
        L = ts_L0(i+1)
        call Firm_Problem

        g_N = ts_g_N0(i+1)
        g = B_N*g_N       
        g_c = (rr-rho)/theta-g 
        c0 = c1/(1d0+g_c*dt)

        I_tilde = ts_I0(i)
        k = ts_k0(i)
        L = ts_L0(i)
        call Firm_Problem

        ts_a_c(i) = c0
        ts_a_dk(i) = (rr-g)*k+w*L+pi-c0

        V0 = (u(c0)*dt+V1)/((rho-(1d0-theta)*g)*dt+1d0)
        ts_V(i) = V0

        ! Update
        c1 = c0
        V1 = V0
    enddo

end subroutine Policy_Function


subroutine Decompose(variable)
    implicit none 
character(*), intent(in) :: variable

    select case (variable)
        case ('rr')
            rr = rr0
            R = rr0+delta
            s_K = (1d0-p)*(R/A/Eta)**(1d0-sigma)
            s_L = (1d0-p)-s_K
            w = A*Gamma*(s_L/(1d0-p))**(1d0/(1d0-sigma))  
            k = w*L/R*s_k/S_L    
        case ('g_hH')
            g_hH = g_hH0
        case ('g_hL')
            g_hL = g_hL0
        case default 
            return 
    end select

end subroutine Decompose 


! ======================================================================
! Steady State
! ======================================================================

subroutine BGP(SPP0)
    implicit none
integer, intent(in) :: SPP0
integer, parameter :: NI = 7
integer :: i
real(8) :: x(2), fvec(2)

    SPP = SPP0
    if (spp == 1) p = 0d0

    if (rho+delta > A) then 
        do i = 1,NI
            x(1) = 0.1d0*real(i,8); x(2) = 1d0/8d0
            call fsolve(BGP_Equation,2,x,fvec,TOL,flag)
            if (flag == 1) exit
        enddo 
        
        N = 1d0
        k = w*L/R*s_k/S_L
        call Firm_Problem
        c = (rr-g)*k+w*L+pi
        Welfare = u(c)/(rho-(1d0-theta)*g)
    elseif (rho+delta .le. A) then 
        I_tilde = 1d0

        N = 1d0 
        k = 1d0 
        call Firm_Problem
        g = (rr-rho)/theta
        c = (rr-g)*k+w*L+pi
        Welfare = u(c)/(rho-(1d0-theta)*g)
    endif 
    
end subroutine BGP


subroutine BGP_Equation(n,x,fvec)
    implicit none
integer, intent(in) :: n
real(8), intent(in) :: x(n)
real(8), intent(out) :: fvec(n)

    I_tilde = x(1)
    epsilon_N = x(2)
    epsilon_I = (mu_I/mu_N)**(1d0/lambda)*epsilon_N
    L = 1d0-epsilon_N-epsilon_I

    if (I_tilde < 0d0 .or. I_tilde > 1d0 .or.           &
        epsilon_N < 0d0 .or. L < 0d0) then
        fvec = -1d6 
        return 
    endif
   
    g_N = 1d0/mu_N*epsilon_N**lambda
    g_I = 1d0/mu_I*epsilon_I**lambda

    g = B_N*g_N
    rr = rho+theta*g
    R = rr+delta
    Eta = I_tilde**(1d0/(sigma-1d0))
    Gamma = ((1d0-exp(-B_N*(1d0-I_tilde)*(sigma-1d0)))/(B_N*(sigma-1d0)))**(1d0/(sigma-1d0))
    s_K = (1d0-p)*(R/A/Eta)**(1d0-sigma)
    s_L = (1d0-p)-s_K
    w = A*Gamma*(s_L/(1d0-p))**(1d0/(1d0-sigma))
    call RD_Problem

    fvec(1) = lambda/mu_N*V_N*epsilon_N**(lambda-1d0)-s_L/L
    fvec(2) = lambda/mu_I*V_I*epsilon_I**(lambda-1d0)-s_L/L
    if (I_tilde == TOL) fvec(2) = max(fvec(2),0d0)

end subroutine BGP_Equation



subroutine RD_Problem
    implicit none
 real(8) :: gamma_I

    gamma_I = exp(B_N*(1d0-I_tilde))
    if (SPP == 0) then
        S_NN = p*(w/A)**(1d0-sigma)
        S_NI = p*(w/A*gamma_I)**(1d0-sigma)
        S_II = p*(R/A)**(1d0-sigma)
    elseif (SPP == 1) then
        S_NN = (w/A)**(1d0-sigma)
        S_NI = (w/A*gamma_I)**(1d0-sigma)
        S_II = (R/A)**(1d0-sigma)
    endif
    V_NN = S_NN/(rr-g+(sigma-1d0)*B_N*g_N)
    V_NI = S_NI/(rr-g+(sigma-1d0)*B_N*g_N)
    V_II = S_II/(rr-g)
    S_N = S_NN-S_II
    S_I = S_II-S_NI 
    V_N = V_NN-V_II
    V_I = V_II-V_NI

end subroutine RD_Problem



subroutine Firm_Problem
    implicit none

    ! Productivity
    Eta = I_tilde**(1d0/(sigma-1d0))
    Gamma = ((1d0-exp(-B_N*(1d0-I_tilde)*(sigma-1d0)))/(B_N*(sigma-1d0)))**(1d0/(sigma-1d0))
            
    if (I_tilde .le. TOL) then
        ! Output
        k = 0d0; y = A/(1d0-p)*Gamma*L
        pi = p*y
        rr = rho+theta*g
        w = A*Gamma
        s_K = R*k/y
        s_L = w*L/y
    elseif (I_tilde .ge. 1d0-TOL) then
        ! Output
        k = 1d0; y = A/(1d0-p)*Eta*k
        pi = p*y
        R = A*Eta
        rr = R-delta
        w = 0d0
        s_K = R*k/y
        s_L = w*L/y
    elseif (TOL < I_tilde .and. I_tilde < 1d0-TOL) then
        ! Output
        y = A/(1d0-p)*((Eta*k)**((sigma-1d0)/sigma)+(Gamma*L)**((sigma-1d0)/sigma))**(sigma/(sigma-1d0))
        pi = p*y
        R = A*Eta*(y*(1d0-p)/(A*Eta*k))**(1d0/sigma)
        rr = R-delta
        w = A*Gamma*(y*(1d0-p)/(A*Gamma*L))**(1d0/sigma)
        s_K = R*k/y
        s_L = w*L/y
    endif 

end subroutine Firm_Problem


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

end module Simulation

