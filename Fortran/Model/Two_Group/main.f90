!
!   main.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
!
!   main.f90
!   Growth
!
!   Created by Rong Fan on 10/10/21.
!   Copyright 2020 FanRong. All rights reserved.
!
program main
use Global_Variable
use Simulation
    implicit none
integer :: j, num
integer :: num_total, iostatus
real(8) :: param(Nparam), param_range(2,Nparam)

    moment_data(1) = 0.04d0                         ! interest rate
    moment_data(2) = 0.028d0                        ! growth rate (RD)
    moment_data(3) = 0.003d0                        ! growth rate (human capital)
    moment_data(4) = 0.625d0                        ! labor share (1980)
    moment_data(5) = 1.4d0                          ! wage premium (1980)
    moment_data(6) = 0.028d0                        ! RD/GDP ratio
    moment_data(7) = 0.605d0                        ! labor share (2005)
    moment_data(8) = 1.6d0                          ! wage premium (2005)
    moment_data(9) = 0.141d0                        ! Change of training time (skilled 1999-2005)
    moment_data(10) = 0.160d0                       ! Change of training time (unskilled 1999-2005)

    moment_weight(1) = 1d0                          ! interest rate
    moment_weight(2) = 1d0                          ! growth rate
    moment_weight(3) = 1d0                          ! growth rate (human capital 1980)
    moment_weight(4) = 1d0                          ! labor share (1980)
    moment_weight(5) = 1d0                          ! wage premium (1980)
    moment_weight(6) = 1d0                          ! RD/GDP ratio (1980)
    moment_weight(7) = 1d0                          ! lasbor share (2005)
    moment_weight(8) = 1d0                          ! wage premium (2005)
    moment_weight(9) = 1d0                          ! Change of training time (skilled)
    moment_weight(10) = 1d0                         ! Change of training time (unskilled)

if (calibration == 1) then 
    call File_Create1
    open(2, file='L_N.txt', status='replace')
        FMT = '(12A10)'
        write(2,FMT) 'number', 'L_N', 'r', 'g_N', 'g_h', 's_L0', 'w0', 'RD',        &
                                      's_L1', 'w1', 'l_H', 'l_L'
        FMT = '(2I10,10F10.4)'
        write(2,FMT) 0, 0, moment_data
    close(2)

    param(1) = 0.1125d0; param_range(1,1) = 0.1d0; param_range(2,1) = 0.15d0            ! p
    param(2) = 0.119d0; param_range(1,2) = 0.1d0; param_range(2,2) = 0.15d0             ! A
    param(3) = 2.2651d0; param_range(1,3) = 2d0; param_range(2,3) = 3d0                 ! B_NH
    param(4) = 2.3619d0; param_range(1,4) = 1d0; param_range(2,4) = 10d0                ! mu_N
    param(5) = 9.0135d0; param_range(1,5) = 1d0; param_range(2,5) = 10d0                ! mu_I
    param(6) = 193.5608d0; param_range(1,6) = 150d0; param_range(2,6) = 200d0           ! mu_h
    param(7) = 0.9026d0; param_range(1,7) = 0.8d0; param_range(2,7) = 0.95d0            ! alhpa_H
    param(8) = 0.2797d0; param_range(1,8) = 0.1d0; param_range(2,8) = 0.4d0             ! alhpa_L
    param(9) = 0.7675d0; param_range(1,9) = 0.7d0; param_range(2,9) = 0.95d0            ! lambda_h
    param(10) = 0.7853d0; param_range(1,10) = 0.75d0; param_range(2,10) = 0.95d0        ! tech shock

    iunit = 1; L_N0 = 1d6
    Model = 3; Np = Nparam; call nlopt(6,param,param_range)

elseif (calibration == 0) then  
    call File_Create1
    open(unit=2, file='BGP', status='old')
        num_total = 0
        do
            num_total = num_total+1
            read(2,*,IOSTAT=iostatus) param1(:,num_total)
            if (iostatus < 0) exit
        end do
    close(2)

    select case (Task)
        case (1) ! Static
            do num = 1,num_total-1
                iunit = num 
                Model = 3 
                call Read_Param1 
                call BGP(0) 
                call Static(0) 
            end do

        case (2) ! BGP
            do num = 1,num_total-1
                iunit = num 
                Model = 3
                call Read_Param1 
                mu_I = mu_I*shock
                call BGP(0)
                s_K0 = s_K; s_LH0 = s_LH; s_LL0 = s_LL
                k0 = k; L0 = (s_LH*L_H+s_LL*L_L)/(s_LH+s_LL)
                w0 = w; h_HL0 = h_HL
                call Save_BGP1(1)

                write(*,*) ''    
                mu_I = shock*mu_I
                call BGP(0)
                s_K1 = s_K; s_LH1 = s_LH; s_LL1 = s_LL
                k1 = k; L1 = (s_LH*L_H+s_LL*L_L)/(s_LH+s_LL)
                w1 = w; h_HL1 = h_HL

                open(2, file='Decompose.txt', status='old', position='append')
                    FMT = '(3I10,7F10.4)'
                    write(2,FMT) iunit, Model, flag,                                                &
                    (B_NH-B_NL)/a_S*a_I,                                                            &            
                    (B_NH-B_NL)/a_S*a_I*(log(S_LH1+S_LL1)-log(S_LH0+S_LL0)),                        & 
                    (B_NH-B_NL)/a_S*a_I*(sigma-1d0)/sigma*s_K/(1d0-p),                              & 
                    (B_NH-B_NL)/a_S*a_I*(sigma-1d0)/sigma*s_K/(1d0-p)*                              & 
                    (log(k1/L1)-log(k0/L0)+s_LL1/(s_LH1+s_LL1)*h_HL1-s_LL0/(s_LH0+s_LL0)*h_HL0),    & 
                    1d0-(B_NH-B_NL)/a_S/sigma,                                                      & 
                    (1d0-(B_NH-B_NL)/a_S/sigma)*(h_HL1-h_HL0),                                      &
                    log(w1/w0)
                close(2)
            end do

        case (3) ! BGP and Transition  
            do num = 1,num_total-1
                iunit = num
                ! Get h_HL       
                Model = 3 
                call Read_Param1
                tau_k = 0d0; tau_N = 0d0; tau_I = 0d0
                tau_hH = 0d0; tau_hL = 0d0
                call BGP(0)

                g_hH0 = g_hH; g_hL0 = g_hL; h_HL0 = h_HL

                do j = 1,3,2
                    Model = j  
                    call Read_Param1
                    tau_k = 0d0; tau_N = 0d0; tau_I = 0d0
                    tau_hH = 0d0; tau_hL = 0d0
                    call BGP(0)
                    call Save_BGP1(1)
                    write(*,*) '==============================================================================' 
                    
                    call Read_Param1   
                    mu_I = shock*mu_I
                    call BGP(0)
                    call Save_BGP1(iter_IR)
                    call Transition
                    call Export_Result1
                enddo 
            end do

    end select

endif 

end program main


