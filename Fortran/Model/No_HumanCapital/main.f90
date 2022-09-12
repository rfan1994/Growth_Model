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
integer :: num, num_total, iostatus
    
    open(unit=2, file='BGP', status='old')
        num_total = 0
        do
            num_total = num_total+1
            read(2,*,IOSTAT=iostatus) param0(:,num_total)
            if (iostatus < 0) exit
        end do
    close(2)

    call File_Create0
    do num = 1,num_total-1
        iunit = num
        call Read_Param0
        call BGP(0)
        call Save_BGP0(1)
     
        mu_I = shock*mu_I
        call BGP(0)
        call Save_BGP0(iter_IR)

        call Transition
        call Export_Result0
    end do

end program main
