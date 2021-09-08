!!====================================================================
!!
!!
!!====================================================================

subroutine CPUTIME_INFO(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM,	only: UNIT_INFO, SOLVE_FLUID
use PARAM_PHYS, only: STEADY
use CHECK_CPU


implicit none


!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
!!- Cycle number
integer, intent(in) :: NCYCLE

!!- Index
integer :: I
!!--------------------------------------------------------------------


I = 1
write(UNIT_INFO(I),10500)
write(UNIT_INFO(I),10500)
write(UNIT_INFO(I),10500)
write(UNIT_INFO(I),10500)'====================================================================='
write(UNIT_INFO(I),10500)'		     CPU TIME INFORMATION'
write(UNIT_INFO(I),10500)'====================================================================='
write(UNIT_INFO(I),10500)
write(UNIT_INFO(I),10505)' Full Elapsed CPU TIME:',CPU_ELAPSED,' s'
write(UNIT_INFO(I),10504)'        Run initiation:',CPU_INITIATION,            ' s',&
                                                   CPU_INITIATION/CPU_ELAPSED,' %'
write(UNIT_INFO(I),10506)'       Mean cycle time:',CPU_CYCLE/NCYCLE,' s/it'
write(UNIT_INFO(I),10500)
write(UNIT_INFO(I),10500)'====================================================================='
write(UNIT_INFO(I),10500)'               |     Time      |    Time/it   | % Phase |  % run '
write(UNIT_INFO(I),10500)'====================================================================='


if(SOLVE_FLUID>0) then
write(UNIT_INFO(I),10502)' Fluid         |   ',CPU_FLUID(1),             &
                                               CPU_FLUID(1)/NCYCLE,      &
                                             100.,                        &
                                             100.*CPU_FLUID(1)/CPU_ELAPSED

write(UNIT_INFO(I),10500)'---------------------------------------------------------------------'
write(UNIT_INFO(I),10502)'FFT            |   ',CPU_FLUID(8),             &
                                               CPU_FLUID(8)/NCYCLE,      &
                                          100.*CPU_FLUID(8)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(8)/CPU_ELAPSED
                                                                                    
write(UNIT_INFO(I),10502)'Non-linear     |   ',CPU_FLUID(2),             &
                                               CPU_FLUID(2)/NCYCLE,      &
                                          100.*CPU_FLUID(2)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(2)/CPU_ELAPSED
write(UNIT_INFO(I),10502)'Adv. Fluid     |   ',CPU_FLUID(3),             &
                                               CPU_FLUID(3)/NCYCLE,      &
                                          100.*CPU_FLUID(3)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(3)/CPU_ELAPSED 
write(UNIT_INFO(I),10502)'Proj div free  |   ',CPU_FLUID(5),             &
                                               CPU_FLUID(5)/NCYCLE,      &
                                          100.*CPU_FLUID(5)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(5)/CPU_ELAPSED
if(STEADY)&
write(UNIT_INFO(I),10502)'Forcing        |   ',CPU_FLUID(4),             &
                                               CPU_FLUID(4)/NCYCLE,      &
                                          100.*CPU_FLUID(4)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(4)/CPU_ELAPSED
write(UNIT_INFO(I),10502)'Statistics     |   ',CPU_FLUID(6),             &
                                               CPU_FLUID(6)/NCYCLE,      &
                                          100.*CPU_FLUID(6)/CPU_FLUID(1),&
                                          100.*CPU_FLUID(6)/CPU_ELAPSED
write(UNIT_INFO(I),10503)'Save           |   ',CPU_FLUID(7),             &
                                               CPU_FLUID(7)/NCYCLE,      &
                                          100.*CPU_FLUID(7)/CPU_ELAPSED
end if !!- if(SOLVE_FLUID>0)

write(UNIT_INFO(I),10500)'====================================================================='


!!--------------------------------------------------------------------
10500 format (2x,A)
10501 format (2x,A,E13.6)
10502 format (2x,A,F9.2,' s |',E11.4,' s |',F6.2,' % |',F6.2,' %  ')
10503 format (2x,A,F9.2,' s |',E11.4,' s |    -    |'  ,F6.2,' %  ')

10504 format (2x,A,F9.2,A,F6.2,A)


10505 format (2x,A,F9.2,A)
10506 format (2x,A,E11.4,A)

end subroutine CPUTIME_INFO
