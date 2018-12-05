 !!====================================================================
!!
!! 
!!> @brief
!!> Routine prescribing the phase shift between unsteady swimmers
!!
!! Date :  08/04/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PHASE_SHIFT


use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE

implicit none

real(kind=8) :: RAND_FACTOR

integer :: IP

! - FCM_SWIMMING==3 : All swimmers synchronized
if ((FCM_NSWIM(3)>0).and.(FCM_NSWIM(4).eq.0)) then

 FCM_PHASE_SHIFT = 0.0
                      
! - FCM_SWIMMING==4 : Phase shift between two halves of the population
else if ((FCM_NSWIM(4)>0).and.(FCM_NSWIM(5).eq.0)) then
 FCM_PHASE_SHIFT(1:int(real(FCM_NSWIM(4))/2.0)) = 0.0
 FCM_PHASE_SHIFT(int(real(FCM_NSWIM(4))/2.0)+1:FCM_NSWIM(4)) = 0.5*TWOPI
 
! - FCM_SWIMMING==5 : Phase shift randomly distributed
else if (FCM_NSWIM(5)>0) then
 do IP = 1, FCM_NSWIM(5)
  call random_number(RAND_FACTOR)
  FCM_PHASE_SHIFT(IP) = RAND_FACTOR*TWOPI
 end do

end if


end subroutine FCM_INITIATION_PHASE_SHIFT
