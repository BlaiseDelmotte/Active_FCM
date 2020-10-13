 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating FBLEVEL
!!
!! Date :  26/02/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_REP_WALL_MAGNITUDE


use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE



implicit none


if (FCM_NSPHERE>0) then

 FCM_WALL_LEVEL =  20.0*6.0*PPI*20.0*maxval(FCM_SPHERE_RADP)
 FCM_WALL_RANGE = 1.1*maxval(FCM_SPHERE_RADP)
 
else

!!!!--- TO BE DONE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end if


end subroutine FCM_INITIATION_REP_WALL_MAGNITUDE
