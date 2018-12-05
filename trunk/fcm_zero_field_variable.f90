 !!====================================================================
!!
!! 
!!> @brief
!!> Routine that resets field variables at the new time step
!!
!! Date :  09/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ZERO_FIELD_VARIABLE


use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none

!---------------------------------------!
!1.FLUID VELOCITY                       !
!---------------------------------------!
UFLU(:,:,:) = 0.0
VFLU(:,:,:) = 0.0
WFLU(:,:,:) = 0.0

!---------------------------------------!
!2.FORCING FIELD                        !
!---------------------------------------!
FCM_FORCING_X(:,:,:) = 0.0
FCM_FORCING_Y(:,:,:) = 0.0
FCM_FORCING_Z(:,:,:) = 0.0


 


   
end subroutine FCM_ZERO_FIELD_VARIABLE
