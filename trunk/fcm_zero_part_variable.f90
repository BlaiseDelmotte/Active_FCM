 !!====================================================================
!!
!! 
!!> @brief
!!> Routine that resets particle variables at the new time step
!!
!! Date :  09/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ZERO_PART_VARIABLE


use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none

!---------------------------------------!
!1.FORCE MONOPOLE  AND TORQUE           !
!---------------------------------------!
FCM_FORCE(:,:) = 0.0
FCM_TORQUE(:,:) = 0.0

!---------------------------------------!
!2.ROTLET                               !
!---------------------------------------!
FCM_AIJ(:,:) = 0.0

!---------------------------------------!
!3.VELOCITY                             !
!---------------------------------------!
FCM_UP(:,1) = 0.0
FCM_VP(:,1) = 0.0
FCM_WP(:,1) = 0.0

!---------------------------------------!
!4.ROTATION                             !
!---------------------------------------!
FCM_OMPX(:) = 0d0
FCM_OMPY(:) = 0d0
FCM_OMPZ(:) = 0d0




   
end subroutine FCM_ZERO_PART_VARIABLE
