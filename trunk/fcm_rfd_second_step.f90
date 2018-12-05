 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Compute RFD second step:
!!  1) Disturb positions
!!  2) Distribute random RFD forces ONLY
!!  3) Stresslet iterations
!!  4) Get a second set of velocity
!!  5) Sum this set with the previous one
!!  CAUTION: ONLY VALID FOR FORWARD EULER
!! Date :  27/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_RFD_SECOND_STEP

!!====================================================================
!!
!!====================================================================

use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use P3DFFT

implicit none

FCM_XP0 = FCM_XP
FCM_YP0 = FCM_YP
FCM_ZP0 = FCM_ZP

FCM_UP0 = FCM_UP(:,1)
FCM_VP0 = FCM_VP(:,1)
FCM_WP0 = FCM_WP(:,1)

FCM_OMPX0 = FCM_OMPX 
FCM_OMPY0 = FCM_OMPY
FCM_OMPZ0 = FCM_OMPZ

FCM_FORCING_X(:,:,:) = 0.0
FCM_FORCING_Y(:,:,:) = 0.0
FCM_FORCING_Z(:,:,:) = 0.0

!- Disturb positions: x = x + Delta_q*epsilon
FCM_XP = FCM_XP0 + FCM_RFD_RAND_FORCE(:,1)*FCM_EPS_RFD**2
FCM_YP = FCM_YP0 + FCM_RFD_RAND_FORCE(:,2)*FCM_EPS_RFD**2
FCM_ZP = FCM_ZP0 + FCM_RFD_RAND_FORCE(:,3)*FCM_EPS_RFD**2

!!!!!!!!!! TO COMMENT !!!!!!!!!!!!!!!!!!!!
!! JUST A TEST TO CHECK THAT MIRROR IMAGE HAS SAME VELOCITY
!~ FCM_XP(2) = LXMAX - FCM_XP(1)
!~ FCM_YP(2) = FCM_YP(1)
!~ FCM_ZP(2) = FCM_ZP(1)
!!!!!!!!!! TO COMMENT !!!!!!!!!!!!!!!!!!!!

if (FCM_NELLIPSOID==0) then
 ! If only spheres, store gaussians
 call FCM_COMPUTE_GAUSSIAN_ENV    
else  
 ! If only ellipsoid, only store gaussian support
 call FCM_COMPUTE_GAUSSIAN_SUPPORT  
end if

FCM_FORCE = FCM_RFD_RAND_FORCE

call FCM_DISTRIB_MONOPOLE

if (FCM_ACTIVATE_STRESSLET==1) then
 ! Don't take shear into account in 2nd RFD step
 FCM_CONSIDER_ROS_SHEAR = 0

 FCM_EIJ = 0.0
 call FCM_DISTRIB_STRESSLET(FCM_SIJ)
 
 if (FCM_BC==2) then
  call FCM_MIRROR_FORCES_X
 end if

 call FCM_FLUID_PREDICTION 
 call FCM_RATE_OF_STRAIN_FILTER

 FCM_FORCING_X(:,:,:) = 0.0
 FCM_FORCING_Y(:,:,:) = 0.0
 FCM_FORCING_Z(:,:,:) = 0.0

 call FCM_CONJUGATE_GRADIENT_STRESSLET

 FCM_FORCING_X(:,:,:) = 0.0
 FCM_FORCING_Y(:,:,:) = 0.0
 FCM_FORCING_Z(:,:,:) = 0.0

 call FCM_DISTRIB_STRESSLET(FCM_SIJ)

 call FCM_DISTRIB_MONOPOLE

end if ! if (FCM_ACTIVATE_STRESSLET==1) then

if (FCM_BC==2) then
  call FCM_MIRROR_FORCES_X
end if
 
call FCM_FLUID_PREDICTION 

call FCM_VELOCITY_FILTER
call FCM_ROTATION_FILTER

!!!! TO COMMENT !!!!!!!!!!!!!!!!
if (FCM_ACTIVATE_STRESSLET==1) then
 call FCM_RATE_OF_STRAIN_FILTER
end if
!!!! TO COMMENT !!!!!!!!!!!!!!!!

FCM_UP(:,1) = FCM_UP(:,1) + FCM_UP0
FCM_VP(:,1) = FCM_VP(:,1) + FCM_VP0
FCM_WP(:,1) = FCM_WP(:,1) + FCM_WP0

FCM_OMPX = FCM_OMPX + FCM_OMPX0
FCM_OMPY = FCM_OMPY + FCM_OMPY0
FCM_OMPZ = FCM_OMPZ + FCM_OMPZ0


! Put the particles back where they were for advection and saves
FCM_XP = FCM_XP0
FCM_YP = FCM_YP0
FCM_ZP = FCM_ZP0

end subroutine FCM_RFD_SECOND_STEP
