 !!====================================================================
!!
!! 
!!> @author 
!!> Blaise Delmotte
!!
!! DESCRIPTION: 
!!> @brief
!!> Compute the first of stresslet from the particle ROS
!!
!!
!!
!! Date :  16/06/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADD_STRESSLET_FORCING

!!====================================================================
!! Here the stresslet forcing is set according to FCM_EIJ
!! using Stokes Law
!!====================================================================
!! Forcing: 
!!------------------------------
!! TO DO : 
!!        1) If ellipse ...
!! !!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use PARAM_PHYS,  only: VISC

implicit none

! Indices for loops
integer :: IP, IND



do IP = 1, FCM_NSPHERE

 ! Compute the first guess of stresslet with the particle ROS
 FCM_SIJ(IP,:) = 20.0/3.0 * VISC * PPI * FCM_SPHERE_RADP(IP)**3 * FCM_EIJ(IP,:)

  
end do

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1

 ! Compute the first guess of stresslet with the particle ROS
 FCM_SIJ(IP,:) = 20.0/3.0 * VISC * PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**3 * FCM_EIJ(IP,:)
  
end do



   
end subroutine FCM_ADD_STRESSLET_FORCING
