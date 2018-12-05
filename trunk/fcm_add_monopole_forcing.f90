 !!====================================================================
!!
!! 
!!> @brief
!!> Routine adding the monopole forcing for each particle
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADD_MONOPOLE_FORCING

!!====================================================================
!! Here the monopole forcing is set according to FCM_EXT_FORCE
!! using Stokes Law
!!====================================================================
!! Forcing: 
!!------------------------------
!! TO DO : 
!!        1) If lub
!!        2) If fb
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use PARAM_PHYS,  only: VISC

implicit none

! Indices for loops
integer :: IP, IND



do IP = 1, FCM_NSPHERE

! Define the force monopole from the external unique force if any
 FCM_FORCE(IP,1) = FCM_FORCE(IP,1) + FCM_EXT_FORCE(1) * 6.0 * PPI * VISC * FCM_SPHERE_RADP(IP) !* (-1d0)**(IP+1)
 FCM_FORCE(IP,2) = FCM_FORCE(IP,2) + FCM_EXT_FORCE(2) * 6.0 * PPI * VISC * FCM_SPHERE_RADP(IP) ! * (-1d0)**(IP+1)
 FCM_FORCE(IP,3) = FCM_FORCE(IP,3) + FCM_EXT_FORCE(3) * 6.0 * PPI * VISC * FCM_SPHERE_RADP(IP) ! * (-1d0)**(IP+1)
  
end do

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL
 IND = IND + 1

! Define the force monopole from the external unique force if any
 FCM_FORCE(IP,1) = FCM_FORCE(IP,1) + FCM_EXT_FORCE(1) * 6.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:)) ! * (-1d0)**(IP+1)
 FCM_FORCE(IP,2) = FCM_FORCE(IP,2) + FCM_EXT_FORCE(2) * 6.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:)) ! * (-1d0)**(IP+1)
 FCM_FORCE(IP,3) = FCM_FORCE(IP,3) + FCM_EXT_FORCE(3) * 6.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:)) ! * (-1d0)**(IP+1)
  
end do

   
end subroutine FCM_ADD_MONOPOLE_FORCING
