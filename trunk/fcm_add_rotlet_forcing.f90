 !!====================================================================
!!
!! 
!!> @author 
!!> Blaise Delmotte
!!
!! DESCRIPTION: 
!!> @brief
!!> Compute the rotlet from the external torque
!!
!!
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADD_ROTLET_FORCING

!!====================================================================
!! Here the rotlet forcing is set according to FCM_EXT_TORQUE
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

 
 ! Define the torque from the external unique torque if any
 ! If torque is "integer" then it follows stokes law, otherwise it 
 ! is prescribed as it is
 if (mod(FCM_EXT_TORQUE(1),1.0).ne.0.0) then 
  FCM_TORQUE(IP,1) = FCM_TORQUE(IP,1) + FCM_EXT_TORQUE(1) 
 else
  FCM_TORQUE(IP,1) = FCM_TORQUE(IP,1) + FCM_EXT_TORQUE(1) * 8.0 * PPI * VISC * FCM_SPHERE_RADP(IP)**3 
 end if
 
 if (mod(FCM_EXT_TORQUE(2),1.0).ne.0.0) then 
  FCM_TORQUE(IP,2) = FCM_TORQUE(IP,2) + FCM_EXT_TORQUE(2) 
 else 
  FCM_TORQUE(IP,2) = FCM_TORQUE(IP,2) + FCM_EXT_TORQUE(2) * 8.0 * PPI * VISC * FCM_SPHERE_RADP(IP)**3 
 end if
 
 if (mod(FCM_EXT_TORQUE(3),1.0).ne.0.0) then 
  FCM_TORQUE(IP,3) = FCM_TORQUE(IP,3) + FCM_EXT_TORQUE(3) 
 else 
  FCM_TORQUE(IP,3) = FCM_TORQUE(IP,3) + FCM_EXT_TORQUE(3) * 8.0 * PPI * VISC * FCM_SPHERE_RADP(IP)**3 
 end if

  
end do


IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND +1

 ! Define the torque from the external unique torque if any
 ! If torque is "integer" then it follows stokes law, otherwise it 
 ! is prescribed as it is
 if (mod(FCM_EXT_TORQUE(1),1.0).ne.0.0) then 
  FCM_TORQUE(IP,1) = FCM_TORQUE(IP,1) + FCM_EXT_TORQUE(1) 
 else
  FCM_TORQUE(IP,1) = FCM_TORQUE(IP,1) + FCM_EXT_TORQUE(1) * 8.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:))**3 
 end if

 if (mod(FCM_EXT_TORQUE(2),1.0).ne.0.0) then 
  FCM_TORQUE(IP,2) = FCM_TORQUE(IP,2) + FCM_EXT_TORQUE(2) 
 else 
  FCM_TORQUE(IP,2) = FCM_TORQUE(IP,2) + FCM_EXT_TORQUE(2) * 8.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:))**3 
 end if
 
 if (mod(FCM_EXT_TORQUE(3),1.0).ne.0.0) then 
  FCM_TORQUE(IP,3) = FCM_TORQUE(IP,3) + FCM_EXT_TORQUE(3) 
 else 
  FCM_TORQUE(IP,3) = FCM_TORQUE(IP,3) + FCM_EXT_TORQUE(3) * 8.0 * PPI * VISC * maxval(FCM_ELLIPSOID_RADP(IND,:))**3  
 end if
 

  
end do



   
end subroutine FCM_ADD_ROTLET_FORCING
