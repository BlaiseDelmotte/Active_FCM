 !!====================================================================
!!
!! 
!!> @brief
!!> Routine removong self induced VEL created by swimming potential dipole
!!
!! Date :  22/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_REMOVE_SELF_VEL

!!====================================================================
!! Here thehe fluid vorticity field is computed with VORTICITY
!! and  filtered with the Gaussian dipole enveloppe to get the particle
!! rotational velocities
!!====================================================================
!! Rotationfiltering: 
!!------------------------------
!!
!! TO DO : 
!!        1) If ellipsoid
!!------------------------------
!! WARNING: We use temporary variables before performing the sum of each
!! proc contribution
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none

!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------

! Indices for loops
integer :: IP, J, I, IND


!~ FCM_UP(:,1) = 0.0
!~ FCM_VP(:,1) = 0.0
!~ FCM_WP(:,1) = 0.0

!!====================================================================
!! 1. SPHERICAL SQUIRMERS
!!====================================================================
do IP = 1, FCM_NSWIM(2)	

 do J = 1, 3
   
   FCM_UP(IP,1) = FCM_UP(IP,1) - FCM_INT_P(IP,1,J)*FCM_HI(IP,J)
   FCM_VP(IP,1) = FCM_VP(IP,1) - FCM_INT_P(IP,2,J)*FCM_HI(IP,J)
   FCM_WP(IP,1) = FCM_WP(IP,1) - FCM_INT_P(IP,3,J)*FCM_HI(IP,J)

  
 end do 
   
end do

!!====================================================================
!! 2. ELLIPSOIDAL SQUIRMERS
!!====================================================================
IND = 0
if (FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1) then
 do IP = FCM_NSPHERE + 1, NPART_FULL
  IND = IND + 1
  FCM_UP(IP,1) = FCM_UP(IP,1) - FCM_P22(IND)*FCM_HI(IP,1)
  FCM_VP(IP,1) = FCM_VP(IP,1) - FCM_P22(IND)*FCM_HI(IP,2)
  FCM_WP(IP,1) = FCM_WP(IP,1) - FCM_P22(IND)*FCM_HI(IP,3)  
 end do
else
 do IP = FCM_NSPHERE + 1, NPART_FULL
  IND = IND + 1
  FCM_UP(IP,1) = FCM_UP(IP,1) - FCM_P11(IND)*FCM_HI(IP,1)
  FCM_VP(IP,1) = FCM_VP(IP,1) - FCM_P11(IND)*FCM_HI(IP,2)
  FCM_WP(IP,1) = FCM_WP(IP,1) - FCM_P11(IND)*FCM_HI(IP,3)  
 end do
end if
   
end subroutine FCM_REMOVE_SELF_VEL
