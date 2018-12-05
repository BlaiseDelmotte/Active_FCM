 !!====================================================================
!!
!! 
!!> @brief
!!> Routine removong self induced ROS created by swimming stresslets
!!
!! Date :  22/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_REMOVE_SELF_ROS

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


! Temporary variables
real(kind=8), dimension(3,3) :: GTEMP
real(kind=8), dimension(3,3) :: ROS
real(kind=8), dimension(3,3) :: Q

!!====================================================================
!! 1. SPHERICAL SQUIRMERS
!!====================================================================
do IP = 1, FCM_NSWIM(1)

 GTEMP(1,1) = FCM_SPIJ(IP,1)
 GTEMP(1,2) = FCM_SPIJ(IP,2)
 GTEMP(1,3) = FCM_SPIJ(IP,3)
 GTEMP(2,1) = FCM_SPIJ(IP,2)
 GTEMP(2,2) = FCM_SPIJ(IP,4)
 GTEMP(2,3) = FCM_SPIJ(IP,5)
 GTEMP(3,1) = FCM_SPIJ(IP,3)
 GTEMP(3,2) = FCM_SPIJ(IP,5)
 GTEMP(3,3) = -FCM_SPIJ(IP,1) - FCM_SPIJ(IP,4)
 
 do I = 1,3
	
  do J = 1, 3
    
    FCM_EIJ(IP,1) = FCM_EIJ(IP,1) - FCM_INT_R11(IP,I,J)*GTEMP(I,J)
    FCM_EIJ(IP,2) = FCM_EIJ(IP,2) - FCM_INT_R12(IP,I,J)*GTEMP(I,J)
    FCM_EIJ(IP,3) = FCM_EIJ(IP,3) - FCM_INT_R13(IP,I,J)*GTEMP(I,J)
    FCM_EIJ(IP,4) = FCM_EIJ(IP,4) - FCM_INT_R22(IP,I,J)*GTEMP(I,J)
    FCM_EIJ(IP,5) = FCM_EIJ(IP,5) - FCM_INT_R23(IP,I,J)*GTEMP(I,J)
    
  end do  

 end do
 
end do



!!====================================================================
!! 2. ELLIPSOIDAL SQUIRMERS
!!====================================================================
IND = 0

GTEMP = 0.0
ROS = 0.0
! Express Stresslet in the ellipsoid frame
if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2).eq.0)) then  
 if (FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1) then
  GTEMP(1,1)  = -1.0/3.0 * FCM_S0
  GTEMP(2,2)  = -1.0/3.0 * FCM_S0
  GTEMP(3,3)  =  2.0/3.0 * FCM_S0
 else
  GTEMP(1,1)  =  2.0/3.0 * FCM_S0
  GTEMP(2,2)  = -1.0/3.0 * FCM_S0
  GTEMP(3,3)  = -1.0/3.0 * FCM_S0
 end if
elseif (FCM_NSWIM(2)>0) then  
 if (FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1) then
  GTEMP(1,1)  = -1.0 * FCM_S0
  GTEMP(2,2)  = -1.0 * FCM_S0
  GTEMP(3,3)  =  2.0 * FCM_S0
 else
  GTEMP(1,1)  =  2.0 * FCM_S0
  GTEMP(2,2)  = -1.0 * FCM_S0
  GTEMP(3,3)  = -1.0 * FCM_S0
 end if
end if

do IP = FCM_NSPHERE + 1, NPART_FULL

  IND = IND + 1
  
  !- Transformation matrix calculated from Nikravesh 1985 Eq 7
  Q = FCM_ROT_MAT(IP,1:3,1:3)

  ! Express ROS in the ellipsoid frame
  ROS(1,1) = FCM_A1111(IND)*GTEMP(1,1) &
           + FCM_A1313(IND)*( GTEMP(2,2) + GTEMP(3,3) )
  ROS(2,2) = FCM_A2222(IND)*GTEMP(2,2) &
           + FCM_A2323(IND)*GTEMP(3,3) &
           + FCM_A1313(IND)*GTEMP(1,1) 
  ROS(3,3) = FCM_A2222(IND)*GTEMP(3,3) &
           + FCM_A2323(IND)*GTEMP(2,2) &
           + FCM_A1313(IND)*GTEMP(1,1) 
               
  ! Express ROS in the reference frame
  ROS = matmul(matmul(transpose(Q),ROS),Q)
  
!~   print*, 'ROS = ', ROS
  
  FCM_EIJ(IP,1) = FCM_EIJ(IP,1) - ROS(1,1)
  FCM_EIJ(IP,2) = FCM_EIJ(IP,2) - ROS(1,2)
  FCM_EIJ(IP,3) = FCM_EIJ(IP,3) - ROS(1,3)
  FCM_EIJ(IP,4) = FCM_EIJ(IP,4) - ROS(2,2)
  FCM_EIJ(IP,5) = FCM_EIJ(IP,5) - ROS(2,3)
  
  
end do
   

end subroutine FCM_REMOVE_SELF_ROS
