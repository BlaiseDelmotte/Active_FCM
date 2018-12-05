 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the gaussian monopole enveloppe for a given particle
!!> position
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_GAUSSIAN_SUPPORT

!!====================================================================
!! Here the gaussian is computed according to the current particle position
!!====================================================================
!! Gaussian monopole computation: 
!!------------------------------
!! TO DO : 
!!        1) If ellipsoid then else endif
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none


! Indices for loops
integer :: IP, ND, IND

!!====================================================================
!! 2. ELLIPSOID GAUSSIAN ENVELOPPE
!!====================================================================

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1

 !- Nearest node of particle center
 FCM_LHNODE(IP,1) = int( FCM_XP(IP)/DX ) 
 FCM_LHNODE(IP,2) = int( FCM_YP(IP)/DY ) 
 FCM_LHNODE(IP,3) = int( FCM_ZP(IP)/DZ ) 
 
  
 do ND = 1, FCM_ELLIPSOID_NGD
 
  ! Gaussian support mapped on the cartesian grid  (very smart formula as periodicity is taken into account)
  FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH + ND 
  FCM_ELLIPSOID_IZP(IND,ND) = (mod(FCM_JZ,NZ) + NZ*(1-isign(1,FCM_JZ))/2)+1
    
  FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH + ND
  FCM_ELLIPSOID_IYP(IND,ND) = (mod(FCM_JY,NY) + NY*(1-isign(1,FCM_JY))/2)+1   
   
  FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH + ND  
  FCM_ELLIPSOID_IXP(IND,ND) = (mod(FCM_JX,NX) + NX*(1-isign(1,FCM_JX))/2)+1
 
 end do !do ND = 1, FCM_ELLIPSOID_NGD

end do !do IP = NSPHERE + 1, NPART_FULL
!~ 
   
end subroutine FCM_COMPUTE_GAUSSIAN_SUPPORT
