 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the gaussian monopole enveloppe for a given particle
!!> position
!!
!! Date :  26/02/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_GAUSSIAN_COEFF_POTDIP_ENV_FUNCTION(IND, &
                                                          XGLOB, &
                                                          Q, &
                                                          FCM_COEFF_POTDIP )

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

!- Arguments
integer, intent(in) :: IND
real(kind=8), dimension(3,1), intent(in) :: XGLOB  !< Global distance to center
real(kind=8), dimension(3,3), intent(in) :: Q  !< Rotation matrix
real(kind=8), intent(out) :: FCM_COEFF_POTDIP
!- Local parameters
real(kind=8) :: XX_DIP, YY_DIP, ZZ_DIP  !< Global distance to center
real(kind=8) :: A11, A22, A33  !< Temporary variables
    


! Coeff of the Pot Dipole Gaussian 2nd derivative  
XX_DIP = -( Q(1,1)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
          + Q(2,1)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
          + Q(3,1)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )
                                           

YY_DIP = -( Q(1,2)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
          + Q(2,2)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
          + Q(3,2)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )  
                          
                        
ZZ_DIP = -( Q(1,3)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
          + Q(2,3)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
          + Q(3,3)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) ) 
                                           

A11 = Q(1,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
    + Q(2,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
    + Q(3,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) 
    
A22 = Q(1,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
    + Q(2,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
    + Q(3,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) 
    
A33 = Q(1,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
    + Q(2,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
    + Q(3,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3)     

FCM_COEFF_POTDIP = XX_DIP**2 + YY_DIP**2 + ZZ_DIP**2 &
                 -(A11 + A22 + A33)
   
end subroutine FCM_COMPUTE_GAUSSIAN_COEFF_POTDIP_ENV_FUNCTION
