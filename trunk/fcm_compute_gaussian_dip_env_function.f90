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

subroutine FCM_COMPUTE_GAUSSIAN_DIP_ENV_FUNCTION(IND, &
                                                 XGLOB, &
                                                 FCM_ELL_DIP_GAUSS )

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
real(kind=8), intent(out) :: FCM_ELL_DIP_GAUSS  

    
!- Gaussian functions centered on yb(ib,:)   cf. Eq.6 and Eq.9 p177 Climent Maxey
!~ FCM_ELL_GAUSS = FCM_ELLIPSOID_ANORM(IND)*dexp( &
!~                                    -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,1) &
!~                                    -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,2) &
!~                                    -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,3) )

!- Dipole Gaussian functions centered on yb(ib,:)
FCM_ELL_DIP_GAUSS = FCM_ELLIPSOID_ANORM_DIP(IND)*dexp( &
                             -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,1) &
                             -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,2) &
                             -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,3) )
                                       

!~ !- Factor of the Dipole_Gaussian derivative
!~ FCM_ELL_GRAD_DIP_GAUSS1 = -( Q(1,1)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                  + Q(2,1)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                  + Q(3,1)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )
!~                                            
!~ 
!~ FCM_ELL_GRAD_DIP_GAUSS2 = -( Q(1,2)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                  + Q(2,2)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                  + Q(3,2)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )  
!~                           
!~                         
!~ FCM_ELL_GRAD_DIP_GAUSS3 = -( Q(1,3)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                  + Q(2,3)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                  + Q(3,3)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) ) 
                                           

    

!~ 
   
end subroutine FCM_COMPUTE_GAUSSIAN_DIP_ENV_FUNCTION
