 !!====================================================================
!!
!! 
!!> @brief
!!> Routine generating random forces for RFD
!!
!! Date :  26/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_RF_SETUP


!!====================================================================
!! CAREFUL : BY NOW CAN ONLY BE USED WITH ONE PROCESSOR !!!!
!! NEED TO INSTALL SPRNG                             
!!====================================================================
use DNS_DIM
use PARAM_PHYS
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!======================================================================

!- Local variables
!- Random var N(0,1)
real(kind=8) :: GAUSS_VAR

integer :: I ,J

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do I = 1, NPART_FULL
 do J = 1,3

  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RFD_RAND_FORCE(I,J) = GAUSS_VAR/FCM_EPS_RFD

 end do
end do

end subroutine FCM_RF_SETUP
