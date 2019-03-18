 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Compute the L2 norm of NPART*5 vector
!!
!! Date :  17/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_L2NORM(VEC, NORM)

!!====================================================================
!!
!!====================================================================

use DNS_DIM
use FCM_FORCING_VARIABLE

implicit none

!- Input argument
real(kind=8), dimension(FCM_ACTIVATE_STRESSLET,5), intent(in) :: VEC

!- Output argument
real(kind=8), intent(out) :: NORM

! Indices for loops
integer :: I, J

NORM = 0.0

do I = 1, FCM_ACTIVATE_STRESSLET


 do J = 1, 5
  
  NORM  = NORM + 2.0*(VEC(I,J))**2 
  
 end do
 
 ! + E33^2 - E11^2 - E22^2
 NORM  = NORM + ( -VEC(I,1)-VEC(I,4) )**2 - VEC(I,1)**2 - VEC(I,4)**2
 
 
end do

NORM = sqrt(NORM)

   
end subroutine FCM_L2NORM
