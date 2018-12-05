 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Compute the scalar product of NPART*5 vectors
!!
!! Date :  17/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_PROD_SCAL(VEC1, VEC2, PROD)

!!====================================================================
!!
!!====================================================================

use DNS_DIM
use FCM_FORCING_VARIABLE

implicit none

!- Input argument
real(kind=8), dimension(FCM_ACTIVATE_STRESSLET,5), intent(in) :: VEC1
real(kind=8), dimension(FCM_ACTIVATE_STRESSLET,5), intent(in) :: VEC2

!- Output argument
real(kind=8), intent(out) :: PROD

! Indices for loops
integer :: I, J

PROD = 0.0

do I = 1, FCM_ACTIVATE_STRESSLET
 
 do J = 1, 5
  
  PROD  =  PROD + VEC1(I,J) * VEC2(I,J)
 
 end do
 
end do
   
end subroutine FCM_PROD_SCAL
