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

subroutine FCM_L2NORM(VEC, TOL, NORM_MAX)

!!====================================================================
!!
!!====================================================================

use DNS_DIM
use FCM_FORCING_VARIABLE

implicit none

!- Input argument
real(kind=8), dimension(FCM_ACTIVATE_STRESSLET,5), intent(in) :: VEC
real(kind=8), intent(in) :: TOL

!- Output argument
real(kind=8), intent(out) :: NORM_MAX

real(kind=8) :: NORM_PREV
! Indices for loops
integer :: I, J

NORM_MAX = 0.0
NORM_PREV = 0.0

do I = 1, FCM_ACTIVATE_STRESSLET


 do J = 1, 5
  
  NORM_PREV  = NORM_PREV + 2.0*(VEC(I,J))**2 
  
 end do
 
 ! + E33^2 - E11^2 - E22^2
 NORM_PREV  = NORM_PREV + ( -VEC(I,1)-VEC(I,4) )**2 - VEC(I,1)**2 - VEC(I,4)**2
 
 if (NORM_PREV>NORM_MAX) then
  NORM_MAX = NORM_PREV
  if (NORM_MAX>TOL**2) then
   exit
  end if 
 end if
 

 
end do

NORM_MAX = sqrt(NORM_MAX)

   
end subroutine FCM_L2NORM
