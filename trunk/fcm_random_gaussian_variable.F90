 !!====================================================================
!!
!! 
!!> @brief
!!> Routine generating random gaussian variable 
!!> Use Marsaglia Polar Method (WIKIPEDIA)
!!
!! Date :  22/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_RANDOM_GAUSSIAN_VARIABLE(VAR)

!!====================================================================
!! Random values are obtained from "init_random_seed" defined at the bottom
!!====================================================================
!! CAREFUL : BY NOW CAN ONLY BE USED WITH ONE PROCESSOR !!!!
!! NEED TO INSTALL SPRNG                             
!!====================================================================
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE

implicit none

#define SIMPLE_SPRNG	
#define USE_MPI
#include "sprng_f.h"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


real(kind=8),          intent(out) :: VAR


!!======================================================================

!- Local variables
real(kind=8) ::  X1, X2, W


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

W = 10.0

do while ( (W>1.0).or.(W==1) )

 X1 = sprng()
 X2 = sprng()
 
 X1 = 2.0*X1 -1.0
 X2 = 2.0*X2 -1.0
 W = X1**2 + X2**2
end do

W = dsqrt( (-2.0*log(W))/W )
VAR = W*X1

end subroutine FCM_RANDOM_GAUSSIAN_VARIABLE
