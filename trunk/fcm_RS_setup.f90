 !!====================================================================
!!
!! 
!!> @brief
!!> Routine generating random stresses for brownian motion
!!
!! Date :  22/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_RS_SETUP


!!====================================================================
!! Random values are obtained from "init_random_seed" defined at the bottom
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
!- RS STD's
real(kind=8) :: RS_STD_DIAG, RS_STD_NODIAG
!- Random var N(0,1)
real(kind=8) :: GAUSS_VAR

real(kind=8) :: MEAN_RS_XX, MEAN_RS_XY, MEAN_RS_XZ, &
                MEAN_RS_YY, MEAN_RS_YZ, MEAN_RS_ZZ

integer :: I ,J, K, IND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


RS_STD_DIAG = 2.0*dsqrt(KBT/(DX*DY*DZ*DTIME))
RS_STD_NODIAG = dsqrt(2.0*KBT/(DX*DY*DZ*DTIME))

do K = ISTART(3),IEND(3)
 do J = ISTART(2),IEND(2)
  do I = ISTART(1),IEND(1)
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_XX(I,J,K) = GAUSS_VAR * RS_STD_DIAG
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_YY(I,J,K) = GAUSS_VAR * RS_STD_DIAG
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_ZZ(I,J,K) = GAUSS_VAR * RS_STD_DIAG
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_XY(I,J,K) = GAUSS_VAR * RS_STD_NODIAG
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_XZ(I,J,K) = GAUSS_VAR * RS_STD_NODIAG
   
   call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
   FCM_RS_YZ(I,J,K) = GAUSS_VAR * RS_STD_NODIAG
   


  end do
 end do
end do


end subroutine FCM_RS_SETUP
