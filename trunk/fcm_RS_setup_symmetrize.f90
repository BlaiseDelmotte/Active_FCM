 !!====================================================================
!!
!! 
!!> @brief
!!> Routine generating random stresses for brownian motion
!!> taking into account symmetries due to image system
!! Date :  04/09/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_RS_SETUP_SYMMETRIZE


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

real(kind=8) :: MEAN_RS_XX, MEAN_RS_XY, MEAN_RS_XZ, &
                MEAN_RS_YY, MEAN_RS_YZ, MEAN_RS_ZZ
!- Random var N(0,1)
real(kind=8) :: GAUSS_VAR


integer :: I ,J, K, IND

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IND = 0
MEAN_RS_XX = 0.0
MEAN_RS_XY = 0.0
MEAN_RS_XZ = 0.0
MEAN_RS_YY = 0.0
MEAN_RS_YZ = 0.0
MEAN_RS_ZZ = 0.0

RS_STD_DIAG = 2.0d0*dsqrt(KBT/(DX*DY*DZ*DTIME))
RS_STD_NODIAG = dsqrt(2.0d0*KBT/(DX*DY*DZ*DTIME))

do K = ISTART(3),IEND(3)
 do J = ISTART(2),IEND(2)
  do I = 2,NX/2
  !do I = NX/2+2,NX
   ! Fill the first half of the domain
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
   
   ! Use symmetries to fill the image part
   FCM_RS_XX(NX+2-I,J,K) = FCM_RS_XX(I,J,K)
   FCM_RS_XY(NX+2-I,J,K) = -FCM_RS_XY(I,J,K)
   FCM_RS_XZ(NX+2-I,J,K) = -FCM_RS_XZ(I,J,K)
   FCM_RS_YY(NX+2-I,J,K) = FCM_RS_YY(I,J,K)
   FCM_RS_YZ(NX+2-I,J,K) = FCM_RS_YZ(I,J,K)
   FCM_RS_ZZ(NX+2-I,J,K) = FCM_RS_ZZ(I,J,K)


  end do
  
  ! Prescribe the values at the boundaries
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_XX(1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_XX(NX/2+1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
  
  FCM_RS_XY(1,J,K) = 0.0
  FCM_RS_XY(NX/2+1,J,K) = 0.0
  
  FCM_RS_XZ(1,J,K) = 0.0
  FCM_RS_XZ(NX/2+1,J,K) = 0.0
  
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_YY(1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_YY(NX/2+1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
   
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_YZ(1,J,K) = GAUSS_VAR * RS_STD_NODIAG *dsqrt(2.0d0)
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_YZ(NX/2+1,J,K) = GAUSS_VAR * RS_STD_NODIAG *dsqrt(2.0d0)
 
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_ZZ(1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
  call FCM_RANDOM_GAUSSIAN_VARIABLE(GAUSS_VAR)
  FCM_RS_ZZ(NX/2+1,J,K) = GAUSS_VAR * RS_STD_DIAG *dsqrt(2.0d0)
  
 end do
end do


!MEAN_RS_XX = MEAN_RS_XX/IND
!MEAN_RS_XY = MEAN_RS_XY/IND
!MEAN_RS_XZ = MEAN_RS_XZ/IND
!MEAN_RS_YY = MEAN_RS_YY/IND
!MEAN_RS_YZ = MEAN_RS_YZ/IND
!MEAN_RS_ZZ = MEAN_RS_ZZ/IND

!MEAN_RS_XX_TIME = MEAN_RS_XX_TIME + MEAN_RS_XX
!MEAN_RS_XY_TIME = MEAN_RS_XY_TIME + MEAN_RS_XY
!MEAN_RS_XZ_TIME = MEAN_RS_XZ_TIME + MEAN_RS_XZ 
!MEAN_RS_YY_TIME = MEAN_RS_YY_TIME + MEAN_RS_YY
!MEAN_RS_YZ_TIME = MEAN_RS_YZ_TIME + MEAN_RS_YZ
!MEAN_RS_ZZ_TIME = MEAN_RS_ZZ_TIME + MEAN_RS_ZZ

!print*,MYID
!print*,'MEAN_RS_XX = ', MEAN_RS_XX 
!print*,'MEAN_RS_XY = ', MEAN_RS_XY 
!print*,'MEAN_RS_XZ = ', MEAN_RS_XZ 
!print*,'MEAN_RS_YY = ', MEAN_RS_YY 
!print*,'MEAN_RS_YZ = ', MEAN_RS_YZ 
!print*,'MEAN_RS_ZZ = ', MEAN_RS_ZZ 

!print*,'MEAN_RS_XX_TIME = ', MEAN_RS_XX_TIME
!print*,'MEAN_RS_XY_TIME = ', MEAN_RS_XY_TIME 
!print*,'MEAN_RS_XZ_TIME = ', MEAN_RS_XZ_TIME 
!print*,'MEAN_RS_YY_TIME = ', MEAN_RS_YY_TIME 
!print*,'MEAN_RS_YZ_TIME = ', MEAN_RS_YZ_TIME 
!print*,'MEAN_RS_ZZ_TIME = ', MEAN_RS_ZZ_TIME 


end subroutine FCM_RS_SETUP_SYMMETRIZE
