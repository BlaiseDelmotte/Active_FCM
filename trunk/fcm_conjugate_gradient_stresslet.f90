 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Compute the stresslet coefficient using CG from Yeo 2010
!!
!! Date :  17/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_CONJUGATE_GRADIENT_STRESSLET

!!====================================================================
!!
!!====================================================================

use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use P3DFFT

implicit none

! Indices for loops
integer :: K


K = 0


FCM_MAX_ITER = 100

! Distribution in physical space
FCM_RES_ROS = -FCM_EIJ
FCM_DIR_MIN = -FCM_EIJ

! Compute max of L2 norm
call FCM_L2NORM(FCM_RES_ROS,FCM_TOL_L2RES_ROS,FCM_L2RES_ROS)

!~ print*, 'FCM_L2RES_ROS= ',FCM_L2RES_ROS
!~ read(*,*)

do while ((FCM_L2RES_ROS.gt.FCM_TOL_L2RES_ROS).and.(K.lt.FCM_MAX_ITER))

 K = K + 1
!~   print*, 'FCM_EIJ = ', FCM_EIJ

 if (K.ge.(FCM_MAX_ITER-1)) then
  if (MYID==0) then
   print*, K, 'th CG iteration'
   print*,'MAX_L2NORM = ', FCM_L2RES_ROS
   print*,'MAX_SIJ = ', maxval(FCM_SIJ)
   print*,'MAX_FORCE = ', maxval(FCM_FORCE)
  end if
 end if


 call FCM_ZERO_FIELD_VARIABLE
 
 ! Zeros ROS
 FCM_EIJ(:,:) = 0.0
 
!~  
!~   print*, 'FCM_DIR_MIN = ', FCM_DIR_MIN

 call FCM_DISTRIB_STRESSLET(FCM_ACTIVATE_STRESSLET,FCM_DIR_MIN)
 
 if (FCM_BC==2) then
  call FCM_MIRROR_FORCES_X
 end if
 
!~  print*, 'maxval(FCM_FORCING_X),maxval(FCM_FORCING_Y),maxval(FCM_FORCING_Z) = ', maxval(FCM_FORCING_X),maxval(FCM_FORCING_Y),maxval(FCM_FORCING_Z) 
 
 call FCM_FLUID_PREDICTION
 
 call FCM_RATE_OF_STRAIN_FILTER

!~   if (MYID==0) then
!~         print*,'K= ', K
!~         print*,'maxval(FCM_EIJ)', maxval(abs(FCM_EIJ))
!~   !     read(*,*)
!~    end if

 
 call FCM_PROD_SCAL(FCM_RES_ROS,FCM_RES_ROS,FCM_INC_DIR_BETA)
 
!~   print*, 'FCM_INC_DIR_BETA= ', FCM_INC_DIR_BETA
 
 call FCM_PROD_SCAL(FCM_EIJ,FCM_DIR_MIN,FCM_INC_DIR_ALPHA)
 
!~   print*, 'FCM_INC_DIR_ALPHA= ', FCM_INC_DIR_ALPHA
 
 FCM_INC_DIR_ALPHA = FCM_INC_DIR_BETA/FCM_INC_DIR_ALPHA
 
!~   print*, 'FCM_INC_DIR_ALPHA= ', FCM_INC_DIR_ALPHA

 
 FCM_SIJ = FCM_SIJ + FCM_INC_DIR_ALPHA*FCM_DIR_MIN
 
    
 
 FCM_RES_ROS = FCM_RES_ROS - FCM_INC_DIR_ALPHA*FCM_EIJ
 
!~   print*, 'FCM_RES_ROS = ', FCM_RES_ROS

! Compute max of L2 norm
 call FCM_L2NORM(FCM_RES_ROS,FCM_TOL_L2RES_ROS,FCM_L2RES_ROS)!~  
 
 
 call FCM_PROD_SCAL(FCM_RES_ROS,FCM_RES_ROS,FCM_INC_DIR_ALPHA)
  
 FCM_INC_DIR_BETA = FCM_INC_DIR_ALPHA/FCM_INC_DIR_BETA
 
  
!~  print*, 'FCM_INC_DIR_BETA = ', FCM_INC_DIR_BETA
 
 FCM_DIR_MIN = FCM_RES_ROS + FCM_INC_DIR_BETA*FCM_DIR_MIN
 
!~  print*, 'FCM_DIR_MIN = ', FCM_DIR_MIN
 
!~    print*, 'FCM_SIJ = ', FCM_SIJ
  
  
!~   read(*,*)

 
 
end do

!if (MYID==0) then
! print*,' NITE_STRESSLET = ', K
!end if

   
end subroutine FCM_CONJUGATE_GRADIENT_STRESSLET
