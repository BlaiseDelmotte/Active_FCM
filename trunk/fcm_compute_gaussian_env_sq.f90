 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the gaussian potential dipole enveloppe
!!> and squirming stresslet enveloppe for a given particle
!!> position
!!
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_GAUSSIAN_ENV_SQ

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

!- Local parameters
real, dimension(3,3) :: Q      !< Transformation matrix from body frame to global coordinate
real, dimension(3,1) :: XLOC   !< Local distance to center
real, dimension(3,1) :: XGLOB  !< Global distance to center

real(kind=8) :: XX, YY, ZZ  !< Global distance to center
real(kind=8) :: XX_DIP, YY_DIP, ZZ_DIP  !< Global distance to center
real(kind=8) :: A11, A22, A33  !< Temporary variables

! Indices for loops
integer :: IP, ND, NDX, NDY, NDZ, IND

!!====================================================================
!! 1. SPHERE SQUIRMING GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSWIM(2)

 
 do ND = 1, FCM_NGD(IP)
 
  FCM_JX = FCM_LHNODE(IP,1) - FCM_NGDH(IP)+ ND
  FCM_JY = FCM_LHNODE(IP,2) - FCM_NGDH(IP)+ ND
  FCM_JZ = FCM_LHNODE(IP,3) - FCM_NGDH(IP)+ ND
  
  
  XX = dble(FCM_JX)*DX-FCM_XP(IP)
  YY = dble(FCM_JY)*DY-FCM_YP(IP)
  ZZ = dble(FCM_JZ)*DZ-FCM_ZP(IP)
  

  ! Factor of the Monopole Gaussian derivative
  FCM_SPHERE_GRAD_GAUSS1_SQ(IP,ND) = -XX /FCM_SPHERE_SIGSQ(IP)
  FCM_SPHERE_GRAD_GAUSS2_SQ(IP,ND) = -YY /FCM_SPHERE_SIGSQ(IP)
  FCM_SPHERE_GRAD_GAUSS3_SQ(IP,ND) = -ZZ /FCM_SPHERE_SIGSQ(IP)
  
  
  ! Coeff of the Pot Dipole Gaussian 2nd derivative
  FCM_COEFF_SPHERE_POTDIP_GAUSS1_SQ(IP,ND) = XX**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_COEFF_SPHERE_POTDIP_GAUSS2_SQ(IP,ND) = YY**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_COEFF_SPHERE_POTDIP_GAUSS3_SQ(IP,ND) = ZZ**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  

 end do 
  

end do


   
end subroutine FCM_COMPUTE_GAUSSIAN_ENV_SQ
