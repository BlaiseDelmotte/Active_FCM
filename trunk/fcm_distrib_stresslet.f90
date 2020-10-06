 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Distribute the stresslet on the gaussian envelope, it computes 
!!>  \f$ \sum_{\alpha=1}^{N_{p}}\left(S_{ij}^{\alpha}\frac{\partial}{\partial x_{j}}\Delta_{D}\left(\mathbf{x}-\mathbf{Y}^{\alpha}(t)\right)\right) \f$. 
!!
!!
!! Date :  17/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_DISTRIB_STRESSLET(DIM1,FCM_STRESS)

!!====================================================================
!! Here the stresslet forcing is distributed on the  
!! dipole Gaussian enveloppe.
!!------------------------------
!!====================================================================

use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use P3DFFT

implicit none
integer, intent(in) :: DIM1
real(kind=8), dimension(DIM1,5), intent(in) :: FCM_STRESS

! Indices for loops
integer :: IP, K, J, I, IND

! Indices for position
integer :: IX, IY, IZ

! Temporary variables
! Rotation matrix
real(kind=8), dimension(3,3) :: Q
! Position from particle center
real(kind=8), dimension(3,1) :: XLOC
real(kind=8), dimension(3,1) :: XGLOB
real(kind=8) :: FCM_STEMP1X, FCM_STEMP2X, FCM_STEMP3X 
real(kind=8) :: FCM_STEMP1Y, FCM_STEMP2Y, FCM_STEMP3Y 
real(kind=8) :: FCM_STEMP1Z, FCM_STEMP2Z, FCM_STEMP3Z 
real(kind=8) :: FCM_COEFF1, FCM_COEFF2, FCM_COEFF3

! Initiation
TMPPHY(:,:,:)  = 0d0
TMPPHY2(:,:,:) = 0d0
TMPPHY3(:,:,:) = 0d0

! Distribution in physical space


!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, DIM1

 do K = 1, FCM_NGD(IP)
 
  IZ = FCM_SPHERE_IZP(IP,K)
  
  ! Only compute the force projection for the piece of gaussian within the block
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then 
  
   FCM_COEFF1 = FCM_SPHERE_DIP_GAUSS3(IP,K) 
   
   FCM_STEMP1X = FCM_STRESS(IP,3) * FCM_SPHERE_GRAD_DIP_GAUSS3(IP,K)   
   FCM_STEMP1Y = FCM_STRESS(IP,5) * FCM_SPHERE_GRAD_DIP_GAUSS3(IP,K)   
   FCM_STEMP1Z = -( FCM_STRESS(IP,4)+FCM_STRESS(IP,1) ) * FCM_SPHERE_GRAD_DIP_GAUSS3(IP,K)
                
      
   do J = 1, FCM_NGD(IP)
   
    IY = FCM_SPHERE_IYP(IP,J)
      
    ! Only compute the force projection for the piece of gaussian within the block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
         
     FCM_COEFF2  = FCM_SPHERE_DIP_GAUSS2(IP,J)*FCM_COEFF1
     
     FCM_STEMP2X = FCM_STEMP1X + FCM_STRESS(IP,2) * FCM_SPHERE_GRAD_DIP_GAUSS2(IP,J)
	    FCM_STEMP2Y = FCM_STEMP1Y + FCM_STRESS(IP,4) * FCM_SPHERE_GRAD_DIP_GAUSS2(IP,J)
     FCM_STEMP2Z = FCM_STEMP1Z + FCM_STRESS(IP,5) * FCM_SPHERE_GRAD_DIP_GAUSS2(IP,J) 
			
     do I = 1, FCM_NGD(IP)
     
      IX = FCM_SPHERE_IXP(IP,I)
      
      FCM_COEFF3 = FCM_SPHERE_DIP_GAUSS1(IP,I)*FCM_COEFF2
      
      FCM_STEMP3X = FCM_STEMP2X + FCM_STRESS(IP,1) * FCM_SPHERE_GRAD_DIP_GAUSS1(IP,I)
	     FCM_STEMP3Y = FCM_STEMP2Y + FCM_STRESS(IP,2) * FCM_SPHERE_GRAD_DIP_GAUSS1(IP,I)
      FCM_STEMP3Z = FCM_STEMP2Z + FCM_STRESS(IP,3) * FCM_SPHERE_GRAD_DIP_GAUSS1(IP,I)
      
           
      TMPPHY(IX,IY,IZ)  = TMPPHY(IX,IY,IZ)  & 
                        + FCM_COEFF3 * FCM_STEMP3X
						
      TMPPHY2(IX,IY,IZ) = TMPPHY2(IX,IY,IZ) & 
                        + FCM_COEFF3 * FCM_STEMP3Y
                        
      TMPPHY3(IX,IY,IZ) = TMPPHY3(IX,IY,IZ) & 
                        + FCM_COEFF3 * FCM_STEMP3Z								
                        
     end do		
     		
    endif			
   end do		
  endif	
 end do  
end do


!!====================================================================
!! 2. ELLIPSOID GAUSSIAN ENVELOPPE
!!====================================================================

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1 
 
 !- Transformation matrix calculated from Nikravesh 1985 Eq 7
 Q = FCM_ROT_MAT(IP,1:3,1:3)

 do K = 1, FCM_ELLIPSOID_NGD
 
  IZ = FCM_ELLIPSOID_IZP(IND,K)
  
  ! Only compute the force projection for the piece of gaussian within the block
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then 
  
   FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH + K      
   XLOC(3,1) = dble(FCM_JZ)*DZ-FCM_ZP(IP)
     
   do J = 1, FCM_ELLIPSOID_NGD
   
    IY = FCM_ELLIPSOID_IYP(IND,J)
      
    ! Only compute the force projection for the piece of gaussian within the block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
    
     FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH + J  
     XLOC(2,1) = dble(FCM_JY)*DY-FCM_YP(IP)
       			
     do I = 1, FCM_ELLIPSOID_NGD
     
      IX = FCM_ELLIPSOID_IXP(IND,I)   
      FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH + I     
      XLOC(1,1) = dble(FCM_JX)*DX-FCM_XP(IP)   
      
      XGLOB = matmul(Q,XLOC)

      call FCM_COMPUTE_GAUSSIAN_DIP_ENV_FUNCTION(IND, &
                                                 XGLOB, &
                                                 FCM_COEFF3 )
                                                 
      call FCM_COMPUTE_GAUSSIAN_GRAD_DIP_ENV_FUNCTION(IND, &
                                                 XGLOB, &
                                                 Q, &                                                
                                                 FCM_STEMP1X, &
                                                 FCM_STEMP1Y, &
                                                 FCM_STEMP1Z  )                                           
                                                 
      
      FCM_STEMP3X = FCM_STRESS(IP,3) * FCM_STEMP1Z &
                  + FCM_STRESS(IP,2) * FCM_STEMP1Y &
                  + FCM_STRESS(IP,1) * FCM_STEMP1X
                  
      FCM_STEMP3Y = FCM_STRESS(IP,5) * FCM_STEMP1Z &
                  + FCM_STRESS(IP,4) * FCM_STEMP1Y &
                  + FCM_STRESS(IP,2) * FCM_STEMP1X
                  
      FCM_STEMP3Z = - ( FCM_STRESS(IP,4)+FCM_STRESS(IP,1) ) * FCM_STEMP1Z &
                    + FCM_STRESS(IP,5) * FCM_STEMP1Y &
                    + FCM_STRESS(IP,3) * FCM_STEMP1X
                 
                 
                  
           
      TMPPHY(IX,IY,IZ)  = TMPPHY(IX,IY,IZ)  & 
                        + FCM_COEFF3 * FCM_STEMP3X 
                        
      TMPPHY2(IX,IY,IZ) = TMPPHY2(IX,IY,IZ) & 
                        + FCM_COEFF3 * FCM_STEMP3Y 
                        
      TMPPHY3(IX,IY,IZ) = TMPPHY3(IX,IY,IZ) & 
                        + FCM_COEFF3 * FCM_STEMP3Z 			
                        
     end do		
     		
    endif			
   end do		
  endif	
 end do  
end do


! Summation to compute total forcing
FCM_FORCING_X = FCM_FORCING_X + TMPPHY
FCM_FORCING_Y = FCM_FORCING_Y + TMPPHY2
FCM_FORCING_Z = FCM_FORCING_Z + TMPPHY3


   
end subroutine FCM_DISTRIB_STRESSLET
