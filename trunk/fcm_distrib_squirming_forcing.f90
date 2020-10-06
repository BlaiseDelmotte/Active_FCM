 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Distribute the squirming stresslet on the gradient monopole gaussian envelope, 
!!> and the squirming potential dipole on the 2nd derivative of the 
!!> dipole gaussian enveloppe
!!
!!
!! Date :  22/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_DISTRIB_SQUIRMING_FORCING

!!====================================================================
!! Here the squirming stresslet forcing is distributed on the  
!! gradient of monopole Gaussian enveloppe, and the squirming potential 
!! dipole on the 2nd derivative of the dipole gaussian enveloppe
!!
!!====================================================================
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
real(kind=8) :: FCM_HTEMP1, FCM_HTEMP2, FCM_HTEMP3 
real(kind=8) :: FCM_SCOEFF1, FCM_SCOEFF2, FCM_SCOEFF3
real(kind=8) :: FCM_HCOEFF1, FCM_HCOEFF2, FCM_HCOEFF3
real(kind=8) :: G11, G12, G13, G21, G22, G23, G31, G32, G33
real(kind=8) :: XX, YY, ZZ
real(kind=8) :: XX2, YY2, ZZ2

! Initiation
TMPPHY(:,:,:)  = 0d0
TMPPHY2(:,:,:) = 0d0
TMPPHY3(:,:,:) = 0d0

! Distribution in physical space


!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSWIM(2)

 G11 = FCM_SPIJ(IP,1)
 G12 = FCM_SPIJ(IP,2)
 G13 = FCM_SPIJ(IP,3)
 G21 = FCM_SPIJ(IP,2)
 G22 = FCM_SPIJ(IP,4)
 G23 = FCM_SPIJ(IP,5)
 G31 = FCM_SPIJ(IP,3)
 G32 = FCM_SPIJ(IP,5)
 G33 = -G11 - G22
 

 do K = 1, FCM_NGD(IP)
 
  IZ = FCM_SPHERE_IZP(IP,K)
  
  
  
  ! Only compute the force projection for the piece of gaussian within the block
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then 
  
   
   FCM_HCOEFF1 = FCM_SPHERE_DIP_GAUSS3(IP,K) / FCM_SPHERE_ANORM_POTDIP(IP)**3
   FCM_SCOEFF1 = FCM_SPHERE_GAUSS3(IP,K) 

   ZZ = FCM_SPHERE_GRAD_GAUSS3_SQ(IP,K)
   ZZ2 = FCM_COEFF_SPHERE_POTDIP_GAUSS3_SQ(IP,K)
   
   FCM_HTEMP1 = ZZ2 - 3.0
   
   FCM_STEMP1X = G13*ZZ
   FCM_STEMP1Y = G23*ZZ
   FCM_STEMP1Z = G33*ZZ
      
   do J = 1, FCM_NGD(IP)
   
    IY = FCM_SPHERE_IYP(IP,J)
    
    
      
    ! Only compute the force projection for the piece of gaussian within the block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
     
     FCM_HCOEFF2 = FCM_SPHERE_DIP_GAUSS2(IP,J)*FCM_HCOEFF1
     FCM_SCOEFF2 = FCM_SPHERE_GAUSS2(IP,J)*FCM_SCOEFF1     
     
     YY = FCM_SPHERE_GRAD_GAUSS2_SQ(IP,J)
     YY2 = FCM_COEFF_SPHERE_POTDIP_GAUSS2_SQ(IP,J)
     
     FCM_HTEMP2 = FCM_HTEMP1 + YY2
     
     FCM_STEMP2X = FCM_STEMP1X + G12*YY
     FCM_STEMP2Y = FCM_STEMP1Y + G22*YY
     FCM_STEMP2Z = FCM_STEMP1Z + G32*YY
         
			
     do I = 1, FCM_NGD(IP)
     
      IX = FCM_SPHERE_IXP(IP,I)
      
      
      FCM_HCOEFF3 = FCM_SPHERE_DIP_GAUSS1(IP,I)*FCM_HCOEFF2
      FCM_SCOEFF3 = FCM_SPHERE_GAUSS1(IP,I)*FCM_SCOEFF2
      
      XX = FCM_SPHERE_GRAD_GAUSS1_SQ(IP,I)
      XX2 = FCM_COEFF_SPHERE_POTDIP_GAUSS1_SQ(IP,I)
      
      FCM_HTEMP3 = FCM_HTEMP2 + XX2
      
      FCM_STEMP3X = FCM_STEMP2X + G11*XX
      FCM_STEMP3Y = FCM_STEMP2Y + G21*XX
      FCM_STEMP3Z = FCM_STEMP2Z + G31*XX
           
      TMPPHY(IX,IY,IZ)  = TMPPHY(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,1) &
                        + FCM_SCOEFF3 * FCM_STEMP3X
						
      TMPPHY2(IX,IY,IZ) = TMPPHY2(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,2) &
                        + FCM_SCOEFF3 * FCM_STEMP3Y
                        
      TMPPHY3(IX,IY,IZ) = TMPPHY3(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,3) &
                        + FCM_SCOEFF3 * FCM_STEMP3Z				
                 
     					
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

 G11 = FCM_SPIJ(IP,1)
 G12 = FCM_SPIJ(IP,2)
 G13 = FCM_SPIJ(IP,3)
 G21 = FCM_SPIJ(IP,2)
 G22 = FCM_SPIJ(IP,4)
 G23 = FCM_SPIJ(IP,5)
 G31 = FCM_SPIJ(IP,3)
 G32 = FCM_SPIJ(IP,5)
 G33 = -G11 - G22

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
      
      call FCM_COMPUTE_GAUSSIAN_GRAD_ENV_FUNCTION(IND, &
                                                  XGLOB, &
                                                  Q, &                                                
                                                  XX, &
                                                  YY, &
                                                  ZZ  )                                                   
                                                  
      call FCM_COMPUTE_GAUSSIAN_ENV_FUNCTION(IND, XGLOB, FCM_SCOEFF3)  
      
      
      call FCM_COMPUTE_GAUSSIAN_DIP_ENV_FUNCTION(IND, XGLOB, FCM_HCOEFF3)

      call FCM_COMPUTE_GAUSSIAN_COEFF_POTDIP_ENV_FUNCTION(IND, &
                                                          XGLOB, &
                                                          Q, &
                                                          FCM_HTEMP3 ) 
  
      FCM_STEMP3X = G11*XX + G12*YY + G13*ZZ 
      FCM_STEMP3Y = G21*XX + G22*YY + G23*ZZ 
      FCM_STEMP3Z = G31*XX + G32*YY + G33*ZZ       
                  
           
      TMPPHY(IX,IY,IZ)  = TMPPHY(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,1) &
                        + FCM_SCOEFF3 * FCM_STEMP3X
						
      TMPPHY2(IX,IY,IZ) = TMPPHY2(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,2) &
                        + FCM_SCOEFF3 * FCM_STEMP3Y
                        
      TMPPHY3(IX,IY,IZ) = TMPPHY3(IX,IY,IZ)  & 
                        + FCM_HCOEFF3 * FCM_HTEMP3 * FCM_HI(IP,3) &
                        + FCM_SCOEFF3 * FCM_STEMP3Z												
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


   
end subroutine FCM_DISTRIB_SQUIRMING_FORCING
