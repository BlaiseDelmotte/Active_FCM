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
  
  
!~   ! Gaussian support mapped on the cartesian grid  (very smart formula as periodicity is taken into account)
!~   FCM_SPHERE_IXP_SQ(IP,ND) = (mod(FCM_JX,NX) + NX*(1-isign(1,FCM_JX))/2)+1
!~   FCM_SPHERE_IYP_SQ(IP,ND) = (mod(FCM_JY,NY) + NY*(1-isign(1,FCM_JY))/2)+1
!~   FCM_SPHERE_IZP_SQ(IP,ND) = (mod(FCM_JZ,NZ) + NZ*(1-isign(1,FCM_JZ))/2)+1	
  
  
  XX = dble(FCM_JX)*DX-FCM_XP(IP)
  YY = dble(FCM_JY)*DY-FCM_YP(IP)
  ZZ = dble(FCM_JZ)*DZ-FCM_ZP(IP)
  
    
!~   ! Gaussian functions centered on yp(ip,:)   cf. Eq.6 and Eq.9 p177 Climent Maxey
!~   FCM_SPHERE_GAUSS1_SQ(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( &
!~                             -XX**2/FCM_SPHERE_SIG2SQ(IP) )
!~   FCM_SPHERE_GAUSS2_SQ(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( & 
!~                             -YY**2/FCM_SPHERE_SIG2SQ(IP) )
!~   FCM_SPHERE_GAUSS3_SQ(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( &
!~                             -ZZ**2/FCM_SPHERE_SIG2SQ(IP) )


  ! Factor of the Monopole Gaussian derivative
  FCM_SPHERE_GRAD_GAUSS1_SQ(IP,ND) = -XX /FCM_SPHERE_SIGSQ(IP)
  FCM_SPHERE_GRAD_GAUSS2_SQ(IP,ND) = -YY /FCM_SPHERE_SIGSQ(IP)
  FCM_SPHERE_GRAD_GAUSS3_SQ(IP,ND) = -ZZ /FCM_SPHERE_SIGSQ(IP)
  
  
  ! Coeff of the Pot Dipole Gaussian 2nd derivative
  FCM_COEFF_SPHERE_POTDIP_GAUSS1_SQ(IP,ND) = XX**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_COEFF_SPHERE_POTDIP_GAUSS2_SQ(IP,ND) = YY**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_COEFF_SPHERE_POTDIP_GAUSS3_SQ(IP,ND) = ZZ**2 /FCM_SPHERE_SIGSQ_DIP(IP)
  

 end do 
  
!~  ! Dipole enveloppe can be computed directly from the monopole enveloppe cf. Thesis Liu P 99
!~  FCM_SPHERE_POTDIP_GAUSS1_SQ(IP,:) = ( & 
!~                                           FCM_SPHERE_GAUSS1_SQ(IP,:) & 
!~                                        * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
!~                                        * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
!~                                      )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP)) &
!~                                      / FCM_SPHERE_ANORM_POTDIP(IP)
!~                                 
!~  FCM_SPHERE_POTDIP_GAUSS2_SQ(IP,:) = ( &
!~                                           FCM_SPHERE_GAUSS2_SQ(IP,:) &
!~                                        * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
!~                                        * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
!~                                      )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP)) &
!~                                      / FCM_SPHERE_ANORM_POTDIP(IP)
!~                                 
!~  FCM_SPHERE_POTDIP_GAUSS3_SQ(IP,:) = ( &
!~                                           FCM_SPHERE_GAUSS3_SQ(IP,:) &
!~                                        * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
!~                                        * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
!~                                      )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP)) &
!~                                      / FCM_SPHERE_ANORM_POTDIP(IP)


end do


!!!!! NOT NEEDED ANYMORE AS IT IS DONE INSIDE THE DISTRIBUTION/FILTER ROUTINES
!~ !!====================================================================
!~ !! 2. ELLIPSOID SQUIRMING GAUSSIAN ENVELOPPE
!~ !!====================================================================
!~ 
!~ IND = 0
!~ 
!~ do IP = FCM_NSPHERE + 1, NPART_FULL
!~ 
!~  IND = IND + 1
!~  
!~  !- Transformation matrix calculated from Nikravesh 1985 Eq 7
!~  Q(1,1) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,2)**2 -0.5 
!~  Q(1,2) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,3) - FCM_QUAT(IP,1)*FCM_QUAT(IP,4) 
!~  Q(1,3) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,4) + FCM_QUAT(IP,1)*FCM_QUAT(IP,3) 
!~  
!~  Q(2,1) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,3) + FCM_QUAT(IP,1)*FCM_QUAT(IP,4)  
!~  Q(2,2) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,3)**2 -0.5 
!~  Q(2,3) =  FCM_QUAT(IP,3)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,2) 
!~  
!~  Q(3,1) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,3) 
!~  Q(3,2) =  FCM_QUAT(IP,3)*FCM_QUAT(IP,4) + FCM_QUAT(IP,1)*FCM_QUAT(IP,2) 
!~  Q(3,3) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,4)**2 -0.5 
!~  
!~  Q = 2.0 * Q
!~   
!~  do NDZ = 1, FCM_ELLIPSOID_NGD
!~  
!~   ! Gaussian support mapped on the cartesian grid  (very smart formula as periodicity is taken into account)
!~   FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH+ NDZ  
 !FCM_ELLIPSOID_IZP_SQ(IND,NDZ) = (mod(FCM_JZ,NZ) + NZ*(1-isign(1,FCM_JZ))/2)+1
!~     
!~   XLOC(3,1) = dble(FCM_JZ)*DZ-FCM_ZP(IP)
!~   
!~   do NDY = 1, FCM_ELLIPSOID_NGD
!~   
!~    FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH+ NDY
  ! FCM_ELLIPSOID_IYP_SQ(IND,NDY) = (mod(FCM_JY,NY) + NY*(1-isign(1,FCM_JY))/2)+1
!~    
!~    XLOC(2,1) = dble(FCM_JY)*DY-FCM_YP(IP)
!~    
!~    do NDX = 1, FCM_ELLIPSOID_NGD
!~    
!~     FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH+ NDX    
   ! FCM_ELLIPSOID_IXP_SQ(IND,NDX) = (mod(FCM_JX,NX) + NX*(1-isign(1,FCM_JX))/2)+1
!~     
!~     XLOC(1,1) = dble(FCM_JX)*DX-FCM_XP(IP)
!~     
!~     XGLOB = matmul(Q,XLOC)
!~     
    !- Gaussian functions centered on yb(ib,:)   cf. Eq.6 and Eq.9 p177 Climent Maxey
   ! FCM_ELLIPSOID_GAUSS_SQ(IND,NDX,NDY,NDZ) = FCM_ELLIPSOID_ANORM(IND)*dexp( &
   !                                          -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,1) &
   !                                          -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,2) &
   !                                          -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,3) )
!~                                        
!~                                        
   ! !- Dipole Gaussian functions centered on yb(ib,:)
   ! FCM_ELLIPSOID_POTDIP_GAUSS_SQ(IND,NDX,NDY,NDZ) = FCM_ELLIPSOID_ANORM_DIP(IND)*dexp( &
   !                                           -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,1)   &
   !                                           -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,2)   & 
   !                                           -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,3) ) 
!~ 
!~ 
!~    ! Factor of the Monopole Gaussian derivative
!~     XX = ( Q(1,1)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ(IND,1) &
!~          + Q(2,1)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ(IND,2) &
!~          + Q(3,1)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ(IND,3) )
!~     FCM_ELLIPSOID_GRAD_GAUSS1_SQ(IND,NDX,NDY,NDZ) = -XX
!~                                                
!~     YY = ( Q(1,2)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ(IND,1) &
!~          + Q(2,2)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ(IND,2) &
!~          + Q(3,2)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ(IND,3) ) 
!~     FCM_ELLIPSOID_GRAD_GAUSS2_SQ(IND,NDX,NDY,NDZ) = -YY
!~                                                
!~     ZZ = ( Q(1,3)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ(IND,1) &
!~          + Q(2,3)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ(IND,2) &
!~          + Q(3,3)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ(IND,3) )                                            
!~     FCM_ELLIPSOID_GRAD_GAUSS3_SQ(IND,NDX,NDY,NDZ) = -ZZ
!~                                                      
!~                                                      
!~     ! Coeff of the Pot Dipole Gaussian 2nd derivative    
!~     A11 = Q(1,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~         + Q(2,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~         + Q(3,1)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) 
!~     A22 = Q(1,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~         + Q(2,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~         + Q(3,2)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) 
!~     A33 = Q(1,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~         + Q(2,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~         + Q(3,3)**2/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) 
!~         
!~     XX_DIP = ( Q(1,1)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~              + Q(2,1)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~              + Q(3,1)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )
!~          
!~          
!~     YY_DIP = ( Q(1,2)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~              + Q(2,2)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~              + Q(3,2)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )
!~          
!~     ZZ_DIP = ( Q(1,3)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~              + Q(2,3)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~              + Q(3,3)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) ) 
!~              
!~     FCM_COEFF_ELLIPSOID_POTDIP_GAUSS_SQ(IND,NDX,NDY,NDZ) = XX_DIP**2 &
!~                                                          + YY_DIP**2 &
!~                                                          + ZZ_DIP**2 &
!~                                                          - (A11 + A22 + A33)
!~ 
!~   
!~    end do !do NDX = 1, FCM_ELLIPSOID_NGD
!~ 
!~   end do !do NDY = 1, FCM_ELLIPSOID_NGD
!~  
!~  end do !do NDZ = 1, FCM_ELLIPSOID_NGD
!~  
!~ 
!~ 
!~ end do !do IP = NSPHERE + 1, NPART_FULL




   
end subroutine FCM_COMPUTE_GAUSSIAN_ENV_SQ
