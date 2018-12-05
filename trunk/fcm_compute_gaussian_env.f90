 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the gaussian monopole enveloppe for a given particle
!!> position
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_GAUSSIAN_ENV

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

! Indices for loops
integer :: IP, ND, NDX, NDY, NDZ, IND


!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSPHERE

! Nearest node of particle center
 FCM_LHNODE(IP,1) = int( FCM_XP(IP)/DX ) 
 FCM_LHNODE(IP,2) = int( FCM_YP(IP)/DY ) 
 FCM_LHNODE(IP,3) = int( FCM_ZP(IP)/DZ ) 
 
 do ND = 1, FCM_NGD(IP)
 
  FCM_JX = FCM_LHNODE(IP,1) - FCM_NGDH(IP) + ND
  FCM_JY = FCM_LHNODE(IP,2) - FCM_NGDH(IP) + ND
  FCM_JZ = FCM_LHNODE(IP,3) - FCM_NGDH(IP) + ND
  
  ! Gaussian support mapped on the cartesian grid  (very smart formula as periodicity is taken into account)
  FCM_SPHERE_IXP(IP,ND) = (mod(FCM_JX,NX) + NX*(1-isign(1,FCM_JX))/2)+1
  FCM_SPHERE_IYP(IP,ND) = (mod(FCM_JY,NY) + NY*(1-isign(1,FCM_JY))/2)+1
  FCM_SPHERE_IZP(IP,ND) = (mod(FCM_JZ,NZ) + NZ*(1-isign(1,FCM_JZ))/2)+1	
  
  XX = dble(FCM_JX)*DX-FCM_XP(IP)
  YY = dble(FCM_JY)*DY-FCM_YP(IP)
  ZZ = dble(FCM_JZ)*DZ-FCM_ZP(IP)
    
  ! Gaussian functions centered on yp(ip,:)   cf. Eq.6 and Eq.9 p177 Climent Maxey
  FCM_SPHERE_GAUSS1(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( &
                            -XX**2/FCM_SPHERE_SIG2SQ(IP) )
  FCM_SPHERE_GAUSS2(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( & 
                            -YY**2/FCM_SPHERE_SIG2SQ(IP) )
  FCM_SPHERE_GAUSS3(IP,ND) = FCM_SPHERE_ANORM(IP)*dexp( &
                            -ZZ**2/FCM_SPHERE_SIG2SQ(IP) )
  
  ! Factor of the Dipole_Gaussian derivative
  FCM_SPHERE_GRAD_DIP_GAUSS1(IP,ND) = -XX /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_SPHERE_GRAD_DIP_GAUSS2(IP,ND) = -YY /FCM_SPHERE_SIGSQ_DIP(IP)
  FCM_SPHERE_GRAD_DIP_GAUSS3(IP,ND) = -ZZ /FCM_SPHERE_SIGSQ_DIP(IP)

 end do 
 
 !- Truncate (=0) the Gaussian enveloppes when going beyond the boundaries.
 !- Gaussians are non-zero on the boundary itself
 if (FCM_BC == 2) then
  if ( (FCM_LHNODE(IP,1)+FCM_NGDH(IP)).gt.(NX/2+1) ) then
    FCM_SPHERE_GAUSS1(IP,FCM_NGDH(IP)-FCM_LHNODE(IP,1)+NX/2+1:ND) = 0.0
  end if
  if ( (FCM_LHNODE(IP,1)-FCM_NGDH(IP)+1).lt.1 ) then
    FCM_SPHERE_GAUSS1(IP,1:FCM_NGDH(IP)-FCM_LHNODE(IP,1)-1) = 0.0
  end if
 end if  
 
!~   
!~  print*,'FCM_LHNODE(IP,1) = ', FCM_LHNODE(IP,1)
!~  print*,'FCM_LHNODE(IP,1)+FCM_NGDH = ', FCM_LHNODE(IP,1)+FCM_NGDH
!~  print*,'FCM_LHNODE(IP,1)-FCM_NGDH+1 = ', FCM_LHNODE(IP,1)-FCM_NGDH-1
!~  print*,'FCM_NGDH-FCM_LHNODE(IP,1)+NX/2+1 = ', FCM_NGDH-FCM_LHNODE(IP,1)+NX/2+1
!~  print*,'FCM_NGDH-FCM_LHNODE(IP,1) = ', FCM_NGDH-FCM_LHNODE(IP,1)
!~  print*,'FCM_SPHERE_IXP(IP,FCM_NGDH-FCM_LHNODE(IP,1)+NX/2+1) = ', FCM_SPHERE_IXP(IP,FCM_NGDH-FCM_LHNODE(IP,1)+NX/2+1)
!~  print*,'FCM_SPHERE_IXP(IP,FCM_NGDH-FCM_LHNODE(IP,1)-1) = ', FCM_SPHERE_IXP(IP,FCM_NGDH-FCM_LHNODE(IP,1)-1)
!~  print*,'FCM_SPHERE_GAUSS1(IP,:) = ', FCM_SPHERE_GAUSS1(IP,:)


 ! Dipole enveloppe can be computed directly from the monopole enveloppe cf. Thesis Liu P 99
 FCM_SPHERE_DIP_GAUSS1(IP,:) = ( & 
                                FCM_SPHERE_GAUSS1(IP,:) & 
                             * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
                             * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
                                )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP))
                                
 FCM_SPHERE_DIP_GAUSS2(IP,:) = ( &
                                FCM_SPHERE_GAUSS2(IP,:) &
                             * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
                             * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
                                )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP))
                                
 FCM_SPHERE_DIP_GAUSS3(IP,:) = ( &
                                FCM_SPHERE_GAUSS3(IP,:) &
                             * (FCM_SPHERE_SIGMA(IP)/FCM_SPHERE_SIGMA_DIP(IP)) &
                             * (1.0/FCM_SPHERE_ANORM_DIP(IP))**(1.0-FCM_SPHERE_SIG2SQ_DIP(IP)/FCM_SPHERE_SIG2SQ(IP)) &
                                )**(FCM_SPHERE_SIG2SQ(IP)/FCM_SPHERE_SIG2SQ_DIP(IP))


!~  print*,'FCM_SPHERE_DIP_GAUSS1(IP,:) = ', FCM_SPHERE_DIP_GAUSS1(IP,:)
!~  read(*,*)
end do




!!!!! NOT NEEDED ANYMORE AS IT IS DONE INSIDE THE DISTRIBUTION/FILTER ROUTINES
!~ !!====================================================================
!~ !! 2. ELLIPSOID GAUSSIAN ENVELOPPE
!~ !!====================================================================
!~ 
!~ IND = 0
!~ 
!~ do IP = FCM_NSPHERE + 1, NPART_FULL
!~ 
!~  IND = IND + 1
!~ 
!~  !- Nearest node of particle center
!~  FCM_LHNODE(IP,1) = int( FCM_XP(IP)/DX ) 
!~  FCM_LHNODE(IP,2) = int( FCM_YP(IP)/DY ) 
!~  FCM_LHNODE(IP,3) = int( FCM_ZP(IP)/DZ ) 
!~  
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
!~   FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH + NDZ  
!~   FCM_ELLIPSOID_IZP(IND,NDZ) = (mod(FCM_JZ,NZ) + NZ*(1-isign(1,FCM_JZ))/2)+1
!~     
!~   XLOC(3,1) = dble(FCM_JZ)*DZ-FCM_ZP(IP)
!~   
!~   do NDY = 1, FCM_ELLIPSOID_NGD
!~   
!~    FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH + NDY
!~    FCM_ELLIPSOID_IYP(IND,NDY) = (mod(FCM_JY,NY) + NY*(1-isign(1,FCM_JY))/2)+1
!~    
!~    XLOC(2,1) = dble(FCM_JY)*DY-FCM_YP(IP)
!~    
!~    do NDX = 1, FCM_ELLIPSOID_NGD
!~    
!~     FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH + NDX    
!~     FCM_ELLIPSOID_IXP(IND,NDX) = (mod(FCM_JX,NX) + NX*(1-isign(1,FCM_JX))/2)+1
!~     
!~     XLOC(1,1) = dble(FCM_JX)*DX-FCM_XP(IP)
!~     
!~     XGLOB = matmul(Q,XLOC)
!~     
!~     !- Gaussian functions centered on yb(ib,:)   cf. Eq.6 and Eq.9 p177 Climent Maxey
!~     FCM_ELLIPSOID_GAUSS(IND,NDX,NDY,NDZ) = FCM_ELLIPSOID_ANORM(IND)*dexp( &
!~                                          -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,1) &
!~                                          -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,2) &
!~                                          -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ(IND,3) )
!~                                        
    !print*, 'NDX, NDY, NDZ = ', NDX, NDY, NDZ       
!~     
   ! print*, 'XGLOB = ', XGLOB                                      
   ! print*, 'FCM_ELLIPSOID_GAUSS(IP,NDX,NDY,NDZ) = ', FCM_ELLIPSOID_GAUSS(IP,NDX,NDY,NDZ)
!~                                        
!~     !- Dipole Gaussian functions centered on yb(ib,:)
!~     FCM_ELLIPSOID_DIP_GAUSS(IND,NDX,NDY,NDZ) = FCM_ELLIPSOID_ANORM_DIP(IND)*dexp( &
!~                                        -(XGLOB(1,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,1) &
!~                                        -(XGLOB(2,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,2) &
!~                                        -(XGLOB(3,1))**2/FCM_ELLIPSOID_SIG2SQ_DIP(IND,3) )
!~                                        
  !  print*, 'FCM_ELLIPSOID_DIP_GAUSS(IP,NDX,NDY,NDZ) = ', FCM_ELLIPSOID_DIP_GAUSS(IP,NDX,NDY,NDZ)
!~   
!~     !- Factor of the Dipole_Gaussian derivative
!~     FCM_ELLIPSOID_GRAD_DIP_GAUSS1(IND,NDX,NDY,NDZ) = -( Q(1,1)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                                       + Q(2,1)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                                       + Q(3,1)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )
!~                                                
   ! print*, 'FCM_ELLIPSOID_GRAD_DIP_GAUSS1(IP,NDX,NDY,NDZ) = ', FCM_ELLIPSOID_GRAD_DIP_GAUSS1(IP,NDX,NDY,NDZ)
!~ 
!~     FCM_ELLIPSOID_GRAD_DIP_GAUSS2(IND,NDX,NDY,NDZ) = -( Q(1,2)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                                       + Q(2,2)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                                       + Q(3,2)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) )  
!~                                                
    !print*, 'FCM_ELLIPSOID_GRAD_DIP_GAUSS2(IP,NDX,NDY,NDZ) = ', FCM_ELLIPSOID_GRAD_DIP_GAUSS2(IP,NDX,NDY,NDZ)
!~                                                
!~     FCM_ELLIPSOID_GRAD_DIP_GAUSS3(IND,NDX,NDY,NDZ) = -( Q(1,3)*XGLOB(1,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,1) &
!~                                                       + Q(2,3)*XGLOB(2,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,2) &
!~                                                       + Q(3,3)*XGLOB(3,1)/FCM_ELLIPSOID_SIGSQ_DIP(IND,3) ) 
!~                                                
   ! print*, 'FCM_ELLIPSOID_GRAD_DIP_GAUSS3(IP,NDX,NDY,NDZ) = ', FCM_ELLIPSOID_GRAD_DIP_GAUSS3(IP,NDX,NDY,NDZ)
    
    
    !read(*,*)
!~     
!~    end do !do NDX = 1, FCM_ELLIPSOID_NGD
!~ 
!~   end do !do NDY = 1, FCM_ELLIPSOID_NGD
!~  
!~  end do !do NDZ = 1, FCM_ELLIPSOID_NGD
!~ 
!~ end do !do IP = NSPHERE + 1, NPART_FULL
!~ 
   
end subroutine FCM_COMPUTE_GAUSSIAN_ENV
