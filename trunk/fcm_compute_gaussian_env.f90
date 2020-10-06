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
    FCM_SPHERE_GAUSS1(IP,FCM_NGDH(IP)-FCM_LHNODE(IP,1)+NX/2+1:FCM_NGD(IP)) = 0.0
  end if
  if ( (FCM_LHNODE(IP,1)-FCM_NGDH(IP)+1).lt.1 ) then
    FCM_SPHERE_GAUSS1(IP,1:FCM_NGDH(IP)-FCM_LHNODE(IP,1)-1) = 0.0
  end if
 end if  
 


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


end do


end subroutine FCM_COMPUTE_GAUSSIAN_ENV
