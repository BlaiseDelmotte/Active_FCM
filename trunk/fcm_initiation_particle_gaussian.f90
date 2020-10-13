 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating the chraracteristic values for Gaussian functions
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PARTICLE_GAUSSIAN

!!====================================================================
!! Here the gaussian width and normalization factors are initiated
!! according to the prescribed size of the particle
!!====================================================================
!! Gaussian initiation: 
!!------------------------------
!! TO DO : 
!!        1) Add polydispersity
!!        2) Add calculation volumic fraction, maxradp etc....
!!        3) If ellipsoid then else endif
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none

! Indices for loops
integer :: IP, IND


FCM_NGDH = FCM_NGD  /2.0

FCM_ELLIPSOID_NGDH = FCM_ELLIPSOID_NGD /2.0

FCM_FRACTION = 0.0

FCM_SPHERE_SIGMA(1:FCM_NSPHERE_1)      = FCM_SIGP * DX * FCM_SPHERE_SIZE(1) 
FCM_SPHERE_SIGMA(FCM_NSPHERE_1+1:FCM_NSPHERE)      = FCM_SIGP * DX * FCM_SPHERE_SIZE(2) 
FCM_SPHERE_RADP(1:FCM_NSPHERE)       = dsqrt(PPI) * FCM_SPHERE_SIGMA
FCM_SPHERE_SIGMA_DIP(1:FCM_NSPHERE)  = FCM_SPHERE_RADP/( (36.0*PPI)**(1.0/6.0) )
FCM_SPHERE_ANORM(1:FCM_NSPHERE)      = 1.0/( dsqrt(TWOPI)*FCM_SPHERE_SIGMA )
FCM_SPHERE_ANORM_DIP(1:FCM_NSPHERE)  = 1.0/( dsqrt(TWOPI)*FCM_SPHERE_SIGMA_DIP )
FCM_SPHERE_SIG2SQ(1:FCM_NSPHERE)     = 2.0  * FCM_SPHERE_SIGMA**2
FCM_SPHERE_SIGSQ_DIP(1:FCM_NSPHERE)  = FCM_SPHERE_SIGMA_DIP**2
FCM_SPHERE_SIG2SQ_DIP(1:FCM_NSPHERE) = 2.0  * FCM_SPHERE_SIGMA_DIP**2


if (FCM_NSWIM(2)>0) then
 FCM_SPHERE_SIGSQ(1:FCM_NSWIM(2))       = FCM_SPHERE_SIGMA(1:FCM_NSWIM(2))**2
 FCM_SPHERE_ANORM_POTDIP(1:FCM_NSWIM(2)) = FCM_SPHERE_SIGMA_DIP(1:FCM_NSWIM(2))**(2.0/3.0) 
end if

do IP = 1, FCM_NSPHERE
 FCM_FRACTION = FCM_FRACTION + 4.0/3.0*PPI*FCM_SPHERE_RADP(IP)**3
end do


IND = 0
do IP = FCM_NSPHERE + 1, NPART_FULL
 IND = IND + 1
 FCM_ELLIPSOID_SIGMA(IND,1:3)      = FCM_SIGP * DX * FCM_ELLIPSOID_SIZE(1:3) ! TO DO for polydispersity *SIZEP(IP) ! cf. Climent Maxey p 186
 FCM_ELLIPSOID_RADP(IND,1:3)       = dsqrt(PPI) * FCM_ELLIPSOID_SIGMA(IND,1:3)
 FCM_ELLIPSOID_SIGMA_DIP(IND,1:3)  = FCM_ELLIPSOID_RADP(IND,1:3)/( (36.0*PPI)**(1.0/6.0) )
 FCM_ELLIPSOID_ANORM(IND)          = 1.0/( dsqrt(TWOPI)**3 &
                                         *FCM_ELLIPSOID_SIGMA(IND,1)&
                                         *FCM_ELLIPSOID_SIGMA(IND,2)&
                                         *FCM_ELLIPSOID_SIGMA(IND,3) )
 FCM_ELLIPSOID_ANORM_DIP(IND)      = 1.0/( dsqrt(TWOPI)**3 &
                                         *FCM_ELLIPSOID_SIGMA_DIP(IND,1)&
                                         *FCM_ELLIPSOID_SIGMA_DIP(IND,2)&
                                         *FCM_ELLIPSOID_SIGMA_DIP(IND,3) )
 FCM_ELLIPSOID_SIG2SQ(IND,1:3)     = 2.0  * FCM_ELLIPSOID_SIGMA(IND,1:3)**2
 FCM_ELLIPSOID_SIGSQ_DIP(IND,1:3)  = FCM_ELLIPSOID_SIGMA_DIP(IND,1:3)**2
 FCM_ELLIPSOID_SIG2SQ_DIP(IND,1:3) = 2.0  * FCM_ELLIPSOID_SIGMA_DIP(IND,1:3)**2
 
 FCM_FRACTION = FCM_FRACTION  &
              + 4.0/3.0*PPI*FCM_ELLIPSOID_RADP(IND,1)&
               *FCM_ELLIPSOID_RADP(IND,2)&
               *FCM_ELLIPSOID_RADP(IND,3)
              
 if ((FCM_NSWIM(2)>0).and.(IND.le.FCM_NSWIM(2))) then
  FCM_ELLIPSOID_SIGSQ(IND,:)        = FCM_ELLIPSOID_SIGMA(IND,:)**2
 end if
 
end do

FCM_FRACTION = FCM_FRACTION/(LXMAX*LYMAX*LZMAX)


if (MYID == 0) then

 print*,'Volumetric fraction = ', FCM_FRACTION
 
 
end if
   
end subroutine FCM_INITIATION_PARTICLE_GAUSSIAN
