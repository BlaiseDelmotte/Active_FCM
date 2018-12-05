 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the integral of the Oseen monopole operator  Pij
!!> for self-induced VEL
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_SELF_VEL_TENSOR

!!====================================================================
!! Here the gaussian is computed independtly of particle positions
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

! Indices for loops
integer :: IP, K, J, I, IND
integer :: L, M

! Indices for position
integer :: IX, IY, IZ

! Temporary variables
real(kind=8), dimension(3,3) :: P_TENSOR
real(kind=8) :: FCM_COEFF1, FCM_COEFF2, FCM_COEFF3
real(kind=8) :: FCM_DXDYDZ
real(kind=8) :: X, Y ,Z

! Ellipsoid variables
real(kind=8) ::  KI, GAMMA, TAU
real(kind=8) ::  C, I0, I1

FCM_INT_P = 0.0

!- dx*dy*dz for numerical Riemann sum
FCM_DXDYDZ = DX*DY*DZ

! Distribution in physical space

!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSWIM(2)

 do K = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
 
  IZ = real(K-1)
  Z = IZ * DZ
  
  FCM_COEFF1 = FCM_SPHERE_ANORM(IP)**3  * FCM_DXDYDZ * dexp(-Z**2/FCM_SPHERE_SIG2SQ(IP))

	
   do J = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
   
    IY = real(J-1)
    Y = IY * DY   
			
    FCM_COEFF2 = FCM_COEFF1 * dexp(-Y**2/FCM_SPHERE_SIG2SQ(IP))
      
     do I = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
     
      IX = real(I-1)
      X = IX * DX   
			
      FCM_COEFF3 = FCM_COEFF2 * dexp(-X**2/FCM_SPHERE_SIG2SQ(IP))
      
      call FCM_OSEEN_QUADRUPOLE(X, Y, Z, FCM_SPHERE_SIGMA_DIP(IP), P_TENSOR)
      
       do L = 1,3
        do M = 1,3
        
         FCM_INT_P(IP,L,M) = FCM_INT_P(IP,L,M) + FCM_COEFF3 *  P_TENSOR(L,M)
        
        end do
       end do
      							
     end do		
		
   end do		

 end do  

end do


!!====================================================================
!! 2. ELLIPSOID SELF VEL TENSOR
!!====================================================================
IND = 0
KI = (PPI/6.0)**(-2.0/3.0) 

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1   


 C = 1.0/(32.0*FCM_ELLIPSOID_SIGMA_DIP(IND,2)**3* ( PPI*(1.0+KI)/2.0 )**(1.5))
 GAMMA = FCM_ELLIPSOID_RADP(IND,1)/FCM_ELLIPSOID_RADP(IND,2)
 TAU = dsqrt(abs(GAMMA**2 - 1.0))
 
 if (GAMMA>1.0) then
  I0 = 2.0/GAMMA
  I1 = -2.0/TAU**3*( TAU/GAMMA - log(GAMMA + TAU) )
 else
  I0 = 2.0/GAMMA
  I1 = 2.0/TAU**3*( TAU/GAMMA - atan(TAU/GAMMA) )
 end if
 
 FCM_P11(IND) = -2.0*C*(I0 - I1)
 FCM_P22(IND) = -C*(I0 + I1)

end do


end subroutine FCM_COMPUTE_SELF_VEL_TENSOR
