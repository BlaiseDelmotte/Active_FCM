 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the integral of the Oseen dipole operator  Rijk
!!> for self-induced ROS
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_SELF_ROS_TENSOR

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
! Indices for loops bounds
integer :: NSTART, NEND

! Indices for position
integer :: IX, IY, IZ

! Temporary variables
real(kind=8), dimension(3,3) :: R1, R2, R3
real(kind=8) :: FCM_COEFF1, FCM_COEFF2, FCM_COEFF3
real(kind=8) :: FCM_DXDYDZ
real(kind=8) :: X, Y ,Z

! Ellipsoid variables
real(kind=8) ::  KI, GAMMA, TAU
real(kind=8) ::  C, I0, I1, I2

FCM_INT_R11 = 0.0
FCM_INT_R12 = 0.0
FCM_INT_R13 = 0.0
FCM_INT_R22 = 0.0
FCM_INT_R23 = 0.0
FCM_INT_R33 = 0.0

!- dx*dy*dz for numerical Riemann sum
FCM_DXDYDZ = DX*DY*DZ

! Distribution in physical space

!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSWIM(1)

 do K = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
 
  IZ = real(K-1)
  Z = IZ * DZ
  
  FCM_COEFF1 = FCM_SPHERE_ANORM_DIP(IP)**3  * FCM_DXDYDZ  &
             * dexp(-Z**2/FCM_SPHERE_SIG2SQ_DIP(IP)) /FCM_SPHERE_SIGSQ_DIP(IP)

	
   do J = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
   
    IY = real(J-1)
    Y = IY * DY   
			
    FCM_COEFF2 = FCM_COEFF1 * dexp(-Y**2/FCM_SPHERE_SIG2SQ_DIP(IP))
      
     do I = -FCM_NGDH(IP)+1, FCM_NGDH(IP)
     
      IX = real(I-1)
      X = IX * DX   
			
      FCM_COEFF3 = FCM_COEFF2 * dexp(-X**2/FCM_SPHERE_SIG2SQ_DIP(IP))
      
      if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2).eq.0)) then
       call FCM_OSEEN_DIPOLE(X, Y, Z, FCM_SPHERE_SIGMA_DIP(IP), R1, R2, R3) 
      else if (FCM_NSWIM(2)>0) then
       call FCM_OSEEN_DIPOLE(X, Y, Z, FCM_SPHERE_SIGMA(IP), R1, R2, R3)       
      end if
           
       do L = 1,3
        do M = 1,3
        
         FCM_INT_R11(IP,L,M) = FCM_INT_R11(IP,L,M) + FCM_COEFF3 * X * R1(L,M)
         FCM_INT_R12(IP,L,M) = FCM_INT_R12(IP,L,M) &
                             + 0.5*FCM_COEFF3 * ( Y*R1(L,M) + X*R2(L,M) )
         FCM_INT_R13(IP,L,M) = FCM_INT_R13(IP,L,M) &
                             + 0.5*FCM_COEFF3 * ( Z*R1(L,M) + X*R3(L,M) )                 
         FCM_INT_R22(IP,L,M) = FCM_INT_R22(IP,L,M) + FCM_COEFF3 * Y * R2(L,M)
         FCM_INT_R23(IP,L,M) = FCM_INT_R23(IP,L,M) &
                             + 0.5*FCM_COEFF3 * ( Y*R3(L,M) + Z*R2(L,M) )
         FCM_INT_R33(IP,L,M) = FCM_INT_R33(IP,L,M) + FCM_COEFF3 * Z * R3(L,M)
        
        
        end do
       end do
      							
     end do		
		
   end do		

 end do  

end do


!!====================================================================
!! 2. ELLIPSOID SELF ROS TENSOR
!!====================================================================
IND = 0

 if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2).eq.0)) then
  KI = 1.0  
 elseif (FCM_NSWIM(2)>0) then
  KI = (PPI/6.0)**(-2.0/3.0) 
 end if

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1
 
 C = 1.0/( 32.0*FCM_ELLIPSOID_SIGMA_DIP(IND,2)**3*(PPI*(1.0+KI)/2.0 )**(1.5) )
 GAMMA = FCM_ELLIPSOID_RADP(IND,1)/FCM_ELLIPSOID_RADP(IND,2)
 TAU = dsqrt(abs(GAMMA**2 - 1.0))
 
 if (GAMMA>1.0) then
  I0 = 2.0/GAMMA
  I1 = -2.0/TAU**3*( TAU/GAMMA - log(GAMMA + TAU) )
  I2 = 1.0/TAU**5*( (GAMMA**2+2.0)*TAU/GAMMA - 3.0*log(GAMMA + TAU) )  
 else
  I0 = 2.0/GAMMA
  I1 = 2.0/TAU**3*( TAU/GAMMA - atan(TAU/GAMMA) )
  I2 = 1.0/TAU**5*( (GAMMA**2+2.0)*TAU/GAMMA - 3.0*atan(TAU/GAMMA) )
 end if
 
 FCM_A1111(IND) = -2.0*C*(                      I1 -         I2)
 FCM_A2222(IND) = -C*    ( 1.0/4.0*I0 + 1.0/2.0*I1 - 3.0/4.0*I2)
 FCM_A1313(IND) =  C*    (                      I1 -         I2)
 FCM_A2323(IND) = -C*    (-1.0/2.0*I0 +         I1 - 1.0/2.0*I2)/2.0
 
end do
!~ print*,'C =', C
!~ print*,'KI =', KI
!~ print*,'GAMMA =', GAMMA 
!~ print*,'FCM_A1111 =', FCM_A1111
!~ print*,'FCM_A2222 =', FCM_A2222
!~ print*,'FCM_A1313 =', FCM_A1313
!~ print*,'FCM_A2323 =', FCM_A2323
!~ print*,'FCM_ELLIPSOID_SIGMA(IND,2) =', FCM_ELLIPSOID_SIGMA(IND,2)
end subroutine FCM_COMPUTE_SELF_ROS_TENSOR
