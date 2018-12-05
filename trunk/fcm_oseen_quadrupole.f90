 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the Oseen FCM tensor Pij giving 
!!> flow disturbances of a quadrupole cf. Maxey 2001 IJMF Eq.46
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_OSEEN_QUADRUPOLE(X, Y, Z, SIGMA, P_TENSOR)


!!====================================================================
!! Rijk computation: 
!!------------------------------
!! TO DO : 
!!        1) If ellipsoid then else endif
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none

!- ---------------------------------------------------------------------
! ROUTINE ARGUMENTS
!- ---------------------------------------------------------------------
real(kind=8), intent(in) :: X, Y, Z
real(kind=8), intent(in) :: SIGMA
real(kind=8), dimension(3,3), intent(out) :: P_TENSOR


! Temporary variables
real(kind=8),dimension(3) :: XVEC
real(kind=8) :: SIGSQ 
real(kind=8) :: ANORM
real(kind=8) :: R, RSQ
real(kind=8) :: AF, BF, CF, DF
real(kind=8) :: GERF 
real(kind=8) :: DELTA

! Indices for loops
integer ::  J, I


!- Variable initiation
XVEC(1) = X
XVEC(2) = Y
XVEC(3) = Z

SIGSQ = SIGMA * SIGMA
ANORM = 1.0/dsqrt(TWOPI*SIGSQ)
RSQ = X*X + Y*Y + Z*Z
R = dsqrt(RSQ)

P_TENSOR = 0.0


!- Compute tensor given by  Maxey & Patel IJMF 2001  Eq. 22b
if (abs(R/SIGMA)<0.01) then
!~     r=r+1.0e-16  !! Commented in Keaveny's code, why ?
    AF = -2.0 * ANORM * ANORM * ANORM
    AF = AF/(3.0)
    do J = 1,3
      P_TENSOR(J,J) = AF + P_TENSOR(J,J)
    end do
else
    DELTA = ANORM*ANORM*ANORM*dexp(-0.5*RSQ/SIGSQ)
    GERF = derf(R / (SIGMA*dsqrt(2.0D0)))
    
    AF = GERF/(4.0*PPI*R*RSQ)
    BF = -3.0*GERF/(4.0*PPI*R*RSQ*RSQ)
    CF = (1.0 + (SIGMA/R)*(SIGMA/R))*DELTA
    DF = (-1.0 - 3.0*(SIGMA/R)*(SIGMA/R))*DELTA/RSQ

    AF = AF - CF
    BF = BF - DF

    do J = 1,3
     do I = 1,3
     
      P_TENSOR(I,J) = BF*XVEC(I)*XVEC(J)
     
     end do
     
     P_TENSOR(J,J) = AF + P_TENSOR(J,J)
     
    end do
end if

   
end subroutine FCM_OSEEN_QUADRUPOLE
