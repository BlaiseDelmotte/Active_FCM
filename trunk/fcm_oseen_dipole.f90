 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the Oseen FCM tensor Rijk giving the fluid velocity induced 
!!> by a force dipole distributed on a monopole gaussian function 
!!> (i.e. Squirming stresslet) CF LOHMOLT MAXEY 2003
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_OSEEN_DIPOLE(X, Y, Z, SIGMA, R1, R2, R3)


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
real(kind=8), dimension(3,3), intent(out) :: R1, R2, R3


! Temporary variables
real(kind=8),dimension(3) :: XVEC
real(kind=8) :: SIGSQ 
real(kind=8) :: ANORM
real(kind=8) :: R, RSQ
real(kind=8) :: DADRTEMP, DBDRTEMP, BTEMP
real(kind=8) :: AF, BF, DADRF, DBDRF
real(kind=8) :: GERF 
real(kind=8) :: DELTA

! Indices for loops
integer ::  J, I


!- Variable initiation
XVEC(1) = X
XVEC(2) = Y
XVEC(3) = Z

SIGSQ = SIGMA * SIGMA
ANORM = 1.0/dsqrt(TWOPI)
RSQ = X*X + Y*Y + Z*Z
R = dsqrt(RSQ)

R1 = 0.0
R2 = 0.0
R3 = 0.0


!- Compute tensor given by Lohmolt & Maxey IJMF 2003 Eq. 22b
if (abs(R/SIGMA)<0.01) then

 AF = ANORM/(PPI * SIGMA)
 DADRF = -2.0*AF/(15.0*SIGSQ)
 DBDRF = -32.0*AF/(35.0*SIGSQ*SIGSQ)
 BF = AF/(30.0*SIGSQ) ! B(0), Eq 39 Maxey Patel
 
 do J = 1,3
 
  do I = 1,3
   
   R1(J,I) = R1(J,I) + DBDRF*XVEC(1)*XVEC(I)*XVEC(J)
   R2(J,I) = R2(J,I) + DBDRF*XVEC(2)*XVEC(I)*XVEC(J)
   R3(J,I) = R3(J,I) + DBDRF*XVEC(3)*XVEC(I)*XVEC(J)
   
   if (j == 1) then
     R1(J,I) = R1(J,I) + DADRF*XVEC(I)
   end if
   if (j == 2) then
     R2(J,I) = R2(J,I) + DADRF*XVEC(I)
   end if
   if (j == 3) then
     R3(J,I) = R3(J,I) + DADRF*XVEC(I)
   end if
   if (i == 1) then
     R1(J,I) = R1(J,I) + BF*XVEC(J)
   end if
   if (i == 2) then
     R2(J,I) = R2(J,I) + BF*XVEC(J)
   end if
   if (i == 3) then
     R3(J,I) = R3(J,I) + BF*XVEC(J)
   end if
   
  end do
   
  R1(J,J) = R1(J,J) + BF*XVEC(1)
  R2(J,J) = R2(J,J) + BF*XVEC(2)
  R3(J,J) = R3(J,J) + BF*XVEC(3)
  
 end do
 
else

 DELTA = ANORM * dexp(-0.5*RSQ/SIGSQ)
 GERF = derf(R/(SIGMA*dsqrt(2.0D0)))
 DADRTEMP = (1.0 + 3.0*SIGSQ/RSQ)*GERF - (4.0*R/SIGMA + 6.0*SIGMA/R)*DELTA
 BTEMP = (1.0 - 3.0*SIGSQ/RSQ)*GERF + (6.0*SIGMA/R)*DELTA
 DBDRTEMP = (3.0 - 15.0*SIGSQ/RSQ)*GERF + (4.0*R/SIGMA + 30.0*SIGMA/R)*DELTA
 DADRF = -DADRTEMP/(8.0*PPI*RSQ)
 BF = BTEMP/(8.0*PPI*R*RSQ)
 DBDRF = -DBDRTEMP/(8.0*PPI*RSQ*RSQ)
 
 do J = 1,3
 
  do I = 1,3
   
   R1(J,I) = R1(J,I) + DBDRF*XVEC(1)*XVEC(I)*XVEC(J)/R
   R2(J,I) = R2(J,I) + DBDRF*XVEC(2)*XVEC(I)*XVEC(J)/R
   R3(J,I) = R3(J,I) + DBDRF*XVEC(3)*XVEC(I)*XVEC(J)/R
   
   if (j == 1) then
     R1(J,I) = R1(J,I) + DADRF*XVEC(I)/R
   end if
   if (j == 2) then
     R2(J,I) = R2(J,I) + DADRF*XVEC(I)/R
   end if
   if (j == 3) then
     R3(J,I) = R3(J,I) + DADRF*XVEC(I)/R
   end if
   if (i == 1) then
     R1(J,I) = R1(J,I) + BF*XVEC(J)
   end if
   if (i == 2) then
     R2(J,I) = R2(J,I) + BF*XVEC(J)
   end if
   if (i == 3) then
     R3(J,I) = R3(J,I) + BF*XVEC(J)
   end if
   
  end do
  
  R1(J,J) = R1(J,J) + BF*XVEC(1)
  R2(J,J) = R2(J,J) + BF*XVEC(2)
  R3(J,J) = R3(J,J) + BF*XVEC(3)
  
 end do
 
end if

   
end subroutine FCM_OSEEN_DIPOLE
