 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the relative center positions bewteen particles
!!> being in a neighbourhood defined by neighbouring buckets
!! Date :  26/06/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_BARRIER_FORCE_NO_BUCKET

!!====================================================================
!! POS_REL is an antisymmetric matrix
!!====================================================================
!! Forcing: 
!!------------------------------
!! TO DO : 
!!        1) check if necessary to save POS_REL as a vector instead
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE



implicit none


! Temp variables
real(kind=8)     :: XI, YI, ZI
real(kind=8)     :: XIJ, YIJ, ZIJ
real(kind=8)     :: RIJSQ
real(kind=8)     :: RREFSQ
real(kind=8)     :: TWOA
real(kind=8)     :: FOURASQ
real(kind=8)     :: TEMP, TEMP2
real(kind=8)     :: FXI, FYI, FZI
real(kind=8)     :: FXIJ, FYIJ, FZIJ




! Indices for loops
integer :: I, J

!- Zeros necessary variables
FXI = 0.0
FYI = 0.0
FZI = 0.0

!- Zeros temp variable
FCM_FORCE_TEMP(:,:) = 0.0



do I = 1, NPART_FULL
 
  XI = FCM_XP(I)
  YI = FCM_YP(I)
  ZI = FCM_ZP(I)
  
  FXI = FCM_FORCE_TEMP(I,1)
  FYI = FCM_FORCE_TEMP(I,2)
  FZI = FCM_FORCE_TEMP(I,3)
  
  
  do J = I+1, NPART_FULL
  
   XIJ = XI - FCM_XP(J)
   YIJ = YI - FCM_YP(J)
   ZIJ = ZI - FCM_ZP(J)
   
   XIJ = XIJ - 2.0*PPI* real(int(XIJ/PPI))
   YIJ = YIJ - 2.0*PPI* real(int(YIJ/PPI))
   ZIJ = ZIJ - 2.0*PPI* real(int(ZIJ/PPI))
   
   RIJSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ;

   TWOA = FCM_SPHERE_RADP(I) + FCM_SPHERE_RADP(J)
   
   RREFSQ = FCM_FBRANGE * TWOA
   
   if (RIJSQ.lt.RREFSQ) then
    
    
    FOURASQ = TWOA**2
    TEMP = RREFSQ - FOURASQ
    TEMP2 = (RREFSQ - RIJSQ)/TEMP
    TEMP2 = TEMP2*TEMP2
    
    FXIJ = FCM_FBLEVEL*TEMP2*TEMP2*XIJ/(TWOA)
    FYIJ = FCM_FBLEVEL*TEMP2*TEMP2*YIJ/(TWOA)
    FZIJ = FCM_FBLEVEL*TEMP2*TEMP2*ZIJ/(TWOA)
    
    FXI = FXI + FXIJ
    FYI = FYI + FYIJ
    FZI = FZI + FZIJ
    
    FCM_FORCE_TEMP(J,1) = FCM_FORCE_TEMP(J,1) - FXIJ
    FCM_FORCE_TEMP(J,2) = FCM_FORCE_TEMP(J,2) - FYIJ
    FCM_FORCE_TEMP(J,3) = FCM_FORCE_TEMP(J,3) - FZIJ

   
   end if
   
  end do

 
  FCM_FORCE_TEMP(I,1) = FXI
  FCM_FORCE_TEMP(I,2) = FYI
  FCM_FORCE_TEMP(I,3) = FZI 
 
end do




! Simple addition of the contribution of each processor to the velocity average
call MPI_ALLREDUCE(FCM_FORCE_TEMP,FCM_FORCE,NPART_FULL*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

   
end subroutine FCM_BARRIER_FORCE_NO_BUCKET
