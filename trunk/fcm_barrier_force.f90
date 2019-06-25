 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the repulsive force to avoid overlapping with 
!!> bucket sorting strategy (cf Keaveny listed list)
!! Date :  10/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_BARRIER_FORCE

!!====================================================================
!! 
!!====================================================================
!! Forcing: 
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use FCM_BUCKET_VARIABLE



implicit none


! Temp variables
real(kind=8)     :: XI, YI, ZI
real(kind=8)     :: XIJ, YIJ, ZIJ
real(kind=8)     :: RIJSQ
real(kind=8)     :: RIJ
real(kind=8)     :: RREFSQ
real(kind=8)     :: TWOA
real(kind=8)     :: FOURASQ
real(kind=8)     :: TEMP, TEMP2
real(kind=8)     :: FXI, FYI, FZI
real(kind=8)     :: FXIJ, FYIJ, FZIJ
real(kind=8)     :: MINRAD
real(kind=8)     :: MAXRAD
real(kind=8)     :: RATIO_RAD


! Indices which indicate current bucket
integer :: JCELL, JCELLO

! Indices for loops
integer :: I, J, ICELL, NABOR
!- MPI Variables
integer :: ERRCODE

!- Zeros necessary variables
FXI = 0.0
FYI = 0.0
FZI = 0.0

!- Zeros temp variable
FCM_FORCE_TEMP(:,:) = 0.0

do ICELL = FCM_LOC_BUCKET_START, FCM_LOC_BUCKET_STOP 
 
 I = FCM_BUCKET_HEAD(ICELL)
 
 do while(I>-1)
 
  XI = FCM_XP(I)
  YI = FCM_YP(I)
  ZI = FCM_ZP(I)
  
  FXI = FCM_FORCE_TEMP(I,1)
  FYI = FCM_FORCE_TEMP(I,2)
  FZI = FCM_FORCE_TEMP(I,3)
  
  J = FCM_BUCKET_PART_LIST(I)
  
  do while(J>-1)
  
   XIJ = XI - FCM_XP(J)
   YIJ = YI - FCM_YP(J)
   ZIJ = ZI - FCM_ZP(J)
     
   XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
   YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
   ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
   
   RIJSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ;
   
   TWOA = FCM_SPHERE_RADP(I) + FCM_SPHERE_RADP(J)
   
   RREFSQ = (FCM_FBRANGE * TWOA)**2
   
   if (RIJSQ.lt.RREFSQ) then    
    
    MINRAD = min(FCM_SPHERE_RADP(I),FCM_SPHERE_RADP(J))
    MAXRAD = max(FCM_SPHERE_RADP(I),FCM_SPHERE_RADP(J))
    RATIO_RAD = MINRAD/MAXRAD
    FOURASQ = TWOA**2

    if (RIJSQ.lt.(0.49*FOURASQ)) then
      print*,MYID,'RIJ/TWOA = ', dsqrt(RIJSQ)/TWOA, 'I,J = ', I,J
      call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR)
    end if

    TEMP = RREFSQ - FOURASQ
    TEMP2 = (RREFSQ - RIJSQ)/TEMP
    TEMP2 = TEMP2*TEMP2
    RIJ = dsqrt(RIJSQ)
    
    FXIJ = FCM_FBLEVEL*TEMP2*TEMP2*XIJ/(TWOA)*RATIO_RAD
    FYIJ = FCM_FBLEVEL*TEMP2*TEMP2*YIJ/(TWOA)*RATIO_RAD
    FZIJ = FCM_FBLEVEL*TEMP2*TEMP2*ZIJ/(TWOA)*RATIO_RAD
    
    FXI = FXI + FXIJ
    FYI = FYI + FYIJ
    FZI = FZI + FZIJ
    
    FCM_FORCE_TEMP(J,1) = FCM_FORCE_TEMP(J,1) - FXIJ
    FCM_FORCE_TEMP(J,2) = FCM_FORCE_TEMP(J,2) - FYIJ
    FCM_FORCE_TEMP(J,3) = FCM_FORCE_TEMP(J,3) - FZIJ

   
   end if
   
   J = FCM_BUCKET_PART_LIST(J)
   
  end do
  
  JCELLO = 13*(ICELL-1)
  
  do NABOR = 1, 13
  
   JCELL = FCM_BUCKET_MAPLIST(JCELLO + NABOR)
   J = FCM_BUCKET_HEAD(JCELL)
   
   do while(J>-1)
  
    XIJ = XI - FCM_XP(J)
    YIJ = YI - FCM_YP(J)
    ZIJ = ZI - FCM_ZP(J)
   
    XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
    YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
    ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
    
    RIJSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ;
!~     print*,'RIJSQ = ', RIJSQ
        
    TWOA = FCM_SPHERE_RADP(I) + FCM_SPHERE_RADP(J)
   
    RREFSQ = (FCM_FBRANGE * TWOA)**2
        
    if (RIJSQ.lt.RREFSQ) then
    
     MINRAD = min(FCM_SPHERE_RADP(I),FCM_SPHERE_RADP(J))
     MAXRAD = max(FCM_SPHERE_RADP(I),FCM_SPHERE_RADP(J))
     RATIO_RAD = MINRAD/MAXRAD

     FOURASQ = TWOA**2
     if (RIJSQ.lt.(0.49*FOURASQ)) then
      print*,MYID,'RIJ/TWOA = ', dsqrt(RIJSQ)/TWOA, 'I,J = ', I,J
      call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR)
     end if

     TEMP = RREFSQ - FOURASQ
     TEMP2 = (RREFSQ - RIJSQ)/TEMP
     TEMP2 = TEMP2*TEMP2
     RIJ = dsqrt(RIJSQ)
     
     FXIJ = FCM_FBLEVEL*TEMP2*TEMP2*XIJ/(TWOA)*RATIO_RAD
     FYIJ = FCM_FBLEVEL*TEMP2*TEMP2*YIJ/(TWOA)*RATIO_RAD
     FZIJ = FCM_FBLEVEL*TEMP2*TEMP2*ZIJ/(TWOA)*RATIO_RAD
     
     FXI = FXI + FXIJ
     FYI = FYI + FYIJ
     FZI = FZI + FZIJ
     
     FCM_FORCE_TEMP(J,1) = FCM_FORCE_TEMP(J,1) - FXIJ
     FCM_FORCE_TEMP(J,2) = FCM_FORCE_TEMP(J,2) - FYIJ
     FCM_FORCE_TEMP(J,3) = FCM_FORCE_TEMP(J,3) - FZIJ
   
    end if
    
!~     print*,'Neighbouring Bucket'
!~     print*, 'I = ', I
!~     print*, 'J = ', J
!~     print*, 'FCM_POS_REL_X(I,J) = ', FCM_POS_REL_X(I,J)
!~     print*, 'FCM_POS_REL_Y(I,J) = ', FCM_POS_REL_Y(I,J)
!~     print*, 'FCM_POS_REL_Z(I,J)  = ',FCM_POS_REL_Z(I,J) 
!~     read(*,*)
       
   
    J = FCM_BUCKET_PART_LIST(J)
   
   end do 
  
  end do
  
  FCM_FORCE_TEMP(I,1) = FXI
  FCM_FORCE_TEMP(I,2) = FYI
  FCM_FORCE_TEMP(I,3) = FZI 
  
  I = FCM_BUCKET_PART_LIST(I)
  
 end do
 
end do

! Simple addition of the contribution of each processor to the force barrier
call MPI_ALLREDUCE(FCM_FORCE_TEMP,FCM_FORCE,NPART_FULL*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)


end subroutine FCM_BARRIER_FORCE
