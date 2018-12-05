 !!====================================================================
!!
!! 
!!> @brief
!!> Compute repulsive force between walls and particles
!! Date :  28/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_REPULSIVE_WALL
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use FCM_BUCKET_VARIABLE



implicit none


! Temp variables
real(kind=8)     :: XI
! Indices for loops
integer :: I, ICELL


!- Zeros temp variable
FCM_FORCE_TEMP(:,:) = 0.0
FCM_FORCE_TEMP2(:,:) = 0.0

do ICELL = FCM_LOC_BUCKET_START, FCM_LOC_BUCKET_STOP 
 
 I = FCM_BUCKET_HEAD(ICELL)
 
 do while(I>-1)
 
  XI = FCM_XP(I)
  
  if ( XI<=FCM_WALL_RANGE ) then
   FCM_FORCE_TEMP(I,1) = FCM_FORCE_TEMP(I,1) &
                       + FCM_WALL_LEVEL*( FCM_WALL_RANGE - XI )
  end if
  
  if ( (LXMAX/2.0-XI)<=FCM_WALL_RANGE ) then
   FCM_FORCE_TEMP(I,1) = FCM_FORCE_TEMP(I,1) &
                       - FCM_WALL_LEVEL*( FCM_WALL_RANGE - (LXMAX/2.0-XI) ) 
  end if
  
  I = FCM_BUCKET_PART_LIST(I)
  
!  print*,'XI = ', XI
!  print*,'LXMAX/2.0-XI = ', LXMAX/2.0-XI
!  print*,'FCM_WALL_RANGE = ', FCM_WALL_RANGE
!  read(*,*)
  
 end do
 
end do


! Simple addition of the contribution of each processor to the velocity average
call MPI_ALLREDUCE(FCM_FORCE_TEMP,FCM_FORCE_TEMP2,NPART_FULL*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

FCM_FORCE = FCM_FORCE + FCM_FORCE_TEMP2

end subroutine FCM_REPULSIVE_WALL
