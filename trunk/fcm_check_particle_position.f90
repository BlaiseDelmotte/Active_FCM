
!!====================================================================
!!
!! 
!!> @brief
!!> Check poistions and velocities have physical values
!!
!! Date :  10/05/2019
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_CHECK_PARTICLE_POSITION

use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use DNS_DIM
use PARAM_PHYS, only: DTIME 

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time control variable
real(kind=8) :: TIME_START, TIME_END
!- Index
integer :: I
!- MPI Variables
integer :: ERRCODE
!------------------------------------------------------------------



!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if

 do I = 1, NPART_FULL
  if ( (isnan(FCM_XP(I))).or.(isnan(FCM_YP(I))).or.(isnan(FCM_ZP(I))) ) then
   if (MYID==0) then
    print*, 'Nan position:'
    print*, I, FCM_XP(I), FCM_YP(I), FCM_ZP(I)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
  if ( (FCM_XP(I).gt.LXMAX).or.(FCM_XP(I).lt.0.0).or.(FCM_YP(I).gt.LYMAX).or.(FCM_YP(I).lt.0.0).or.(FCM_ZP(I).gt.LZMAX).or.(FCM_ZP(I).lt.0.0) ) then
   if (MYID==0) then
    print*, 'Particle crossed more than one computational domain in one time iteration:'
    print*, I, FCM_XP(I), FCM_YP(I), FCM_ZP(I)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
  if ( ((FCM_XP(I)+1).eq.(FCM_XP(I))).or.((FCM_YP(I)+1).eq.(FCM_YP(I))).or.((FCM_ZP(I)+1).eq.(FCM_ZP(I))) ) then
   if (MYID==0) then
    print*, 'Infinity position:'
    print*, I, FCM_XP(I), FCM_YP(I), FCM_ZP(I)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
  if ( (isnan(FCM_UP(I,1))).or.(isnan(FCM_VP(I,1))).or.(isnan(FCM_WP(I,1))) ) then
   if (MYID==0) then
    print*, 'Nan velocity:'
    print*, I, FCM_UP(I,1), FCM_VP(I,1), FCM_WP(I,1)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
  if ( ((FCM_UP(I,1)+1).eq.(FCM_UP(I,1))).or.((FCM_VP(I,1)+1).eq.(FCM_VP(I,1))).or.((FCM_WP(I,1)+1).eq.(FCM_WP(I,1))) ) then
   if (MYID==0) then
    print*, 'Infinity velocity:'
    print*, I, FCM_UP(I,1), FCM_VP(I,1), FCM_WP(I,1)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
  if ( (FCM_BC.eq.2).and.(FCM_XP(I).gt.(LXMAX/2.0)) ) then
   if (MYID==0) then
    print*, 'Particle went through boundary:'
    print*, I, 0, LXMAX/2.0, FCM_XP(I)
    print*, 'Stop simulation:'
   end if
   call MPI_Abort(MPI_COMM_WORLD, ERRCODE, IERR) 
  end if
 end do


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
end if




end subroutine FCM_CHECK_PARTICLE_POSITION
