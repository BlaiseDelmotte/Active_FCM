
!!====================================================================
!!
!! 
!!> @brief
!!> Define and Assign the local number of buckets per proc
!! Date :  04/12/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine BUCKET_DISTRIBUTION(LXMAX, &
                               LYMAX, &
                               LZMAX, &
                               BOXSIZE, &
                               MYID, &
                               NPROC )
!!=====================================================================
!
!!---------------------------------------------------------------------
!! Warning: 
!!          
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
!!
!!=====================================================================


use BUCKET_VARIABLE


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Maximal radius  
real(kind=8), intent(in) :: LXMAX
real(kind=8), intent(in) :: LYMAX
real(kind=8), intent(in) :: LZMAX
real(kind=8), intent(in) :: BOXSIZE
integer, intent(in) :: MYID
integer, intent(in) :: NPROC

!- Index
integer :: IP
!------------------------------------------------------------------


!- Defines bucket size in terms of maximal radius
FCM_BUCKET_SIZE_X =  BOXSIZE
FCM_BUCKET_SIZE_Y =  FCM_BUCKET_SIZE_X 
FCM_BUCKET_SIZE_Z =  FCM_BUCKET_SIZE_X 


!~ !- Match bucket size to have an integer number of buckets
!~ if(FCM_BUCKET_SIZE.gt.LXMAX/2.0)then
!~  FCM_BUCKET_SIZE = LXMAX/2.0
!~ end if
 

if(mod(LXMAX,FCM_BUCKET_SIZE_X).ne.0)then
 FCM_BUCKET_SIZE_X = LXMAX/real(floor(LXMAX/FCM_BUCKET_SIZE_X))
end if

if(mod(LYMAX,FCM_BUCKET_SIZE_Y).ne.0)then
 FCM_BUCKET_SIZE_Y = LYMAX/real(floor(LYMAX/FCM_BUCKET_SIZE_Y))
end if

if(mod(LZMAX,FCM_BUCKET_SIZE_Z).ne.0)then
 FCM_BUCKET_SIZE_Z = LZMAX/real(floor(LZMAX/FCM_BUCKET_SIZE_Z))
end if




!- Number of bucket in each direction
FCM_BUCKET_NB_DIR_X = floor(LXMAX/FCM_BUCKET_SIZE_X)
FCM_BUCKET_NB_DIR_Y = floor(LYMAX/FCM_BUCKET_SIZE_Y)
FCM_BUCKET_NB_DIR_Z = floor(LZMAX/FCM_BUCKET_SIZE_Z)

!!!!! POTENTIAL OPTIMIZATION FOR X-WALL TO IMPROVE LOAD BALANCING
!!! NEED TO BE THOUGHT 
!if (FCM_BC==1) then
! !- if 3D periodic domain we mesh the full length in each directions
! !- with buckets
! FCM_BUCKET_NB_DIR_X = floor(LXMAX/FCM_BUCKET_SIZE)
!else if (FCM_BC==2) then
! !- if slip at x=LXMAX/2, only mesh that part 
! FCM_BUCKET_NB_DIR_X = floor(LXMAX/2.0/FCM_BUCKET_SIZE)
!end if
!!!!! POTENTIAL OPTIMIZATION FOR X-WALL TO IMPROVE LOAD BALANCING

!- Number of bucket in the whole box
FCM_BUCKET_NB_TOT = FCM_BUCKET_NB_DIR_X &
                   *FCM_BUCKET_NB_DIR_Y &
                   *FCM_BUCKET_NB_DIR_Z

!- Size of the bucket map: 13 neighbours per bucket
FCM_MAPSIZE = 13*FCM_BUCKET_NB_TOT


FCM_LOC_BUCKET_START = (MYID) * (floor( real(FCM_BUCKET_NB_TOT) / real(NPROC) ))+1;
FCM_LOC_BUCKET_STOP = (MYID+1) * (floor( real(FCM_BUCKET_NB_TOT) / real(NPROC) ));
    
if(MYID == NPROC-1)then
  FCM_LOC_BUCKET_STOP = FCM_BUCKET_NB_TOT
end if


if (FCM_BUCKET_NB_TOT.lt.27) then
 if (MYID == 0) then
  print*, '!!!!!!!!!!!!!!!!!!!!!!!! ERROR !!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  print*, '!Bucket number smaller than possible periodic images of 27!'
  print*, '!Chose a smaller bucket size                              !'
  print*,  '!!!!!!!!!!!!!! EXIT PROGRAM !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  stop
 end if
end if

 Print*,'FCM_BUCKET_SIZE_X  = ',FCM_BUCKET_SIZE_X
 Print*,'FCM_BUCKET_SIZE_Y  = ',FCM_BUCKET_SIZE_Y
 Print*,'FCM_BUCKET_SIZE_Z  = ',FCM_BUCKET_SIZE_Z
 Print*,'FCM_BUCKET_NB_DIR_X  = ',FCM_BUCKET_NB_DIR_X 
 Print*,'FCM_BUCKET_NB_DIR_Y  = ',FCM_BUCKET_NB_DIR_Y
 Print*,'FCM_BUCKET_NB_DIR_Z  = ',FCM_BUCKET_NB_DIR_Z
 Print*,'FCM_LOC_BUCKET_START  = ',FCM_LOC_BUCKET_START 
 Print*,'FCM_LOC_BUCKET_STOP  = ',FCM_LOC_BUCKET_STOP
 Print*,'mod(LXMAX,FCM_BUCKET_SIZE_X)  = ',mod(LXMAX,FCM_BUCKET_SIZE_X)
 Print*,'mod(LYMAX,FCM_BUCKET_SIZE_Y)  = ',mod(LYMAX,FCM_BUCKET_SIZE_Y)
 Print*,'mod(LZMAX,FCM_BUCKET_SIZE_Z)  = ',mod(LZMAX,FCM_BUCKET_SIZE_Z)
!~  read(*,*)

end subroutine BUCKET_DISTRIBUTION
