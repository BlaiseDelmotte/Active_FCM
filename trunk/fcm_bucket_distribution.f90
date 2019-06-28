
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

subroutine FCM_BUCKET_DISTRIBUTION
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

use GEOMETRIC_VARIABLE
use FCM_FORCING_VARIABLE
use FCM_PART_VARIABLE
use FCM_BUCKET_VARIABLE
use DNS_DIM

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Maximal radius  
real(kind=8) :: FCM_MAXRAD


!- Index
integer :: IP
!------------------------------------------------------------------


!- Defines bucket size in terms of maximal radius
if(FCM_NELLIPSOID == 0)then
 FCM_MAXRAD = maxval(FCM_SPHERE_RADP)
else
 FCM_MAXRAD = maxval(FCM_ELLIPSOID_RADP)
end if

FCM_BUCKET_SIZE_X =  2.0*FCM_FBRANGE*FCM_MAXRAD
FCM_BUCKET_SIZE_Y =  FCM_BUCKET_SIZE_X 
FCM_BUCKET_SIZE_Z =  FCM_BUCKET_SIZE_X 

! If we have a boundary at LX/2 then check we don't overlap the boundary
if (FCM_BC.eq.2) then
 !- Match bucket size to have an integer number of buckets
 if(FCM_BUCKET_SIZE_X.gt.LXMAX/2.0)then
   FCM_BUCKET_SIZE_X = LXMAX/2.0
 end if
 if(mod(LXMAX/2.0,FCM_BUCKET_SIZE_X).ne.0)then
  FCM_BUCKET_SIZE_X = LXMAX/2.0/real(floor(LXMAX/2.0/FCM_BUCKET_SIZE_X))
 end if
 FCM_BUCKET_NB_DIR_X = floor(LXMAX/2.0/FCM_BUCKET_SIZE_X)
else ! If periodic boundary at LX
 if(mod(LXMAX,FCM_BUCKET_SIZE_X).ne.0)then
  FCM_BUCKET_SIZE_X = LXMAX/real(floor(LXMAX/FCM_BUCKET_SIZE_X))
 end if
 !- Number of bucket along direction
 FCM_BUCKET_NB_DIR_X = floor(LXMAX/FCM_BUCKET_SIZE_X)
end if 

if (MYID.eq.0) then
 print*, 'FCM_BUCKET_SIZE_X/(LXMAX/2.0), FCM_BUCKET_NB_DIR_X = ',FCM_BUCKET_SIZE_X/(LXMAX/2.0), FCM_BUCKET_NB_DIR_X
end if


if(mod(LYMAX,FCM_BUCKET_SIZE_Y).ne.0)then
 FCM_BUCKET_SIZE_Y = LYMAX/real(floor(LYMAX/FCM_BUCKET_SIZE_Y))
end if
FCM_BUCKET_NB_DIR_Y = floor(LYMAX/FCM_BUCKET_SIZE_Y)

if(mod(LZMAX,FCM_BUCKET_SIZE_Z).ne.0)then
 FCM_BUCKET_SIZE_Z = LZMAX/real(floor(LZMAX/FCM_BUCKET_SIZE_Z))
end if
FCM_BUCKET_NB_DIR_Z = floor(LZMAX/FCM_BUCKET_SIZE_Z)


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

!~  Print*,'FCM_MAXRAD = ',FCM_MAXRAD
!~  Print*,'FCM_BUCKET_SIZE_X  = ',FCM_BUCKET_SIZE_X
!~  Print*,'FCM_BUCKET_SIZE_X/FCM_MAXRAD = ',FCM_BUCKET_SIZE_X / FCM_MAXRAD 
!~  Print*,'FCM_BUCKET_SIZE_Y  = ',FCM_BUCKET_SIZE_Y
!~  Print*,'FCM_BUCKET_SIZE_Y/FCM_MAXRAD = ',FCM_BUCKET_SIZE_Y / FCM_MAXRAD
!~  Print*,'FCM_BUCKET_SIZE_Z  = ',FCM_BUCKET_SIZE_Z
!~  Print*,'FCM_BUCKET_SIZE_Z/FCM_MAXRAD = ',FCM_BUCKET_SIZE_Z / FCM_MAXRAD
!~  Print*,'FCM_BUCKET_NB_DIR_X  = ',FCM_BUCKET_NB_DIR_X 
!~  Print*,'FCM_BUCKET_NB_DIR_Y  = ',FCM_BUCKET_NB_DIR_Y
!~  Print*,'FCM_BUCKET_NB_DIR_Z  = ',FCM_BUCKET_NB_DIR_Z
!~  Print*,'FCM_LOC_BUCKET_START  = ',FCM_LOC_BUCKET_START 
!~  Print*,'FCM_LOC_BUCKET_STOP  = ',FCM_LOC_BUCKET_STOP
!~  Print*,'mod(LXMAX,FCM_BUCKET_SIZE_X)  = ',mod(LXMAX,FCM_BUCKET_SIZE_X)
!~  Print*,'mod(LYMAX,FCM_BUCKET_SIZE_Y)  = ',mod(LYMAX,FCM_BUCKET_SIZE_Y)
!~  Print*,'mod(LZMAX,FCM_BUCKET_SIZE_Z)  = ',mod(LZMAX,FCM_BUCKET_SIZE_Z)
!~  read(*,*)

end subroutine FCM_BUCKET_DISTRIBUTION
