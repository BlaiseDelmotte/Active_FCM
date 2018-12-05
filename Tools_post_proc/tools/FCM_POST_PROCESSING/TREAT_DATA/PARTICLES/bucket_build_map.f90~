
!!====================================================================
!!
!! 
!!> @brief
!!> Define the map of neighbouring buckets
!! Date :  04/12/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine BUCKET_BUILD_MAP
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

!- Temporary variables  
integer :: FCM_IMAP
integer :: FCM_TEMPMAP

!- Index
integer :: IX, IY, IZ
!------------------------------------------------------------------


FCM_IMAP = 0
FCM_TEMPMAP = 0


do IZ = 1, FCM_BUCKET_NB_DIR_Z
 do IY = 1, FCM_BUCKET_NB_DIR_Y
  do IX = 1, FCM_BUCKET_NB_DIR_X
   
   !- Gives the number of current bucket
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX, IY, IZ)
   
   !- Index in the array of the map neighbours
   FCM_IMAP = (FCM_TEMPMAP-1)*13
  
   !-Store the 13 right-upper neighbours
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY, IZ)
   FCM_BUCKET_MAPLIST(FCM_IMAP+1) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY+1, IZ)
   FCM_BUCKET_MAPLIST(FCM_IMAP+2) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX, IY+1, IZ)
   FCM_BUCKET_MAPLIST(FCM_IMAP+3) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX-1, IY+1, IZ)
   FCM_BUCKET_MAPLIST(FCM_IMAP+4) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY, IZ-1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+5) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY+1, IZ-1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+6) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX, IY+1, IZ-1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+7) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX-1, IY+1, IZ-1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+8) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY, IZ+1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+9) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX+1, IY+1, IZ+1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+10) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX, IY+1, IZ+1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+11) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX-1, IY+1, IZ+1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+12) = FCM_TEMPMAP
   
   call BUCKET_INUMBER(FCM_TEMPMAP, &
                           FCM_BUCKET_NB_DIR_X, &
                           FCM_BUCKET_NB_DIR_Y, &
                           FCM_BUCKET_NB_DIR_Z, &
                           IX, IY, IZ+1)
   FCM_BUCKET_MAPLIST(FCM_IMAP+13) = FCM_TEMPMAP
   

   
  end do
 end do
end do

end subroutine BUCKET_BUILD_MAP
