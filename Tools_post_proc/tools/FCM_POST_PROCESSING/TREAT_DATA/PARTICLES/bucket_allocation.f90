
!!====================================================================
!!
!! 
!!> @brief
!!> Allocates the arrays for bucket sorting algorithm
!!>Set the map of buckets, define their 13 upper-right neighbors and 
!!> assign the local number of buckets per proc
!! Date :  04/12/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine BUCKET_ALLOCATION(NPART)
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

integer, intent(in) :: NPART

!------------------------------------------------------------------


allocate(FCM_BUCKET_MAPLIST(FCM_MAPSIZE))
allocate(FCM_BUCKET_HEAD(FCM_BUCKET_NB_TOT))
allocate(FCM_BUCKET_PART_LIST(NPART))


end subroutine BUCKET_ALLOCATION
