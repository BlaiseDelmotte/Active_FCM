
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

subroutine FCM_BUCKET_ALLOCATIONS
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
use FCM_BUCKET_VARIABLE
use DNS_DIM

implicit none


!------------------------------------------------------------------


allocate(FCM_BUCKET_MAPLIST(FCM_MAPSIZE))
allocate(FCM_BUCKET_HEAD(FCM_BUCKET_NB_TOT))
allocate(FCM_BUCKET_PART_LIST(NPART_FULL))


end subroutine FCM_BUCKET_ALLOCATIONS
