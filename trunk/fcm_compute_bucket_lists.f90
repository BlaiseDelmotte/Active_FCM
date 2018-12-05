
!!====================================================================
!!
!! 
!!> @brief
!!> Establish linked-list for nf and collsion interactions
!! Date :  04/12/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_COMPUTE_BUCKET_LISTS
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
use FCM_PART_VARIABLE
use FCM_BUCKET_VARIABLE
use DNS_DIM

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

real(kind=8) :: TEMP_M

!- Index which tells which bucket the particle belongs to
integer :: FCM_IND

!- Index
integer :: IP, IB
!------------------------------------------------------------------

do IB = 1, FCM_BUCKET_NB_TOT
 FCM_BUCKET_HEAD(IB) = -1
end do


!TEMP_M = real(FCM_BUCKET_NB_DIR)



do IP = 1, NPART_FULL
 FCM_IND = (floor( FCM_XP(IP)/FCM_BUCKET_SIZE_X)) &
         + (floor( FCM_YP(IP)/FCM_BUCKET_SIZE_Y )) * FCM_BUCKET_NB_DIR_X &
         + (floor( FCM_ZP(IP)/FCM_BUCKET_SIZE_Z )) * FCM_BUCKET_NB_DIR_X * FCM_BUCKET_NB_DIR_Y +1 

 FCM_BUCKET_PART_LIST(IP) = FCM_BUCKET_HEAD(FCM_IND)
 FCM_BUCKET_HEAD(FCM_IND) = IP
 
!~  if (MYID==0) then
!~   Print*,'IP = ',IP
!~   Print*,'FCM_IND = ',FCM_IND
!~   Print*, 'FCM_BUCKET_PART_LIST(IP) = ', FCM_BUCKET_PART_LIST(IP)
!~  end if
 
end do




end subroutine FCM_COMPUTE_BUCKET_LISTS
