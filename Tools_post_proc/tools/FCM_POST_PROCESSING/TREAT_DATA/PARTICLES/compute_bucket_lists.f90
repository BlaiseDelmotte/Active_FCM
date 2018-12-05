
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

subroutine COMPUTE_BUCKET_LISTS(NPART, &
                                POSI )
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
real(kind=8), dimension(NPART,3), intent(in) :: POSI
integer, intent(in) :: NPART


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



do IP = 1, NPART 
 FCM_IND = (floor( POSI(IP,1)/FCM_BUCKET_SIZE_X)) &
         + (floor( POSI(IP,2)/FCM_BUCKET_SIZE_Y )) * FCM_BUCKET_NB_DIR_X &
         + (floor( POSI(IP,3)/FCM_BUCKET_SIZE_Z )) * FCM_BUCKET_NB_DIR_X * FCM_BUCKET_NB_DIR_Y +1 

 FCM_BUCKET_PART_LIST(IP) = FCM_BUCKET_HEAD(FCM_IND)
 FCM_BUCKET_HEAD(FCM_IND) = IP
 
!~  if (MYID==0) then
!~   Print*,'IP = ',IP
!~   Print*,'FCM_IND = ',FCM_IND
!~   Print*, 'FCM_BUCKET_PART_LIST(IP) = ', FCM_BUCKET_PART_LIST(IP)
!~  end if
 
end do




end subroutine COMPUTE_BUCKET_LISTS
