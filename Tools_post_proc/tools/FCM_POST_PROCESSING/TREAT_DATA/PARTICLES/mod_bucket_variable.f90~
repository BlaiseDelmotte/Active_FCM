!!=====================================================================
!!
!! FCM Buckets variables
!!
!! Description :
!!> @brief
!!> Module containing all variable useful for bucket sorting algorithm
!!
!!
!! Date :  21/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!
!!=====================================================================
module BUCKET_VARIABLE

implicit none

!!====================================================================
!! 1. BUCKET SORTING VARIABLES
!!====================================================================

integer :: FCM_BUCKET_NB_DIR_X      !< Bucket number along x-direction
integer :: FCM_BUCKET_NB_DIR_Y      !< Bucket number along y-direction
integer :: FCM_BUCKET_NB_DIR_Z      !< Bucket number along z-direction

integer :: FCM_BUCKET_NB_TOT        !< Total Bucket number

integer :: FCM_MAPSIZE              !< Size of the bucket map

integer :: FCM_LOC_BUCKET_START     !- Bucket Interval for each proc (local)
integer :: FCM_LOC_BUCKET_STOP      !- Bucket Interval for each proc (local)

real(kind=8) :: FCM_BUCKET_SIZE_X, FCM_BUCKET_SIZE_Y,  FCM_BUCKET_SIZE_Z    !< Bucket size

integer, dimension(:), allocatable :: FCM_BUCKET_MAPLIST      !< Map of neighboring buckets

integer, dimension(:), allocatable :: FCM_BUCKET_HEAD      !< Tells which is the first part in each bucket

integer, dimension(:), allocatable :: FCM_BUCKET_PART_LIST      !< Define the following part in the list


end module BUCKET_VARIABLE
