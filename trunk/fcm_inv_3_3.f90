 !!====================================================================
!!
!! 
!!> @brief
!!> Compute inverse of 3*3 matrix
!! Date :  2014/02/07
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INV_3_3(MAT_H, MAT_HINV)

!!====================================================================
!! 
!!====================================================================
!! Forcing: 
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use FCM_BUCKET_VARIABLE



implicit none


real(kind=8), intent(in), dimension(3,3) :: MAT_H
real(kind=8), intent(out), dimension(3,3) :: MAT_HINV

MAT_HINV(1, 1) = (MAT_H(2, 2)*MAT_H(3, 3)-MAT_H(2, 3)*MAT_H(3, 2)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(1, 2) = -(MAT_H(1, 2)*MAT_H(3, 3)-MAT_H(1, 3)*MAT_H(3, 2)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(1, 3) = (MAT_H(1, 2)*MAT_H(2, 3)-MAT_H(1, 3)*MAT_H(2, 2)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(2, 1) = -(MAT_H(2, 1)*MAT_H(3, 3)-MAT_H(2, 3)*MAT_H(3, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(2, 2) = (MAT_H(1, 1)*MAT_H(3, 3)-MAT_H(1, 3)*MAT_H(3, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(2, 3) = -(MAT_H(1, 1)*MAT_H(2, 3)-MAT_H(1, 3)*MAT_H(2, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(3, 1) = (MAT_H(2, 1)*MAT_H(3, 2)-MAT_H(2, 2)*MAT_H(3, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(3, 2) = -(MAT_H(1, 1)*MAT_H(3, 2)-MAT_H(1, 2)*MAT_H(3, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))
MAT_HINV(3, 3) = (MAT_H(1, 1)*MAT_H(2, 2)-MAT_H(1, 2)*MAT_H(2, 1)) &
               / (MAT_H(1, 1)*MAT_H(2, 2)*MAT_H(3, 3)&
                 -MAT_H(1, 1)*MAT_H(2, 3)*MAT_H(3, 2)&
                 -MAT_H(1, 2)*MAT_H(2, 1)*MAT_H(3, 3)&
                 +MAT_H(1, 2)*MAT_H(2, 3)*MAT_H(3, 1)&
                 +MAT_H(1, 3)*MAT_H(2, 1)*MAT_H(3, 2)&
                 -MAT_H(1, 3)*MAT_H(2, 2)*MAT_H(3, 1))

end subroutine FCM_INV_3_3
