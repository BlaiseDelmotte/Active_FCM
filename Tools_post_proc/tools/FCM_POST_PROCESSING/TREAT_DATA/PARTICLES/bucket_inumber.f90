
!!====================================================================
!!
!! 
!!> @brief
!!>  Assign the bucket number in the bucket list
!! Date :  04/12/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine BUCKET_INUMBER(INUMBER, M_X, M_Y, M_Z, IX, IY, IZ)
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


implicit none

!---------------------------------------------------------------------
! ARGUMENTS
!---------------------------------------------------------------------
integer, intent(out) :: INUMBER
integer, intent(in) :: M_X
integer, intent(in) :: M_Y
integer, intent(in) :: M_Z
integer, intent(in) :: IX
integer, intent(in) :: IY
integer, intent(in) :: IZ


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Temporary variables  
real(kind=8) :: TEMP_INUMBER
real(kind=8) :: TEMP_M_X
real(kind=8) :: TEMP_M_Y
!------------------------------------------------------------------

TEMP_M_X = real(M_X)
TEMP_M_Y = real(M_Y)

TEMP_INUMBER = real(mod((IX-1+M_X),M_X)) &
             + real(mod((IY-1+M_Y),M_Y))*TEMP_M_X & 
             + real(mod((IZ-1+M_Z),M_Z))*TEMP_M_X*TEMP_M_Y &
             +1.0

INUMBER = floor(TEMP_INUMBER)


end subroutine BUCKET_INUMBER
