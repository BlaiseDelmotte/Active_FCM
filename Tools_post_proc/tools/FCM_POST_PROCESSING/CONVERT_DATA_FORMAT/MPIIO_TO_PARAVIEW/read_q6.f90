!!====================================================================
!!
!!
!!====================================================================


subroutine READ_AND_WRITE_Q6(NPART, &
                             FILENAME, &
                             Q)

!!====================================================================
!!
!!
!!====================================================================

implicit none


!--------------------------------------------------------------------
! ARGUMENTS AND ARRAY
!--------------------------------------------------------------------

integer, intent(in) :: NPART
character(len=40) :: FILENAME
real(kind=8), dimension(NPART) :: Q
integer :: i


open (unit = 302, file = trim(FILENAME),status = 'old', action = 'read')

do i = 1, NPART

   read(302,*) Q(i)
 
end do
return
end subroutine READ_AND_WRITE_Q6
