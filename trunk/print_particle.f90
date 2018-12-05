!!=====================================================================
!!
!!
!!  Print particle position and velocity
!!
!!
!!=====================================================================

subroutine PRINT_PARTICLE(TIME,NFILEOUT)

!!- Commenting -PF - 29/02/2012
!!subroutine PRINT_PARTICLE(TIME,NFILEOUT)

!!=====================================================================
!!
!!
!!=====================================================================

use PARTICLE_PARALLEL

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time
real(kind=8), intent(in) :: TIME

!!- Commenting -PF - 29/02/2012
!!integer, intent(in) :: NFILEOUT

!- File name
character(len=30) :: FILENAME


!- Index
integer :: I, J, K
!---------------------------------------------------------------------



do J = 1, NIG


!-Print filename
write(FILENAME,10101)'p',J,'_t',NFILEOUT,'.p',MYID,'.dat'


!- ASCII
open(unit=300,file=trim(FILENAME))

!- Ecriture ASCII
write(300,2000)
write(300,2001)NPART_LOC(J),1,1,TIME

do I = 1,NPART_LOC(J)
 write(300,'(8(e17.7))')PART(I,J)%XP, PART(I,J)%YP, PART(I,J)%ZP, &
                        PART(I,J)%UP, PART(I,J)%VP, PART(I,J)%WP !, &
!                        PART(I,J)%UFAP, PART(I,J)%VFAP, PART(I,J)%WFAP, &
!                        REAL(PART(I,J)%ID), REAL(PART(I,J)%PROC_ID)                                   
end do


!- close file
close(300)

end do


!===========================================================================
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2001 format ('ZONE F=POINT I=',i6,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

10101 format (A,I2.2,A,I2.2,A,I2.2,A)
10205 format(A,I4.4)

end subroutine PRINT_PARTICLE
