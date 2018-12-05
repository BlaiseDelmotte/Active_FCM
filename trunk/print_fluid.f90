!!=====================================================================
!!
!!   Print_fluid
!!
!!=====================================================================

subroutine PRINT_FLUID(TIME)

!!=====================================================================
!!
!!
!!=====================================================================

use dns_dim
use fluid_variable
use scalar_variable
use geometric_variable
use param_phys        
use ensight_var        

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time
real(kind=8), intent(in) :: TIME

!- File name
character(len=30) :: FILENAME

integer :: TOTO

!- Index
integer :: I, J, K
!---------------------------------------------------------------------


!-Print filename
write(FILENAME,10101)'uf_t',NFILEOUT,trim(FILE_EXT),'.dat'


!- ASCII
open(unit=300,file=trim(FILENAME))

!- Ecriture ASCII
if(SOLVE_SCALAR) then
 write(300,2001)
else
 write(300,2000)
end if
write(300,2010)ISIZE(1),ISIZE(2),ISIZE(3),TIME

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
if(SOLVE_SCALAR) then
   write(300,10000)XMESH(I), YMESH(J), ZMESH(K), &
                            UFLU(I,J,K)-UREF(1), &
                            VFLU(I,J,K)-UREF(2), &
                            WFLU(I,J,K)-UREF(3), &
                           THETA(I,J,K)
else
   write(300,10000)XMESH(I), YMESH(J), ZMESH(K), &
                            UFLU(I,J,K)-UREF(1), &
                            VFLU(I,J,K)-UREF(2), &
                            WFLU(I,J,K)-UREF(3)
end if
  end do
 end do
end do


!- close file
close(300)


NFILEOUT = NFILEOUT + 1

TOTO = 10

!===========================================================================
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2001 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "T"')
2010 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

10000 format (10(e17.7))
10101 format (A,I2.2,A,A)

end subroutine PRINT_FLUID
