subroutine GET_DIMS_NEW(MYID,NX, NY, NZ, IPROC, JPROC, IS,IE,ISZE)

implicit none

integer, intent(in) :: MYID
integer, intent(in) :: NX, NY, NZ
integer, intent(in) :: IPROC, JPROC

integer, dimension(3), intent(out) :: IS
integer, dimension(3), intent(out) :: IE
integer, dimension(3), intent(out) :: ISZE

!!------------------------------------------------------------------------------

integer :: DXCPU, DYCPU
!!==============================================================================


DXCPU = NX/IPROC
DYCPU = NY/JPROC





IS(1) = 1
IE(1) = NX


IS(3) = int(MYID/IPROC)*DXCPU+1
IE(3) = IS(3) + DXCPU -1

IS(2) = 1 + MYID*DXCPU - IPROC*(IS(3)-1)
IE(2) = IS(2) + DXCPU - 1


ISZE(1) = NX
ISZE(2) = IE(2) - IS(2) + 1
ISZE(3) = IE(3) - IS(3) + 1 




end subroutine GET_DIMS_NEW
