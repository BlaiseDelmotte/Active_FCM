!!====================================================================
!!
!!
!!====================================================================

subroutine SAVE_SCALAR(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use SCALAR_VARIABLE
use PARAM_PHYS
use MPI_STRUCTURES
use RHS_VARIABLES

!!====================================================================
!!
!!
!!====================================================================

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Cycle number
integer, intent(in) :: NCYCLE

!- File name 
character(len=40) :: FILENAME

!- File descriptor
integer :: DESCRIPTEUR

!- Size of record
integer :: RECSIZE

!- Statistics of the droped field
real(kind=8), dimension(4) :: STATSCL

!- Index
integer :: I, J ,K
!---------------------------------------------------------------------




!!====================================================================
!! 1. Fluid velocity
!!====================================================================
if(ISAVEFLUID == 1) then


!!--------------------------------------------------------------------
!! 1.1. Save scalar field
!!--------------------------------------------------------------------

!- Define file name
if(NCYCLE>0) then
 write(FILENAME,10200)'scl',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10100)'scl',trim(FILE_EXT),'.end'
end if 


!- Open file containing the last fluid velocity field
open(unit = 150, file = trim(FILENAME),form='unformatted')

write(150)ISIZE(1),ISTART(1),IEND(1)
write(150)ISIZE(2),ISTART(2),IEND(2)
write(150)ISIZE(3),ISTART(3),IEND(3)

write(150)(((THETA(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(150)

if(MYID==0) write(*,*) 'Save fluid solution --> Multiple binary files'


!!--------------------------------------------------------------------
!! 1.2.  Save the right-hand-side of scalar field
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10200)'rhs_scl',trim(FILE_EXT),'.end'

!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')


write(120)TN, TNM1, TNM2
write(120)FSIZE(1),FSTART(1),FEND(1)
write(120)FSIZE(2),FSTART(2),FEND(2)
write(120)FSIZE(3),FSTART(3),FEND(3)

!- Drop rhs
write(120)(((RHS_SCL(I,J,K,TN  ),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_SCL(I,J,K,TNM1),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_SCL(I,J,K,TNM2),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))

!- Close file
close(120)


if(MYID==0) write(*,*) 'Save rhs of scalar --> Multiple Binary Files'



!!====================================================================
!! 2. MPI I/O 
!!====================================================================
elseif(ISAVEFLUID == 2) then

!!--------------------------------------------------------------------
!! 2.1. Save scalar field
!!--------------------------------------------------------------------
!!- Print file name
 if(NCYCLE>0) then
  write(FILENAME,10200)'scl',trim(FILE_EXT),'_t',NCYCLE,'.bin'
 else
  FILENAME='scl.end'
 end if 


 call SAVE_MPIIO(THETA,FILENAME)


 if(MYID==0) write(*,*) 'Save scalar solution --> MPI I/O'



!!--------------------------------------------------------------------
!! 2.2.  Save the right-hand-side of scalar field
!!--------------------------------------------------------------------
  FILENAME='rhs_scl.end'
  call SAVE_MPIIO_RHS(RHS_SCL,FILENAME)

  if(MYID==0) write(*,*) 'Save rhs of scalar --> MPI I/O'

end if


STATSCL(:) = ZERO

!!- 1st order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))   ),STATSCL(1))

!!- 2nd order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATSCL(2))

!!- 3rd order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATSCL(3))

!!- 4th order moment
call RSUMCPU(SUM(THETA(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATSCL(4))


if(MYID==0) write(*,*) '        <T> = ',STATSCL(1) / NGLOB
if(MYID==0) write(*,*) '      <T^2> = ',STATSCL(2) / NGLOB
if(MYID==0) write(*,*) '      <T^3> = ',STATSCL(3) / NGLOB
if(MYID==0) write(*,*) '      <T^4> = ',STATSCL(4) / NGLOB
if(MYID==0) write(*,*)

!I = (IEND(1)-ISTART(1))/2
!J = (IEND(2)-ISTART(2))/2
!K = (IEND(3)-ISTART(3))/2
!write(*,10500)'Id=',MYID,' T(',I,',',J,',',K,')=',THETA(I,J,K)



if(MYID==0) write(*,*) 'Final scalar field --> Saved'


!!====================================================================

10100 format (A,A,A)
10200 format (A,A,A,I8.8,A)
10500 format (A,I3,A,I3,A,I3,A,I3,A,E13.6)

end subroutine SAVE_SCALAR
