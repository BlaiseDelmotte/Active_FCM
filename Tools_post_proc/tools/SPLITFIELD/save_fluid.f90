!!====================================================================
!!
!!                   Save solution
!!
!!====================================================================

subroutine SAVE_FLUID

!!====================================================================
!!
!! ISAVEM = 0: Direct access file
!!
!!        = 1: Using of MPI_FILE_WRITE_ORDERED 
!!             This method do not allow a restart with a different 
!!             number of CPU
!!
!!        = 2: Using of MPI_FILE_WRITE_AT_ALL
!!             Build a binary allowing a restart with a different
!!             number of CPU 
!!
!!====================================================================

use dns_dim
use fluid_variable
use param_phys
use forcing
use MPI_structures

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- File name 
character(len=40) :: FILENAME

!- File descriptor
integer :: DESCRIPTEUR

!- Record size
integer :: RECSIZE


integer(kind=MPI_OFFSET_KIND) :: POS_FILE
integer                       :: SIZE_REAL,SIZE_INT

!- Index
integer :: I, J ,K
!---------------------------------------------------------------------



!!====================================================================
!! 1. Multiple binary file
!!====================================================================
if(ISAVEM ==1) then


!- Define file name
write(FILENAME,10101)'FLUID',trim(FILE_EXT),'.end'

!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')

write(120)ISIZE(1),ISTART(1),IEND(1)
write(120)ISIZE(2),ISTART(2),IEND(2)
write(120)ISIZE(3),ISTART(3),IEND(3)

write(120)(((UFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
write(120)(((VFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
write(120)(((WFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)

if(MYID==0) write(*,*) 'Save fluid solution --> OK'
if(MYID==0) write(*,*) '                    --> Multiple Binary Files'



!!====================================================================
!! 2. Direct access method
!!====================================================================
elseif(ISAVEM == 2) then

 RECSIZE = 8*ISIZE(1)*ISIZE(2)*ISIZE(3)

!!--------------------------------------------------------------------
!! 2.1 x-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='uf.end'

 open(unit=120,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 write(unit=120,rec=MYID+1)UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(120)

!!--------------------------------------------------------------------
!! 2.2 y-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='vf.end'

 open(unit=121,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 write(unit=121,rec=MYID+1)VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(121)

!!--------------------------------------------------------------------
!! 2.3 z-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='wf.end'

 open(unit=122,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 write(unit=122,rec=MYID+1)WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(122)


if(MYID==0) write(*,*) 'Save fluid solution --> OK'
if(MYID==0) write(*,*) '                    --> Direct Access File'






!!====================================================================
!! 3. MPI I/O MPI_FILE_WRITE_ORDERED
!!====================================================================
!!
!!- The following is usefull for collective print in a single file.
!!  However it can not allow to restart with a different CPU
!!  decomposition.
!!--------------------------------------------------------------------
elseif(ISAVEM == 3) then

!!--------------------------------------------------------------------
!! 3.1 x-component
!!--------------------------------------------------------------------
!!- Print file name
FILENAME='uf.end'


call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                      trim(FILENAME), &
   MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, &
                         DESCRIPTEUR, &
                                IERR  )


RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)

!call MPI_FILE_WRITE_ORDERED(DESCRIPTEUR, & 
call MPI_FILE_WRITE_ALL(DESCRIPTEUR, & 
UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                 RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

!call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)


!!--------------------------------------------------------------------
!! 3.2 y-component
!!--------------------------------------------------------------------
!!- Print file name
FILENAME='vf.end'


call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                      trim(FILENAME), &
   MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, &
                         DESCRIPTEUR, &
                                IERR  )


RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)

!call MPI_FILE_WRITE_ORDERED(DESCRIPTEUR, & 
call MPI_FILE_WRITE_ALL(DESCRIPTEUR, & 
VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                 RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

!call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)


!!--------------------------------------------------------------------
!! 3.3 x-component
!!--------------------------------------------------------------------
!!- Print file name
FILENAME='wf.end'


call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                      trim(FILENAME), &
   MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, &
                         DESCRIPTEUR, &
                                IERR  )


RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)


!call MPI_FILE_WRITE_ORDERED(DESCRIPTEUR, & 
call MPI_FILE_WRITE_ALL(DESCRIPTEUR, & 
WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                 RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

!call MPI_BARRIER(MPI_COMM_WORLD,IERR)
call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)






if(MYID==0) write(*,*) 'Save fluid solution --> OK'
if(MYID==0) write(*,*) '                    --> MPI_FILE_WRITE_ORDERED'


elseif(ISAVEM ==4) then

!!- Print file name
FILENAME = 'uf.end'
call SAVE_MPIIO(UFLU,FILENAME)
FILENAME = 'vf.end'
call SAVE_MPIIO(VFLU,FILENAME)
FILENAME = 'wf.end'
call SAVE_MPIIO(WFLU,FILENAME)

if(MYID==0) write(*,*) 'Save fluid solution --> OK'
if(MYID==0) write(*,*) '                    --> MPI_FILE_WRITE_ALL'

end if





!!====================================================================
!! 2.Forcing
!!====================================================================
if(STEADY.and.MYID==0) then

 !- Define file name
 FILENAME='FORCING.end'

 !- Open file containing the forcing
 open(unit = 120, file = trim(FILENAME),form='unformatted')

 !- Number of forced waves
 write(120) NFORCE_FULL

 !- Random seed
 write(120) IDFORCE

 !- Forcing coefficients
 write(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
 write(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
 write(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)

 !- Index of hermit coefficient
! write(120)(NHERM(I),I=1,NFORCE_FULL)


 !- Close file
 close(120)

 write(*,*) 'Final forcing coefficients --> Saved'

end if

!!====================================================================
10101 format (A,A,A)

end subroutine SAVE_FLUID
