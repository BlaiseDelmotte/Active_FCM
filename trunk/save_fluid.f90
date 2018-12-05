!!====================================================================
!!
!!                   Save solution
!!
!!====================================================================

subroutine SAVE_FLUID(NCYCLE)

!!====================================================================
!!
!! ISAVEFLUID = 1 : Multiple binary files
!!            = 2 : MPI I/O
!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use PARAM_PHYS
use MPI_STRUCTURES
use RHS_VARIABLES


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

!- Record size
integer :: RECSIZE


integer(kind=MPI_OFFSET_KIND) :: POS_FILE
integer                       :: SIZE_REAL,SIZE_INT

!- Statistics of the droped field
real(kind=8), dimension(4) :: STATX, STATY, STATZ

!- Index
integer :: I, J ,K
!---------------------------------------------------------------------

!!====================================================================
!! 1. Multiple binary file
!!====================================================================
if(ISAVEFLUID ==1) then

!!--------------------------------------------------------------------
!! 1.1. Save the fluid velocity
!!--------------------------------------------------------------------

!!--------------------------------------------------------------------
!! 1.1.1. Uf
!!--------------------------------------------------------------------
!- Define file name
if(NCYCLE>0) then
 write(FILENAME,10200)'uf',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10100)'uf',trim(FILE_EXT),'.end'
end if 




!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')

write(120)ISIZE(1),ISTART(1),IEND(1)
write(120)ISIZE(2),ISTART(2),IEND(2)
write(120)ISIZE(3),ISTART(3),IEND(3)

write(120)(((UFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)


!!--------------------------------------------------------------------
!! 1.1.2. Vf
!!--------------------------------------------------------------------
!- Define file name
if(NCYCLE>0) then
 write(FILENAME,10200)'vf',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10100)'vf',trim(FILE_EXT),'.end'
end if 


!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')

write(120)ISIZE(1),ISTART(1),IEND(1)
write(120)ISIZE(2),ISTART(2),IEND(2)
write(120)ISIZE(3),ISTART(3),IEND(3)

write(120)(((VFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)

!!--------------------------------------------------------------------
!! 1.1.3. Wf
!!--------------------------------------------------------------------
!- Define file name
if(NCYCLE>0) then
 write(FILENAME,10200)'wf',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10100)'wf',trim(FILE_EXT),'.end'
end if 


!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')

write(120)ISIZE(1),ISTART(1),IEND(1)
write(120)ISIZE(2),ISTART(2),IEND(2)
write(120)ISIZE(3),ISTART(3),IEND(3)

write(120)(((WFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

!- Close file
close(120)


if(MYID==0) write(*,*) 'Save fluid solution --> Multiple Binary Files'



!!--------------------------------------------------------------------
!! 1.2. Save the right-hand-side of fluid momentum
!!--------------------------------------------------------------------
if(NCYCLE<0.and.SOLVE_FLUID==1) then


!!--------------------------------------------------------------------
!! 1.2.1. Uf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10200)'rhs_uf',trim(FILE_EXT),'.end'

!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')


write(120)TN, TNM1, TNM2
write(120)FSIZE(1),FSTART(1),FEND(1)
write(120)FSIZE(2),FSTART(2),FEND(2)
write(120)FSIZE(3),FSTART(3),FEND(3)

!- Drop rhs
write(120)(((RHS_UFOU(I,J,K,TN  ),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_UFOU(I,J,K,TNM1),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_UFOU(I,J,K,TNM2),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))

!- Close file
close(120)



!!--------------------------------------------------------------------
!! 1.2.2. Vf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10200)'rhs_vf',trim(FILE_EXT),'.end'

!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')


write(120)TN, TNM1, TNM2
write(120)FSIZE(1),FSTART(1),FEND(1)
write(120)FSIZE(2),FSTART(2),FEND(2)
write(120)FSIZE(3),FSTART(3),FEND(3)

!- Drop rhs
write(120)(((RHS_VFOU(I,J,K,TN  ),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_VFOU(I,J,K,TNM1),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_VFOU(I,J,K,TNM2),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))

!- Close file
close(120)



!!--------------------------------------------------------------------
!! 1.2.3. Wf
!!--------------------------------------------------------------------
!- Define file name
write(FILENAME,10200)'rhs_wf',trim(FILE_EXT),'.end'

!- Open file containing the last fluid velocity field
open(unit = 120, file = trim(FILENAME),form='unformatted')


write(120)TN, TNM1, TNM2
write(120)FSIZE(1),FSTART(1),FEND(1)
write(120)FSIZE(2),FSTART(2),FEND(2)
write(120)FSIZE(3),FSTART(3),FEND(3)

!- Drop rhs
write(120)(((RHS_WFOU(I,J,K,TN  ),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_WFOU(I,J,K,TNM1),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))
write(120)(((RHS_WFOU(I,J,K,TNM2),I=FSTART(1),FEND(1)),J=FSTART(2),FEND(2)),K=FSTART(3),FEND(3))

!- Close file
close(120)



if(MYID==0) write(*,*) 'Save rhs of fluid momentum --> Multiple Binary Files'
if(MYID==0) write(*,*) '  + TN   = ',TN
if(MYID==0) write(*,*) '  + TNM1 = ',TNM1
if(MYID==0) write(*,*) '  + TNM2 = ',TNM2

end if




!!====================================================================
!! 2. MPI I/O 
!!====================================================================
elseif(ISAVEFLUID == 2) then

!!--------------------------------------------------------------------
!! 2.1. Save the fluid velocity
!!--------------------------------------------------------------------


!!--------------------------------------------------------------------
!! 2.1.1. Uf
!!--------------------------------------------------------------------
!!- Print file name
if(NCYCLE>0) then
 write(FILENAME,10201)'uf_t',NCYCLE,'.bin'
else
 FILENAME='uf.end'
end if 


call SAVE_MPIIO(UFLU,FILENAME)



!!--------------------------------------------------------------------
!! 2.1.2. Vf
!!--------------------------------------------------------------------
if(NCYCLE>0) then
 write(FILENAME,10201)'vf_t',NCYCLE,'.bin'
else
 FILENAME='vf.end'
end if 


call SAVE_MPIIO(VFLU,FILENAME)



!!--------------------------------------------------------------------
!! 2.1.3. Wf
!!--------------------------------------------------------------------
if(NCYCLE>0) then
 write(FILENAME,10201)'wf_t',NCYCLE,'.bin'
else
 FILENAME='wf.end'
end if 

call SAVE_MPIIO(WFLU,FILENAME)



if(MYID==0) write(*,*) 'Save fluid solution --> MPI I/O'

!!--------------------------------------------------------------------
!! 2.2. Save the right-hand-side of fluid momentum
!!--------------------------------------------------------------------

if (SOLVE_FLUID==1) then
!!--------------------------------------------------------------------
!! 2.2.1. Uf
!!--------------------------------------------------------------------
  FILENAME='rhs_uf.end'
  call SAVE_MPIIO_RHS(RHS_UFOU,FILENAME)


!!--------------------------------------------------------------------
!! 2.2.2. Vf
!!--------------------------------------------------------------------
  FILENAME='rhs_vf.end'
  call SAVE_MPIIO_RHS(RHS_VFOU,FILENAME)


!!--------------------------------------------------------------------
!! 2.2.3. Wf
!!--------------------------------------------------------------------
  FILENAME='rhs_wf.end'
  call SAVE_MPIIO_RHS(RHS_wFOU,FILENAME)

  if(MYID==0) write(*,*) 'Save RHS --> MPI I/O'
end if

end if



!!====================================================================
!! 3. Statistics of the saved file
!!====================================================================

!~ STATX(:) = ZERO
!~ STATY(:) = ZERO
!~ STATZ(:) = ZERO

!~ !!- 1st order moment
!~ call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATX(1))
!~ call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATY(1))
!~ call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),STATZ(1))

!~ !!- 2nd order moment
!~ call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATX(2))
!~ call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATY(2))
!~ call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),STATZ(2))
!~ 
!~ !!- 3rd order moment
!~ call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATX(3))
!~ call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATY(3))
!~ call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**3),STATZ(3))
!~ 
!~ !!- 4th order moment
!~ call RSUMCPU(SUM(UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATX(4))
!~ call RSUMCPU(SUM(VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATY(4))
!~ call RSUMCPU(SUM(WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**4),STATZ(4))


!~ if(MYID==0) write(*,*) '       <Ux> = ',STATX(1) / NGLOB
!~ if(MYID==0) write(*,*) '       <Uy> = ',STATY(1) / NGLOB
!~ if(MYID==0) write(*,*) '       <Uz> = ',STATZ(1) / NGLOB
!~ if(MYID==0) write(*,*) '     <Ux^2> = ',STATX(2) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uy^2> = ',STATY(2) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uz^2> = ',STATZ(2) / NGLOB
!~ if(MYID==0) write(*,*) '     <Ux^3> = ',STATX(3) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uy^3> = ',STATY(3) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uz^3> = ',STATZ(3) / NGLOB
!~ if(MYID==0) write(*,*) '     <Ux^4> = ',STATX(4) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uy^4> = ',STATY(4) / NGLOB
!~ if(MYID==0) write(*,*) '     <Uz^4> = ',STATZ(4) / NGLOB
if(MYID==0) write(*,*)


!!====================================================================
10100 format (A,A,A)
10200 format (A,A,A,I8.8,A)
10201 format (A,I8.8,A)
10300 format (A,A)
10400 format (A,F8.2,A)


end subroutine SAVE_FLUID
