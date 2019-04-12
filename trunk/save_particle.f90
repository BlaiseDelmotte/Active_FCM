!!====================================================================
!!
!!
!!====================================================================

subroutine SAVE_PARTICLE(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use CHECK_CPU
use MPI_STRUCTURES


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Cycle number
integer, intent(in) :: NCYCLE

!- File name 
character(len=40) :: FILENAME

integer :: IDUMMY
integer :: RECSIZE, SIZE_INT, SIZE_REAL, NVARIABLE

!- Variable for MPI I/O
integer(kind=MPI_OFFSET_KIND) :: POS_FILE, POS_OFFSET
integer :: FILETYPE
integer :: LOC_SIZE
integer :: DESCRIPTEUR
real (kind=8), dimension(:), allocatable :: TEMP_VAR


!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J, NP
!---------------------------------------------------------------------


!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if



!!====================================================================
!! 1. Multiple binary files
!!====================================================================
if(ISAVEPART == 1) then 

!- Define file name
if(NCYCLE>=0) then
 write(FILENAME,10403)'PART',trim(FILE_EXT),'_t',NCYCLE,'.bin'
else
 write(FILENAME,10303)'PART',trim(FILE_EXT),'.end'
end if 



!- Open file containing the last particle position and velocity
open(unit = 150, file = trim(FILENAME),form='unformatted')

do J = 1, NIG

 write(150)NPART_LOC(J)
 write(150)(PART(I,J)%XP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%YP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%ZP,I=1,NPART_LOC(J))

 write(150)(PART(I,J)%UP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WP,I=1,NPART_LOC(J))
 
 write(150)(PART(I,J)%UFAP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VFAP,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WFAP,I=1,NPART_LOC(J))

 if(SOLVE_SCALAR) then
   write(150)(PART(I,J)%TP,I=1,NPART_LOC(J))
   write(150)(PART(I,J)%TFAP,I=1,NPART_LOC(J))   
 end if

end do !!- Loop:J = 1, NIG


!- Close file
close(150)



!!===============================================================
!! 1.2 Drop the right-hand-side
!!===============================================================
if(NCYCLE<0) then  

!! As the RHS is computed at each time step, up(tn-1), up(tn-2),
!! ufap(t-1) and ufap(tn-2) are droppped

write(FILENAME,10303)'rhs_part',trim(FILE_EXT),'.end'
!- Open file containing the
open(unit = 150, file = trim(FILENAME),form='unformatted')

do J = 1, NIG

 write(150)NPART_LOC(J)

 write(150)(PART(I,J)%UP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%UP_NM2,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VP_NM2,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WP_NM2,I=1,NPART_LOC(J))
 
 write(150)(PART(I,J)%UFAP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%UFAP_NM2,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VFAP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%VFAP_NM2,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WFAP_NM1,I=1,NPART_LOC(J))
 write(150)(PART(I,J)%WFAP_NM2,I=1,NPART_LOC(J))

 if(SOLVE_SCALAR) then
   write(150)(PART(I,J)%TP_NM1,I=1,NPART_LOC(J))
   write(150)(PART(I,J)%TP_NM2,I=1,NPART_LOC(J))
   write(150)(PART(I,J)%TFAP_NM1,I=1,NPART_LOC(J))
   write(150)(PART(I,J)%TFAP_NM2,I=1,NPART_LOC(J))
 end if
  
end do !!- Loop:J = 1, NIG

!- Close file
close(150)


 if(MYID==0)write(*,*) 'Final particle position and velocity --> dropped'
 if(MYID==0)write(*,*) '     + Multiple Binary Files'
 if(MYID==0)write(*,*) 'Rhs particle position and velocity --> dropped'
 if(MYID==0)write(*,*) '     + Multiple Binary Files'

end if


if(NCYCLE>0) then 
 if(MYID==0)write(*,*) 'Particle position and velocity file dropped --> ',trim(FILENAME)
end if

!!====================================================================
!! 2. Direct access file
!!====================================================================
!!- IMPORTANT
!!  When using the standard FORTRAN binary file declaration
!!  "unformatted" we found binaries 4 tiles larger than the
!!  expected size. Up to now the solution is to use "binary" 
!!  that probably write C formatted binary ... 
!!--------------------------------------------------------------------
elseif(ISAVEPART == 2) then


!! 
call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,SIZE_REAL,IERR)
call MPI_TYPE_SIZE(MPI_INTEGER,SIZE_INT,IERR)


!!- Number of saved variables
!! Basically 6 variables are stored: xp, yp, zp, up, vp, wp
 NVARIABLE = 6

!!- in case of add temperature
if(SOLVE_SCALAR) NVARIABLE = NVARIABLE + 1


!!- Compute the length of the record 
RECSIZE = (SIZE_INT + NPMAX_LOC*SIZE_REAL*NVARIABLE)*NIG


!- Define file name
 FILENAME = 'traj.end'
 
open(unit=150,          &
      file=trim(FILENAME),&
      status='replace',   &
      access='direct',    &
!      action='write',      &
      recl=RECSIZE,       &
      form='unformatted') ! Had to change 'binary' to 'unformatted' to make it compile with mpifort

if(SOLVE_SCALAR) then
write(unit=150,rec=1+MYID)  &
            (NPART_LOC(J), &
	     (PART(I,J)%XP, &
              PART(I,J)%YP, &
              PART(I,J)%ZP, &
              PART(I,J)%UP, &
              PART(I,J)%VP, &
              PART(I,J)%WP, &
	      PART(I,J)%TP, &
	      I=1,NPMAX_LOC), J=1,NIG)
else
write(unit=150,rec=1+MYID)  &
            (NPART_LOC(J),  &
	     (PART(I,J)%XP, &
              PART(I,J)%YP, &
              PART(I,J)%ZP, &
              PART(I,J)%UP, &
              PART(I,J)%VP, &
              PART(I,J)%WP, &
	      I=1,NPMAX_LOC), J=1,NIG)
end if

call MPI_BARRIER(MPI_COMM_WORLD,IERR)
close(150)


if(NCYCLE<0) then 
 if(MYID==0)write(*,*) 'Final particle position and velocity --> Saved'
 if(MYID==0)write(*,*) '     + Direct Access Binary File'
else
 if(MYID==0)write(*,*) 'Particle position and velocity file dropped --> ',trim(FILENAME)
end if


!!---------------------------------
!!- Number of saved variables
!! Basically  1 variable is stored: Tp
! NVARIABLE = 1
!
!!- Compute the length of the record 
!RECSIZE = NPMAX_LOC*SIZE_REAL*NVARIABLE*NIG
!
!
!- Define file name
! FILENAME = 'traj_scl.end'
! 
!open(unit=150,          &
!      file=trim(FILENAME),&
!      status='replace',   &
!      access='direct',    &
!!      action='write',      &
!      recl=RECSIZE,       &
!      form='binary')
!
!write(unit=150,rec=1+MYID)  &
!            ((PART(I,J)%TP,I=1,NPMAX_LOC), J=1,NIG)


!call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!close(150)

!!---------------------------------

if(DEBUG) then
do J=1, NIG
call ISUMCPU(NPART_LOC(J),IDUMMY)
if (MYID==0) write(*,*) 'Save Part -> Class:',J,' Full number of particles ',IDUMMY
end do
!do J=1, NIG
!write(600+MYID,*)MYID,'IG=',J,'  NPART_LOC=',NPART_LOC(J)
!end do
end if



if(DEBUG) then
 if(MYID==0)write(*,*) 'Save Part -> NVARIABLE = ',NVARIABLE
 if(MYID==0)write(*,*) 'Save Part ->   RECSIZE = ',RECSIZE
 if(MYID==0)write(*,*) 'Save Part -> NPMAX_LOC = ',NPMAX_LOC
 if(MYID==0)write(*,*) 'Save Part -> SIZE_REAL = ',SIZE_REAL
 if(MYID==0)write(*,*) 'Save Part ->  SIZE_INT = ',SIZE_INT
end if
if(MYID==0)write(*,*) '    + Direct access file'
if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)
if(MYID==0)write(*,10200) '     + File size = ',real(RECSIZE*NPROC/(2**20)),' Mo'





!!====================================================================
!! 3. MPI I/O
!!====================================================================
elseif(ISAVEPART == 3) then

FILENAME = 'traj.end'

call SAVE_PART_MPIIO(PART,FILENAME)


if(MYID==0)write(*,*) 'Final particle position and velocity --> Saved'
if(MYID==0)write(*,*) '    + MPI I/O'
if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)


end if !!- If: ISAVEPART



!!--------------------------------------------------------------------
!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(7) = CPU_PART(7) + TIME_END - TIME_START
end if

!!--------------------------------------------------------------------
10000 format (15(e17.7))
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10200 format (A,F8.2,A)
10300 format (A,A)

10303 format (A,A,A)
10403 format (A,A,A,I8.8,A)
20000 format (I2,2x,I2,2x,I2,10(e17.7))

end subroutine SAVE_PARTICLE
