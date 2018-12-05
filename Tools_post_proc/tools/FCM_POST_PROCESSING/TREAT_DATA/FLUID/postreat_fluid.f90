!!====================================================================
!!
!!
!!====================================================================

Program POSTREAT_FLUID

!!====================================================================
!!
!!
!!====================================================================

use MPI

implicit none


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

!- File name 
character(len=40) :: FILENAME


!!- Vraiables from direct access files
real(kind=8), allocatable, dimension(:,:,:) :: UFLU
real(kind=8), allocatable, dimension(:,:,:) :: VFLU
real(kind=8), allocatable, dimension(:,:,:) :: WFLU


!!- Variables from post-treatment
real(kind=8), allocatable, dimension(:,:,:) :: NORM_VEL

!!- Dicrestization points in each direction
integer :: NX, NY, NZ

!- Record size
integer :: RECSIZE

!- Index
integer :: I, J, K, IJK

!-  MPI Variables
integer :: IERR,NPROC,MYID



!!--------------------------------------------------------------------
!! 0.0 MPI World initiation
!!--------------------------------------------------------------------
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)




!---------------------------------------------------------------------
!=====================================================================
! 0. DISCRETIZATION PARAMETERS AND ARRAY ALLOCATION
!=====================================================================

NX = 128
NY = NX
NZ = NX


allocate( UFLU(NX,NY,NZ), VFLU(NX,NY,NZ), WFLU(NX,NY,NZ) )

allocate( NORM_VEL(NX,NY,NZ) )

!=====================================================================
! 1. FLUID VELOCITY
!=====================================================================

!!--------------------------------------------------------------------
!! 1.1 Read Direct access file written with MPI/IO
!!--------------------------------------------------------------------
!!--------------------------------------------------------------------
!!- x-component
!!--------------------------------------------------------------------
!!- Print file name
FILENAME='uf.end'
call READ_MPIIO(NX,NY,NZ,UFLU,FILENAME)

!!--------------------------------------------------------------------
!! y-component
!!--------------------------------------------------------------------
!!- Print file name
FILENAME='vf.end'
call READ_MPIIO(NX,NY,NZ,VFLU,FILENAME)

!!--------------------------------------------------------------------
!! z-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='wf.end'
call READ_MPIIO(NX,NY,NZ,WFLU,FILENAME)


!!--------------------------------------------------------------------
!! 1.2 Compute velocity norm
!!--------------------------------------------------------------------
 NORM_VEL = dsqrt( UFLU**2 + VFLU**2 + WFLU**2 ) 

 
 print*,'maxval(UFLU) = ', maxval(UFLU)
 print*,'maxval(VFLU) = ', maxval(VFLU)
 print*,'maxval(WFLU) = ', maxval(WFLU)
 
 print*,'maxval(NORM_VEL) = ', maxval(NORM_VEL)
 
 
 !- Free MPI environement
 call MPI_FINALIZE (ierr)

!~ 
!~ !!--------------------------------------------------------------------
!~ !! 3.3 MPI I/O Ordered
!~ !!--------------------------------------------------------------------
!~ elseif(ISAVEFLUID == 3) then 
!~ 
!~ !!--------------------------------------------------------------------
!~ !!- x-component
!~ !!--------------------------------------------------------------------
!~ !!- Print file name
!~  FILENAME='uf.ini'
!~ 
!~ 
!~  RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)
!~ 
!~  call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
!~                        trim(FILENAME), &
!~                       MPI_MODE_RDONLY, &
!~                         MPI_INFO_NULL, &
!~                           DESCRIPTEUR, &
!~                                  IERR  )
!~ 
!~ ! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
!~  call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
!~  UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
!~                   RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)
!~ 
!~ ! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!~  call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)
!~ 
!~ !!--------------------------------------------------------------------
!~ !!- y-component
!~ !!--------------------------------------------------------------------
!~ !!- Print file name
!~  FILENAME='vf.ini'
!~ 
!~ 
!~  RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)
!~ 
!~  call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
!~                        trim(FILENAME), &
!~                       MPI_MODE_RDONLY, &
!~                         MPI_INFO_NULL, &
!~                           DESCRIPTEUR, &
!~                                  IERR  )
!~ 
!~ ! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
!~  call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
!~  VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
!~                   RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)
!~ 
!~ ! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!~  call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)
!~ 
!~ !!--------------------------------------------------------------------
!~ !!- y-component
!~ !!--------------------------------------------------------------------
!~ !!- Print file name
!~  FILENAME='wf.ini'
!~ 
!~ 
!~  RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)
!~ 
!~  call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
!~                        trim(FILENAME), &
!~                       MPI_MODE_RDONLY, &
!~                         MPI_INFO_NULL, &
!~                           DESCRIPTEUR, &
!~                                  IERR  )
!~ 
!~ ! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
!~  call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
!~  WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
!~                   RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)
!~ 
!~ ! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
!~  call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)
!~ 
!~ 
!~ 
!~ 
!~ !!--------------------------------------------------------------------
!~ !! 3.4 MPI I/O 
!~ !!--------------------------------------------------------------------
!~ elseif(ISAVEFLUID == 4) then
!~ 
!~  FILENAME='uf.ini'
!~  call READ_MPIIO(UFLU,FILENAME)
!~  FILENAME='vf.ini'
!~  call READ_MPIIO(VFLU,FILENAME)
!~  FILENAME='wf.ini'
!~  call READ_MPIIO(WFLU,FILENAME)
!~ 
!~ 
!~ end if
!~ 
!~ 
!~ 
!~ if(MYID==0) write(*,*)'Fluid velocity initiation: Read from file --> OK'
!~ 
!~ 
!~ !!--------------------------------------------------------------------
!~ !! 3.4. Read Stochastic force
!~ !!--------------------------------------------------------------------
!~ if(STEADY.and. MYID==0) then
!~ 
!~  !- Define file name
!~  FILENAME='forcing.ini'
!~ 
!~  !- Open file containing the forcing coefficients
!~  open(unit =120, file=trim(FILENAME), form='unformatted')
!~ 
!~  !- Number of forced waves
!~  read(120) NFORCE_FULL
!~ 
!~  !- Random seed
!~  read(120) IDFORCE
!~ !!  read(120) IDUMMY
!~ 
!~  !- Forcing coefficients
!~  read(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
!~  read(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
!~  read(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)
!~ 
!~   !- Index of hermit coefficient
!~ !  read(120)(NHERM(I),I=1,NFORCE_FULL)
!~ 
!~  !- Close file
!~  close(120)
!~ 
!~  write(*,*)'Forcing coefficients initiation: Read from file --> OK'
!~ 
!~ end if
!~ 
!~ if(STEADY.and. NPROC>1) then
!~  call MPI_BCAST(FORCING_UFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
!~  call MPI_BCAST(FORCING_VFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
!~  call MPI_BCAST(FORCING_WFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
!~ end if
!~ 
!~ 
!~ 
!~ 
!~ 
!~ !!====================================================================
!~ !! 4. Restart from NTMIX
!~ !!====================================================================
!~ elseif(INIT_FLUID_VELOCITY == 4) then
!~ 
!~  LNTMIX =    LXMAX / (2.*PPI)
!~  UNTMIX = 608*VISC / LNTMIX  !!- 608 is the accoustic Reynolds number
!~  TNTMIX =   LNTMIX / UNTMIX
!~ 
!~  if(MYID==0) then 
!~   write(*,*)'Initiation from NTMIX files: UREF = ',UNTMIX
!~   write(*,*)'                             LREF = ',LNTMIX
!~   write(*,*)'                             TREF = ',TNTMIX
!~  end if
!~ 
!~  !- Define file name
!~  FILENAME='FLUID.NTMIX'
!~ 
!~  !- Open file containing the initial fluid velocity field
!~  open(unit = 120, file =trim(FILENAME), form='unformatted')
!~ 
!~  read(120) FX,FY,FZ
!~ !! if(MYID==0) write(*,*) 'FX,FY,FZ',FX,FY,FZ
!~ 
!~  read(120) IX,IY,IZ,FX,FY,FZ
!~ !! if(MYID==0) write(*,*) 'IX,IY,IZ,FX,FY,FZ',IX,IY,IZ,FX,FY,FZ
!~ 
!~  read(120) SPEC,REACT
!~ !! if(MYID==0) write(*,*) 'spec,react',SPEC,REACT
!~ 
!~  read(120) RDUMMY
!~ !! if(MYID==0) write(*,*) 'Time',RDUMMY
!~ 
!~  read(120) RDUMMY,RDUMMY,(RDUMMY,i=1,NSPEC)
!~  read(120) RDUMMY,RDUMMY,RDUMMY,RDUMMY
!~  read(120) RDUMMY,RDUMMY,RDUMMY,RDUMMY,RDUMMY,RDUMMY
!~ 
!~  read(120) (RDUMMY,jx=IX,FX)
!~  read(120) (RDUMMY,jy=IY,FY)
!~  read(120) (RDUMMY,jz=IZ,FZ)
!~ 
!~ 
!~  allocate(UREAD(FX,FY,FZ,6))
!~ 
!~  do IVAR = 1, 6
!~   do JZ = IZ,FZ
!~    read(120)((UREAD(JX,JY,JZ,IVAR),JX=IX,FX),JY=IY,FY)
!~   end do
!~  end do
!~ 
!~  !- Close file
!~  close(120)
!~ 
!~ 
!~  do K = ISTART(3), IEND(3)
!~   do J = ISTART(2), IEND(2)
!~    do I = ISTART(1), IEND(1)
!~     UFLU(I,J,K) = UREAD(I,J,K,2)*UNTMIX
!~     VFLU(I,J,K) = UREAD(I,J,K,3)*UNTMIX
!~     WFLU(I,J,K) = UREAD(I,J,K,4)*UNTMIX
!~    end do
!~   end do
!~  end do
!~ 
!~  deallocate(UREAD)
!~ 
!~ 
!~ 
!~ end if
!~ 
!~ 
!~ 
!~ if(MYID==0) write(*,*) 'Fluid initiation --> OK'


!!====================================================================
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10101 format (A,A,A)

end program POSTREAT_FLUID
