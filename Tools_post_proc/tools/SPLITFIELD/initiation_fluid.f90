!!====================================================================
!!
!!
!!====================================================================

subroutine INITIATION_FLUID

!!====================================================================
!!
!!
!!====================================================================

use dns_dim
use geometric_variable
use fluid_variable
use param_phys
use forcing
use mpi_structures

implicit none


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
real(kind=8) :: X0, Y0, Z0
real(kind=8) :: RADIUS, U0, V0, W0, OMEGA, PSI

!- random number
real(kind=8) :: XRAND

!- File name 
character(len=40) :: FILENAME
!- 
integer, dimension(3) :: ISIZE_READ, ISTART_READ, IEND_READ

!!- Read from NTMIX
real(kind=8), allocatable, dimension(:,:,:,:) :: UREAD
real(kind=8) :: UNTMIX, LNTMIX, TNTMIX
real(kind=8) :: RDUMMY

integer :: NSPEC,NREAC,SPEC,REACT
integer :: JX,JY,JZ,IVAR,ILFIN,JLFIN,KLFIN,KCUT
integer :: IX,IY,IZ,FX,FY,FZ,IDUMMY,ILDEB,JLDEB,KLDEB

!- File descriptor
integer :: DESCRIPTEUR

integer :: RECSIZE

!- Index
integer :: I, J, K, IJK

!---------------------------------------------------------------------

UREF(:) = ZERO

!!====================================================================
!! 0. FLUID VELOCITY EQUAL TO ZERO
!!====================================================================
if(INIT_FLUID_VELOCITY == 0) then

 UFLU = ZERO 
 VFLU = ZERO
 WFLU = ZERO

 !- Init forcing
 if(STEADY) then 
  FORCING_UFOU(:) = cmplx(ZERO,ZERO)
  FORCING_VFOU(:) = cmplx(ZERO,ZERO)
  FORCING_WFOU(:) = cmplx(ZERO,ZERO)

  if(MYID==0) write(*,*)'Forcing coefficients initiation: F=0 --> OK'

 end if

 if(MYID==0) write(*,*)'Fluid velocity initiation: Uf(t=0)=0 --> OK'

!!====================================================================
!! 1. SINGLE EDDY
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 1) then

 RADIUS = LXMAX/10.
 OMEGA = 0.1 !- 1/s

 U0 = 1.
 V0 = ZERO
 W0 = ZERO

 X0 = LXMAX/2.
 Y0 = LYMAX/2.
 Z0 = ZERO

 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)

    PSI = OMEGA*exp(-((XMESH(I)-X0)**2+(YMESH(J)-Y0)**2)/(2.*RADIUS**2))

    UFLU(I,J,K) = U0 - (YMESH(J)-Y0)/RADIUS**2*PSI
    VFLU(I,J,K) = V0 + (XMESH(I)-X0)/RADIUS**2*PSI
    WFLU(I,J,K) = ZERO

!    RADIUS = sqrt((XMESH(I)-LXMAX/2)**2 + (YMESH(J)-LYMAX/2)**2)
!    call random_number(XRAND)
!    if(RADIUS < LXMAX/5) then
!      UFLU(I,J,K) = 10. + XRAND
!    else
!      UFLU(I,J,K) = 0.
!    end if     


   end do
  end do
 end do

 !- Init forcing
 if(STEADY) then
  FORCING_UFOU(:) = cmplx(ZERO,ZERO)
  FORCING_VFOU(:) = cmplx(ZERO,ZERO)
  FORCING_WFOU(:) = cmplx(ZERO,ZERO)

  if(MYID==0) write(*,*)'Forcing coefficients initiation: F=0 --> OK'

 end if



 UREF(1) = U0
 UREF(2) = V0
 UREF(3) = W0


 if(MYID==0) write(*,*)'Fluid velocity initiation: 2d eddy --> OK'



!!====================================================================
!! 2. RANDOM FLUID VELOCITY
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 2) then

 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)

    call random_number(XRAND)
    UFLU(I,J,K) = 0.5 - XRAND

    call random_number(XRAND)
    VFLU(I,J,K) = 0.5 - XRAND

    call random_number(XRAND)
    WFLU(I,J,K) = 0.5 - XRAND

   end do
  end do
 end do

 !- Init forcing
 if(STEADY) then
  FORCING_UFOU(:) = cmplx(ZERO,ZERO)
  FORCING_VFOU(:) = cmplx(ZERO,ZERO)
  FORCING_WFOU(:) = cmplx(ZERO,ZERO)

  if(MYID==0) write(*,*)'Forcing coefficients initiation: F=0 --> OK'

 end if



 if(MYID==0) write(*,*)'Fluid velocity initiation: Random --> OK'




!!====================================================================
!! 3. Restart file
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 3) then

!!--------------------------------------------------------------------
!! 3.1 Multiple Binary Files
!!--------------------------------------------------------------------
if(ISAVEM == 1) then 

 !- Define file name
 write(FILENAME,10101)'FLUID',trim(FILE_EXT),'.ini'

 !- Open file containing the initial fluid velocity field
 open(unit = 120, file=trim(FILENAME), form='unformatted')

 !- Read size of stored file
 read(120)ISIZE_READ(1),ISTART_READ(1),IEND_READ(1)
 read(120)ISIZE_READ(2),ISTART_READ(2),IEND_READ(2)
 read(120)ISIZE_READ(3),ISTART_READ(3),IEND_READ(3)

 !- Check size of the stored field 
 if(	(ISIZE_READ(1)/=ISIZE(1)) &
    .or.(ISIZE_READ(2)/=ISIZE(2)) &
    .or.(ISIZE_READ(3)/=ISIZE(3)) ) then
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!! 		    ERROR		     !!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!'
  write(*,*)'!! file	 : init_fluid.f90'
  write(*,*)'!!'
  write(*,*)'!! PROBLEM DESCRIPTION: '
  write(*,*)'!! initial fluid velocity field has'
  write(*,*)'!! wrong dimensions' 
  write(*,*)'!!'
  write(*,*)'!! NX_READ=',ISIZE_READ(1),'  local NX=',ISIZE(1)
  write(*,*)'!! NY_READ=',ISIZE_READ(2),'  local NY=',ISIZE(2)
  write(*,*)'!! NZ_READ=',ISIZE_READ(3),'  local NZ=',ISIZE(3)
  write(*,*)'!!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  stop
 end if  

 read(120)(((UFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
 read(120)(((VFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))
 read(120)(((WFLU(I,J,K),I=ISTART(1),IEND(1)),J=ISTART(2),IEND(2)),K=ISTART(3),IEND(3))

 !- Close file
 close(120)



!!--------------------------------------------------------------------
!! 3.2 Direct access file
!!--------------------------------------------------------------------
elseif(ISAVEM == 2) then 

 RECSIZE = 8*ISIZE(1)*ISIZE(2)*ISIZE(3)

!!--------------------------------------------------------------------
!!- x-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='uf.ini'

 open(unit=120,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 read(unit=120,rec=MYID+1)UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(120)

!!--------------------------------------------------------------------
!! y-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='vf.ini'

 open(unit=121,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 read(unit=121,rec=MYID+1)VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(121)

!!--------------------------------------------------------------------
!! z-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='wf.ini'

 open(unit=122,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

 read(unit=122,rec=MYID+1)WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

 call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 close(122)



!!--------------------------------------------------------------------
!! 3.3 MPI I/O Ordered
!!--------------------------------------------------------------------
elseif(ISAVEM == 3) then 

!!--------------------------------------------------------------------
!!- x-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='uf.ini'


 RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)

 call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                       trim(FILENAME), &
                      MPI_MODE_RDONLY, &
                        MPI_INFO_NULL, &
                          DESCRIPTEUR, &
                                 IERR  )

! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
 call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
 UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                  RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)

!!--------------------------------------------------------------------
!!- y-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='vf.ini'


 RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)

 call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                       trim(FILENAME), &
                      MPI_MODE_RDONLY, &
                        MPI_INFO_NULL, &
                          DESCRIPTEUR, &
                                 IERR  )

! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
 call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
 VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                  RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)

!!--------------------------------------------------------------------
!!- y-component
!!--------------------------------------------------------------------
!!- Print file name
 FILENAME='wf.ini'


 RECSIZE = ISIZE(1)*ISIZE(2)*ISIZE(3)

 call MPI_FILE_OPEN(   MPI_COMM_WORLD, &
                       trim(FILENAME), &
                      MPI_MODE_RDONLY, &
                        MPI_INFO_NULL, &
                          DESCRIPTEUR, &
                                 IERR  )

! call MPI_FILE_READ_ORDERED(DESCRIPTEUR, & 
 call MPI_FILE_READ_ALL(DESCRIPTEUR, & 
 WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)), &
                  RECSIZE, MPI_DOUBLE_PRECISION, statut, IERR)

! call MPI_BARRIER(MPI_COMM_WORLD,IERR)
 call MPI_FILE_CLOSE(DESCRIPTEUR, IERR)




!!--------------------------------------------------------------------
!! 3.3 MPI I/O READ_AT_ALL
!!--------------------------------------------------------------------
elseif(ISAVEM == 4) then

 FILENAME='uf.ini'
 call READ_MPIIO(UFLU,FILENAME)
 FILENAME='vf.ini'
 call READ_MPIIO(VFLU,FILENAME)
 FILENAME='wf.ini'
 call READ_MPIIO(WFLU,FILENAME)

end if

if(MYID==0) write(*,*)'Fluid velocity initiation: Read from file --> OK'


!!--------------------------------------------------------------------
!! 3.4. Read Stochastic force
!!--------------------------------------------------------------------
if(STEADY.and. MYID==0) then

 !- Define file name
 FILENAME='FORCING.ini'

 !- Open file containing the forcing coefficients
 open(unit =120, file=trim(FILENAME), form='unformatted')

 !- Number of forced waves
 read(120) NFORCE_FULL

 !- Random seed
 read(120) IDFORCE
!!  read(120) IDUMMY

 !- Forcing coefficients
 read(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
 read(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
 read(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)

  !- Index of hermit coefficient
!  read(120)(NHERM(I),I=1,NFORCE_FULL)

 !- Close file
 close(120)

 write(*,*)'Forcing coefficients initiation: Read from file --> OK'

end if

if(STEADY.and. NPROC>1) then
 call MPI_BCAST(FORCING_UFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_VFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_WFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
end if





!!====================================================================
!! 4. Restart from NTMIX
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 4) then

 LNTMIX =    LXMAX / (2.*PPI)
 UNTMIX = 608*VISC / LNTMIX  !!- 608 is the accoustic Reynolds number
 TNTMIX =   LNTMIX / UNTMIX

 if(MYID==0) then 
  write(*,*)'Initiation from NTMIX files: UREF = ',UNTMIX
  write(*,*)'                             LREF = ',LNTMIX
  write(*,*)'                             TREF = ',TNTMIX
 end if

 !- Define file name
 FILENAME='FLUID.NTMIX'

 !- Open file containing the initial fluid velocity field
 open(unit = 120, file =trim(FILENAME), form='unformatted')

 read(120) FX,FY,FZ
!! if(MYID==0) write(*,*) 'FX,FY,FZ',FX,FY,FZ

 read(120) IX,IY,IZ,FX,FY,FZ
!! if(MYID==0) write(*,*) 'IX,IY,IZ,FX,FY,FZ',IX,IY,IZ,FX,FY,FZ

 read(120) SPEC,REACT
!! if(MYID==0) write(*,*) 'spec,react',SPEC,REACT

 read(120) RDUMMY
!! if(MYID==0) write(*,*) 'Time',RDUMMY

 read(120) RDUMMY,RDUMMY,(RDUMMY,i=1,NSPEC)
 read(120) RDUMMY,RDUMMY,RDUMMY,RDUMMY
 read(120) RDUMMY,RDUMMY,RDUMMY,RDUMMY,RDUMMY,RDUMMY

 read(120) (RDUMMY,jx=IX,FX)
 read(120) (RDUMMY,jy=IY,FY)
 read(120) (RDUMMY,jz=IZ,FZ)


 allocate(UREAD(FX,FY,FZ,6))

 do IVAR = 1, 6
  do JZ = IZ,FZ
   read(120)((UREAD(JX,JY,JZ,IVAR),JX=IX,FX),JY=IY,FY)
  end do
 end do

 !- Close file
 close(120)


 do K = ISTART(3), IEND(3)
  do J = ISTART(2), IEND(2)
   do I = ISTART(1), IEND(1)
    UFLU(I,J,K) = UREAD(I,J,K,2)*UNTMIX
    VFLU(I,J,K) = UREAD(I,J,K,3)*UNTMIX
    WFLU(I,J,K) = UREAD(I,J,K,4)*UNTMIX
   end do
  end do
 end do

 deallocate(UREAD)



end if



if(MYID==0) write(*,*) 'Fluid initiation --> OK'


!!====================================================================
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10101 format (A,A,A)

end subroutine INITIATION_FLUID
