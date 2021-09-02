!!====================================================================
!! 
!!  
!!====================================================================

program ns3d

!!====================================================================
!! Numerical code solving the 3-dimensionnal Navier-Stokes equations
!! using a spectral decomposition. The code is coupled with a 
!! Lagrangian particle tracking module.
!!
!! The numerical features are:
!!  * Fluid: 
!!       - 3d Navier-Stokes equation
!!       - Incompressible flow
!!       - 3rd order Adam-Bashforth (time-advancing)
!!       - Integrating factor (viscous terms)
!!
!!  * Scalar: 
!!       - Transport equations
!!       - 3rd order Adam-Bashforth (time-advancing)
!!       - Integrating factor (diffusive terms)
!!
!!  * Particles: 
!!       - fixed points, fluid elements, inertial particles
!!       - 3rd order Adam-Bashforth (time-advancing)
!!       - Integrating factor (diffusive terms)
!!====================================================================

use DNS_DIM            !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use FORCING            !- Forcing
use FLUID_VARIABLE     !- Fluid velocity
use SCALAR_VARIABLE
use GEOMETRIC_VARIABLE !- 
use RHS_VARIABLES      !- Right-Hand-Side
use STATISTICS         !- Statistics
use WORK_ARRAYS
use CHECK_CPU          !- Variable for cpu time checking
use CPUTIME_CONTROL
use FCM_PART_VARIABLE  !- FCM particles variables
use FCM_FORCING_VARIABLE  !- FCM forcing variables



use MPI_STRUCTURES

use P3DFFT

implicit none


!---------------------------------------------------------------------
! ARGUMENT STATEMENT
!---------------------------------------------------------------------
integer :: IARGC,NB_ARGS
character(len=10) :: ARG

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Current time
real(kind=8) :: TIME

!- Integer flag for subroutine argument
integer      :: IFLAG1, IFLAG2
logical      :: LFLAG1
real(kind=8) :: RDUMMY
integer      :: IDUMMY, IDUMMY2

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Time measure variable
real(kind=8) :: MEASURE_START, MEASURE_END

real etime          ! Declare the type of etime()
real elapsed(2)     ! For receiving user and system time
real elapsed2(2)     ! For receiving user and system time
real total          ! For receiving total time
real total2          ! For receiving total time

!- Stop flag
logical :: CONT, CONT_CPU

!- Index
integer :: I, J, K, NCYCLE

real(kind=8) :: TIME0, TIME1

!-
integer :: NOUT2, NOUT3

!- Temporary index of switching for Time integration scheme
integer :: ISAVE
!---------------------------------------------------------------------





!!====================================================================
!! Get walltime
!!====================================================================

! lecture eventuelle du decoupage sur la ligne de commande
NB_ARGS = IARGC()


if ( NB_ARGS /= 0 ) then
 call getarg(1,ARG)
 read(ARG,*) WALLTIME
else 
 WALLTIME = INFINITY
end if

DEBUG = .false.


!!====================================================================
!! 1. INITIATION
!!====================================================================

!!--------------------------------------------------------------------
!! 1.1 MPI World initiation
!!--------------------------------------------------------------------
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)

if (MYID==0) write(*,*) 'WallTime =', WALLTIME


!- CPU
CPU_FLUID(:)=0.
CPU_ELAPSED = 0.
CPU_CYCLE = 0.
CPU_INITIATION = 0.

!!- Check cpu time
TIME_START = MPI_WTIME()

!!--------------------------------------------------------------------
!! 1.2 Read parameter
!!--------------------------------------------------------------------
call READPARAM

if (SOLVE_FLUID ==2)  call FCM_READPARAM


!!--------------------------------------------------------------------
!! 1.3 Run initiation
!!--------------------------------------------------------------------
call INIT_RUN




NCYCLE = 1 !- Initiation of time cycle
TIME = 0.  !- Initiation of time


!!--------------------------------------------------------------------
!! 1.4 Domain splitting
!!--------------------------------------------------------------------
!- Dimentionality of the cpu decomposition (=1: slab, =2: squared)
NDIM = 2

if(NDIM == 1) then
 DIMS(1) = 1
 DIMS(2) = NPROC
else if(NDIM == 2) then
 if (MYID==0) print *, 'Creating proc. grid with mpi_dims_create'
 DIMS(1) = 0
 DIMS(2) = 0
 call MPI_DIMS_CREATE(NPROC,2,DIMS,IERR)
 if(DIMS(1) > DIMS(2)) then
  DIMS(1) = DIMS(2)
  DIMS(2) = NPROC / DIMS(1)
 endif
endif
IPROC = DIMS(1)
JPROC = DIMS(2)

if(MYID == 0)write(*,*)'Using processor grid ',iproc,' x ',jproc


!!- FFT flag 
FFTFLAG = 'fft'

CPU_INIT = MPI_WTIME()

!!- FFt Initiation
! Call for EOS
  call P3DFFT_SETUP(DIMS,NX,NY,NZ,MPI_COMM_WORLD,NX,NY,NZ,.TRUE.)
! Call for hyperion and IMFT
! call P3DFFT_SETUP(DIMS,NX,NY,NZ,.TRUE.)


!!- Split the geometry
call P3DFFT_GET_DIMS(ISTART,IEND,ISIZE,1)
call P3DFFT_GET_DIMS(FSTART,FEND,FSIZE,2)

!!- Dimension initiation
NTOT = FSIZE(1)*FSIZE(2)*FSIZE(3)
NGLOB = NX * NY * NZ
FACTOR = 1.0D0/real(NGLOB)


!!- Create ghost cell if needed
if(NGHTCELL>0) then

! Create MPI structures for ghost cells MPI-exchange
 call CREATE_MPI_VECTOR
! Find neighbouring for each processor
 call NEIGHBOURING

end if




!!--------------------------------------------------------------------
!! 1.5. Allocate arrays
!!--------------------------------------------------------------------
call ALLOCATE_ARRAYS


!!--------------------------------------------------------------------
!! 1.6. Mesh, wavenumbers, integrating factor and aliasing control
!!--------------------------------------------------------------------
call MESHING



!!--------------------------------------------------------------------
!! 1.7. Initiation of Physics
!!--------------------------------------------------------------------


!!- FCM Initiation
if  (SOLVE_FLUID ==2) then

call FCM_INITIATION

end if

 

!!- Fluid Initiation
if(SOLVE_FLUID == 1) call INITIATION_FLUID

!- Update Ghost cells
if(NGHTCELL>0) then
 call FLUIDCOMM(UFLU)
 call FLUIDCOMM(VFLU)
 call FLUIDCOMM(WFLU)
end if



!!- Initiation forcing
if(STEADY) call INITIATION_FORCING



!!- Scalar initiation
if(SOLVE_SCALAR) then

 call INITIATION_SCALAR

 !- Update Ghost cells
 if(NGHTCELL>0) call FLUIDCOMM(THETA)

end if




!!- Time-averaging initiation
if(STAT_TIME) then
 if(LEVEL0_STFLU) MEAN_TIME_FLUID = ZERO
 if(LEVEL0_STSCL) MEAN_TIME_SCL = ZERO
 NEVEN = 0
end if


!!--------------------------------------------------------------------
!! 1.8. Opening files
!!--------------------------------------------------------------------
LFLAG1 = .true.
call OPENCLOSE(LFLAG1)


!!--------------------------------------------------------------------
!! 1.9. Print info about the numerical simulation
!!--------------------------------------------------------------------
if(MYID==0) call INFORUN




!!- End initiation
TIME_END = MPI_WTIME()
CPU_INITIATION = TIME_END - TIME_START

CPU_ELAPSED = CPU_ELAPSED + CPU_INITIATION

if (MYID == 0) then
 write(*,10700)100.*NCYCLE/NCYCLEMAX,CPU_ELAPSED
end if


CONT     = .true.
CONT_CPU = .true.

NOUT2 = 1
NOUT3 = 1


!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!
!== Time loop ==!== Time loop ==!== Time loop ==!== Time loop ==!

!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)


if(MYID==0)write(*,*)
if(MYID==0)write(*,*) 'Start the time-loop !!'
if(MYID==0)write(*,*)



do while(CONT)

 !!- CPU check
 if(MYID == 0) TIME_START = MPI_WTIME()


!!====================================================================
!! 2. STATISTICS
!!====================================================================
!! The statistic are compute at the begining of the time loop
!! because all variables are at the same cycle.
!!--------------------------------------------------------------------
 if(mod(NCYCLE,FOUT2) == 0) then


!!--------------------------------------------------------------------
!! 2.1 Fluid Statistics 
!!--------------------------------------------------------------------
 if(LEVEL0_STFLU) call STAT_FLUID(NCYCLE,TIME)


!!--------------------------------------------------------------------
!! 2.2 Scalar Statistics 
!!--------------------------------------------------------------------
 if(LEVEL0_STSCL) call STAT_SCALAR(NCYCLE,TIME)


!!- Event count for time-averaged statistics
 if((LEVEL0_STFLU.or.LEVEL0_STSCL.or.LEVEL0_STPAR).and.STAT_TIME)  NEVEN = NEVEN + 1

 end if




!!====================================================================
!! 4. Fluid flow prediction
!!====================================================================

!!--------------------------------------------------------------------
!! 4.1 Direct numerical simulation
!!--------------------------------------------------------------------
if(SOLVE_FLUID == 1) then

 call FLUID_PREDICTION(NCYCLE)
 
!- Update Ghost cell
 if(NGHTCELL>0) then
  call FLUIDCOMM(UFLU)
  call FLUIDCOMM(VFLU)
  call FLUIDCOMM(WFLU)
 end if



!!--------------------------------------------------------------------
!! 4.2 Solve Stokes equation --> Blaise your job !
!!--------------------------------------------------------------------
elseif (SOLVE_FLUID ==2)   then
 if (FCM_FLOW_TYPE==1) then
  call FCM_ALGORITHM(NCYCLE)
 else if (FCM_FLOW_TYPE==2) then
  call FCM_ALGORITHM_BROWNIAN(NCYCLE)
 end if

! if (MYID==0) then
!  if (mod(NCYCLE,FOUT3)==0) then
!   do I = 1, NPART_FULL
!    write(666,'(8e17.7)') real(NCYCLE),  real(I), &
!                         FCM_XP_NOPER(I), FCM_YP_NOPER(I), FCM_ZP_NOPER(I), &
!                         FCM_PSWIM(I,1), FCM_PSWIM(I,2), FCM_PSWIM(I,3)
!   end do
!  end if
! end if
end if


!!--------------------------------------------------------------------
!! 4.3 Channel flow --> Pascal your job !
!!--------------------------------------------------------------------



!!====================================================================
!! X. Manage time stepping
!!====================================================================
!! Two possibilities:
!! + First, the run reaches the user defined maximum step number
!! + The elapsed time reaches the limit imposed on the computer.
!!--------------------------------------------------------------------
 NCYCLE = NCYCLE + 1
 TIME = TIME + DTIME

!- Switch time integration pointers for RHS of fluid
 if(NCYCLE <= NCYCLEMAX) then
  ISAVE = TNM2
  TNM2 = TNM1
  TNM1 = TN
  TN = ISAVE
 end if !-> (NCYCLE <= NCYCLEMAX)


 if (SOLVE_FLUID/=2) then
 !!- Switch to 2nd order time-integration scheme
  if(NCYCLE == 2.and.INIT_FLUID_VELOCITY < 3) FLAG_TSCHEME = 2

 !!- Switch to 3rd order time-integration scheme
  if(NCYCLE == 3.and.INIT_FLUID_VELOCITY < 3) FLAG_TSCHEME = 3
  
 else
  if ((FCM_FLOW_TYPE==1).and.(FCM_RUNTUMBLE==0)) then
   if ((NCYCLE == 2).and.(SOLVE_COLLISION<4)) FLAG_TSCHEME = 2

   if ((NCYCLE == 3).and.(SOLVE_COLLISION<4)) FLAG_TSCHEME = 3

   ! test with AB4
   if ((NCYCLE == 4).and.(SOLVE_COLLISION<4)) FLAG_TSCHEME = 4
   
  !- If Brownian Motion-->> only use Forward Euler Scheme
  else if ((FCM_FLOW_TYPE==2).or.(FCM_RUNTUMBLE==1)) then
   
   FLAG_TSCHEME = 1
 
  end if
   

 end if

 if(NCYCLE == NCYCLEMAX) CONT = .false.


!!- Check Wall time 
 call TREMAIN(CONT_CPU)
 if(CONT_CPU) CONT = .false.
!!--------------------------------------------------------------------





!- CPU check
TIME_END = MPI_WTIME()
!!- Full CPU time elapsed
CPU_ELAPSED = CPU_ELAPSED + TIME_END - TIME_START

!!- Averaged CPU cycle time
CPU_CYCLE = CPU_CYCLE + TIME_END - TIME_START


!!- Print percentage of simulation accomplished
if((NCYCLEMAX>=FOUT0).and.(mod(NCYCLE,FOUT0)==0).and.(MYID==0)) then
 write(*,10700)100.*NCYCLE/NCYCLEMAX,CPU_ELAPSED
end if


end do




if(MYID==0)write(*,*)
if(MYID==0)write(*,*) ' Time loop ended !!!'
 


!!====================================================================
!! 5. STATISTICS FINALIZING
!!====================================================================
!- Synchronize all the process
!!call MPI_BARRIER(MPI_COMM_WORLD,IERR)


!!--------------------------------------------------------------------
!! 5.1 Last cycle
!!--------------------------------------------------------------------
 if(LEVEL0_STFLU) call STAT_FLUID(NCYCLEMAX,TIME)

 if(LEVEL0_STSCL) call STAT_SCALAR(NCYCLEMAX,TIME)


!!- Event count for time-averaged statistics
 if((LEVEL0_STFLU.or.LEVEL0_STSCL.or.LEVEL0_STPAR).and.STAT_TIME)  NEVEN = NEVEN + 1

 if(LEVEL0_STFLU.or.LEVEL0_STPAR.or.LEVEL0_STSCL) then
   if(MYID==0) call PRINT_LASTCYCLE(NCYCLE)
 end if


!!--------------------------------------------------------------------
!! 5.2. Print time-averaged statistics
!!--------------------------------------------------------------------
if(STAT_TIME) call PRINT_TIMESTAT(NCYCLE)


if(LEVEL2_STPAR) call PRINT_LAGFUNCTION



!!--------------------------------------------------------------------
!! 5.3. Compute and print last spectrum
!!--------------------------------------------------------------------
IFLAG1 = 99
if (SOLVE_FLUID == 1) call SPEC3D(IFLAG1)


IFLAG1 = 99
if(SOLVE_SCALAR) call SPEC3D_SCALAR(IFLAG1)



!!--------------------------------------------------------------------
!! 5.4. File closing
!!--------------------------------------------------------------------
LFLAG1 = .false.
call OPENCLOSE(LFLAG1)


!!====================================================================
!! 6. SAVE SOLUTION FOR RESTART
!!====================================================================

!!--------------------------------------------------------------------
!! 6.1. Fluid velocity field 
!!--------------------------------------------------------------------
!!- Print last fluid solution for restart
IDUMMY = -99
if(SOLVE_FLUID>0) call SAVE_FLUID(IDUMMY)


if(STEADY) call SAVE_FORCING


!!--------------------------------------------------------------------
!! 6.2. Scalar field
!!--------------------------------------------------------------------
!!- Print last scalar solution for restart
IDUMMY = -99
!!ISAVEFLUID = 1
if(SOLVE_SCALAR) call SAVE_SCALAR(IDUMMY)
!!ISAVEFLUID = 4



!!--------------------------------------------------------------------
!! 6.4. Particle velocities and orientations
!!--------------------------------------------------------------------
!!- Print final particle kinematics for restart
if (SOLVE_FLUID ==2)  then
 call FCM_SAVE_PARTICLE_CHARACTERISTICS
 
 IDUMMY = -99
 call FCM_SAVE_PARTICLE_KINEMATICS(IDUMMY)
 !if (FCM_ACTIVATE_STRESSLET>0) then
 ! call FCM_SAVE_PARTICLE_DYNAMICS(IDUMMY) 
 !end if
end if




if(MYID==0) then
write(*,*)
write(*,*) '**************************************'
write(*,*) '       END COMPUTATION FOR NCYCLE'
write(*,*) '**************************************'
write(*,*) '           Cycle =', NCYCLE
write(*,*) '      Ending time=', TIME
write(*,*) '**************************************'
write(*,*)
end if



!!- runtime CPU control
if(MYID==0) call CPUTIME_INFO(NCYCLE)



!- Synchronize all the process
call MPI_BARRIER(MPI_COMM_WORLD,IERR)



if(MYID==0) then
 close(UNIT_INFO(1))
 close(UNIT_INFO(2))
 close(UNIT_INFO(3))
 close(UNIT_INFO(4))
end if



!- Clean P3DFFT
 call P3DFFT_CLEAN

!- Free MPI environement
 call MPI_FINALIZE (ierr)


!!- runtime CPU control
if(MYID==0) write(*,*)' This is the END ... '

!!----------------------------------------------------------------------

10601 format (2x,A,  1x,E13.6 )
10602 format (2x,A,2(1x,E13.6))
10603 format (2x,A,3(3x,E13.6))


10604 format (6(A,E13.6,2x))

10000 format (30(e17.7))

!10700 format (2x,' Computation at ',f6.2,' %')
10700 format (2x,' Computation at ',f6.2,' %, Elapsed time:',f12.3,' s')
10701 format (2x,' Class ',i3,' has ',i7)
10702 format (2x,' C',i2.2,' has exchanged ',i5,' part. and occupied ',i6,' in one cpu')

end program ns3d

