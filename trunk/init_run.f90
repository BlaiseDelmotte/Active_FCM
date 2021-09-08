!!====================================================================
!!
!! This routine defines the constants for the simulation
!!
!!====================================================================

subroutine INIT_RUN

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use FORCING
use GEOMETRIC_VARIABLE
use ENSIGHT_VAR
use STATISTICS
use RANDOM_VAR


implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

integer :: I, N
!---------------------------------------------------------------------

!!====================================================================
!! 1. Basics
!!====================================================================


!!- P. Fede - 12/06/2013
!! Spatial correlation disconnected
!!LEVEL4_STFLU = .false.
!!- P. Fede - 12/06/2013




!! In case of FCM (i.e. Stokes flow) scalar and DPS are switch off
if(SOLVE_FLUID == 2) then
 SOLVE_SCALAR = .false.

!!- no forcing because no turbulence !
 STEADY = .false.
end if


if(SOLVE_FLUID == 0) then

 SOLVE_SCALAR = .false.

!!- No-fluid -> no-stat
 LEVEL0_STFLU = .false.
 LEVEL1_STFLU = .false.
 LEVEL2_STFLU = .false.
 LEVEL3_STFLU = .false.
 LEVEL4_STFLU = .false.

!!- no forcing because no turbulence !
 STEADY = .false.
end if




!!- Ghost cells
NGHTCELL = 0
if(SOLVE_FLUID>0) then
 NGHTCELL = 2
end if



!- Pointers of time-integration scheme
TN = 1
TNM1 = 2
TNM2 = 3

!- First time step treated as Euler scheme
FLAG_TSCHEME = 1


!- File containing all informations
if(MYID==0) then
 UNIT_INFO(1) = 151
 UNIT_INFO(2) = 152
 UNIT_INFO(3) = 153
 UNIT_INFO(4) = 154
 UNIT_INFO(5) = 155

 !- Opening file
 open(unit=UNIT_INFO(1), file='run.info', action='write', status='replace')

 !- Opening file
 open(unit=UNIT_INFO(2), file='stat.info', action='write', status='replace')

 !- Opening file
 open(unit=UNIT_INFO(3), file='cpu.info', action='write', status='replace')



 !- Opening file
 open(unit=UNIT_INFO(5), file='numerics.info', action='write', status='replace')


end if



!!- Method for saving restart file
!!  This variable is filled in subroutine INIT_RUN
!! 
!! ISAVEFLUID = 1: Multiple binary files
!!            = 2: Using of MPI I/O
ISAVEFLUID = 2

!! 

!! ISAVEFORCE= 1: Multiple binary files
!!           = 2: Direct access file
ISAVEFORCE = 2


!!====================================================================
!! 2. Physics
!!====================================================================
PPI = 4.*atan(1.D0) !- Pi
TWOPI = 2.*PPI      !- 2.*Pi

ICMPL = CMPLX(0.D0,1.D0) !!- complex number

KFIRST = TWOPI/LXMAX   !!- First wavenumber
KLAST = PPI/(LXMAX/NX) !!- Last wavenumber


!!- Forced wavenumber
if(STEADY) then
 KFORCE_MIN = KFORCE_MIN*KFIRST
 KFORCE_MAX = KFORCE_MAX*KFIRST
end if



!!- Seed for random number generator
!! Must be the same for all process for the forcing
IDFORCE = -2
ISET = 0
IY = 0
IV(:) = 0



FILTERING = .false.
KCUT = 2*KLAST
if(FILTERING) then
 KCUT = 8*KFIRST
 if(MYID==0) write(*,*) ' Filtering activated'
 if(MYID==0) write(*,*) '     + Kcut/K0 =',KCUT/KFIRST
end if



!!====================================================================
!! 3. Outputing
!!====================================================================
ENSIGHT_BIN = .false.
ENSIGHT_DIM = 3
SLCNUM = NZ/2


NFILEOUT = 0

write(FILE_EXT,10205)'.p',MYID




if(ENSIGHT_OUT > 0) then
CASEFLU = 'fluid'
NBNODE = 2

allocate(NODELIST(NBNODE))
allocate(NODETYPE(NBNODE))

end if



!!====================================================================
!! 4. Statistics
!!====================================================================

if(.not.SOLVE_SCALAR)  then
 LEVEL0_STSCL = .false.
 LEVEL1_STSCL = .false.
end if


if(LEVEL1_STFLU.or.LEVEL2_STFLU.or.LEVEL3_STFLU.or.LEVEL4_STFLU) then
 LEVEL0_STFLU=.true.
else
 LEVEL0_STFLU=.false.
end if

if(LEVEL1_STSCL.or.LEVEL2_STSCL) then
 LEVEL0_STSCL=.true.
else
 LEVEL0_STSCL=.false.
end if


if(SOLVE_FLUID==0) then
 LEVEL0_STFLU = .false.
 LEVEL0_STSCL = .false.
 
end if


if(LEVEL0_STFLU.or.LEVEL0_STSCL) then
 NSTAT = 100
else
 NSTAT = 1
end if



if(MYID==0 .and. DEBUG) then
 write(*,*)
 write(*,*)'-------------------'
 write(*,*)'--> DEBUG MODE <---'
 write(*,*)'-------------------'
 write(*,*)
end if



!!--------------------------------------------------------------------
10201 format(A,I1)
10202 format(A,I2)
10203 format(A,I3)
10204 format(A,I4)

10205 format(A,I4.4)

10600 format (A,I3,A,I4,A,E13.6)


end subroutine INIT_RUN
