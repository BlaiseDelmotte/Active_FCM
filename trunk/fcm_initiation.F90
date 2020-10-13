!!====================================================================
!!
!! 
!!> @brief
!!> Routine containing all the FCM initiations
!!
!! Date :  10/10/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION

!!====================================================================
!! 
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


use PARTICLE_PARALLEL

use MPI_STRUCTURES

use P3DFFT

implicit none

#define SIMPLE_SPRNG	
#define USE_MPI 
#include "sprng_f.h"

!- SPRNG Variables
integer :: SEED
SPRNG_POINTER :: MY_STREAM

! Initialize the streams for the random number generationwith SPRNG
SEED = make_sprng_seed()
MY_STREAM = init_sprng(SEED,SPRNG_DEFAULT,4)

call FCM_INITIATION_PARTICLE_GAUSSIAN 

! If more than 4 part: random position seeding (orientation can be prescribed)
if ((NPART_FULL.gt.0).and.((FCM_INIT_PART_POS==1).or.(FCM_INIT_PART_POS==4))) then

 ! Particle initiation must be called with only 1 proc because 
 ! random numbers are generated with local clock time
 if (MYID==0) then
  call FCM_INITIATION_PARTICLE_POS_ORIENT
 end if
 
 call MPI_BCAST(FCM_XP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_YP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_ZP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_XP_NOPER,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_YP_NOPER,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_ZP_NOPER,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_PSWIM,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
  if (FCM_USE_QUAT == 1) then
   call MPI_BCAST(FCM_QUAT,4*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR) 
  else
   call MPI_BCAST(FCM_P2,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
   call MPI_BCAST(FCM_P3,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)   
  end if   
 if (FCM_NELLIPSOID>0) then
  call MPI_BCAST(FCM_ROT_MAT,NPART_FULL*3*3,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 end if
 
else 
 call FCM_INITIATION_PARTICLE_POSITION
 call FCM_INITIATION_PARTICLE_ORIENTATION

end if



if (FCM_NSWIM(1)>0) then

! Compute Self induced tensors 
 call FCM_COMPUTE_SELF_ROS_TENSOR

 if (FCM_NSWIM(2)>0) then
  call FCM_COMPUTE_SELF_VEL_TENSOR
 end if

 ! Initialize phase shift if time-dependent stroke 
 if (FCM_NSWIM(3)>0) then
  if (MYID==0) then
   call FCM_INITIATION_PHASE_SHIFT
  end if   
  call MPI_BCAST(FCM_PHASE_SHIFT,FCM_NSWIM(3),MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 end if

 ! Initialize swimming velocities
 call FCM_INITIATION_PARTICLE_VELOCITY
end if

if (SOLVE_COLLISION>=3) then
 ! Set the bucket map for bucket sorting
 call FCM_BUCKET_DISTRIBUTION
 call FCM_BUCKET_ALLOCATIONS
 call FCM_BUCKET_BUILD_MAP
 
 ! Initialize Rep force magnitude
 call FCM_INITIATION_BARRIER_MAGNITUDE
 
 if (FCM_BC==2) then
  ! Initialize Rep wall force magnitude
  call FCM_INITIATION_REP_WALL_MAGNITUDE
 end if
 
end if !if (SOLVE_COLLISION>=3) then



!- Zeros necessary variables
call FCM_ZERO_PART_VARIABLE
call FCM_ZERO_FIELD_VARIABLE
FCM_SIJ = 0.0
if (FCM_NSWIM(1)>0) then
 FCM_SPIJ = 0.0
end if


if  (SOLVE_FLUID ==2)  then
 call FCM_SAVE_PARTICLE_CHARACTERISTICS
end if

if ((MYID==0).and.(SOLVE_FLUID ==2)) then
call FCM_WRITE_INFO
end if
 
end subroutine FCM_INITIATION
