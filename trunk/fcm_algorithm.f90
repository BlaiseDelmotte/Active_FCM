!!====================================================================
!!
!! 
!!> @brief
!!> Routine containing all the FCM Algorithm 
!!
!! Date :  22/09/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ALGORITHM(NCYCLE)

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



! Iteration number
integer, intent(in) :: NCYCLE

integer :: I
real(kind=8) :: seconds  


 !- Zeros necessary variables
 call FCM_ZERO_PART_VARIABLE
 call FCM_ZERO_FIELD_VARIABLE

 
 if (FCM_NELLIPSOID==0) then
 ! If only spheres, store gaussians
  call FCM_COMPUTE_GAUSSIAN_ENV  
  
  if (FCM_NSWIM(2)>0) then
   call FCM_COMPUTE_GAUSSIAN_ENV_SQ
  end if   
 else  
 ! If only ellipsoid, only store gaussian support
  call FCM_COMPUTE_GAUSSIAN_SUPPORT  
 end if
 
 if (SOLVE_COLLISION >= 3) then
  !- Establish linked-list for nf and collision interactions
  call FCM_COMPUTE_BUCKET_LISTS
  
  !- Use forces/torques for collisions
  if (SOLVE_COLLISION == 3) then
  
   if (FCM_NELLIPSOID>0) then
    !- GB for ellipsoids
    call FCM_GB_SHIFTED_UNIAXIAL
   else 
   
    !- isotropic barrier for spheres
    call FCM_BARRIER_FORCE
!~ 
    !- Add Repulsive wall force if walls along x
    if (FCM_BC==2) then
     call FCM_REPULSIVE_WALL
    end if
    
   end if !if (FCM_NELLIPSOID>0) then
  end if !if (SOLVE_COLLISION == 3) then
 end if !if (SOLVE_COLLISION >= 3) then
 
 if (maxval(abs(FCM_EXT_FORCE))>0) then
  call FCM_ADD_MONOPOLE_FORCING
 end if
 

 if (maxval(abs(FCM_EXT_TORQUE))>0) then
  call FCM_ADD_ROTLET_FORCING
 end if

 if (FCM_NSWIM(1)>0) then
 
  if (FCM_NSWIM(3)>0) then
   ! Add time dependent stroke 
   call FCM_TIME_STROKE(NCYCLE)

  end if 
  
  call FCM_ADD_SWIMMING_FORCING
  
 end if 


 if (maxval(abs(FCM_FORCE))>0) then 
  call FCM_DISTRIB_MONOPOLE
 end if

 if (maxval(abs(FCM_TORQUE))>0) then 
  call FCM_DISTRIB_ROTLET
 end if

 if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2)==0)) then 
  call FCM_DISTRIB_STRESSLET(FCM_NSWIM(1),FCM_SPIJ)
 else if (FCM_NSWIM(2)>0) then
  call FCM_DISTRIB_SQUIRMING_FORCING
 end if

 if (FCM_ACTIVATE_STRESSLET>0) then
   
  FCM_EIJ = 0.0  
  
  if (NCYCLE==1) then
   FCM_SIJ = 0.0
  end if

  call FCM_DISTRIB_STRESSLET(FCM_ACTIVATE_STRESSLET,FCM_SIJ)
   
  ! Take shear into account for initiation (cf. Yeo)
  FCM_CONSIDER_ROS_SHEAR = 1
     
  !- If mirror x-wall
  if (FCM_BC==2) then
   call FCM_MIRROR_FORCES_X
  end if
  seconds = MPI_Wtime ( );
  call FCM_FLUID_PREDICTION  
  seconds = MPI_Wtime ( )-seconds;
  if (MYID.eq.0) then
   print*, 'time_fluid_prediction = ', seconds
  end if
  call FCM_RATE_OF_STRAIN_FILTER  
  
      
  if (FCM_NSWIM(1)>0) then
   call FCM_REMOVE_SELF_ROS
  end if
  
    
  FCM_FORCING_X(:,:,:) = 0.0
  FCM_FORCING_Y(:,:,:) = 0.0
  FCM_FORCING_Z(:,:,:) = 0.0

  ! Don't take shear into account in the PCG
  FCM_CONSIDER_ROS_SHEAR = 0

  seconds = MPI_Wtime ( );
  call FCM_CONJUGATE_GRADIENT_STRESSLET
  seconds = MPI_Wtime ( )-seconds;
  if (MYID.eq.0) then
   print*, 'time_stresslet = ', seconds
  end if

  FCM_FORCING_X(:,:,:) = 0.0
  FCM_FORCING_Y(:,:,:) = 0.0
  FCM_FORCING_Z(:,:,:) = 0.0
  
  call FCM_DISTRIB_STRESSLET(FCM_ACTIVATE_STRESSLET,FCM_SIJ)
  
  if (maxval(abs(FCM_FORCE))>0) then 
   call FCM_DISTRIB_MONOPOLE
  end if
  
  if (maxval(abs(FCM_TORQUE))>0) then
   call FCM_DISTRIB_ROTLET
  end if
 
  if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2)==0)) then 
   call FCM_DISTRIB_STRESSLET(FCM_NSWIM(1),FCM_SPIJ)
  else if (FCM_NSWIM(2)>0) then
   call FCM_DISTRIB_SQUIRMING_FORCING
  end if

 end if  !if (FCM_ACTIVATE_STRESSLET==1) then


  
 !- If mirror x-wall
 if (FCM_BC==2) then
  call FCM_MIRROR_FORCES_X
 end if

 call FCM_FLUID_PREDICTION  
   
 call FCM_VELOCITY_FILTER

 
 if (FCM_NSWIM(2)>0) then
  call FCM_REMOVE_SELF_VEL   
 end if
  
 call FCM_ROTATION_FILTER
 if (mod(NCYCLE,FOUT1)==1) then
  call SAVE_FLUID(NCYCLE)
 end if 

 if (mod(NCYCLE,FOUT3)==1) then
  call FCM_SAVE_PARTICLE_KINEMATICS(NCYCLE) 
!    if (FCM_ACTIVATE_STRESSLET==1) then
!     call FCM_SAVE_PARTICLE_DYNAMICS(NCYCLE)
!    end if
 end if 
 
 call FCM_ADV_PARTICLE_POSITION
 
 if (FCM_USE_QUAT == 1) then
  call FCM_ADV_PARTICLE_ORIENTATION
 else
  call FCM_ADV_PARTICLE_ORIENTATION_PSWIM_ONLY(DTIME)
 end if

 call FCM_CHECK_PARTICLE_POSITION

end subroutine FCM_ALGORITHM
