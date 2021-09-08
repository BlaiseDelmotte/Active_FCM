 !!====================================================================
!!
!! 
!!
!! DESCRIPTION: 
!!
!!
!
!!> @brief
!!> Compute ERIC second step:
!
!!  CAUTION: ONLY VALID FOR FORWARD EULER
!! Date :  21/10/2014
!!   1) Advect particles with previous fluctuating velocity (with dt/2)
!!   2) Distribute fluctutating stress (the same as previous step) + Monopoles at current positions
!!   3) Provides a new velocity which will be used to advect particles from previous step
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_KS_SECOND_STEP(NCYCLE)

!!====================================================================
!!
!!====================================================================

use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use P3DFFT
use PARAM_PHYS         !- Physical & numerical parameters

implicit none

integer, intent(in):: NCYCLE 
integer :: IP, IC


FCM_XP0 = FCM_XP
FCM_YP0 = FCM_YP
FCM_ZP0 = FCM_ZP

! X(k+1/2)
FCM_XP = FCM_XP + DTIME/2.0*FCM_UP(:,1)
FCM_YP = FCM_YP + DTIME/2.0*FCM_VP(:,1)
FCM_ZP = FCM_ZP + DTIME/2.0*FCM_WP(:,1)

FCM_XP = FCM_XP - LXMAX*real(floor(FCM_XP/LXMAX))
FCM_YP = FCM_YP - LYMAX*real(floor(FCM_YP/LYMAX))
FCM_ZP = FCM_ZP - LZMAX*real(floor(FCM_ZP/LZMAX))

FCM_FORCING_X(:,:,:) = 0.0
FCM_FORCING_Y(:,:,:) = 0.0
FCM_FORCING_Z(:,:,:) = 0.0

FCM_FORCE = 0.0
FCM_TORQUE = 0.0

FCM_UP = 0.0
FCM_VP = 0.0
FCM_WP = 0.0

FCM_OMPX = 0.0
FCM_OMPY = 0.0
FCM_OMPZ = 0.0

! If only spheres, store gaussians
call FCM_COMPUTE_GAUSSIAN_ENV  



if (FCM_NSWIM(2)>0) then
 call FCM_COMPUTE_GAUSSIAN_ENV_SQ
end if   


if (SOLVE_COLLISION > 0) then
 !- Establish linked-list for nf and collision interactions
 call FCM_COMPUTE_BUCKET_LISTS
 !- isotropic barrier for spheres
 call FCM_BARRIER_FORCE
 !- Add Repulsive wall force if walls along x
 if (FCM_BC==2) then
  call FCM_REPULSIVE_WALL
 end if
end if !if (SOLVE_COLLISION > 0) then


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
 FCM_SIJ = 0.0

 ! Take shear into account for initiation (cf. Yeo)
 FCM_CONSIDER_ROS_SHEAR = 1

 !- If mirror x-wall
 if (FCM_BC==2) then
  call FCM_MIRROR_FORCES_X
 end if

 call FCM_FLUID_PREDICTION_BROWNIAN

 call FCM_RATE_OF_STRAIN_FILTER

 if (FCM_NSWIM(1)>0) then
  call FCM_REMOVE_SELF_ROS
 end if

 FCM_FORCING_X(:,:,:) = 0.0
 FCM_FORCING_Y(:,:,:) = 0.0
 FCM_FORCING_Z(:,:,:) = 0.0

 ! Don't take shear into account in the PCG
 FCM_CONSIDER_ROS_SHEAR = 0

 call FCM_CONJUGATE_GRADIENT_STRESSLET

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

call FCM_FLUID_PREDICTION_BROWNIAN


!  if (FCM_ACTIVATE_STRESSLET==1) then
!    FCM_CONSIDER_ROS_SHEAR = 1
!  call FCM_RATE_OF_STRAIN_FILTER  
!    if (FCM_SWIMMING>=1) then
!      call FCM_REMOVE_SELF_ROS
!    end if
! 
!    if (MYID==0) then
!         print*,'FCM_SIJ = ', FCM_SIJ
!         print*,'FCM_EIJ = ', FCM_EIJ
!        read(*,*)
!    end if
! end if !if (FCM_ACTIVATE_STRESSLET==1) then


call FCM_VELOCITY_FILTER


if (FCM_NSWIM(2)>0) then
 call FCM_REMOVE_SELF_VEL   
end if


call FCM_ROTATION_FILTER


! Put the particles back where they were for advection and saves
 FCM_XP = FCM_XP0
 FCM_YP = FCM_YP0
 FCM_ZP = FCM_ZP0
 
end subroutine FCM_KS_SECOND_STEP
