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

subroutine FCM_ALGORITHM_BROWNIAN(NCYCLE)

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
  

!!~  if (MYID==0) then
!~ if (NCYCLE==1) then
!~  FCM_XP(1) = LXMAX/2.0 - FCM_SPHERE_RADP(1) 
! FCM_XP(1) =  FCM_SPHERE_RADP(1) 
!~ !!~   FCM_XP(1) = 01.57
!~ 
!~  if( NPART_FULL==2) then
!~   FCM_XP(2) = LXMAX/2.0 + FCM_SPHERE_RADP(1)
!~   
!~   FCM_YP(2) = FCM_YP(1)
!~   FCM_ZP(2) = FCM_ZP(1)
!~  end if
 
!~ end if
!~ if (MYID==0) then
!~    print*,NCYCLE
!~    print*,FCM_XP/FCM_SPHERE_RADP(2)
!~    print*,FCM_YP/FCM_SPHERE_RADP(2)
!~    print*,FCM_ZP/FCM_SPHERE_RADP(2)
!~    print*,'FCM_PSWIM(1,:) = ',FCM_PSWIM(1,:)
!~    print*,'FCM_PSWIM(2,:) = ',FCM_PSWIM(2,:)   
!~    print*,'distance /(FCM_SPHERE_RADP(1)+FCM_SPHERE_RADP(2)) = '
!~    print*,dsqrt( (FCM_XP(2)-FCM_XP(1))**2 +(FCM_YP(2)-FCM_YP(1))**2 + (FCM_ZP(2)-FCM_ZP(1))**2 ) &
!~                /(FCM_SPHERE_RADP(1)+FCM_SPHERE_RADP(2))
!~  end if  
!~     read(*,*)
!~  end if
!~    print*,'FCM_YP/FCM_SPHERE_RADP(1) = ', FCM_YP/FCM_SPHERE_RADP(1)
!~    print*,'FCM_ZP/FCM_SPHERE_RADP(1) = ', FCM_ZP/FCM_SPHERE_RADP(1)
!~    read(*,*)
!~    print*,'FCM_PSWIM = ', FCM_PSWIM
!~ !  print*,'FCM_VSW = ', FCM_VSW
!~   read(*,*)


 !- Zeros necessary variables
 call FCM_ZERO_PART_VARIABLE
 call FCM_ZERO_FIELD_VARIABLE


! If only spheres, store gaussians
 call FCM_COMPUTE_GAUSSIAN_ENV  


 if (FCM_BC==1) then
  call FCM_RS_SETUP
 else
  call FCM_RS_SETUP_SYMMETRIZE
 end if

 call FCM_FLUID_PREDICTION_BROWNIAN

 


 ! if (FCM_ACTIVATE_STRESSLET==1) then
 !  FCM_CONSIDER_ROS_SHEAR = 1
 !  call FCM_RATE_OF_STRAIN_FILTER  
 !  if (FCM_SWIMMING>=1) then
 !    call FCM_REMOVE_SELF_ROS
 !  end if

 !  if (MYID==0) then
 !       print*,'maxval(FCM_EIJ)', maxval(abs(FCM_EIJ))
  !     read(*,*)
 !  end if
 ! end if !if (FCM_ACTIVATE_STRESSLET==1) then


   
 call FCM_VELOCITY_FILTER
 
 if (FCM_FLOW_TYPE == 2) then
  !! - SECOND STEP OF RFD SCHEME
  call FCM_KS_SECOND_STEP(NCYCLE) 
 end if ! if (FCM_FLOW_TYPE == 2) then

 
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

  
!!!!!!!!!! TO COMMENT !!!!!!!!!!!!!!!!!!!!
 !! JUST A TEST TO CHECK THAT MIRROR IMAGE HAS SAME VELOCITY
!~  FCM_XP(2) = LXMAX - FCM_XP(1)
!~  FCM_YP(2) = FCM_YP(1)
!~  FCM_ZP(2) = FCM_ZP(1)

!~  FCM_XP(1) = FCM_XP(1) - (LXMAX/2.0 - 2.0*FCM_SPHERE_RADP(1))/real(NCYCLEMAX - 2)
!~  if( NPART_FULL==2) then
!~   FCM_XP(2) = LXMAX - FCM_XP(1)
!~  end if
!!!!!!!!!! TO COMMENT !!!!!!!!!!!!!!!!!!!!

!~ if (MYID==0) then
!~    print*,'FCM_FORCE = ', FCM_FORCE
!~    print*,'FCM_VSW = ', FCM_VSW
!~    print*,'FCM_VSW/rad = ', FCM_VSW/FCM_SPHERE_RADP
!~   print*,'FCM_SIJ= ', FCM_SIJ
!~    print*,'FCM_VEL= '
!~   print*,FCM_UP(:,1)
!~   print*,FCM_VP(:,1)
!~   print*,FCM_WP(:,1)
!~    print*,'FCM_ROT= '
!~   print*,FCM_OMPX
!~   print*,FCM_OMPY
!~   print*,FCM_OMPZ
!~  end if

!~   read(*,*)
!~   read(*,*) !/(6.0*PPI*FCM_SPHERE_RADP(1))
!~   print*,FCM_VP(2,1) !/(6.0*PPI*FCM_SPHERE_RADP(1))
!~   print*,FCM_VP(1,1)
!~   print*,FCM_WP(1,1)

!~   read(*,*)

!~   print*,FCM_WP(:,1)
!~   read(*,*)
!   print*,'maxval(FCM_VP(:,1))= ',maxval(FCM_VP(:,1))
!   print*,'maxval(FCM_WP(:,1)) = ',maxval(FCM_WP(:,1))
!~     read(*,*)
! end if 
!~  if (MYID==0) then
!~     print*, 'FCM_XP = ',FCM_XP
!~     read(*,*)
!~  end if
  
!~  if (MYID==0) then
!~   print*,'NCYCLE = ', NCYCLE
!~   print*,'FCM_TORQUE = ', FCM_TORQUE
!~   print*,'FCM_OMPY = ', FCM_OMPY
!~   print*,'FCM_PSWIM = ', FCM_PSWIM 
!~   read(*,*)
!~  end if
!~   if ( (maxval(abs(FCM_TORQUE))>0).or.(maxval(abs(FCM_FORCE))>0) ) then
!~    print*,'FCM_TORQUE = ', FCM_TORQUE
!~   end if
!~   
!~   if ((maxval(FCM_UP(:,1))/FCM_VSW>10.0).or.(maxval(FCM_VP(:,1))/FCM_VSW>10.0).or.(maxval(FCM_WP(:,1))/FCM_VSW>10.0)) then

!~   end if
       
!~  end if
!~  
!~  if (MYID==0) then
  !if (mod(NCYCLE,FOUT3)==0) then
!~     print*,'NCYCLE ', NCYCLE
!~    print*,'FCM_XP = ', FCM_XP
!~    print*,'FCM_YP = ', FCM_YP
!~    print*,'FCM_ZP = ', FCM_ZP
 !  print*, FCM_WP(:,1)
!~    print*,'FCM_VP(:,1) = ', FCM_VP(:,1)-FCM_VSW*FCM_PSWIM(:,2)
!~    print*,'FCM_WP(:,1) = ', FCM_WP(:,1)-FCM_VSW*FCM_PSWIM(:,3)
!~    read(*,*)
!~    print*,'FCM_PSWIM = ', FCM_PSWIM   
!~    read(*,*)
!~    print*,'FCM_OMPX = ', FCM_OMPX
!~    print*,'FCM_OMPY = ', FCM_OMPY
!~    print*,'FCM_OMPZ = ', FCM_OMPZ
!~    read(*,*)
!~    print*,'FCM_SIJ = ', FCM_SIJ
!~    print*,'FCM_EIJ = ', FCM_EIJ
!~    read(*,*)
  !end if
!~  end if
!~ 

end subroutine FCM_ALGORITHM_BROWNIAN
