!!====================================================================
!!
!! 
!!> @brief
!!> Routine reading the fcm parameters from the ASCII file  'fcm_param.in'
!!
!! Date :  14/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_READPARAM

!!====================================================================
!! 
!!====================================================================

use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use DNS_DIM

implicit none


! Process information variables
integer :: IERROR

!- Open file
open(unit=200,file='fcm_param.in',status='old', iostat=IERROR)

if(IERROR /= 0) then

 print *, "CPU -- ", MYID, ":: ERROR: Input data file fcm_param.in open error!"
 call MPI_FINALIZE(IERROR)
 call abort()
 
else

 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line

!!====================================================================
!! 1. INITIATION PARAMETERS
!!====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_INIT_PART_POS    !>  Particle position initiation mode
 read(200,*) ! comment line
 read(200,*) FCM_INIT_PART_ORIENT !>  Particle orientation initiation mode
 
!!====================================================================
!! 2. GAUSSIAN PARAMETERS
!!====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line

 read(200,*) FCM_NGD_REF             !> REf Gaussian support
 
 read(200,*) ! comment line

 read(200,*) FCM_SIGP            !>  Gaussian width
  
!!====================================================================
!! 3. PARTICLE PROPERTIES
!!====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line

 read(200,*) FCM_NSPHERE, FCM_NSPHERE_1       !>  Number of spheres, number of spheres of type 1 (for bidispersity)
 
 read(200,*) ! comment line
 
 read(200,*) FCM_SPHERE_SIZE     !>  Default sphere size
 
 read(200,*) ! comment line

 read(200,*) FCM_NELLIPSOID      !>  Number of ellipsoids
 
 read(200,*) ! comment line
 
 read(200,*) FCM_ELLIPSOID_SIZE  !>  Default ellipsoid size
 
!!====================================================================
!! 4. FORCING PARAMETERS
!!====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 
 read(200,*) FCM_EXT_FORCE      !>  External unique force
 
 read(200,*) ! comment line
 
 read(200,*) FCM_EXT_TORQUE     !>  External unique torque
 
 !!====================================================================
!! 5. STRESSLET PARAMETERS
!!====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 
 read(200,*) FCM_ACTIVATE_STRESSLET  !> rigidity stresslet. 0: not activated, 1: activated
 
 read(200,*) ! comment line
 
 read(200,*) FCM_TOL_L2RES_ROS       !>  Tolerance on the L2 norm of the particle ROS
!====================================================================
! 6. FLOW PARAMETERS
!=====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_SHEAR                 !> shear ux = shear*z
 !====================================================================
! 7. SWIMMING PARAMETERS
!=====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_SWIMMING               !> activate swimming stresslet.0: not activated, 1: activated, 2: squirmer 3: time-dep sync, 4: time-dep dphi=pi/2, 5: time-dep dphi=rand
 read(200,*) ! comment line
 read(200,*) FCM_VSW_REF                    !< Ref Intrinsic swimming velocity
 read(200,*) ! comment line
 read(200,*) FCM_BETA                   !< Ratio between swim vel and swim stresslet magnitude
 !====================================================================
! 8. CONTACT HANDLING PARAMETERS
!=====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_FBLEVEL                    !< Repulsive Force Magnitude 
 read(200,*) ! comment line
 read(200,*) FCM_FBRANGE                    !< Repulsive Force Magnitude 
!====================================================================
! 9. BOUNDARY CONDITION AND FLOW PARAMETERS
!=====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_FLOW_TYPE             !< Chose if Brownian motion(=2) or not(=1)
 read(200,*) ! comment line
 read(200,*) FCM_BC                    !< BC=1: 3D Periodic, BC=2: stress free surface along x
 read(200,*) ! comment line
 read(200,*) KBT                       !< Thermal Energy
!====================================================================
! 10. RUN AND TUMBLE PARAMETERS
!=====================================================================
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) ! comment line
 read(200,*) FCM_RUNTUMBLE             !< 0: NO RUN AND TUMBLE, 1: YES
 read(200,*) ! comment line
 read(200,*) FCM_TAU_RUN                    !< MEAN  RUN TIME IN SECOND


 close(200)

 
 ! Exceptionally allocate NGD here 
 allocate(FCM_NGD(NPART_FULL))
 allocate(FCM_NGDH(NPART_FULL))
 
 if (FCM_NELLIPSOID>0) then
  !- Ellipsoid Gaussian support, required for allocation
  FCM_ELLIPSOID_NGD = maxval(FCM_ELLIPSOID_SIZE) * FCM_NGD_REF
 else
  FCM_NGD(1:FCM_NSPHERE_1) = nint(FCM_SPHERE_SIZE(1) * real(FCM_NGD_REF))
  FCM_NGD(FCM_NSPHERE_1+1:NPART_FULL) = nint(FCM_SPHERE_SIZE(2) * real(FCM_NGD_REF))
  FCM_NGD_MAX = maxval(FCM_NGD)
  
 ! print*,'FCM_NGD = ',  FCM_NGD
 ! print*,'FCM_NGD_MAX = ', FCM_NGD_MAX 

 end if
 
 FCM_NSWIM = 0
 if ( (FCM_SWIMMING(1).ge.1).and.(FCM_SWIMMING(2).ge.1) ) then
  FCM_NSWIM(1) = FCM_NSPHERE
 else if ( (FCM_SWIMMING(1).ge.1).and.(FCM_SWIMMING(2).eq.0) )then
  FCM_NSWIM(1) = FCM_NSPHERE_1
 end if
 
 if ( (FCM_SWIMMING(1).ge.2).and.(FCM_SWIMMING(2).ge.2) ) then
  FCM_NSWIM(2) = FCM_NSPHERE
 else if ( (FCM_SWIMMING(1).ge.2).and.(FCM_SWIMMING(2)<2) )then
  FCM_NSWIM(2) = FCM_NSPHERE_1
 end if
 
 if ( (FCM_SWIMMING(1).ge.3).and.(FCM_SWIMMING(2).ge.3) ) then
  FCM_NSWIM(3) = FCM_NSPHERE
 else if ( (FCM_SWIMMING(1).ge.3).and.(FCM_SWIMMING(2)<3) )then
  FCM_NSWIM(3) = FCM_NSPHERE_1
 end if
 
 if ( (FCM_SWIMMING(1).ge.4).and.(FCM_SWIMMING(2).ge.4) ) then
  FCM_NSWIM(4) = FCM_NSPHERE
 else if ( (FCM_SWIMMING(1).ge.4).and.(FCM_SWIMMING(2)<4) )then
  FCM_NSWIM(4) = FCM_NSPHERE_1
 end if
 
 if ( (FCM_SWIMMING(1).ge.5).and.(FCM_SWIMMING(2).ge.5) ) then
  FCM_NSWIM(5) = FCM_NSPHERE
 else if ( (FCM_SWIMMING(1).ge.5).and.(FCM_SWIMMING(2)<5) )then
  FCM_NSWIM(5) = FCM_NSPHERE_1
 end if 
 
 
!~  read(*,*)
 
 FCM_USE_QUAT = 0
 
 

 if(MYID==0) write(*,*) 'FCM parameters reading --> OK'

end if

end subroutine FCM_READPARAM
