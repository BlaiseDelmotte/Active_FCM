
!!=====================================================================
!!
!! FCM Forcing Variables
!!
!! Description :
!!> @brief
!!> Module containing all forcing terms of FCM particles
!!
!!
!! Date :  15/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!
!!=====================================================================
module FCM_FORCING_VARIABLE


implicit none

!!====================================================================
!! 0.0. FLAGS FOR FLOW AND BC
!!====================================================================
integer :: FCM_FLOW_TYPE                 !< Chose if Brownian motion(=2) or not(=1)
integer :: FCM_BC                        !< BC=1: 3D Periodic, BC=2: stress free surface along x
real(kind=8) :: KBT                      !< Thermal energy
!!====================================================================
!! 0. FORCE FIELD 
!!====================================================================
!! Forces   
real(kind=8), dimension(:,:,:), allocatable :: FCM_FORCING_X
real(kind=8), dimension(:,:,:), allocatable :: FCM_FORCING_Y
real(kind=8), dimension(:,:,:), allocatable :: FCM_FORCING_Z

!! Random Stresses if Brownian Motion in Real space
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_XX
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_YY
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_ZZ
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_XY
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_XZ
real(kind=8), dimension(:,:,:), allocatable :: FCM_RS_YZ

!! Random Stresses if Brownian Motion in Fourier space
double complex, dimension(:,:,:), allocatable :: FCM_RS_XX_FOU
double complex, dimension(:,:,:), allocatable :: FCM_RS_YY_FOU
double complex, dimension(:,:,:), allocatable :: FCM_RS_ZZ_FOU
double complex, dimension(:,:,:), allocatable :: FCM_RS_XY_FOU
double complex, dimension(:,:,:), allocatable :: FCM_RS_XZ_FOU
double complex, dimension(:,:,:), allocatable :: FCM_RS_YZ_FOU


!!====================================================================
!! 1. MONOPOLE PARAMETERS
!!====================================================================

integer :: FCM_ACTIVATE_FLUB                                     !< activate Flub. 0: not activated, 1: activated
integer :: FCM_PLACE_FLUB                                        !< place Flub. 0: before swimming dipole, 1: after swimming dipole

real(kind=8) :: FCM_LUB_THR                                      !< Threshold for lubrication activation

real(kind=8) :: FCM_FBRANGE                                      !<  K. Yeo repulsive force range
real(kind=8) :: FCM_FBLEVEL                                      !<  K. Yeo repulsive force intensity

real(kind=8) :: FCM_EPS_RFD                                      !< RFD for Brownian motion

real(kind=8) :: FCM_WALL_RANGE                                   !< Range of repulsive wall force
real(kind=8) :: FCM_WALL_LEVEL                                   !< Magnitude of repulsive wall force

real(kind=8), dimension(3) :: FCM_EXT_FORCE                      !< External force monopole, same value prescribed for each part

real(kind=8), dimension(:,:),allocatable :: FCM_FORCE_TEMP       !<  Temporary Force monopole for coll barrier
real(kind=8), dimension(:,:),allocatable :: FCM_FORCE_TEMP2      !<  Temporary Force monopole for wall barrier


real(kind=8), dimension(:,:),allocatable :: FCM_FORCE            !< Force monopole

real(kind=8), dimension(:,:),allocatable :: FCM_RFD_RAND_FORCE   !< Random force for  RFD integration


!!====================================================================
!! 2. TORQUE PARAMETERS
!!====================================================================

real(kind=8), dimension(3) :: FCM_EXT_TORQUE             !< External torque, same value prescribed for each part

real(kind=8), dimension(:,:),allocatable :: FCM_TORQUE   !< Force monopole

real(kind=8), dimension(:,:),allocatable :: FCM_AIJ      !< Rotlet


!!====================================================================
!! 3. STRESSLET PARAMETERS
!!====================================================================

integer :: FCM_ACTIVATE_STRESSLET                        !< activate rigidity stresslet. 0: not activated, N: N particles stresslet activated

real(kind=8), dimension(:,:),allocatable :: FCM_SIJ      !< Stresslet


!!====================================================================
!! 4. PARTICLE RATE OF STRAIN PCG MINIMIZATION PARAMETERS
!!====================================================================

integer :: FCM_MAX_ITER                                  !< Max number of iteration for minimization

integer :: FCM_CONSIDER_ROS_SHEAR                        !< Flag to chose or not to integrate the shear rate in the rate of strain computation


real(kind=8) :: FCM_INC_DIR_ALPHA                        !< Increment in the direction DIR_MIN for Stresslet coeffs
real(kind=8) :: FCM_INC_DIR_BETA                         !< Increment in the direction DIR_MIN to get the new DIR_MIN

real(kind=8) :: FCM_TOL_L2RES_ROS                 	     !< Tolerance on the L2 norm of the particle ROS

real(kind=8) :: FCM_L2RES_ROS                            !< L2 norm of the sum of the particle ROS

real(kind=8), dimension(:,:),allocatable :: FCM_DIR_MIN  !< Increment and direction of minimization 

real(kind=8), dimension(:,:),allocatable :: FCM_RES_ROS   !< Residual Rate of strain applied to particles


!!====================================================================
!! 5. FLOW FORCING
!!====================================================================

real(kind=8) :: FCM_SHEAR                                 !< Shear ux = shear*z


!!====================================================================
!! 6. SWIMMING PARAMETERS
!!====================================================================
integer :: FCM_RUNTUMBLE                   !< 0: NO RUN AND TUMBLE, 1: YES          
  
integer, dimension(5) :: FCM_NSWIM                      !< count type of swimmers. 1: stresslets, 2: squirmer,>3: time-dependent

integer, dimension(2) :: FCM_SWIMMING                                  !< activate swimming stresslet. 0: not activated, 1: activated only dipole, 2: squirming


real(kind=8) :: FCM_VSW_REF                                 !<  Reference value for intrinsic swimming velocity

real(kind=8) :: FCM_BETA                                 !< Ratio FCM_S0/FCM_VSW

real(kind=8) :: FCM_S0                                   !< Swimming stresslet intensity

real(kind=8) :: FCM_H0                                   !< Squirming potential dipole intensity

real(kind=8) :: FCM_TAU_RUN                    !< MEAN  RUN TIME IN SECOND


real(kind=8), dimension(:),allocatable  :: FCM_VSW  !< Swimming velocity depend on the swimmer

real(kind=8), dimension(:),allocatable  :: FCM_VSW_TIME  !< Time dependent swimming velocity depend on the swimmer

real(kind=8), dimension(:),allocatable  :: FCM_B2_TIME  !< Time dependent stresslet magnitude

real(kind=8), dimension(:),allocatable  :: FCM_PHASE_SHIFT !< Time dependent phase shift

real(kind=8), dimension(:,:),allocatable :: FCM_HI     !< Swimming potential dipole

real(kind=8), dimension(:,:),allocatable :: FCM_SPIJ     !< Swimming stresslet

real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_P  !< Integral of tensor Pij to get self induced VEL because of squirming quadrupole

real(kind=8), dimension(:),allocatable :: FCM_A1111  !< Component of self ROS tensor for ellipsoid
real(kind=8), dimension(:),allocatable :: FCM_A2222  !< Component of self ROS tensor for ellipsoid
real(kind=8), dimension(:),allocatable :: FCM_A1313  !< Component of self ROS tensor for ellipsoid
real(kind=8), dimension(:),allocatable :: FCM_A2323  !< Component of self ROS tensor for ellipsoid

real(kind=8), dimension(:),allocatable :: FCM_P11 !< Component of self vel tensor for ellipsoid
real(kind=8), dimension(:),allocatable :: FCM_P22  !< Component of self vel tensor for ellipsoid

real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R11  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet
real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R12  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet
real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R13  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet
real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R22  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet
real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R23  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet
real(kind=8), dimension(:,:,:),allocatable :: FCM_INT_R33  !< Integral of tensor Rijk to get self induced ROS because of squirming stresslet


!!====================================================================
!! 7. REPULSIVE POTENTIAL PARAMETERS
!!====================================================================
real(kind=8)     :: POW_RHO    !< 1/Power of threshold distance
real(kind=8)     :: THRES_RHO  !< Threshold surface distance
real(kind=8)     :: POW_REP    !< Power of dimensionless surface distance
real(kind=8)     :: SIGMA_MIN  !< Minimal distance between centers

!!====================================================================
!! 7. FLAGELLAR SWIMMING PARAMETERS...TO ADD IF EVERYTHING OK
!!====================================================================

end module FCM_FORCING_VARIABLE
