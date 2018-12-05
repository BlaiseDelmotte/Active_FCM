!!=====================================================================
!!
!! FCM Particles Variables
!!
!! Description :
!!> @brief
!!> Module containing all geometrical and cinematical characteristics of FCM particles
!!
!!
!! Date :  14/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!
!!=====================================================================
module FCM_PART_VARIABLE

implicit none



!!====================================================================
!! 0. PARTICLE POLYDISPERSITY
!!====================================================================

integer :: FCM_NSPHERE, 	FCM_NSPHERE_1      !< Number of spheres, number of spheres of type 1 (for bidispersity)
integer :: FCM_NELLIPSOID	                 !< Number of ellipsoids

real(kind=8) :: FCM_FRACTION                 !< Volumic fraction

!!====================================================================
!! 1. COMMON GAUSSIAN PARAMETERS
!!====================================================================

integer :: FCM_JX, FCM_JY, FCM_JZ	                !< Local indices in the Gaussian support
integer :: FCM_NGD_REF                            !< Gaussian Support for particle of size 1.0
integer :: FCM_NGD_MAX                            !<  Max Size of Gaussian Support
integer, dimension(:), allocatable :: FCM_NGD 	                                !< Gaussian support
integer, dimension(:), allocatable :: FCM_NGDH	                                !< Half Gaussian support
real(kind=8) :: FCM_SIGP                            !< Gaussian width

integer, dimension(:,:),allocatable :: FCM_LHNODE                !< Central node of gaussian  


!!====================================================================
!! 1a. SPHERE GAUSSIAN PARAMETERS
!!====================================================================
  
real(kind=8), dimension(2) :: FCM_SPHERE_SIZE  !< Default sphere size

integer, dimension(:,:),allocatable :: FCM_SPHERE_IXP, FCM_SPHERE_IYP, FCM_SPHERE_IZP !< Gaussian nodes on the mesh

real(kind=8), dimension(:), allocatable :: FCM_SPHERE_ANORM           !< Gaussian prefactor
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_ANORM_DIP       !< Dipole Gaussian prefactor
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIG2SQ          !< Denominator in the exp of the Gaussian
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIGSQ_DIP       !< Half of denominator in the exp of the dipole Gaussian
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIG2SQ_DIP      !< Denominator in the exp of the dipole Gaussian

real(kind=8), dimension(:), allocatable :: FCM_SPHERE_RADP            !< Particle radius n
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIGMA           !< Gaussian width for each particle
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIGMA_DIP       !< Dipole Gaussian width for each particle


real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_GAUSS1, &
                                            FCM_SPHERE_GAUSS2, &
                                            FCM_SPHERE_GAUSS3          !< Monopole Gaussian distribution                                            
                                            
real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_DIP_GAUSS1, &
                                            FCM_SPHERE_DIP_GAUSS2, &
                                            FCM_SPHERE_DIP_GAUSS3      !< Dipole Gaussian distribution
                                            
real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_GRAD_DIP_GAUSS1, &
                                            FCM_SPHERE_GRAD_DIP_GAUSS2, &
                                            FCM_SPHERE_GRAD_DIP_GAUSS3 !< Grad dipole Gaussian distribution


!!====================================================================
!! 1b. ELLIPSOID GAUSSIAN PARAMETERS
!!====================================================================
 
integer :: FCM_ELLIPSOID_NGD 	             !< Ellipsoid Gaussian support
integer :: FCM_ELLIPSOID_NGDH 	             !< Ellipsoid half Gaussian support
  
real(kind=8), dimension(3) :: FCM_ELLIPSOID_SIZE  !< Default ellipsoid size

integer, dimension(:,:),allocatable :: FCM_ELLIPSOID_IXP, &
                                       FCM_ELLIPSOID_IYP, &
                                       FCM_ELLIPSOID_IZP !< Gaussian nodes on the mesh

real(kind=8), dimension(:), allocatable :: FCM_ELLIPSOID_ANORM           !< Gaussian prefactor
real(kind=8), dimension(:), allocatable :: FCM_ELLIPSOID_ANORM_DIP       !< Dipole Gaussian prefactor

real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIG2SQ          !< Denominator in the exp of the Gaussian
real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIGSQ_DIP       !< Half of denominator in the exp of the dipole Gaussian
real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIG2SQ_DIP      !< Denominator in the exp of the dipole Gaussian

real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_RADP            !< Particle radius n
real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIGMA           !< Gaussian width for each particle
real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIGMA_DIP       !< Dipole Gaussian width for each particle


!~ real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_GAUSS        !< Monopole Gaussian distribution
!~ real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_DIP_GAUSS    !< Dipole Gaussian distribution
!~ real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_GRAD_DIP_GAUSS1, &
!~                                                 FCM_ELLIPSOID_GRAD_DIP_GAUSS2, &
!~                                                 FCM_ELLIPSOID_GRAD_DIP_GAUSS3  !< Grad dipole Gaussian distribution

!!====================================================================
!! 2. SWIMMING AND SQUIRMING GAUSSIAN PARAMETERS
!!====================================================================


!!====================================================================
!! 2a. SPHERE SQUIRMING GAUSSIAN PARAMETERS
!!====================================================================
  


! Useful only for squirming gaussians
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_SIGSQ      !< Half of denominator in the exp of the Gaussian
real(kind=8), dimension(:), allocatable :: FCM_SPHERE_ANORM_POTDIP   !< Factor for potential dipole enveloppe



!~! integer, dimension(:,:),allocatable :: FCM_SPHERE_IXP_SQ, &
!~!                                        FCM_SPHERE_IYP_SQ, &
!~!                                        FCM_SPHERE_IZP_SQ !< Gaussian nodes on the mesh

!~! real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_GAUSS1_SQ, &
!~!                                             FCM_SPHERE_GAUSS2_SQ, &
!~!                                             FCM_SPHERE_GAUSS3_SQ        !< Monopole Gaussian distribution
                                            
real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_GRAD_GAUSS1_SQ, &
                                            FCM_SPHERE_GRAD_GAUSS2_SQ, &
                                            FCM_SPHERE_GRAD_GAUSS3_SQ    !< Dipole Gaussian distribution
                                            
!~! real(kind=8), dimension(:,:),allocatable :: FCM_SPHERE_POTDIP_GAUSS1_SQ, &
!~!                                             FCM_SPHERE_POTDIP_GAUSS2_SQ, &
!~!                                             FCM_SPHERE_POTDIP_GAUSS3_SQ !< Potential dipole Gaussian distribution
                                            
real(kind=8), dimension(:,:),allocatable :: FCM_COEFF_SPHERE_POTDIP_GAUSS1_SQ, &
                                            FCM_COEFF_SPHERE_POTDIP_GAUSS2_SQ, &
                                            FCM_COEFF_SPHERE_POTDIP_GAUSS3_SQ !< Coeff potential dipole Gaussian distribution (to compute 2nd derivative)     



!~ !!====================================================================
!~ !! 1b. ELLIPSOID GAUSSIAN PARAMETERS
!~ !!====================================================================


real(kind=8), dimension(:), allocatable :: FCM_ELLIPSOID_ANORM_POTDIP       !< Pot Dipole Gaussian prefactor
real(kind=8), dimension(:,:), allocatable :: FCM_ELLIPSOID_SIGSQ          !< Half of the Denominator in the exp of the Gaussian

!~ !integer, dimension(:,:),allocatable :: FCM_ELLIPSOID_IXP_SQ, &
!~!                                        FCM_ELLIPSOID_IYP_SQ, &
!~ !                                       FCM_ELLIPSOID_IZP_SQ !< Gaussian nodes on the mesh

!~! real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_GAUSS_SQ        !< Monopole Gaussian distribution

!~! real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_POTDIP_GAUSS_SQ    !< < Potential dipole Gaussian distribution

!~ real(kind=8), dimension(:,:,:,:),allocatable :: FCM_ELLIPSOID_GRAD_GAUSS1_SQ, &
!~                                                 FCM_ELLIPSOID_GRAD_GAUSS2_SQ, &
!~                                                 FCM_ELLIPSOID_GRAD_GAUSS3_SQ    !< Dipole Gaussian distribution
!~                                             
!~ real(kind=8), dimension(:,:,:,:),allocatable :: FCM_COEFF_ELLIPSOID_POTDIP_GAUSS_SQ


!!====================================================================
!! 3. PARTICLE LOCATION 
!!====================================================================

integer :: FCM_INIT_PART_POS                          !< Initialization mode. 0 : in hard, 1 : random, 2 : from a file, 3 : from a routine (place_spheres_aligned_x)

real(kind=8), dimension(:), allocatable :: FCM_XP     !< Abscissa
real(kind=8), dimension(:), allocatable :: FCM_YP     !< Ordinate
real(kind=8), dimension(:), allocatable :: FCM_ZP     !< Height

real(kind=8), dimension(:), allocatable :: FCM_XP0     !< Abscissa, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_YP0     !< Ordinate, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_ZP0     !< Height, used for RFD

real(kind=8), dimension(:), allocatable :: FCM_XP_NOPER     !< Abscissa with no periodicity
real(kind=8), dimension(:), allocatable :: FCM_YP_NOPER     !< Ordinate with no periodicity
real(kind=8), dimension(:), allocatable :: FCM_ZP_NOPER     !< Height with no periodicity



!!====================================================================
!! 4. PARTICLE ORIENTATION 
!!====================================================================

integer :: FCM_INIT_PART_ORIENT                             !< Initialization mode. 0 : in hard, 1 : random, 2 : from a file, 3 : from a routine 
integer :: FCM_USE_QUAT                     !< Integer which speicfies if we use quaternions of 3 orientation vectors
!- Using quaternions
real(kind=8) :: FCM_QUATNORM 								!< Norm of quaternion vector for renormalization

real(kind=8), dimension(:,:), allocatable :: FCM_PSWIM     	!< Swimming direction for swimmers
real(kind=8), dimension(:,:), allocatable :: FCM_PSWIM0     	!< Temporary Swimming direction for swimmers


real(kind=8), dimension(:,:), allocatable :: FCM_P2     	!< Orthogonal to swimming direction

real(kind=8), dimension(:,:), allocatable :: FCM_P3     	!< Third vector of local frame

real(kind=8), dimension(:,:), allocatable :: FCM_QUAT     	!< Quaternion defining the local body frame



real(kind=8), dimension(:,:,:), allocatable :: FCM_ROT_MAT   !< Rotation matrix between local and global frame

real(kind=8), dimension(:,:,:), allocatable :: FCM_OM_VEC_Q !< RHS for quaternion integration : Omega vec Quat

real(kind=8), dimension(:,:,:), allocatable :: FCM_OM_VEC_PSWIM !< RHS for quaternion integration : Omega vec Quat

real(kind=8), dimension(:,:,:), allocatable :: FCM_OM_VEC_P2 !< RHS for quaternion integration : Omega vec Quat

real(kind=8), dimension(:,:,:), allocatable :: FCM_OM_VEC_P3 !< RHS for quaternion integration : Omega vec Quat


!!====================================================================
!! 5. PARTICLE VELOCITIES
!!====================================================================

real(kind=8), dimension(:,:), allocatable :: FCM_UP        !< x-translationnal velocities, 3 saves in time
real(kind=8), dimension(:,:), allocatable :: FCM_VP        !< y-translationnal velocities, 3 saves in time
real(kind=8), dimension(:,:), allocatable :: FCM_WP        !< z-translationnal velocities, 3 saves in time

real(kind=8), dimension(:), allocatable :: FCM_UP0        !< x-translationnal velocities, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_VP0        !< y-translationnal velocities, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_WP0        !< z-translationnal velocities, used for RFD

real(kind=8), dimension(:), allocatable :: FCM_OMPX      !< x-rotationnal velocities
real(kind=8), dimension(:), allocatable :: FCM_OMPY      !< y-rotationnal  velocities 
real(kind=8), dimension(:), allocatable :: FCM_OMPZ      !< z-rotationnal  velocities

real(kind=8), dimension(:), allocatable :: FCM_OMPX0      !< x-rotationnal velocities, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_OMPY0     !< y-rotationnal  velocities, used for RFD
real(kind=8), dimension(:), allocatable :: FCM_OMPZ0     !< z-rotationnal  velocities, used for RFD

! - Temporary work variables before applying MPI_ALLREDUCE
real(kind=8), dimension(:), allocatable :: FCM_UP_TEMP        !< x-translationnal velocities
real(kind=8), dimension(:), allocatable :: FCM_VP_TEMP        !< y-translationnal velocities 
real(kind=8), dimension(:), allocatable :: FCM_WP_TEMP        !< z-translationnal velocities

real(kind=8), dimension(:), allocatable :: FCM_OMPX_TEMP      !< x-rotationnal velocities
real(kind=8), dimension(:), allocatable :: FCM_OMPY_TEMP      !< y-rotationnal  velocities 
real(kind=8), dimension(:), allocatable :: FCM_OMPZ_TEMP      !< z-rotationnal  velocities

!!====================================================================
!! 6. PARTICLE RATE OF STRAIN
!!====================================================================

real(kind=8), dimension(:,:), allocatable :: FCM_EIJ        !< Rate of strain of particles

! - Temporary work variables before applying MPI_ALLREDUCE
real(kind=8), dimension(:,:), allocatable :: FCM_EIJ_TEMP     !< Rate of strain of particles




end module FCM_PART_VARIABLE
