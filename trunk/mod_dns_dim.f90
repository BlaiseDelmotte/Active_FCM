!!====================================================================
!!
!! Dimension of arrays 
!!
!!====================================================================
!!
!!    2     4     6     8    10    12    14    16     18
!!   20    22    24    28    30    32    36    40     42
!!   44    48    56    60    64    66    70    72     80
!!   84    88    90    96   110   112   120   126    128
!!  132   140   144   154   160   168   176   180    192
!!  198   210   220   224   240   252   256   264    280
!!  288   308   320   330   336   352   360   384    396
!!  420   440   448   462   480   504   512   528    560
!!  576   616   630   640   660   672   704   720    768
!!  770   792   840   880   896   924   960   990   1008
!! 1024
!!
!!====================================================================
module DNS_DIM

use MPI

implicit none


!!--------------------------------------------------------------------
!! Local dimension and mpi stuffs
!!--------------------------------------------------------------------
integer               :: IERR,NU,NDIM,NPROC,MYID
integer, dimension(2) :: DIMS
integer, dimension(3) :: ISTART, IEND, ISIZE
integer, dimension(3) :: FSTART, FEND, FSIZE
integer               :: IPROC,JPROC

integer, dimension(MPI_STATUS_SIZE) :: STATUT

integer      :: NTOT
real(kind=8) :: FACTOR, NGLOB


character(len=3) :: FFTFLAG


!!====================================================================
!!- Array size
!!====================================================================
integer :: NX 
integer :: NY
integer :: NZ

!- Number of forced wavenumber
integer :: NWAVEFORCE


!---------------------------------------------------------------------
!- Size of Ghost cell for interpolation
!---------------------------------------------------------------------
integer ::  NGHTCELL

!- This variable will be replaced by the FLAG_TSCHEME read in the 
!  param.in file. Normaly, the ghost cells will be only defined with 
!  DPS and the size be defined directly by the order of the interpolation
!  scheme.
!
!  Up to now, we define NGHTCELL for checking the numerical development 
!  of the ghost cells.
!---------------------------------------------------------------------


!!--------------------------------------------------------------------
!! Physical constant
!!--------------------------------------------------------------------
real(kind=8) :: PPI
real(kind=8) :: TWOPI

double complex :: ICMPL


!!- Particle's arrays
integer :: NIG       !- Particle class number

!!- Number of particles per CPU for uniforme distribution
integer :: NPCPU_UNIF

!!- Number of particles
integer :: NPART_FULL     

!!- maximum exchanged particles per CPU     
integer :: NPEXCH_MAX

!!- maximum number of particle in a cpu   
integer :: NPMAX_LOC
integer, dimension(:), allocatable :: NPMAX_CPU
integer, dimension(:), allocatable :: NPMIN_CPU


!!- Local number of particle
integer, dimension(:), allocatable :: NPART_LOC


!- Size of array containing statistics
integer :: NSTAT

!- Number of lagrangian correlation
integer, parameter :: NBLGRMAX = 3

integer,dimension(NBLGRMAX) :: NT0

integer :: NTLGRMAX !!- ca sert probablement a rien


!- Length of Lagrangian time-correlation
integer :: DIMLGR




!- Max Number of particle used to compute the lagrangian correlation
integer :: MAXPARTLAG

integer :: DIMLGRMAX


!- Flag for time integration
integer :: FLAG_TSCHEME


integer, dimension(5) :: UNIT_INFO !- unit number 



!!- Method for saving restart file
!!  This variable is filled in subroutine INIT_RUN
!! 
!! ISAVEFLUID = 1: Multiple binary files (one per CPU)
!!            = 2: Direct access file
!!            = 3: MPI_FILE_WRITE_ORDERED
!!            = 4: MPI I/O (recommended)
integer :: ISAVEFLUID

!! 
!! ISAVEPART = 1: Multiple binary files (one per CPU)
!!           = 2: Direct access file
!!           = 3: MPI_FILE_WRITE_ORDERED
!!           = 4: MPI I/O (recommended)
integer :: ISAVEPART


!!====================================================================
!!- time 
!!====================================================================

!- cycle at t
integer :: TN

!- cycle at t-dt
integer :: TNM1

!- cycle at t-2*dt
integer :: TNM2



!!====================================================================
!! Numerical constant
!!====================================================================
real(kind=8), parameter :: ZERO = 1D-33
real(kind=8), parameter :: INFINITY = 1D+33



!- File extension containing the CPU Id
character (len=10) :: FILE_EXT

!- Cycle of file printing
integer :: NFILEOUT


!- Flag for solving Navier-Stokes
!  SOLVE_FLUID = 0 : Frozen flow
!              = 1 : DNS Homogeneous Isotropic
!              = 2 : Stokes flow
integer :: SOLVE_FLUID


!- Flag for Lagrangian particle tracking
logical :: SOLVE_PART

!- Flag for Lagrangian particle collision 
integer :: SOLVE_COLLISION

!- Flag for passive scalar equation
logical :: SOLVE_SCALAR

!!--------------------------------------------------------------------
!! Statistics
!!--------------------------------------------------------------------
logical :: READSTAT     !- Flag to read from previous run
logical :: STAT_TIME    !- Flag to perform temporal statistics

logical :: LEVEL0_STFLU !- Statistics or not
logical :: LEVEL1_STFLU !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
logical :: LEVEL2_STFLU !- Divergence, dissipation
logical :: LEVEL3_STFLU !- Skewness, Flatness, PDF of fluid velocity & gradients
logical :: LEVEL4_STFLU !- Two points statistics (Spatial correlation)
logical :: LEVEL5_STFLU !- Eulerian time-correlation

logical :: LEVEL0_STPAR !- Statistics or not
logical :: LEVEL1_STPAR !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
logical :: LEVEL2_STPAR !- Lagrangian functiond
logical :: LEVEL3_STPAR !- Spatial distribution


logical :: LEVEL0_STSCL !- Statistics or not
logical :: LEVEL1_STSCL !- One point statistics (mean, variance, ...)
logical :: LEVEL2_STSCL !- One point statistics (mean, variance, ...)


!!====================================================================
!! Filtering
!!====================================================================
logical :: FILTERING
real(kind=8) :: KCUT


!!--------------------------------------------------------------------
!! Debugging
!!--------------------------------------------------------------------
!- Debug flag
logical :: DEBUG

end module DNS_DIM
