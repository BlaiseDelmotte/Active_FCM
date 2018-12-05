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


!!--------------------------------------------------------------------
integer, dimension(3) :: ISTARTI, IENDI, ISIZEI
integer, dimension(3) :: FSTARTI, FENDI, FSIZEI
!!--------------------------------------------------------------------


integer, dimension(MPI_STATUS_SIZE) :: STATUT

integer      :: NTOT
real(kind=8) :: FACTOR, NGLOB



!!====================================================================
!!- Array size
!!====================================================================
integer :: NX 
integer :: NY
integer :: NZ

!!- Size of the interpolated solution
integer :: NXI, NYI, NZI


!!- Method for saving restart file
!!  This variable is filled in subroutine INIT_RUN
!! 
!! ISAVEM = 0: Using of MPI_FILE_WRITE_ORDERED 
!!             This method do not allow a restart with a different 
!!             number of CPU
!!        = 1: Using of MPI_FILE_WRITE_AT_ALL
!!             Build a binary allowing a restart with a different
!!             number of CPU 
integer :: ISAVEM

!---------------------------------------------------------------------
!- Size of Ghost cell for interpolation
!---------------------------------------------------------------------
integer, parameter ::  NGHTCELL = 2

!- This variable will be replaced by the FLAG_TSCHEME read in the 
!  param.in file. Normaly, the ghost cells will be only defined with 
!  DPS and the size be defined directly by the order of the interpolation
!  scheme.
!
!  Up to now, we define NGHTCELL for checking the numerical development 
!  of the ghost cells.
!---------------------------------------------------------------------


end module DNS_DIM

module GEOMETRIC

implicit none

!- Box length
real(kind=8) :: LXMAX 
real(kind=8) :: LYMAX
real(kind=8) :: LZMAX

!- Mesh 
real(kind=8), dimension(:), allocatable :: XMESH
real(kind=8), dimension(:), allocatable :: YMESH
real(kind=8), dimension(:), allocatable :: ZMESH

!- Mesh 
real(kind=8), dimension(:), allocatable :: XMESHI
real(kind=8), dimension(:), allocatable :: YMESHI
real(kind=8), dimension(:), allocatable :: ZMESHI


end module GEOMETRIC


!!=====================================================================
!!
!! Fluid variable
!!
!!=====================================================================
module fluid_variable

implicit none

!- 
real(kind=8), dimension(:,:,:), allocatable :: TEMP
!- 
real(kind=8), dimension(:,:,:), allocatable :: TEMPI

end module fluid_variable

