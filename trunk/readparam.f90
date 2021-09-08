!!====================================================================
!!
!! Routine reading the simulation parameters from the ASCII file
!!                             'param.in'
!!
!!====================================================================

subroutine READPARAM

!!====================================================================
!! 
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use FORCING
use GEOMETRIC_VARIABLE

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- index
integer :: I
!---------------------------------------------------------------------


!- Open file
open(unit=200,file='param.in',status='old')

read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


!!====================================================================
!! 1. NUMERICS
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) SOLVE_FLUID
read(200,*) ! comment line

read(200,*) NX        !- Mesh size
read(200,*) NY        !- Mesh size
read(200,*) NZ        !- Mesh size
read(200,*)           ! comment line
read(200,*) DTIME     !- Time step
read(200,*) NCYCLEMAX !- Maximum of cycle


!!====================================================================
!! 2. Fluid numerics
!!====================================================================
read(200,*)                     ! comment line
read(200,*) INIT_FLUID_VELOCITY !- Fluid initiation
read(200,*) STEADY              !- Flag for forcing
read(200,*) SIGMA_FORCE         !- Intensity
read(200,*) TIME_FORCE          !- Time of forcing
read(200,*) KFORCE_MIN          !- Wavenumber
read(200,*) KFORCE_MAX          !- Wavenumber

!!====================================================================
!! 3. Scalar numerics
!!====================================================================
read(200,*) ! comment line
read(200,*) SOLVE_SCALAR
read(200,*) INIT_SCALAR


!!====================================================================
!! 4. Particle numerics
!!====================================================================
read(200,*)            ! comment line
read(200,*) SOLVE_COLLISION !- Contact Handling
read(200,*) NPART_FULL !- Number of particle

!- Interpolation scheme
! INTERP_SCHEME = 1 : Lagrangian polynomial 1st order
!               = 2 : Lagrangian polynomial 2nd order
!               = 3 : Lagrangian polynomial 3rd order
!               = 4 : SFM



!!====================================================================
!! 5. Outputing
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) READSTAT     !- read statistics from previous run
read(200,*) STAT_TIME    !- Flag to perform averaging over the time
read(200,*)              ! comment line
                           !- Level of statistics for the fluid
read(200,*) LEVEL1_STFLU !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
read(200,*) LEVEL2_STFLU !- Divergence, dissipation
read(200,*) LEVEL3_STFLU !- Skewness, Flatness, PDF of fluid velocity & gradients
read(200,*) LEVEL4_STFLU !- Two points statistics (Spatial correlation)
read(200,*)              ! comment line
read(200,*) LEVEL1_STSCL !- One point statistics (mean, Reynolds stress, 3rd & 4th-order)
read(200,*) LEVEL2_STSCL !- Dissipation and gradients
read(200,*)              ! comment line

!- Outputing parameters
!---------------------
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) FOUT0       !- Screen printing
read(200,*) FOUT1       !- Fluid Solution printing
read(200,*) FOUT2       !- Stat printing
read(200,*) FOUT3       !- Part Solution printing
read(200,*) FOUT4       !- To be defined
read(200,*) FOUT5       !- To be defined
read(200,*) ENSIGHT_OUT !- Ensight's format velocity field


!!====================================================================
!! 6. Fluid properties
!!====================================================================
read(200,*) ! comment line
read(200,*) ! comment line
read(200,*) ! comment line


read(200,*) LXMAX !- Computational domain size
read(200,*) LYMAX !- Computational domain size
read(200,*) LZMAX !- Computational domain size
read(200,*) VISC  !- Molecular viscosity
read(200,*) RHOF  !-  Fluid density


!!====================================================================
!! 7. Scalar field properties
!!====================================================================
read(200,*)          ! comment line
read(200,*)          ! comment line
read(200,*)          ! comment line
read(200,*) DIFF_SCL !- Scalar diffusivity
read(200,*) CP_SCL   !- Specific heat of scalar
read(200,*) GRAD_SCL !- Mean gradient of scalar


close(200)


!- Solution printing
if(FOUT1>NCYCLEMAX) FOUT1 = NCYCLEMAX

!- Statistic printing
if(FOUT2>NCYCLEMAX) FOUT2 = NCYCLEMAX

!- Statistic printing
if(FOUT3>NCYCLEMAX) FOUT3 = NCYCLEMAX


if(MYID==0) write(*,*) 'Parameters reading --> OK'


end subroutine READPARAM
