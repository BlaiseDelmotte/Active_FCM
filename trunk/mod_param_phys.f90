!!====================================================================
!!
!! Physical parameters
!!
!!====================================================================
module PARAM_PHYS

implicit none


!!--------------------------------------------------------------------
!!- Fluid physical parameter
!!--------------------------------------------------------------------
real(kind=8) :: VISC      !- Molecular viscosity
real(kind=8) :: RHOF      !- Fluid density

integer :: INIT_FLUID_VELOCITY
!- Fluid Initiation
! = 0 : Uf = 0
! = 1 : Random value
! = 2 : Uniform fluid velocity (defined in fortran file)
! = 3 : Read fluid velocity field read from thi generating code
! = 4 : Read fluid velocity field read from stored files

real(kind=8), dimension(3) :: UREF !- Initial mean velocity field


!!--------------------------------------------------------------------
!!- Scalar parameter
!!--------------------------------------------------------------------
real(kind=8) :: DIFF_SCL  !- Diffusivity of Passive SCalar
real(kind=8) :: CP_SCL    !- "Specific heat" of scalar
real(kind=8) :: GRAD_SCL  !- Gradient of scalar

integer :: INIT_SCALAR

!!--------------------------------------------------------------------
!! Numerical parameter
!!--------------------------------------------------------------------
real(kind=8)  :: DTIME     !- Time step
integer       :: NCYCLEMAX !- Maximum of cycle number

!- Forced turbulence
logical :: STEADY

!!--------------------------------------------------------------------
!! Output parameters
!!--------------------------------------------------------------------
integer :: FOUT0       !- Screen printing
integer :: FOUT1       !- Fluid Solution printing
integer :: FOUT2       !- Stat printing
integer :: FOUT3       !- Particles Solution printing
integer :: FOUT4       !- To define
integer :: FOUT5       !- To define

integer :: ENSIGHT_OUT !- Ensight's format





end module param_phys
