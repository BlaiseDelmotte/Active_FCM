!!=====================================================================
!!
!! Fluid variable
!!
!!=====================================================================
module fluid_variable

implicit none

!- Fluid velocity at tn
real(kind=8), dimension(:,:,:), allocatable :: UFLU
real(kind=8), dimension(:,:,:), allocatable :: VFLU
real(kind=8), dimension(:,:,:), allocatable :: WFLU

!- Complex velocity
double complex, dimension(:,:,:), allocatable :: UFOU
double complex, dimension(:,:,:), allocatable :: VFOU
double complex, dimension(:,:,:), allocatable :: WFOU

end module fluid_variable
