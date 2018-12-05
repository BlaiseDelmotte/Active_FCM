!!=====================================================================
!!
!!=====================================================================
module rhs_variables

implicit none

!- Right-Hand-Side term of Navier-Stokes momentum equation
double complex, dimension(:,:,:,:), allocatable :: RHS_UFOU
double complex, dimension(:,:,:,:), allocatable :: RHS_VFOU
double complex, dimension(:,:,:,:), allocatable :: RHS_WFOU


!- Vorticity
!double complex, allocatable, dimension(:,:,:) :: VORTX, VORTY, VORTZ
!real(kind=8),   allocatable, dimension(:,:,:) :: ZETAX, ZETAY, ZETAZ
!real(kind=8),   allocatable, dimension(:,:,:) :: TX, TY, TZ

end module rhs_variables
