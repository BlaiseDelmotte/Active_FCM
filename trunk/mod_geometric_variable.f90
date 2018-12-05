!!=====================================================================
!!
!! Geometrics (mesh, wavenumber and integrating factor)
!!
!!=====================================================================
module GEOMETRIC_VARIABLE

implicit none

!- mesh step
real(kind=8) :: DX, DY, DZ

!- Box length
real(kind=8) :: LXMAX 
real(kind=8) :: LYMAX
real(kind=8) :: LZMAX

!- Mesh 
real(kind=8), dimension(:), allocatable :: XMESH
real(kind=8), dimension(:), allocatable :: YMESH
real(kind=8), dimension(:), allocatable :: ZMESH

real(kind=8) :: KFIRST !- First solved wavenumber
real(kind=8) :: KLAST  !- Last solved wavenumber

!- Wavenumber
real(kind=8), dimension(:), allocatable :: KX
real(kind=8), dimension(:), allocatable :: KY
real(kind=8), dimension(:), allocatable :: KZ




!- Integrating factor (for time integration)
real(kind=8), dimension(:,:,:), allocatable :: FILTER
 


end module GEOMETRIC_VARIABLE
