!!=====================================================================
!!
!! Statistics
!!
!!=====================================================================
module statistics

implicit none


!- Stored number of event
integer :: NEVEN

integer :: NCYCLELAST

!- Global array containing statistics on the fluid
real(kind=8), dimension(:), allocatable :: MEAN_FLUID

!- Time-averaged statistics on the fluid
real(kind=8), dimension(:), allocatable :: MEAN_TIME_FLUID

!- Global array containing statistics on the scalar
real(kind=8), dimension(:), allocatable :: MEAN_SCL

!- Time-averaged statistics on the scalar
real(kind=8), dimension(:), allocatable :: MEAN_TIME_SCL


!- Dissipation
real(kind=8) :: EPS_FLU

!- Dissipation computed from spectrum
real(kind=8) :: EPS_FLU_SPEC

!- Mean gradients
real(kind=8), dimension(3,3,4) :: DUIDXJ



!- Global array containing statistics on the particles
real(kind=8), dimension(:,:), allocatable :: MEAN_PART_LOC
real(kind=8), dimension(:,:), allocatable :: MEAN_PART

!- Time-averaged statistics on the particles
real(kind=8), dimension(:,:), allocatable :: MEAN_TIME_PART


!integer :: NBLGRMAX

!- 
!real(kind=8), dimension(:), allocatable :: UPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UUPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VVPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WWPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WFAPT0_MEAN
!
!real(kind=8), dimension(:), allocatable :: UUFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: VVFAPT0_MEAN
!real(kind=8), dimension(:), allocatable :: WWFAPT0_MEAN
!
real(kind=8), dimension(:,:,:), allocatable :: RPX_LOC
real(kind=8), dimension(:,:,:), allocatable :: RPY_LOC
real(kind=8), dimension(:,:,:), allocatable :: RPZ_LOC

real(kind=8), dimension(:,:,:), allocatable :: RFAPX_LOC
real(kind=8), dimension(:,:,:), allocatable :: RFAPY_LOC
real(kind=8), dimension(:,:,:), allocatable :: RFAPZ_LOC


real(kind=8), dimension(:,:,:), allocatable :: RTP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RTFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RTFAPVFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RVFAPTFAP_LOC


!!- Fluid velocity spatial correlation
integer :: DIMSCOR
real(kind=8), dimension(:), allocatable :: MEAN_RUXLOC
real(kind=8), dimension(:), allocatable :: MEAN_RVXLOC
real(kind=8), dimension(:), allocatable :: MEAN_RWXLOC



!!- lagrangian correlation function of subgrid turbulence
real(kind=8), dimension(:,:,:), allocatable :: RDUFAPDUFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDVFAPDVFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDWFAPDWFAP_LOC

real(kind=8), dimension(:,:,:), allocatable :: RDUFAPUFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDVFAPVFAP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDWFAPWFAP_LOC

real(kind=8), dimension(:,:,:), allocatable :: RDUFAPUP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDVFAPVP_LOC
real(kind=8), dimension(:,:,:), allocatable :: RDWFAPWP_LOC





!!-
!- size of array for particle concentration
integer :: NXCP
integer :: NYCP
integer :: NZCP
real(kind=8) :: RATIOCP
real(kind=8), dimension(:,:,:), allocatable :: CONCP
integer :: NPDFCP
real(kind=8), dimension(:,:), allocatable :: PDFCP_LOC



end module STATISTICS
