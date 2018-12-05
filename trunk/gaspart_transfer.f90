!!=====================================================================
!!
!!  Model for drag force
!!
!!=====================================================================

subroutine GASPART_TRANSFER(                    &
                         IG,                    &
                      UPART, VPART, WPART,      &
                       UFAP,  VFAP,  WFAP,      &
                    INVTAUP, VRNRM, CDRAG, REP, &
              INVTAUP_THETA,   NUP              )  

!!=====================================================================
!! The following fortran file performs time advancing of particle
!! position. The forces acting on the particles are:
!!   + Drag
!!   + Gravity
!!
!! So the momentum particle governing equation is:
!!
!!    dup,i     up,i - uf@p,i
!!    ----- = - ------------- - gi
!!     dt           tau_p 
!!
!! For time integration, this equation is rewriten as
!!
!!    dup,i     
!!    ----- = -COEFA*up,i + COEFB
!!     dt          
!!
!! with
!!
!!             1              uf@p,i
!!   COEFA = -----   COEFB = -------- - gi
!!           tau_p            tau_p
!!
!!--------------------------------------------------------------------
!! with tau_p the particle response time to drag force computed
!! in terms of particle Reynolds number (Rep) and drag coefficient
!! according to the Schiller & Nauman (1935) correlation:
!!
!!        dp.|up - uf@p|          24
!! Rep = ----------------     Cd=----.(1+0.15*Rep^0.687)
!!            nu_f                Rep
!!
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
!!
!!=====================================================================

use PARAM_PHYS 
use DNS_DIM, only: SOLVE_SCALAR

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
integer     , intent(in) :: IG
real(kind=8), intent(in) :: UPART
real(kind=8), intent(in) :: VPART
real(kind=8), intent(in) :: WPART
real(kind=8), intent(in) :: UFAP
real(kind=8), intent(in) :: VFAP
real(kind=8), intent(in) :: WFAP

!- Inverse of particle response time
real(kind=8), intent(out) :: INVTAUP

!- Fluid-Particle relative velocity
real(kind=8), intent(out) :: VRNRM

!- Particle Reynolds number
real(kind=8), intent(out) :: REP

!- Drag coefficient
real(kind=8), intent(out) :: CDRAG

!- Inverse of particle temperature response time
real(kind=8), intent(out) :: INVTAUP_THETA

!- Nusselt number
real(kind=8), intent(out) :: NUP


!---------------------------------------------------------------------
! LOCAL ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Fluid-Particle relative velocity
real(kind=8) :: VRX
real(kind=8) :: VRY
real(kind=8) :: VRZ

!- Prandtl number
real(kind=8) :: PRANDTL

!------------------------------------------------------------------


!!------------------------------------------------------------------- 
!! Particle Reynolds number 
!!-------------------------------------------------------------------
!- Fluid-Particle Relative velocity
VRX = UPART - UFAP
VRY = VPART - VFAP
VRZ = WPART - WFAP

VRNRM = sqrt(VRX*VRX + VRY*VRY + VRZ*VRZ)


!- Particle Reynolds number
REP = DPART(IG)*VRNRM/VISC

!- Stokes's law
CDRAG = 24./REP

INVTAUP = 18.*RHOF*VISC/RHOP(IG)/DPART(IG)**2

!- Prandtl number
PRANDTL = VISC/DIFF_SCL


!- Stokes's for heat
NUP = 2.0

if(SOLVE_SCALAR) INVTAUP_THETA = NUP*CP_SCL/(3.*CP_PART(IG)*PRANDTL)*INVTAUP



if(REP>1E-5) then
!!-------------------------------------------------------------------
!! Drag response time --> Schiller & Nauman (1935) 
!!-------------------------------------------------------------------
!- Schiller & Nauman (1935) 
 CDRAG = 24./REP*(1. + 0.15*REP**0.687)

!- Inverse of particle response time
 INVTAUP = 3./4.*RHOF/RHOP(IG)*CDRAG/DPART(IG)*VRNRM 


!!-------------------------------------------------------------------
!! Heat transfer --> Ranz-Marshal
!!-------------------------------------------------------------------
 if(SOLVE_SCALAR) then
 
  NUP = 2.0 + 0.55*REP**0.5*PRANDTL**(1./3.)

  INVTAUP_THETA = 6.*NUP*VISC*RHOF*CP_SCL/CP_PART(IG)/DPART(IG)**2/RHOP(IG)/PRANDTL
 
 end if

end if


end subroutine GASPART_TRANSFER
