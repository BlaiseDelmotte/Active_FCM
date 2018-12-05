!!=====================================================================
!!
!!   Build the Right-Hand-Side of Lagrangian particle trajectory
!!
!!=====================================================================

subroutine RHS_PARTICLE

!!=====================================================================
!! The following fortran file performs time advancing of particle
!! position. The considered forces acting on the particle are:
!!   + Drag
!!   + Gravity
!!
!! So the momentum particle governing equation is:
!!
!!    dup,i     up,i - uf@p,i
!!    ----- = - ------------- - gi
!!     dt           tau_p 
!!
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

use dns_dim
use particle_variable
use param_phys 
use check_cpu
use statistics,  only: MEAN_PART

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Fluid-Particle relative velocity
real(kind=8) :: VRX
real(kind=8) :: VRY
real(kind=8) :: VRZ
real(kind=8) :: VRNRM

!- Particle Reynolds number
real(kind=8) :: REP

!- Drag coefficient
real(kind=8) :: CDRAG

!- Inverse of particle response time INVTAUP = 1/TAUP
real(kind=8) :: INVTAUP

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J
!------------------------------------------------------------------


!!- CPU check




!- Loop over particle kind
do J = 1, NIG

!!===================================================================
!! 1. Motionless particle
!!===================================================================
!if((PARTDEF(J) == 0).or.(PARTDEF(J) == 1)) then
!if(PARTDEF(J) == 0) then

! UPART(:,J) = ZERO
! VPART(:,J) = ZERO
! WPART(:,J) = ZERO

!!===================================================================
!! 2. Fluid Element
!!===================================================================
!elseif(PARTDEF(J) == 1) then

! UPART(:,J) = UFAP(:,J)
! VPART(:,J) = VFAP(:,J)
! WPART(:,J) = WFAP(:,J)


!!===================================================================
!! 3. Solid inertial particle
!!===================================================================
!elseif(PARTDEF(J)==2) then
if(PARTDEF(J)==2) then

 do I = 1, NPMAX 

!!------------------------------------------------------------------- 
!! 3.1. Particle Reynolds number 
!!-------------------------------------------------------------------

 !- Fluid-Particle Relative velocity
 VRX = UFAP(I,J) - UPART(I,J)
 VRY = VFAP(I,J) - VPART(I,J)
 VRZ = WFAP(I,J) - WPART(I,J)

 VRNRM = sqrt(VRX*VRX + VRY*VRY + VRZ*VRZ)

 !- Particle Reynolds number
 REP = DPART(J)*VRNRM/VISC

!!-------------------------------------------------------------------
!! 3.1. Drag response time
!!-------------------------------------------------------------------
!!- Schiller & Nauman (1935) 
 CDRAG = 24./REP*(1. + 0.15*REP**0.687)

 !- Particle response time
 INVTAUP = 3./4.*RHOF/RHOP(J)*CDRAG/DPART(J)*VRNRM 


!!-------------------------------------------------------------------
!! 3.2. Particle velocity
!!-------------------------------------------------------------------
 !- Integrating factor
 INTFACTOR_PART(I,J)  = exp(- DTIME * INVTAUP)

 !- Right-Hand-Side
 RHS_UPART(I,J,TN) = UFAP(I,J) * INVTAUP
 RHS_VPART(I,J,TN) = VFAP(I,J) * INVTAUP
 RHS_WPART(I,J,TN) = WFAP(I,J) * INVTAUP - GRAVITY(J)



!!-------------------------------------------------------------------
!! 3.3. Particle statistics
!!-------------------------------------------------------------------

 end do

 end if


 if(LEVEL0_STPAR.and.(PARTDEF(J)==2)) then
 MEAN_PART(28,J) = MEAN_PART(28,J) + VRNRM /real(NPMAX)
 MEAN_PART(29,J) = MEAN_PART(29,J) + REP /real(NPMAX)
 MEAN_PART(30,J) = MEAN_PART(30,J) + CDRAG /real(NPMAX)
 MEAN_PART(31,J) = MEAN_PART(31,J) + INVTAUP /real(NPMAX)
 else
 MEAN_PART(28,J) = ZERO
 MEAN_PART(29,J) = ZERO
 MEAN_PART(30,J) = ZERO
 MEAN_PART(31,J) = ZERO
 end if


end do 


!!===================================================================
!! 4. Right-Hand-Side of particle position
!!===================================================================
do J=1, NIG
if(PARTDEF(J)<2) then
RHS_XPART(:,J,TN) = UFAP(:,J)
RHS_YPART(:,J,TN) = VFAP(:,J)
RHS_ZPART(:,J,TN) = WFAP(:,J)
else
RHS_XPART(:,J,TN) = UPART(:,J)
RHS_YPART(:,J,TN) = VPART(:,J)
RHS_ZPART(:,J,TN) = WPART(:,J)
end if
end do




end subroutine RHS_PARTICLE
