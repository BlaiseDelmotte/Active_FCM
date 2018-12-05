!!=====================================================================
!!
!!   Time-advancing particle
!!
!!=====================================================================

subroutine ADV_PARTICLE_VELOCITY

!!=====================================================================
!! The numerical scheme used for time-advancing is a 2nd order
!! Adams-Bashforth coupled with integrating factor.
!! 
!!  df(t)                      d[f(t).exp(at)]
!! ------ + af(t) + b = 0  <=> --------------- + b.exp(at) = 0
!!   dt                               dt
!!             
!! 3rd order Adams-Bashforth is written as:
!!
!!  df(t)                f^n+1 - f^n     3            1
!! ------ = phi(t)  ==> ------------- = ---.phi^n  -.---.phi^n-1 
!!   dt                       Dt         2            2
!!
!!---------------------------------------------------------------------
!! Warning: The variable FLAG_TSCHEME is modified in main.f90 file 
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
!!
!! The strongest assumption in the is that the COEFA is assumed
!! nearly the same between all time step --> as OCEFA is 1/taup such
!! an assumption is quiet legitime.
!!
!!
!!=====================================================================

use DNS_DIM
use PARTICLE_PARALLEL
use PARAM_PHYS 
use CHECK_CPU

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Coefficient for time integration (3rd order)
real(kind=8), dimension(3) :: COEFA
real(kind=8), dimension(3) :: COEFBX	   
real(kind=8), dimension(3) :: COEFBY	   
real(kind=8), dimension(3) :: COEFBZ   


!- Inverse of particle response time INVTAUP = 1/TAUP
real(kind=8) :: INVTAUP

!- Fluid-Particle relative velocity
real(kind=8) :: VRNRM

!- Particle Reynolds number
real(kind=8) :: REP

!- Drag coefficient
real(kind=8) :: CDRAG

!- Time step
real(kind=8) ::DTIME2, DTIME12


!!-----------------------------------------------------
!! Variable for temperature
!!-----------------------------------------------------
!- Inverse of temperature particle response time INVTAUP = 1/TAUP
real(kind=8) :: INVTAUP_THETA

!- Nusselt
real(kind=8) :: NUP

!- Prandtl
real(kind=8) :: PRANDTL

real(kind=8), dimension(3) :: COEFA_THETA
real(kind=8), dimension(3) :: COEFB_THETA
!!-----------------------------------------------------




!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J
!------------------------------------------------------------------


!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if



DTIME12 = DTIME/12.
DTIME2 = DTIME/2.


!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
if(FLAG_TSCHEME == 1) then

do J =1, NIG

if(PARTDEF(J)==2) then


do I = 1, NPART_LOC(J)

COEFA(:)  = ZERO
COEFBX(:) = ZERO
COEFBY(:) = ZERO
COEFBZ(:) = ZERO


if(SOLVE_SCALAR) then
 COEFA_THETA(:) = ZERO
 COEFB_THETA(:) = ZERO
end if


!!--------------------------------------------------------------------
!! 1.1. Fluid-particle transfer
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(           J,  & ! <- Particle class
                       PART(I,J)%UP,  & ! <- Part x-vel
		       PART(I,J)%VP,  & ! <- Part y-vel
		       PART(I,J)%WP,  & ! <- Part z-vel
		     PART(I,J)%UFAP,  & ! <- Fluid x-vel
		     PART(I,J)%VFAP,  & ! <- Fluid y-vel
		     PART(I,J)%WFAP,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

!!
  COEFA(1) = INVTAUP
 COEFBX(1) = INVTAUP*PART(I,J)%UFAP
 COEFBY(1) = INVTAUP*PART(I,J)%VFAP
 COEFBZ(1) = INVTAUP*PART(I,J)%WFAP - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(1) = INVTAUP_THETA
  COEFB_THETA(1) = INVTAUP_THETA*PART(I,J)%TFAP
 end if


!!--------------------------------------------------------------------
!! 1.2. Solve gas-particle heat tansfer
!!--------------------------------------------------------------------
 if(SOLVE_SCALAR) then
  PART(I,J)%TP = exp(-COEFA_THETA(1)*DTIME)*(PART(I,J)%TP + DTIME*COEFB_THETA(1))
 end if


!!--------------------------------------------------------------------
!! 1.3. Solve particle momentum equation
!!--------------------------------------------------------------------

 PART(I,J)%UP = exp(-COEFA(1)*DTIME)*(PART(I,J)%UP + DTIME*COEFBX(1) )
 PART(I,J)%VP = exp(-COEFA(1)*DTIME)*(PART(I,J)%VP + DTIME*COEFBY(1) )
 PART(I,J)%WP = exp(-COEFA(1)*DTIME)*(PART(I,J)%WP + DTIME*COEFBZ(1) )


end do !!- Loop: NPART_LOC(J)

end if !!- If: PARTDEF(J)==2

end do !!- Loop: NIG



!!====================================================================
!! 2. 2nd Order Adam-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 2) then

do J =1, NIG

if(PARTDEF(J)==2) then


do I = 1, NPART_LOC(J)

COEFA(:)  = ZERO
COEFBX(:) = ZERO
COEFBY(:) = ZERO
COEFBZ(:) = ZERO


!!--------------------------------------------------------------------
!! 2.1. Fluid-particle transfer
!!--------------------------------------------------------------------

!!--------------------------------------------------------------------
!!- Step n
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(           J,  & ! <- Particle class
                       PART(I,J)%UP,  & ! <- Part x-vel
		       PART(I,J)%VP,  & ! <- Part y-vel
		       PART(I,J)%WP,  & ! <- Part z-vel
		     PART(I,J)%UFAP,  & ! <- Fluid x-vel
		     PART(I,J)%VFAP,  & ! <- Fluid y-vel
		     PART(I,J)%WFAP,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number


 COEFA(1) = INVTAUP
 COEFBX(1) = INVTAUP*PART(I,J)%UFAP
 COEFBY(1) = INVTAUP*PART(I,J)%VFAP
 COEFBZ(1) = INVTAUP*PART(I,J)%WFAP - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(1) = INVTAUP_THETA
  COEFB_THETA(1) = INVTAUP_THETA*PART(I,J)%TFAP
 end if


!!--------------------------------------------------------------------
!!- Step n-1
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(           J,  & ! <- Particle class
                   PART(I,J)%UP_NM1,  & ! <- Part x-vel
	           PART(I,J)%VP_NM1,  & ! <- Part y-vel
		   PART(I,J)%WP_NM1,  & ! <- Part z-vel
		 PART(I,J)%UFAP_NM1,  & ! <- Fluid x-vel
		 PART(I,J)%VFAP_NM1,  & ! <- Fluid y-vel
		 PART(I,J)%WFAP_NM1,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

 COEFA(2) = INVTAUP
 COEFBX(2) = INVTAUP*PART(I,J)%UFAP_NM1
 COEFBY(2) = INVTAUP*PART(I,J)%VFAP_NM1
 COEFBZ(2) = INVTAUP*PART(I,J)%WFAP_NM1 - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(2) = INVTAUP_THETA
  COEFB_THETA(2) = INVTAUP_THETA*PART(I,J)%TFAP_NM1
 end if
 


!!--------------------------------------------------------------------
!! 2.2. Solve gas-particle heat tansfer
!!--------------------------------------------------------------------
 if(SOLVE_SCALAR) then
 
 PART(I,J)%TP = exp(-COEFA_THETA(1)*DTIME)                               &
       *(PART(I,J)%TP + DTIME2*(3.*COEFB_THETA(1)                        &
                                 - COEFB_THETA(2)*exp(-COEFA_THETA(2)*DTIME)) )

end if


!!--------------------------------------------------------------------
!! 2.3. Solve particle momentum equation
!!--------------------------------------------------------------------

 PART(I,J)%UP = exp(-COEFA(1)*DTIME)                               &
       *(PART(I,J)%UP + DTIME2*(3.*COEFBX(1)                       &
                                 - COEFBX(2)*exp(-COEFA(2)*DTIME)) )

 PART(I,J)%VP = exp(-COEFA(1)*DTIME)                               &
       *(PART(I,J)%VP + DTIME2*(3.*COEFBY(1)                       &
                                 - COEFBY(2)*exp(-COEFA(2)*DTIME)) )

 PART(I,J)%WP = exp(-COEFA(1)*DTIME)                               &
       *(PART(I,J)%WP + DTIME2*(3.*COEFBZ(1)                       &
                                 - COEFBZ(2)*exp(-COEFA(2)*DTIME)) )


end do !!- Loop: NPART_LOC(J)

end if !!- If: PARTDEF(J)==2
end do !!- Loop: NIG



!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 3 ) then

do J =1, NIG

if(PARTDEF(J)==2) then


do I = 1, NPART_LOC(J)

COEFA(:)  = ZERO
COEFBX(:) = ZERO
COEFBY(:) = ZERO
COEFBZ(:) = ZERO

if(SOLVE_SCALAR) then
 COEFA_THETA(:) = ZERO
 COEFB_THETA(:) = ZERO
end if


!!--------------------------------------------------------------------
!! 3.1. Fluid-particle transfer
!!--------------------------------------------------------------------

!!--------------------------------------------------------------------
!!- Step n
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(         J,  & ! <- Particle class
                       PART(I,J)%UP,  & ! <- Part x-vel
		       PART(I,J)%VP,  & ! <- Part y-vel
		       PART(I,J)%WP,  & ! <- Part z-vel
		     PART(I,J)%UFAP,  & ! <- Fluid x-vel
		     PART(I,J)%VFAP,  & ! <- Fluid y-vel
		     PART(I,J)%WFAP,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number


 COEFA(1) = INVTAUP
 COEFBX(1) = INVTAUP*PART(I,J)%UFAP
 COEFBY(1) = INVTAUP*PART(I,J)%VFAP
 COEFBZ(1) = INVTAUP*PART(I,J)%WFAP - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(1) = INVTAUP_THETA
  COEFB_THETA(1) = INVTAUP_THETA*PART(I,J)%TFAP
 end if


!!--------------------------------------------------------------------
!!- Step n-1
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(         J,  & ! <- Particle class
                   PART(I,J)%UP_NM1,  & ! <- Part x-vel
	           PART(I,J)%VP_NM1,  & ! <- Part y-vel
		   PART(I,J)%WP_NM1,  & ! <- Part z-vel
		 PART(I,J)%UFAP_NM1,  & ! <- Fluid x-vel
		 PART(I,J)%VFAP_NM1,  & ! <- Fluid y-vel
		 PART(I,J)%WFAP_NM1,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

 COEFA(2) = INVTAUP
 COEFBX(2) = INVTAUP*PART(I,J)%UFAP_NM1
 COEFBY(2) = INVTAUP*PART(I,J)%VFAP_NM1
 COEFBZ(2) = INVTAUP*PART(I,J)%WFAP_NM1 - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(2) = INVTAUP_THETA
  COEFB_THETA(2) = INVTAUP_THETA*PART(I,J)%TFAP_NM1
 end if


!!--------------------------------------------------------------------
!!- Step n-2
!!--------------------------------------------------------------------
 call GASPART_TRANSFER(         J,  & ! <- Particle class
                   PART(I,J)%UP_NM2,  & ! <- Part x-vel
	           PART(I,J)%VP_NM2,  & ! <- Part y-vel
		   PART(I,J)%WP_NM2,  & ! <- Part z-vel
		 PART(I,J)%UFAP_NM2,  & ! <- Fluid x-vel
		 PART(I,J)%VFAP_NM2,  & ! <- Fluid y-vel
		 PART(I,J)%WFAP_NM2,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

 COEFA(3) = INVTAUP
 COEFBX(3) = INVTAUP*PART(I,J)%UFAP_NM2
 COEFBY(3) = INVTAUP*PART(I,J)%VFAP_NM2
 COEFBZ(3) = INVTAUP*PART(I,J)%WFAP_NM2 - GRAVITY(J)

 if(SOLVE_SCALAR) then
  COEFA_THETA(3) = INVTAUP_THETA
  COEFB_THETA(3) = INVTAUP_THETA*PART(I,J)%TFAP_NM2
 end if




!!--------------------------------------------------------------------
!! 3.2. Solve gas-particle heat tansfer
!!--------------------------------------------------------------------
if(SOLVE_SCALAR) then

 PART(I,J)%TP = exp(-COEFA_THETA(1)*DTIME)                 &
      *(PART(I,J)%TP                                       &
        + DTIME12*( 23.*COEFB_THETA(1)                          &
                   -16.*COEFB_THETA(2)*exp(-COEFA_THETA(2)   *DTIME)  &
                   + 5.*COEFB_THETA(3)*exp(-COEFA_THETA(3)*2.*DTIME) ) )
 
end if



!!--------------------------------------------------------------------
!! 3.3. Solve particle momentum equation
!!--------------------------------------------------------------------

 PART(I,J)%UP = exp(-COEFA(1)*DTIME)                        &
      *(PART(I,J)%UP                                        &
        + DTIME12*( 23.*COEFBX(1)                           &
                   -16.*COEFBX(2)*exp(-COEFA(2)   *DTIME)   &
                   + 5.*COEFBX(3)*exp(-COEFA(3)*2.*DTIME) ) )

 PART(I,J)%VP = exp(-COEFA(1)*DTIME)                       &
      *(PART(I,J)%VP                                       &
        + DTIME12*( 23.*COEFBY(1)                          &
                   -16.*COEFBY(2)*exp(-COEFA(2)   *DTIME)  &
                   + 5.*COEFBY(3)*exp(-COEFA(3)*2.*DTIME) ) )

 PART(I,J)%WP = exp(-COEFA(1)*DTIME)                       &
      *(PART(I,J)%WP                                       &
        + DTIME12*( 23.*COEFBZ(1)                          &
                   -16.*COEFBZ(2)*exp(-COEFA(2)   *DTIME)  &
                   + 5.*COEFBZ(3)*exp(-COEFA(3)*2.*DTIME) ) )



end do !!- Loop: NPART_LOC(J)

end if !!- If: PARTDEF(J)==2
end do !!- Loop: NIG

end if !!- If FLAG_TSCHEME



!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(3) = CPU_PART(3) + TIME_END - TIME_START
end if


end subroutine ADV_PARTICLE_VELOCITY
