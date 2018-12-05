!!====================================================================
!!
!!          Particle TRACKING
!!
!!====================================================================

subroutine PARTICLE_TRACKING(NCYCLE)

!!====================================================================
!! Full explicit Adam-Bashforth is used for the time-integration of
!! the fluid momentum equation.
!! Note that the first step is performed with a first order Euler
!! scheme.
!!--------------------------------------------------------------------
!! Algorithm
!!----------
!!
!!  1. Fluid velocity at the particle position
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!
!!           uf@p[t(n)] = interp{uf[x[t(n)],t(n)]}
!!
!!
!!  2. Time-advancing the particle velocities
!! ++++++++++++++++++++++++++++++++++++++++++
!!
!!  up[x[t(n)],t(n)] --> up[x[t(n+1)],t(n+1)]
!!
!!        d up(t)    up(t) - uf@p(t)
!!        ------ = - ---------------- + gi
!!          dt            tau_p
!!
!!  3. Time-advancing the particles position:
!! ++++++++++++++++++++++++++++++++++++++++++
!!
!!  x[t(n)] --> x[t(n+1)]
!!
!!        d x(t)
!!        ------ = up(t)
!!          dt
!!
!!
!!  4. Apply periodic boundary conditions
!! ++++++++++++++++++++++++++++++++++++++
!!
!!
!!====================================================================

use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use FLUID_VARIABLE, only: UFLU, VFLU, WFLU
use SCALAR_VARIABLE, only: THETA
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!!- cycle number
integer, intent(in) :: NCYCLE


!- CPU loading
real(kind=8), dimension(NIG) :: NPMAXCPU
real(kind=8), dimension(NIG) :: NPMINCPU
real(kind=8), dimension(NIG) :: VARNPCPU


!!
real(kind=8) :: INVTAUP, REP
real(kind=8) :: VRNRM, VRX, VRY, VRZ


!!- dummy variable for subroutine argument
integer      :: IDUMMY

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!!- Index
integer :: I, J ,K, NP
!!====================================================================
real(kind=8) :: QP2,UP,VP,WP, UPUP,VPVP,WPWP,UPVP,UPWP,VPWP
!!====================================================================

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if



!!====================================================================
!! 0. INITIATION
!!====================================================================
!!- Initiation statistics
if(LEVEL0_STPAR) MEAN_PART(:,:) = ZERO


!!====================================================================
!! 1. Interpolation of fluid velocity at the particle positions
!!====================================================================

do J = 1, NIG


if(SOLVE_FLUID>0) then 
!!--------------------------------------------------------------------
!! 1.1. Save fluid velocity at particle position for 
!!      Time integration
!!--------------------------------------------------------------------
 PART(1:NPART_LOC(J),J)%UFAP_NM2 = PART(1:NPART_LOC(J),J)%UFAP_NM1
 PART(1:NPART_LOC(J),J)%UFAP_NM1 = PART(1:NPART_LOC(J),J)%UFAP

 PART(1:NPART_LOC(J),J)%VFAP_NM2 = PART(1:NPART_LOC(J),J)%VFAP_NM1
 PART(1:NPART_LOC(J),J)%VFAP_NM1 = PART(1:NPART_LOC(J),J)%VFAP

 PART(1:NPART_LOC(J),J)%WFAP_NM2 = PART(1:NPART_LOC(J),J)%WFAP_NM1
 PART(1:NPART_LOC(J),J)%WFAP_NM1 = PART(1:NPART_LOC(J),J)%WFAP

 if(SOLVE_SCALAR) then
  PART(1:NPART_LOC(J),J)%TFAP_NM2 = PART(1:NPART_LOC(J),J)%TFAP_NM1
  PART(1:NPART_LOC(J),J)%TFAP_NM1 = PART(1:NPART_LOC(J),J)%TFAP
 end if

!!=====================================================================
!!- Filtering the fluid velocity
!!=====================================================================
!! This is a specific procedure for investigating the effect of the 
!! Small scales on the dispersion of particle.
!! The cutoff wavenumber is defined in the fortran file FILTERING_FLUID
!!-------------------------------------------------------------------- 
 if(FILTERING) then
  call FILTERING_FLUID(J)
 else

!!--------------------------------------------------------------------
!! 1.2 x-velocity Interpolation
!!--------------------------------------------------------------------
 call INTERPH(INTERP_SCHEME,                &
		      XMESH, YMESH, ZMESH,  &
                       UFLU,                &
                NPART_LOC(J),               &
 	        PART(1:NPART_LOC(J),J)%XP,  &
                PART(1:NPART_LOC(J),J)%YP,  &
                PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%UFAP )

!!--------------------------------------------------------------------
!! 1.3 y-velocity Interpolation
!!--------------------------------------------------------------------
 call INTERPH(INTERP_SCHEME,                &
 		      XMESH, YMESH, ZMESH,  &
                       VFLU,                &
               NPART_LOC(J),                &
 	       PART(1:NPART_LOC(J),J)%XP,   &
               PART(1:NPART_LOC(J),J)%YP,   &
               PART(1:NPART_LOC(J),J)%ZP,   &
 	       PART(1:NPART_LOC(J),J)%VFAP  )


!!--------------------------------------------------------------------
!! 1.3 z-velocity Interpolation
!!--------------------------------------------------------------------
 call INTERPH(INTERP_SCHEME,                &
 	              XMESH, YMESH, ZMESH,  &
                       WFLU,                &
               NPART_LOC(J),                &
 	        PART(1:NPART_LOC(J),J)%XP,  &
                PART(1:NPART_LOC(J),J)%YP,  &
                PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%WFAP )

!!--------------------------------------------------------------------
!! 1.4 Scalar Interpolation
!!--------------------------------------------------------------------
 if(SOLVE_SCALAR) then
 
 call INTERPH(INTERP_SCHEME,                &
                      XMESH, YMESH, ZMESH,  &
                      THETA,                &
               NPART_LOC(J),                &
 	        PART(1:NPART_LOC(J),J)%XP,  &
		PART(1:NPART_LOC(J),J)%YP,  &
		PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%TFAP )
 end if




 end if !!- if(FILTERING)

else

  PART(1:NPART_LOC(J),J)%UFAP = ZERO
  PART(1:NPART_LOC(J),J)%VFAP = ZERO
  PART(1:NPART_LOC(J),J)%WFAP = ZERO

end if !- end if SOLVE_FLUID>1


end do !- end loop  J = 1, NIG




!!====================================================================
!! 2. Advance particle velocities
!!====================================================================

!!--------------------------------------------------------------------
!! 2.1 Save variables
!!--------------------------------------------------------------------
do J = 1, NIG

!!- Particle velocity
!!- n-1 --> n-2
 PART(1:NPART_LOC(J),J)%UP_NM2 = PART(1:NPART_LOC(J),J)%UP_NM1
 PART(1:NPART_LOC(J),J)%VP_NM2 = PART(1:NPART_LOC(J),J)%VP_NM1
 PART(1:NPART_LOC(J),J)%WP_NM2 = PART(1:NPART_LOC(J),J)%WP_NM1

!!- n --> n-1
 PART(1:NPART_LOC(J),J)%UP_NM1 = PART(1:NPART_LOC(J),J)%UP
 PART(1:NPART_LOC(J),J)%VP_NM1 = PART(1:NPART_LOC(J),J)%VP
 PART(1:NPART_LOC(J),J)%WP_NM1 = PART(1:NPART_LOC(J),J)%WP

!!- Particle scalar
 if(SOLVE_SCALAR) then
  PART(1:NPART_LOC(J),J)%TP_NM2 = PART(1:NPART_LOC(J),J)%TP_NM1
  PART(1:NPART_LOC(J),J)%TP_NM1 = PART(1:NPART_LOC(J),J)%TP
 end if


!!- Fluid elements
 if(PARTDEF(J)<2) then
  PART(1:NPART_LOC(J),J)%UP = PART(1:NPART_LOC(J),J)%UFAP
  PART(1:NPART_LOC(J),J)%VP = PART(1:NPART_LOC(J),J)%VFAP
  PART(1:NPART_LOC(J),J)%WP = PART(1:NPART_LOC(J),J)%WFAP
  
 if(SOLVE_SCALAR) PART(1:NPART_LOC(J),J)%TP = PART(1:NPART_LOC(J),J)%TFAP
  
 end if

end do !!- end loop: do J = 1, NIG


!!--------------------------------------------------------------------
!! 2.2. Update particle velocity
!!--------------------------------------------------------------------
!!- up^n --> up^n+1 = f(up^n, ufap^n)
if(SOLVE_FLUID>0) call ADV_PARTICLE_VELOCITY




!!====================================================================
!! 3. Advance particle positions
!!====================================================================
!!- xp^n --> xp^n+1 = f(up^n)
call ADV_PARTICLE_POSITION




!do J = 1, NIG
!UP = ZERO
!VP = ZERO
!WP = ZERO
!UPUP = ZERO
!VPVP = ZERO
!WPWP = ZERO
!UPVP = ZERO
!UPWP = ZERO
!VPWP = ZERO
!do I = 1,NPART_LOC(J)
!UP = UP +  PART(I,J)%UP
!VP = VP +  PART(I,J)%VP
!WP = WP +  PART(I,J)%WP
!end do
!UP = UP / NPART_LOC(J)
!VP = VP / NPART_LOC(J)
!WP = WP / NPART_LOC(J)
!do I = 1,NPART_LOC(J)
!UPUP = UPUP +  PART(I,J)%UP*PART(I,J)%UP
!VPVP = VPVP +  PART(I,J)%VP*PART(I,J)%VP
!WPWP = WPWP +  PART(I,J)%WP*PART(I,J)%WP
!UPVP = UPVP +  PART(I,J)%UP*PART(I,J)%VP
!UPWP = UPWP +  PART(I,J)%UP*PART(I,J)%WP
!VPWP = VPWP +  PART(I,J)%VP*PART(I,J)%WP
!end do
!UPUP = UPUP / NPART_LOC(J)
!VPVP = VPVP / NPART_LOC(J)
!WPWP = WPWP / NPART_LOC(J)
!UPVP = UPVP / NPART_LOC(J)
!UPWP = UPWP / NPART_LOC(J)
!VPWP = VPWP / NPART_LOC(J)
!write(*,*)'PART:',J
!write(*,*)'UP=',UP
!write(*,*)'VP=',VP
!write(*,*)'WP=',WP
!write(*,*)'UPUP=',UPUP
!write(*,*)'VPVP=',VPVP
!write(*,*)'WPWP=',WPWP
!write(*,*)'UPVP=',UPVP
!write(*,*)'UPWP=',UPWP
!write(*,*)'VPWP=',VPWP
!end do



!!====================================================================
!! 4. Periodic boundary conditions
!!====================================================================
call BOUNDARY_PARTICLE





if(mod(NCYCLE,NCYCLEMAX/20)==0) then

 if(SOLVE_PART) then
!!- check the number of particles
 do J = 1, NIG

  call IMAXCPU(NPART_LOC(J),IDUMMY)
  NPMAXCPU(J) = real(IDUMMY)

  call IMINCPU(NPART_LOC(J),IDUMMY)
  NPMINCPU(J) = real(IDUMMY)

  call ISUMCPU((NPART_LOC(J)-NPCPU_UNIF)**2,IDUMMY)
  VARNPCPU(J) = (IDUMMY/NPROC)**0.5

!!  if(MYID==0) write(*,10703)J,int(NPMAXCPU(J)),int(NPMINCPU(J))
 end do

 if (MYID==0) write(UNIT_INFO(4),10000)NCYCLE*DTIME,&
      (NPMAXCPU(J),NPMINCPU(J),VARNPCPU(J),J=1,NIG)

 end if
end if




!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(1) = CPU_PART(1) + TIME_END - TIME_START
end if




!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine PARTICLE_TRACKING

