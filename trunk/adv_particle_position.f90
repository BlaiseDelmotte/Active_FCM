!!=====================================================================
!!
!!   Time-advancing particle position
!!
!!=====================================================================

subroutine ADV_PARTICLE_POSITION

!!=====================================================================
!! The numerical scheme used for time-advancing is a 2nd order
!! Adams-Bashforth coupled with integrating factor.
!! 
!!  df(t)                      d[f(t).exp(at)]
!! ------ + af(t) + b = 0  <=> --------------- + b.exp(at) = 0
!!   dt                               dt
!!             
!! 2nd order Adams-Bashforth is written as:
!!
!!  df(t)                f^n+1 - f^n     3            1
!! ------ = phi(t)  ==> ------------- = ---.phi^n  -.---.phi^n-1 
!!   dt                       Dt         2            2
!!
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
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
!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Time step
real(kind=8) ::DTIME2, DTIME12

!- Index
integer :: I, J
!------------------------------------------------------------------



!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if




DTIME2 = DTIME / 2.
DTIME12 = DTIME / 12.

do J = 1, NIG

!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
if(FLAG_TSCHEME == 1) then


!!- Fluid elements or inertial particle
if(PARTDEF(J)>0) then

 do I = 1, NPART_LOC(J)
  PART(I,J)%XP = PART(I,J)%XP + DTIME*PART(I,J)%UP 
  PART(I,J)%YP = PART(I,J)%YP + DTIME*PART(I,J)%VP
  PART(I,J)%ZP = PART(I,J)%ZP + DTIME*PART(I,J)%WP
 end do

end if!!- If  PARTDEF(J)>0




!!====================================================================
!! 2. 2nd Order Adam-Bashford time integration
!!====================================================================
elseif(FLAG_TSCHEME == 2) then

!!- Fluid elements and inertial particles
if(PARTDEF(J)>0) then

 do I = 1, NPART_LOC(J)
  PART(I,J)%XP = PART(I,J)%XP &
               + DTIME2*(3.*PART(I,J)%UP-PART(I,J)%UP_NM1)

  PART(I,J)%YP = PART(I,J)%YP &
               + DTIME2*(3.*PART(I,J)%VP-PART(I,J)%VP_NM1)

  PART(I,J)%ZP = PART(I,J)%ZP &
               + DTIME2*(3.*PART(I,J)%WP-PART(I,J)%WP_NM1)
 end do

end if!!- If  PARTDEF(J)>0




!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================

elseif(FLAG_TSCHEME == 3 ) then

!!- Fluid elements and inertial particles
if(PARTDEF(J)>0) then


 do I = 1, NPART_LOC(J)

  PART(I,J)%XP = PART(I,J)%XP               &
         +  DTIME12*( 23.* PART(I,J)%UP	  &
                     -16.* PART(I,J)%UP_NM1 &
                     + 5.* PART(I,J)%UP_NM2 )

  PART(I,J)%YP = PART(I,J)%YP               &
         +  DTIME12*( 23.* PART(I,J)%VP	  &
                     -16.* PART(I,J)%VP_NM1 &
                     + 5.* PART(I,J)%VP_NM2 )

  PART(I,J)%ZP = PART(I,J)%ZP               &
         +  DTIME12*( 23.* PART(I,J)%WP	  &
                     -16.* PART(I,J)%WP_NM1 &
                     + 5.* PART(I,J)%WP_NM2 )

 end do

end if !!- If  PARTDEF(J)>0

end if !!- If: FLAG_TSCHEME



end do !!- Loop: 1, NIG




!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(2) = CPU_PART(2) + TIME_END - TIME_START
end if




end subroutine ADV_PARTICLE_POSITION
