!!=====================================================================
!!
!!   Time-advancing fluid solution in Fourier Space
!!
!!=====================================================================

subroutine ADV_FLUID

!!=====================================================================
!! The numerical scheme used for time-advancing is a 2nd order
!! Adams-Bashforth. As the scheme is explicit user has to take
!! care about the timestep.
!!
!!  df(t)                f^n+1 - f^n     3            1
!! ------ = phi(t)  ==> ------------- = ---.phi^n  -.---.phi^n-1 
!!   dt                       Dt         2            2
!!
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
!! The size of INTINTCOEF is NX/P1*NY/P2*NZ this array is filled one time
!! in meshing.f90 fortran file. To reduce the memory it is possible to
!! remove this array and to compute the integrating INTCOEF at each time 
!! step. The choice must be the balance between memory and computational
!! cost ...
!!
!!=====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: KX, KY, KZ
use RHS_VARIABLES
use PARAM_PHYS, only: VISC, DTIME 
use CHECK_CPU

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Squared wavenumber
real(kind=8) :: K2, DTIME2, DTIME12, INTCOEF

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J ,K
!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if

!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adams-Bashforth 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
if(FLAG_TSCHEME == 1) then


do K = FSTART(3),FEND(3)
do J = FSTART(2),FEND(2)
do I = FSTART(1),FEND(1)

 K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2

 if (K2 <= ZERO) K2 = INFINITY

 INTCOEF = exp(-DTIME*VISC*K2)

 UFOU(I,J,K) = INTCOEF*( UFOU(I,J,K) + DTIME*RHS_UFOU(I,J,K,TN) )
 VFOU(I,J,K) = INTCOEF*( VFOU(I,J,K) + DTIME*RHS_VFOU(I,J,K,TN) )
 WFOU(I,J,K) = INTCOEF*( WFOU(I,J,K) + DTIME*RHS_WFOU(I,J,K,TN) )



enddo
enddo
enddo


!!====================================================================
!! 2. 2nd Order Adams-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 2) then

DTIME2 = DTIME/2.

do K = FSTART(3),FEND(3)
do J = FSTART(2),FEND(2)
do I = FSTART(1),FEND(1)

 K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2

 if (K2 <= ZERO) K2 = INFINITY

 INTCOEF = exp(-DTIME*VISC*K2)

 UFOU(I,J,K)=INTCOEF*(UFOU(I,J,K) &
        + DTIME2*(3.*RHS_UFOU(I,J,K,TN)-INTCOEF*RHS_UFOU(I,J,K,TNM1)))

 VFOU(I,J,K)=INTCOEF*(VFOU(I,J,K) &
        + DTIME2*(3.*RHS_VFOU(I,J,K,TN)-INTCOEF*RHS_VFOU(I,J,K,TNM1)))

 WFOU(I,J,K)=INTCOEF*(WFOU(I,J,K) &
        + DTIME2*(3.*RHS_WFOU(I,J,K,TN)-INTCOEF*RHS_WFOU(I,J,K,TNM1)))

enddo
enddo
enddo



!!====================================================================
!! 3. 3rd Order Adams-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 3 ) then

DTIME12 = DTIME/12.

do K = FSTART(3),FEND(3)
do J = FSTART(2),FEND(2)
do I = FSTART(1),FEND(1)

 K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2

 if (K2 <= ZERO) K2 = INFINITY

 INTCOEF = exp(-DTIME*VISC*K2)


UFOU(I,J,K) = INTCOEF*( UFOU(I,J,K) +                     &
          DTIME12*( 23.*             RHS_UFOU(I,J,K,TN)   &
                   -16.* INTCOEF    *RHS_UFOU(I,J,K,TNM1) &
                   + 5.* INTCOEF**2 *RHS_UFOU(I,J,K,TNM2) )  )


VFOU(I,J,K) = INTCOEF*( VFOU(I,J,K) +                     &
          DTIME12*( 23.*             RHS_VFOU(I,J,K,TN)   &
                   -16.* INTCOEF    *RHS_VFOU(I,J,K,TNM1) &
                   + 5.* INTCOEF**2 *RHS_VFOU(I,J,K,TNM2) )  )


WFOU(I,J,K) = INTCOEF*( WFOU(I,J,K) +                     &
          DTIME12*( 23.*             RHS_WFOU(I,J,K,TN)   &
                   -16.* INTCOEF    *RHS_WFOU(I,J,K,TNM1) &
                   + 5.* INTCOEF**2 *RHS_WFOU(I,J,K,TNM2) )  )

enddo
enddo
enddo


end if


!!- CPU check
if(MYID == 0) then
 TIME_END=MPI_WTIME()
 CPU_FLUID(3) = CPU_FLUID(3) + TIME_END - TIME_START
end if


end subroutine ADV_FLUID
