!!=====================================================================
!!
!!   Time-advancing passive scalar in Fourier Space
!!
!!=====================================================================

subroutine ADV_SCALAR

!!=====================================================================
!! The numerical scheme used for time-advancing is a 2nd order
!! Adams-Bashforth. As the scheme is explicit user has to take
!! care about the timestep.
!!
!!  df(t)                f^n+1 - f^n     3            1
!! ------ = phi(t)  ==> ------------- = ---.phi^n  -.---.phi^n-1 
!!   dt                       Dt         2            2
!!
!!
!!=====================================================================

use DNS_DIM
use FLUID_VARIABLE
use SCALAR_VARIABLE
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use WORK_ARRAYS
use CHECK_CPU

use P3DFFT

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Squared wavenumber
real(kind=8) :: K2, DTIME2, DTIME12

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

real(kind=8) :: KAPPA2, INTFCT_SCL

integer :: I, J , K
!---------------------------------------------------------------------

!!- CPU check
!if(MYID == 0) then
! TIME_START=MPI_WTIME()
!end if


!!====================================================================
!! 1. Build Right-Hand-Side of transport equation
!!====================================================================
RHS_SCL(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),TN) = cmplx(ZERO,ZERO)



!!--------------------------------------------------------------------
!!- Term: theta*u
!!--------------------------------------------------------------------
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = THETA(I,J,K)*UFLU(I,J,K)
  end do
 end do
end do


!- Synchronize all the process
! call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Fourier space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  RHS_SCL(I,J,K,TN) = RHS_SCL(I,J,K,TN) -ICMPL*KX(I)*TMPFOU(I,J,K) !!* FILTER(I,J,K)
  end do
 end do
end do

!!--------------------------------------------------------------------
!!- Term: theta*v
!!--------------------------------------------------------------------
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = THETA(I,J,K)*VFLU(I,J,K)
  end do
 end do
end do


!- Synchronize all the process
! call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Fourier space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  RHS_SCL(I,J,K,TN) = RHS_SCL(I,J,K,TN) -ICMPL*KY(J)*TMPFOU(I,J,K) !!* FILTER(I,J,K)
  end do
 end do
end do

!!--------------------------------------------------------------------
!!- Term: theta*w
!!--------------------------------------------------------------------
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = THETA(I,J,K)*WFLU(I,J,K)
  end do
 end do
end do


!- Synchronize all the process
! call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Fourier space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  RHS_SCL(I,J,K,TN) = RHS_SCL(I,J,K,TN) -ICMPL*KZ(K)*TMPFOU(I,J,K)!  !*FILTER(I,J,K)
  end do
 end do
end do


!!--------------------------------------------------------------------
!!- Forcing by mean gradient
!!--------------------------------------------------------------------
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  RHS_SCL(I,J,K,TN) = RHS_SCL(I,J,K,TN)      &
		    - GRAD_SCL*VFOU(I,J,K)
  end do
 end do
end do



!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
if(FLAG_TSCHEME == 1) then

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2.D0 + KY(J)**2.D0 + KZ(K)**2.D0
   
   if (KAPPA2 <= ZERO) KAPPA2 = INFINITY

   INTFCT_SCL = exp(-DIFF_SCL*KAPPA2*DTIME)

   THETAFOU(I,J,K) = INTFCT_SCL*( THETAFOU(I,J,K) + DTIME*RHS_SCL(I,J,K,TN) )

  end do
 end do
end do


!!====================================================================
!! 2. 2nd Order Adam-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 2) then

DTIME2 = DTIME/2.D0

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2.D0 + KY(J)**2.D0 + KZ(K)**2.D0
   
   if (KAPPA2 <= ZERO) KAPPA2 = INFINITY

   INTFCT_SCL   = exp(-DIFF_SCL*KAPPA2*DTIME)


   THETAFOU(I,J,K)=INTFCT_SCL*(THETAFOU(I,J,K) &
            + DTIME2*(3.D0*RHS_SCL(I,J,K,TN)-INTFCT_SCL*RHS_SCL(I,J,K,TNM1)))

  end do
 end do
end do

!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================
elseif(FLAG_TSCHEME == 3 ) then

DTIME12 = DTIME/12.

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KAPPA2 = KX(I)**2.D0 + KY(J)**2.D0 + KZ(K)**2.D0
   
   if (KAPPA2 <= ZERO) KAPPA2 = INFINITY

   INTFCT_SCL   = exp(-DIFF_SCL*KAPPA2*DTIME)
  
   THETAFOU(I,J,K) = INTFCT_SCL*( THETAFOU(I,J,K) +  &
          DTIME12*( 23.D0*                   RHS_SCL(I,J,K,TN  ) &
                   -16.D0* INTFCT_SCL       *RHS_SCL(I,J,K,TNM1) &
                   + 5.D0*(INTFCT_SCL**2.D0)*RHS_SCL(I,J,K,TNM2) )  )

  end do
 end do
end do


end if



!!- Zeroing first wavenumber
!!THETAFOU(1,1,1) = cmplx(ZERO,ZERO)


!!- CPU check
!if(MYID == 0) then
! TIME_END=MPI_WTIME()
! CPU_FLUID(2) = TIME_END - TIME_START
!end if


end subroutine ADV_SCALAR
