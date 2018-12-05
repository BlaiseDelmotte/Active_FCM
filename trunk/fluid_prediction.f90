!!====================================================================
!!
!! Direct Numerical Simulation of turbulent fluid flow
!!
!!====================================================================

subroutine FLUID_PREDICTION(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM             !- Dimension
use PARAM_PHYS          !- Physical & numerical parameters
use FORCING             !- Forcing
use FLUID_VARIABLE      !- Fluid velocity
use SCALAR_VARIABLE    
use GEOMETRIC_VARIABLE 
use RHS_VARIABLES       !- Variable for cpu time checking
use STATISTICS         
use WORK_ARRAYS
use CHECK_CPU       

use MPI_STRUCTURES

use P3DFFT

implicit none

!!====================================================================
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
!- cycle number
integer, intent(in) :: NCYCLE

!!- Shifted coefficient
!real(kind=8) :: DELTA
!- Complex velocity shifted
!double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: USHIFT
!double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VSHIFT
!double complex, dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: WSHIFT


!!
!- Time control variable
real(kind=8) :: TIME_START, TIME_END
real(kind=8) :: TIME_START2

!!real(kind=8) :: SHIFT
double complex :: SHIFT

integer :: I, J, K
!!
!!====================================================================

!!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
 TIME_START2 = TIME_START
end if







TMPPHY(:,:,:) = UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,UFOU,FFTFLAG) !- x-component of fluid velocity

TMPPHY(:,:,:) = VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,VFOU,FFTFLAG) !- y-component of fluid velocity

TMPPHY(:,:,:) = WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,WFOU,FFTFLAG) !- z-component of fluid velocity





!- FFT Normalization
UFOU(:,:,:) = UFOU(:,:,:)*FACTOR
VFOU(:,:,:) = VFOU(:,:,:)*FACTOR
WFOU(:,:,:) = WFOU(:,:,:)*FACTOR



RHS_UFOU(:,:,:,TN) = cmplx(ZERO,ZERO)
RHS_VFOU(:,:,:,TN) = cmplx(ZERO,ZERO)
RHS_WFOU(:,:,:,TN) = cmplx(ZERO,ZERO)


if(SOLVE_SCALAR) then
 TMPPHY(:,:,:) = THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
 call P3DFFT_FTRAN_R2C(TMPPHY,THETAFOU,FFTFLAG)
 
 THETAFOU(:,:,:) = THETAFOU(:,:,:)*FACTOR
end if



if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(8) = CPU_FLUID(8) +  TIME_END - TIME_START
end if







!!- Build Right-Hand-Side of Navier-Stokes equations
!!call BUILD_RHS_ARG(UFOU, VFOU, WFOU, RHS_UFOU, RHS_VFOU, RHS_WFOU)
!!
!!DELTA = DX*0.5 !- Shift
!!
!!do K = FSTART(3), FEND(3)
!! do J = FSTART(2), FEND(2)
!!  do I = FSTART(1), FEND(1)
!!  
!! SHIFT = exp(ICMPL*KX(I)*DELTA)*exp(ICMPL*KY(J)*DELTA)*exp(ICMPL*KZ(K)*DELTA)
!!
!!! USHIFT(I,J,K) = UFOU(I,J,K)*SHIFT
!!! VSHIFT(I,J,K) = VFOU(I,J,K)*SHIFT
!!! WSHIFT(I,J,K) = WFOU(I,J,K)*SHIFT
!! USHIFT(I,J,K) = UFOU(I,J,K)*exp(ICMPL*KX(I)*DELTA)*exp(ICMPL*KY(J)*DELTA)*exp(ICMPL*KZ(K)*DELTA)
!! VSHIFT(I,J,K) = VFOU(I,J,K)*exp(ICMPL*KX(I)*DELTA)*exp(ICMPL*KY(J)*DELTA)*exp(ICMPL*KZ(K)*DELTA)
!! WSHIFT(I,J,K) = WFOU(I,J,K)*exp(ICMPL*KX(I)*DELTA)*exp(ICMPL*KY(J)*DELTA)*exp(ICMPL*KZ(K)*DELTA)
!!
!!  end do 
!! end do 
!!end do 
!!
!!call BUILD_RHS_ARG(USHIFT, VSHIFT, WSHIFT, RHS_UFOU, RHS_VFOU, RHS_WFOU)
!!
!!RHS_UFOU = RHS_UFOU*0.5
!!RHS_VFOU = RHS_VFOU*0.5
!!RHS_WFOU = RHS_WFOU*0.5


!!- Add Forcing
!!if(STEADY) call ADD_FORCING

if(STEADY) call ADD_FORCING2


!!- Build Right-Hand-Side of Navier-Stokes equations
call BUILD_RHS_VORT
!!call BUILD_RHS_VORT_OMEM





!!====================================================================
!! X. Time-integration of fluid momentum
!!====================================================================
!! Full explicit Adams-Bashforth is used for the time-integration of
!! the fluid momentum equation.
!! Note that the first step is performed with a first order Euler
!! scheme.
!!--------------------------------------------------------------------
call ADV_FLUID



!- Projection on soloneidal basis
call PROJ_DIVFREE




!!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
end if





!- Back in physical space
TMPFOU = UFOU
call P3DFFT_BTRAN_C2R(TMPFOU,UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       

TMPFOU = VFOU 
call P3DFFT_BTRAN_C2R(TMPFOU,VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)

TMPFOU = WFOU
call P3DFFT_BTRAN_C2R(TMPFOU,WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       




!!- Check CFL limitation
!if(mod(NCYCLE,FOUT0) == 0) call CHECK_CFL
call CHECK_CFL(NCYCLE)




!!====================================================================
!! X. Passive scalar equation
!!====================================================================
if(SOLVE_SCALAR) then

 call ADV_SCALAR

 TMPFOU = THETAFOU
 call P3DFFT_BTRAN_C2R(TMPFOU,THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       

end if







!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(8) = CPU_FLUID(8) + TIME_END - TIME_START

 CPU_FLUID(1) = CPU_FLUID(1) + TIME_END - TIME_START2
end if




end subroutine FLUID_PREDICTION
