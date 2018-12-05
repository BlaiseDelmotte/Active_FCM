!!====================================================================
!!
!! 
!! This routine build the Right Hand Side of Navier-Stokes equations
!!
!!
!!====================================================================

subroutine VORTICITY(ZETAX,ZETAY,ZETAZ)

!!====================================================================
!!
!!  1- compute the vorticity in Fourier space (VORTX,VORTY,VORTZ)
!!   and FFT-1 to get vorticity in physical space (ZETAX,ZETAY,ZETAZ)
!!
!!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: KX, KY, KZ, FILTER
use RHS_VARIABLES
use PARAM_PHYS,  only: VISC
use CHECK_CPU


use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- Vorticity

!- Physical Space
real(kind=8),intent(out),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAX
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAY
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAZ

!- Fourier Space
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTX
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTY
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: VORTZ


!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Index
integer :: I, J, K
!------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!!=================================================================
!! 1. Vorticity
!!=================================================================

!!-----------------------------------------------------------------
!! 1.1. Vorticity in Fourier space
!!-----------------------------------------------------------------
!! [vort_x]f = i*ky*[wf]fou - i*kz*[vf]fou

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   VORTX(I,J,K) = ICMPL*WFOU(I,J,K)*KY(J) - ICMPL*VFOU(I,J,K)*KZ(K)
   VORTY(I,J,K) = ICMPL*UFOU(I,J,K)*KZ(K) - ICMPL*WFOU(I,J,K)*KX(I)
   VORTZ(I,J,K) = ICMPL*VFOU(I,J,K)*KX(I) - ICMPL*UFOU(I,J,K)*KY(J)
  end do
 end do
end do



!!-----------------------------------------------------------------
!! 1.2. Vorticity back in physical space
!!-----------------------------------------------------------------
!- Synchronize all the process
 !call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
call P3DFFT_BTRAN_C2R(VORTX,ZETAX(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       
call P3DFFT_BTRAN_C2R(VORTY,ZETAY(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)
call P3DFFT_BTRAN_C2R(VORTZ,ZETAZ(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)


!!- CPU check
if(MYID == 0) then
 TIME_END=MPI_WTIME()
 CPU_FLUID(2) = CPU_FLUID(2) + TIME_END - TIME_START
end if

end subroutine VORTICITY
