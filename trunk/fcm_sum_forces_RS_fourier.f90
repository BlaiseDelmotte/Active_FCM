 !!====================================================================
!!
!! 
!!> @brief
!!> Add RF to forcing field in Fourier space
!!
!! Date :  28/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_SUM_FORCES_RS_FOURIER

!!====================================================================
!! Here the monopole forcing is distributed on the monopole 
!! Gaussian enveloppe.
!!====================================================================
!! Forcing: 
!!------------------------------
!!
!! TO DO : 
!!        1) If ellipsoid
!!------------------------------
!!------------------------------
!!====================================================================

use DNS_DIM	        !- Dimension
use PARAM_PHYS          !- Physical & numerical parameters
use FLUID_VARIABLE      !- Fluid velocity
use FCM_FORCING_VARIABLE
use FCM_PART_VARIABLE
use GEOMETRIC_VARIABLE 
use RHS_VARIABLES       !- Variable for cpu time checking
use WORK_ARRAYS
use CHECK_CPU	       

use MPI_STRUCTURES

use P3DFFT

implicit none

!!====================================================================
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------



real(kind=8) :: TEMP

!- Time control variable
real(kind=8) :: TIME_START, TIME_END
real(kind=8) :: TIME_START2

real(kind=8) :: MEASURE_START, MEASURE_END

integer :: I, J, K
!!
!!====================================================================

!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
 TIME_START2 = TIME_START
end if


!~ call cpu_time(MEASURE_START)
!~ call P3DFFT_FTRAN_R2C(FCM_FORCING_X(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),RHS_UFOU(:,:,:,1),FFTFLAG)    !- x-component of fluid velocity
!~ call cpu_time(MEASURE_END)
!~ print*, '("Time = ",f6.9," seconds.")',MEASURE_END - MEASURE_START
!~ 
!~ read(*,*)

!- Projection of forcing terms in Fourier space
TMPPHY(:,:,:) = FCM_FORCING_X(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,RHS_UFOU(:,:,:,1),FFTFLAG)    !- x-component of fluid velocity

TMPPHY(:,:,:) = FCM_FORCING_Y(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,RHS_VFOU(:,:,:,1),FFTFLAG)    !- y-component of fluid velocity

TMPPHY(:,:,:) = FCM_FORCING_Z(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,RHS_WFOU(:,:,:,1),FFTFLAG)    !- z-component of fluid velocity


!- Projection of RS in Fourier space
TMPPHY(:,:,:) = FCM_RS_XX(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_XX_FOU(:,:,:),FFTFLAG)    !- xx-component of fluid velocity

TMPPHY(:,:,:) = FCM_RS_YY(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_YY_FOU(:,:,:),FFTFLAG)    !- yy-component of fluid velocity

TMPPHY(:,:,:) = FCM_RS_ZZ(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_ZZ_FOU(:,:,:),FFTFLAG)    !- zz-component of fluid velocity

TMPPHY(:,:,:) = FCM_RS_XY(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_XY_FOU(:,:,:),FFTFLAG)    !- xy-component of fluid velocity

TMPPHY(:,:,:) = FCM_RS_XZ(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_XZ_FOU(:,:,:),FFTFLAG)    !- xz-component of fluid velocity

TMPPHY(:,:,:) = FCM_RS_YZ(:,ISTART(2):IEND(2),ISTART(3):IEND(3))
call P3DFFT_FTRAN_R2C(TMPPHY,FCM_RS_YZ_FOU(:,:,:),FFTFLAG)    !- yz-component of fluid velocity

!~ call cpu_time(MEASURE_START)


!- Remove continuous forcing values
if(MYID == 0) then
 RHS_UFOU(1,1,1,1) = cmplx(ZERO,ZERO)
 RHS_VFOU(1,1,1,1) = cmplx(ZERO,ZERO)
 RHS_WFOU(1,1,1,1) = cmplx(ZERO,ZERO)

 FCM_RS_XX_FOU(1,1,1) = cmplx(ZERO,ZERO)
 FCM_RS_YY_FOU(1,1,1) = cmplx(ZERO,ZERO)
 FCM_RS_ZZ_FOU(1,1,1) = cmplx(ZERO,ZERO)
 FCM_RS_XY_FOU(1,1,1) = cmplx(ZERO,ZERO)
 FCM_RS_XZ_FOU(1,1,1) = cmplx(ZERO,ZERO)
 FCM_RS_YZ_FOU(1,1,1) = cmplx(ZERO,ZERO)
end if


!- Solve Stokes equation in Fourier Space
do K = FSTART(3),FEND(3)
 do J = FSTART(2),FEND(2)
  do I = FSTART(1),FEND(1)
 

    !- Add the div(RS) to the forcing terms
    RHS_UFOU(I,J,K,1) = RHS_UFOU(I,J,K,1) &
                      + cmplx(0.0,KX(I))*FCM_RS_XX_FOU(I,J,K) &
                      + cmplx(0.0,KY(J))*FCM_RS_XY_FOU(I,J,K) &
                      + cmplx(0.0,KZ(K))*FCM_RS_XZ_FOU(I,J,K)

    RHS_VFOU(I,J,K,1) = RHS_VFOU(I,J,K,1) &
                      + cmplx(0.0,KX(I))*FCM_RS_XY_FOU(I,J,K) &
                      + cmplx(0.0,KY(J))*FCM_RS_YY_FOU(I,J,K) &
                      + cmplx(0.0,KZ(K))*FCM_RS_YZ_FOU(I,J,K)

    RHS_WFOU(I,J,K,1) = RHS_WFOU(I,J,K,1) &
                      + cmplx(0.0,KX(I))*FCM_RS_XZ_FOU(I,J,K) &
                      + cmplx(0.0,KY(J))*FCM_RS_YZ_FOU(I,J,K) &
                      + cmplx(0.0,KZ(K))*FCM_RS_ZZ_FOU(I,J,K)

  end do
 end do
end do

!- FFT Normalization
RHS_UFOU(:,:,:,1) = RHS_UFOU(:,:,:,1)*FACTOR
RHS_VFOU(:,:,:,1) = RHS_VFOU(:,:,:,1)*FACTOR
RHS_WFOU(:,:,:,1) = RHS_WFOU(:,:,:,1)*FACTOR

!- Remove continuous velocity values
if(MYID == 0) then
 RHS_UFOU(1,1,1,1) =cmplx(ZERO,ZERO)
 RHS_VFOU(1,1,1,1) =cmplx(ZERO,ZERO)
 RHS_WFOU(1,1,1,1) =cmplx(ZERO,ZERO)
end if

!~ call cpu_time(MEASURE_END)
!~  print*, '("Time = ",f6.9," seconds.")',MEASURE_END - MEASURE_START 
!~ read(*,*)

!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
end if

!- Back in physical space
TMPFOU = RHS_UFOU(:,:,:,1)
call P3DFFT_BTRAN_C2R(TMPFOU,FCM_FORCING_X(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)   

TMPFOU = RHS_VFOU(:,:,:,1)
call P3DFFT_BTRAN_C2R(TMPFOU,FCM_FORCING_Y(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)

TMPFOU = RHS_WFOU(:,:,:,1)
call P3DFFT_BTRAN_C2R(TMPFOU,FCM_FORCING_Z(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG) 

!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()

 CPU_FLUID(1) = CPU_FLUID(1) + TIME_END - TIME_START2
end if
   
end subroutine FCM_SUM_FORCES_RS_FOURIER
