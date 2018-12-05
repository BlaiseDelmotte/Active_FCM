!!====================================================================
!!
!! 
!! 
!!> @brief
!!> This routine computes the rate of strain of fluid phase
!!
!! Date :  18/06/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_FLUID_RATE_OF_STRAIN( FCM_FLUID_E11, &
									 FCM_FLUID_E12, &
									 FCM_FLUID_E13, &
									 FCM_FLUID_E22, &
									 FCM_FLUID_E23 )

!!====================================================================
!!
!!  1- compute the ROS in Fourier space 
!!   and FFT-1 to get ROS in physical space 
!!
!!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: KX, KY, KZ, FILTER
use FCM_FORCING_VARIABLE, only: FCM_SHEAR, FCM_CONSIDER_ROS_SHEAR
use RHS_VARIABLES
use PARAM_PHYS,  only: VISC
use CHECK_CPU


use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- ROS

!- Physical Space
real(kind=8),intent(out),   &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: FCM_FLUID_E11
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: FCM_FLUID_E12
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: FCM_FLUID_E13
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: FCM_FLUID_E22
real(kind=8),intent(out),    &
   dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: FCM_FLUID_E23

!- Fourier Space
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: FCM_FLUID_E11FOU
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: FCM_FLUID_E12FOU
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: FCM_FLUID_E13FOU
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: FCM_FLUID_E22FOU
double complex,  &
   dimension(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)) :: FCM_FLUID_E23FOU


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
!! 1. Fluid ROS
!!=================================================================

!!-----------------------------------------------------------------
!! 1.1. ROS in Fourier space
!!-----------------------------------------------------------------
!! [E_xx]f = i*ky*[wf]fou - i*kz*[vf]fou

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   FCM_FLUID_E11FOU(I,J,K) = ICMPL*UFOU(I,J,K)*KX(I)
   FCM_FLUID_E12FOU(I,J,K) = 0.5*( ICMPL*UFOU(I,J,K)*KY(J) + ICMPL*VFOU(I,J,K)*KX(I) )
   FCM_FLUID_E13FOU(I,J,K) = 0.5*( ICMPL*UFOU(I,J,K)*KZ(K) + ICMPL*WFOU(I,J,K)*KX(I) )
   FCM_FLUID_E22FOU(I,J,K) = ICMPL*VFOU(I,J,K)*KY(J)
   FCM_FLUID_E23FOU(I,J,K) = 0.5*( ICMPL*VFOU(I,J,K)*KZ(K) + ICMPL*WFOU(I,J,K)*KY(J) )
  end do
 end do
end do




!!-----------------------------------------------------------------
!! 1.2. Vorticity back in physical space
!!-----------------------------------------------------------------
!- Synchronize all the process
 !call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
call P3DFFT_BTRAN_C2R(FCM_FLUID_E11FOU,FCM_FLUID_E11(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)       
call P3DFFT_BTRAN_C2R(FCM_FLUID_E12FOU,FCM_FLUID_E12(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)     
call P3DFFT_BTRAN_C2R(FCM_FLUID_E13FOU,FCM_FLUID_E13(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)     
call P3DFFT_BTRAN_C2R(FCM_FLUID_E22FOU,FCM_FLUID_E22(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)     
call P3DFFT_BTRAN_C2R(FCM_FLUID_E23FOU,FCM_FLUID_E23(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)     


if (FCM_CONSIDER_ROS_SHEAR==1) then
 !- Add shear intensity 
 FCM_FLUID_E13 = FCM_FLUID_E13 + 0.5 * FCM_SHEAR 
end if

!!- CPU check
if(MYID == 0) then
 TIME_END=MPI_WTIME()
 CPU_FLUID(2) = CPU_FLUID(2) + TIME_END - TIME_START
end if

end subroutine FCM_FLUID_RATE_OF_STRAIN
