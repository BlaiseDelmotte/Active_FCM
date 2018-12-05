 !!====================================================================
!!
!! 
!!> @brief
!!> Routine prescribing the time-dependent swimming gait 
!!
!! Date :  08/04/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_TIME_STROKE(NCYCLE)


use DNS_DIM
use PARAM_PHYS
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE

implicit none

!- Arguments
integer, intent(in) :: NCYCLE

!- Rad max
real(kind=8) :: RADMAX

!- Stroke frequency and pulsation
real(kind=8) :: FREQ, PULSATION

!- Amplitudes Fourier Modes Vsw
real(kind=8), dimension(FCM_NSWIM(3)) :: AMP0_VSW,  AMP1C_VSW, AMP1S_VSW, &
                           AMP2C_VSW, AMP2S_VSW


!- Amplitudes Fourier Modes B2
real(kind=8), dimension(FCM_NSWIM(3)) :: AMP0_B2,  AMP1C_B2, AMP1S_B2, &
                          AMP2C_B2, AMP2S_B2, &
                          AMP3C_B2, AMP3S_B2
                          
!- Phases Fourier Modes B2  
real(kind=8) :: PHASE1C_B2, PHASE1S_B2, &
                PHASE2C_B2, PHASE2S_B2, &
                PHASE3C_B2, PHASE3S_B2                    

if (FCM_NELLIPSOID>0) then
 RADMAX = maxval(FCM_ELLIPSOID_RADP)
endif

!- Swimming frequency
FREQ = 1.0/0.0189
PULSATION = FREQ*TWOPI*real(NCYCLE)*DTIME

!- Swimming modes (Guasto Data/rad_Guasto*rad_FCM)
AMP0_VSW = 49.5400* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP1S_VSW = 122.24* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP1C_VSW = - 34.7240* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP2S_VSW = - 8.4400* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP2C_VSW = - 12.7600* FCM_SPHERE_RADP(1:FCM_NSWIM(3))

!- Stresslet modes (Guasto Data/rad_Guasto*rad_FCM)
AMP0_B2 = 4.5347* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP1S_B2 = 91.4529* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP1C_B2 = 64.0530* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP2S_B2 = - 92.6420* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP2C_B2 = - 84.7192* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP3S_B2 = - 5.9849* FCM_SPHERE_RADP(1:FCM_NSWIM(3))
AMP3C_B2 = - 10.4545* FCM_SPHERE_RADP(1:FCM_NSWIM(3))

PHASE1S_B2 = 0.1666
PHASE1C_B2 = 1.7373
PHASE2S_B2 = 2.0054
PHASE2C_B2 = 3.5761
PHASE3S_B2 = - 1.7125
PHASE3C_B2 = - 0.9154


FCM_VSW_TIME   = AMP1S_VSW * dsin(PULSATION + FCM_PHASE_SHIFT) &
               + AMP1C_VSW * dcos(PULSATION + FCM_PHASE_SHIFT) &
               + AMP2S_VSW * dsin(2.0*PULSATION + FCM_PHASE_SHIFT) &
               + AMP2C_VSW * dcos(2.0*PULSATION + FCM_PHASE_SHIFT) &
               + AMP0_VSW
            
FCM_B2_TIME  = AMP1S_B2 * dsin(PULSATION + PHASE1S_B2 + FCM_PHASE_SHIFT) &
             + AMP1C_B2 * dcos(PULSATION + PHASE1C_B2 + FCM_PHASE_SHIFT) &
             + AMP2S_B2 * dsin(2.0*PULSATION + PHASE2S_B2 + FCM_PHASE_SHIFT) &
             + AMP2C_B2 * dcos(2.0*PULSATION + PHASE2C_B2 + FCM_PHASE_SHIFT) &
             + AMP3S_B2 * dsin(3.0*PULSATION + PHASE3S_B2 + FCM_PHASE_SHIFT) &
             + AMP3C_B2 * dcos(3.0*PULSATION + PHASE3C_B2 + FCM_PHASE_SHIFT) &
             + AMP0_B2  


end subroutine FCM_TIME_STROKE
