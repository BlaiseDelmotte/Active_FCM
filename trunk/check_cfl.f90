!!====================================================================
!! 
!! Check CFL condition for time-integration stability
!!
!!====================================================================

subroutine CHECK_CFL(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: DX, DY, DZ
use PARAM_PHYS,         only: DTIME, VISC, FOUT0

implicit none



!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
integer, intent(in) :: NCYCLE

!- Maximum velocity
real(kind=8) :: UMAX, VMAX, WMAX

!- Time step according to CFL condition
real(kind=8) :: DTIME_FOURIER
real(kind=8) :: DTIME_FOURIER2
real(kind=8) :: COEF_CFL
real(kind=8) :: CFL
real(kind=8) :: DTIME_COURANT

integer :: IFLAG1
!-----------------------------------------------------------------

CFL = 0.35

if(NPROC>1) then

 call RMAXCPU(abs(UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),UMAX)
 call RMAXCPU(abs(VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),VMAX)
 call RMAXCPU(abs(WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),WMAX)

else

 UMAX=maxval(UFLU)
 VMAX=maxval(VFLU)
 WMAX=maxval(WFLU)

end if


COEF_CFL = ZERO
if(UMAX> ZERO) COEF_CFL=max(COEF_CFL,(DX/UMAX)**(4./3.))
if(VMAX> ZERO) COEF_CFL=max(COEF_CFL,(DY/VMAX)**(4./3.))
if(WMAX> ZERO) COEF_CFL=max(COEF_CFL,(DZ/WMAX)**(4./3.))


DTIME_FOURIER = 2**(2./3.)*CFL**(1./3.)*COEF_CFL

DTIME_COURANT = min(0.25*DX/UMAX,0.25*DY/VMAX)
DTIME_COURANT = min(DTIME_COURANT,0.25*DZ/WMAX)

DTIME_FOURIER2 = min(VISC/DX**2,VISC/DY**2)
DTIME_FOURIER2 = min(DTIME_FOURIER2,VISC/DZ**2)



!!======================================================================
!! Print in run.info
!!======================================================================
if((mod(NCYCLE,FOUT0) == 0).and.(MYID==0)) then
!if(MYID==0)) then
 IFLAG1 = 5

 write(UNIT_INFO(IFLAG1),*)'-----------------------'
 write(UNIT_INFO(IFLAG1),10602)'max(|uf|) =',UMAX,' Courant (<0.25) =',UMAX*DTIME/DX
 write(UNIT_INFO(IFLAG1),10602)'max(|vf|) =',VMAX,' Courant (<0.25) =',VMAX*DTIME/DY
 write(UNIT_INFO(IFLAG1),10602)'max(|wf|) =',WMAX,' Courant (<0.25) =',WMAX*DTIME/DZ
 write(UNIT_INFO(IFLAG1),10601)'Dt = ',DTIME
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(Courant) =',DTIME/DTIME_COURANT
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(Fourier) =',DTIME/DTIME_FOURIER2
 write(UNIT_INFO(IFLAG1),10601)'Dt/Dt(Stab   ) =',DTIME/DTIME_FOURIER
 write(UNIT_INFO(IFLAG1),*)


end if

!-----------------------------------------------------------------------
10601 format (2x,A,1x,E13.6)
10602 format (2x,A,1x,E13.6,A,2x,E13.6)

end subroutine CHECK_CFL
