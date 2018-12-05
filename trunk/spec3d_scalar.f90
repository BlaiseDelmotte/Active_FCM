!!====================================================================
!!
!!     3-dimensional turbulent spectrum and statistics
!!                       of scalar field
!!
!!====================================================================

subroutine SPEC3D_SCALAR(NOUT)

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM	       !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use SCALAR_VARIABLE    !- Fluid velocity
use GEOMETRIC_VARIABLE !- Mesh & Wavenumbers
use STATISTICS, only:EPS_FLU_SPEC

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------------------------------------------------------------

!- 
integer, intent(in) :: NOUT

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- 3d spectrum
real(kind=8), dimension(:), allocatable :: SPECTRUM, SPECLOC

!- File name
character (len=30) :: FILESPEC

!- 
real(kind=8) :: KXGLOB, KAPPA2, KAPPA

!- Kinetic energy
real(kind=8) :: KF_SCL

!- Dissipation
real(kind=8) :: EPS_SCL

!- Eddy Lifetime
real(kind=8) :: TMIX_SCL

!- Integral length scale
real(kind=8) :: LF_SCL

!- Kolmogorov length scale
real(kind=8) :: ETA_SCL

!- Kolmogorov Time scale
real(kind=8) :: TIME_SCL

!- CFL For scalar resolution
real(kind=8) :: CFL_SCL


!- Index
integer :: I, J, K, N, KMAX, IK
!------------------------------------------------------------------


!!====================================================================
!! 1. initiation
!!====================================================================
KMAX = NX/2 + 1

allocate(SPECTRUM(KMAX))
allocate(SPECLOC(KMAX))


SPECTRUM(:) = 0.
SPECLOC(:) = 0.


!!====================================================================
!! 2. 3d spectrum
!!====================================================================
do K = FSTART(3),FEND(3)
 do J = FSTART(2),FEND(2)
  do I = FSTART(1),FEND(1)

    KAPPA2 = (KX(I)/KFIRST)**2 + (KY(J)/KFIRST)**2 + (KZ(K)/KFIRST)**2

    IK = sqrt(KAPPA2) + 1
    if(IK>KMAX) IK = NX-IK

    SPECLOC(IK) = SPECLOC(IK) + real(THETAFOU(I,J,K)*conjg(THETAFOU(I,J,K)))

  enddo
 enddo
enddo


call MPI_ALLREDUCE(SPECLOC,SPECTRUM, &
              KMAX,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!--------------------------------------------------------------------
!! 1.2. Normalisation and printing 
!!--------------------------------------------------------------------
do N = 1, KMAX
 SPECTRUM(N) = SPECTRUM(N)/KFIRST
end do


!!====================================================================
!! 3. Statistics
!!====================================================================

!!--------------------------------------------------------------------
!! 3.1. Spectrum statistics
!!--------------------------------------------------------------------
if(MYID==0) then

 KF_SCL = 0.
EPS_SCL = 0.
 LF_SCL = 0.


!!--------------------------------------------------------------------
!!- Note: KFIRST takes place for dk in the integral
!!--------------------------------------------------------------------
do N = 2, KMAX
   KXGLOB = (N-1)*KFIRST
  KF_SCL =  KF_SCL +                       SPECTRUM(N)       *KFIRST
 EPS_SCL = EPS_SCL + 2.*DIFF_SCL*KXGLOB**2*SPECTRUM(N)       *KFIRST
  LF_SCL =  LF_SCL +                       SPECTRUM(N)/KXGLOB*KFIRST
end do


 ETA_SCL = (DIFF_SCL**3/EPS_FLU_SPEC)**0.25
TIME_SCL = (DIFF_SCL   /EPS_FLU_SPEC)**0.5

LF_SCL = LF_SCL*PPI/(4.*KF_SCL/3.)

!!- Mixing scalar time
TMIX_SCL = 2.*KF_SCL/EPS_SCL

!!--------------------------------------------------------------------
!! 3.2. Print in file
!!--------------------------------------------------------------------
if(NOUT<0) then
 !- Case for interpolation checking
 FILESPEC = 'spec_theta_int.dat'
else
 write(FILESPEC,4000)NOUT
end if
open(unit=200, file=trim(FILESPEC))

!!- File header
write(200,9000)

do N=2, KMAX-1
 KXGLOB = (N-1)*KFIRST
 write (200,10000)KXGLOB, KXGLOB*ETA_SCL,                            &
                           SPECTRUM(N),                              &
                           SPECTRUM(N)/(EPS_SCL*EPS_FLU_SPEC**(-3./4)*DIFF_SCL**1.25), &
                           SPECTRUM(N)*(1./EPS_SCL*EPS_FLU_SPEC**(1./3.)*KXGLOB**(5./3.))
end do
close (200)



deallocate(SPECTRUM)
deallocate(SPECLOC)

!!- Compute the optimal time step 
CFL_SCL = DIFF_SCL*DTIME/DX**2.D0


!!--------------------------------------------------------------------
!! 2.2. Print in "stat.info" the last spectrum statistics
!!--------------------------------------------------------------------
if((NOUT == FOUT1).or.(NOUT == 99)) then
 I=1
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),*)'====================================================================='
 write(UNIT_INFO(I),*)'SCALAR STATISTICS FROM SPECTRUM'
 write(UNIT_INFO(I),*)'====================================================================='
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)'Scalar variance:          <theta^2> = ',KF_SCL,' K2'
 write(UNIT_INFO(I),10603)'Longitudinal length scale: Lf_theta = ',LF_SCL,' m'
 write(UNIT_INFO(I),10604)'                        Lx/Lf_theta = ',LXMAX/LF_SCL
 write(UNIT_INFO(I),10603)'Mixing timescale:             T_mix = ',TMIX_SCL,' s'
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)'Dissipation:                 eps_theta = ',EPS_SCL,' K2/s'
 write(UNIT_INFO(I),10603)'Dissipative length scale:    eta_theta = ',ETA_SCL,' m'
 write(UNIT_INFO(I),10603)'Dissipative   time scale:    tau_theta = ',TIME_SCL,' s'
 write(UNIT_INFO(I),10603)'                        eta_theta.kmax = ',ETA_SCL*KLAST
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),*)
 write(UNIT_INFO(I),10603)' CFL stability = ',CFL_SCL,' [-] '
 write(UNIT_INFO(I),10603)' Time step for CFL=1: ',DX**2.D0/DIFF_SCL,' s'
 write(UNIT_INFO(I),*)
end if

end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
4000 format ('spec_theta_',I2.2,'.stat')

9000 format('# k,k*eta_t,E(k),E(k)/(eps_t*eps_f^-0.75*Kf^1.25),E(k)*(eps_t^-1*eps_f^1/3*k^5/3)')

10000 format (10(e17.7))

10601 format (2x,A,E13.6,A,E13.6)
10602 format (2x,A,E13.6)
10603 format (2x,A,E13.6,A)
10604 format (2x,A,f6.2)

end subroutine SPEC3D_SCALAR
