!!====================================================================
!!
!!  Forcing scheme for statistically steady turbulent flows
!!
!!====================================================================

subroutine ADD_FORCING2

!!====================================================================
!!
!! The forcing scheme, proposed by Eswaran & Pope, J. Comp. & Fluids
!! (1988) is based on a stochastic force added at low-wavenumbers.
!!
!! The stochastic force is given by a Ornstein-Uhlenbeck process,
!! parameterized by a timescale (TIME_FORCE) and a variancce
!! (SIGMA_FORCE). The range of modified wavenumber is controled by
!! KFORCING_MIN and KFORCING_MAX.
!!
!! 
!!   KFORCING_MIN: minimum forced wavenumber
!!   KFORCING_MAX: maximum forced wavenumber
!!   TIME_FORCE: timescale of forcing
!!   SIGMA_FORCE: variance of forcing
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use RHS_VARIABLES
use FLUID_VARIABLE
use PARAM_PHYS
use FORCING
use CHECK_CPU
use RANDOM_VAR

use P3DFFT

implicit none



!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
real(kind=8) :: XAMPF
real(kind=8) :: ALPHA, BETA

!- Random number generator
real(kind=8) :: GAUSS

real(kind=8) :: XAMPL, XE

!- number of forced waves
integer :: NF

!- Wavenumber
real(kind=8) :: KF

!- Imaginary and real part of forcing
real(kind=8) :: FORCE_REAL, FORCE_IMAG

!- Divergence
double complex :: DIVC


!- Statistics of the droped field
real(kind=8) :: STATX, STATY, STATZ

!- index
integer :: I, J, K, MJ, MK

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!!====================================================================
!! 1. Compute stochastic force
!!====================================================================

ALPHA =                      DTIME / TIME_FORCE
BETA  = (2.D0*SIGMA_FORCE**2*DTIME / TIME_FORCE)**0.5


NF = 0

do K = FSTART(3), FEND(3)
do J = FSTART(2), FEND(2)
do I = FSTART(1), FEND(1)


!! KF = (KX(I)**2 + KY_FULL(J)**2 + KZ_FULL(K)**2)**0.5
 KF = (KX(I)**2 + KY(J)**2 + KZ(K)**2)**0.5

 if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then

  NF = NF + 1

  XAMPF = (KFORCE_MAX - KF) / 0.2D0 / KFORCE_MAX
  XAMPL = (exp(2.D0*XAMPF) - 1.D0) / (exp(2.D0*XAMPF) + 1.D0)

!!-------------------------------------------------------------------- 
!! 1.1. x-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(XRANDF)

  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_UFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(XRANDF)

  XE = XE * XAMPL
  FORCE_IMAG = aimag(FORCING_UFOU(NF))*(1.D0-ALPHA) + XE*BETA
 
 
!- Stochastic acceleration
  FORCING_UFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


!!-------------------------------------------------------------------- 
!! 1.2. y-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(YRANDF)

  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_VFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(YRANDF)

  XE = XE * XAMPL
  FORCE_IMAG = aimag(FORCING_VFOU(NF))*(1.D0-ALPHA) + XE*BETA

 
!- Stochastic acceleration
  FORCING_VFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


!!-------------------------------------------------------------------- 
!! 1.3. z-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(ZRANDF)

  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_WFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
!  XE = GAUSS2(ZRANDF)

  XE = XE * XAMPL
  FORCE_IMAG = aimag(FORCING_WFOU(NF))*(1.D0-ALPHA) + XE*BETA

 
!- Stochastic acceleration
  FORCING_WFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 



!!-----------------------------
!!- Project on soloneidal basis
!!-----------------------------
!- Compute the divergence
! DIVC = KX(I)*FORCING_UFOU(NF) + KY(J)*FORCING_VFOU(NF) + KZ(K)*FORCING_WFOU(NF)

! FORCING_UFOU(NF) = FORCING_UFOU(NF) - KX(I)*DIVC/KF**2.
! FORCING_VFOU(NF) = FORCING_VFOU(NF) - KY(J)*DIVC/KF**2.
! FORCING_WFOU(NF) = FORCING_WFOU(NF) - KZ(K)*DIVC/KF**2.


 RHS_UFOU(I,J,K,TN) = RHS_UFOU(I,J,K,TN) + FORCING_UFOU(NF)
 RHS_VFOU(I,J,K,TN) = RHS_VFOU(I,J,K,TN) + FORCING_VFOU(NF)
 RHS_WFOU(I,J,K,TN) = RHS_WFOU(I,J,K,TN) + FORCING_WFOU(NF)

end if

end do
end do
end do 



!STATX = ZERO
!STATY = ZERO
!STATZ = ZERO
!do I=1,NFORCE_FULL
! STATX = STATX + real(FORCING_UFOU(I)*conjg(FORCING_UFOU(I)))
! STATY = STATY + real(FORCING_VFOU(I)*conjg(FORCING_VFOU(I)))
! STATZ = STATZ + real(FORCING_VFOU(I)*conjg(FORCING_WFOU(I)))
!end do
!
!call RSUMCPU(STATX,STATX)
!call RSUMCPU(STATY,STATY)
!call RSUMCPU(STATZ,STATZ)
!
!if(MYID==0) write(*,*) '       <fx> = ',STATX / NFORCE_FULL
!if(MYID==0) write(*,*) '       <fy> = ',STATY / NFORCE_FULL
!if(MYID==0) write(*,*) '       <fz> = ',STATZ / NFORCE_FULL


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(4) = CPU_FLUID(4) + TIME_END - TIME_START
end if


end subroutine ADD_FORCING2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function Gauss(ID)
!
! use RANDOM_VAR
!
!  implicit none
!
!  integer, intent(inout) :: ID
!
!  real(kind=8) :: Gauss
!  real(kind=8) :: FAC,RSQ,V1,V2
!  real(kind=8) :: XRAND
!  real(kind=8), external :: RAN1
!
!  intrinsic LOG, SQRT
!
!  if (ISET /= 0) then
!    Gauss = GSET
!    ISET = 0
!    return
!  end if
!  do
!    ! Generateur aleatoire tire de Numerical receipes (routine interne)
!    call random_number(XRAND)
!    V1 = 2.*XRAND - 1.
!    call random_number(XRAND)
!    V2 = 2.*XRAND - 1.
!
!!    V1 = 2.*RAN1(ID) - 1.
!!    V2 = 2.*RAN1(ID) - 1.
!    !
!    RSQ = V1**2 + V2**2
!    if (RSQ<1. .and. RSQ/=0.) exit
!  end do
!  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
!  GSET = V1 * FAC
!  Gauss = V2 * FAC
!  ISET = 1
!end function Gauss
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!function Gauss2(XRAND)
!
! use RANDOM_VAR
!
!  implicit none
!
!!!  integer, intent(inout) :: ID
!
!  real(kind=8) :: Gauss2
!  real(kind=8) :: FAC,RSQ,V1,V2
!  real(kind=8),intent(inout) :: XRAND
!
!!!  real(kind=8), external :: RAN1
!
!  intrinsic LOG, SQRT
!
!  if (ISET /= 0) then
!    Gauss2 = GSET
!    ISET = 0
!    return
!  end if
!  do
!    ! Generateur aleatoire tire de Numerical receipes (routine interne)
!    call random_number(XRAND)
!    V1 = 2.*XRAND - 1.
!
!    call random_number(XRAND)
!    V2 = 2.*XRAND - 1.
!
!!    V1 = 2.*RAN1(ID) - 1.
!!    V2 = 2.*RAN1(ID) - 1.
!    !
!    RSQ = V1**2 + V2**2
!    if (RSQ<1. .and. RSQ/=0.) exit
!  end do
!  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
!  GSET = V1 * FAC
!  Gauss2 = V2 * FAC
!  ISET = 1
!end function Gauss2
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
function RAN1(IDUM)

  use RANDOM_VAR

  implicit none

  integer,      parameter :: IA = 16807
  integer,      parameter :: IM = 2147483647
  real(kind=8), parameter :: AM = 1./IM
  integer,      parameter :: IQ = 127773
  integer,      parameter :: IR = 2836
  integer,      parameter :: NTAB = 32
  integer,      parameter :: NDIV = 1+(IM-1)/NTAB
  real(kind=8), parameter :: EPS = 1.2e-7
  real(kind=8), parameter :: RNMX = 1.-EPS

  integer, intent(inout) :: IDUM

  real(kind=8) :: RAN1

  integer :: J,K

!  intrinsic MAX, MIN

!  common /GENE2/ IV,IY
!  integer, dimension(32) :: IV
!  integer :: IY = 0

!!  data IV/NTAB*0/


  if (IDUM<=0 .or. IY==0) then
    IDUM = MAX(-IDUM,1)
    do J = NTAB+8,1,-1
      K = IDUM / IQ
      IDUM = IA*(IDUM-K*IQ) - IR*K
      if (IDUM < 0) then
        IDUM = IDUM + IM
      end if
      if (J <= NTAB) then
        IV(J) = IDUM
      end if
    end do
    IY = IV(1)
  end if
  K = IDUM / IQ
  IDUM = IA*(IDUM-K*IQ) - IR*K
  if (IDUM < 0) then
    IDUM = IDUM + IM
  end if
  J = 1 + IY/NDIV
  IY = IV(J)
  IV(J) = IDUM
  ran1 = MIN(AM*IY,RNMX)

end function RAN1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Gauss(ID)

 use RANDOM_VAR

 implicit none

 integer, intent(inout) :: ID

 real(kind=8) :: Gauss
 real(kind=8) :: FAC,RSQ,V1,V2

 real(kind=8), external :: RAN1

! intrinsic LOG, SQRT

! common /GENE1/ ISET,GSET

! integer :: ISET = 0
! real(kind=8) :: GSET

  if (ISET /= 0) then
    Gauss = GSET
    ISET = 0
    return
  end if
  do
    ! Generateur aleatoire tire de Numerical receipes (routine interne)
    V1 = 2.*RAN1(ID) - 1.
    V2 = 2.*RAN1(ID) - 1.
    !
    RSQ = V1**2 + V2**2
    if (RSQ<1. .and. RSQ/=0.) exit
  end do
  FAC = SQRT(-2.*LOG(RSQ)/RSQ)
  GSET = V1 * FAC
  Gauss = V2 * FAC
  ISET = 1
end function Gauss

!
