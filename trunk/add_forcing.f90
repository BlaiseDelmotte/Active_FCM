!!====================================================================
!!
!!  Forcing scheme for statistically steady turbulent flows
!!
!!====================================================================

subroutine ADD_FORCING

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

use P3DFFT

implicit none



!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
real(kind=8) :: XAMPF
real(kind=8) :: ALPHA, BETA

!- Random number generator
real(kind=8) :: GAUSS, XAMPL, XE

!- number of forced waves
integer :: NF

!- Wavenumber
real(kind=8) :: KF

!- Imaginary and real part of forcing
real(kind=8) :: FORCE_REAL, FORCE_IMAG

!- Divergence
double complex :: DIVC

!- index
integer :: I, J, K, MJ, MK

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


!! Tout le boulot est effectue par P0
!! Tout le boulot est effectue par P0
!! Tout le boulot est effectue par P0


if(MYID == 0) then 


!!====================================================================
!! 1. Compute stochastic force
!!====================================================================

ALPHA =                      DTIME / TIME_FORCE
BETA  = (2.D0*SIGMA_FORCE**2*DTIME / TIME_FORCE)**0.5


 do NF = 1, NFORCE_FULL

  I = IFORCING(NF)
  J = JFORCING(NF)
  K = KFORCING(NF)


  KF = (KX(I)**2 + KY_FULL(J)**2 + KZ_FULL(K)**2)**0.5

  XAMPF = (KFORCE_MAX - KF) / 0.2D0 / KFORCE_MAX
  XAMPL = (exp(2.D0*XAMPF) - 1.D0) / (exp(2.D0*XAMPF) + 1.D0)

!!-------------------------------------------------------------------- 
!! 1.1. x-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_UFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL

  FORCE_IMAG = aimag(FORCING_UFOU(NF))*(1.D0-ALPHA) + XE*BETA
 
 
!- Stochastic acceleration
  FORCING_UFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


!!-------------------------------------------------------------------- 
!! 1.2. y-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_VFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL
  FORCE_IMAG = aimag(FORCING_VFOU(NF))*(1.D0-ALPHA) + XE*BETA

 
!- Stochastic acceleration
  FORCING_VFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


!!-------------------------------------------------------------------- 
!! 1.3. z-velocity
!!-------------------------------------------------------------------- 

!- real part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL
  FORCE_REAL = real(FORCING_WFOU(NF))*(1.D0-ALPHA) + XE*BETA

!- complex part
  XE = GAUSS(IDFORCE)
  XE = XE * XAMPL
  FORCE_IMAG = aimag(FORCING_WFOU(NF))*(1.D0-ALPHA) + XE*BETA

 
!- Stochastic acceleration
  FORCING_WFOU(NF) = cmplx(FORCE_REAL,FORCE_IMAG) 


end do 




!!====================================================================
!! 2. Enforce divergence free
!!====================================================================
 do NF = 1, NFORCE_FULL

  I = IFORCING(NF)
  J = JFORCING(NF)
  K = KFORCING(NF)

  KF = KX(I)**2 + KY_FULL(J)**2 + KZ_FULL(K)**2

  if(KF <=ZERO) KF = INFINITY


  !- Compute the divergence
  DIVC = KX(I)     *FORCING_UFOU(NF) & 
       + KY_FULL(J)*FORCING_VFOU(NF) &
       + KZ_FULL(K)*FORCING_WFOU(NF)


  FORCING_UFOU(NF) = FORCING_UFOU(NF)  - KX(I)     /KF*DIVC
  FORCING_VFOU(NF) = FORCING_VFOU(NF)  - KY_FULL(J)/KF*DIVC
  FORCING_WFOU(NF) = FORCING_WFOU(NF)  - KZ_FULL(K)/KF*DIVC
 

 end do



end if !!- End loop MYID==0



!!====================================================================
!! 3. Enforce Hermitian symetry
!!====================================================================
!! We apply conjugating relations in order that the stochastic force
!! follows the Hermitian symetry properties. This will ensure that the
!! backward Fourier transform of the stochastic force is REAL.
!!====================================================================
do NF = 1, NFORCE_FULL
 if(NHERM(NF) > 0) then 
  FORCING_UFOU(NHERM(NF)) = conjg(FORCING_UFOU(NF))
  FORCING_VFOU(NHERM(NF)) = conjg(FORCING_VFOU(NF))
  FORCING_WFOU(NHERM(NF)) = conjg(FORCING_WFOU(NF))
 end if 
end do






!!====================================================================
!! 3. Envois a tout le monde !!!
!!====================================================================
if(NPROC>1) then
 call MPI_BCAST(FORCING_UFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_VFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FORCING_WFOU,NFORCE_FULL,MPI_DOUBLE_COMPLEX,0,MPI_COMM_WORLD,IERR)
end if




!!====================================================================
!! 4. Add forcing to the fluid momentum equation
!!====================================================================
do NF = 1, NFORCE_FULL

 I = IFORCING(NF)
 J = JFORCING(NF)
 K = KFORCING(NF)

 if(I>=FSTART(1).and. I<=FEND(1)) then
 if(J>=FSTART(2).and. J<=FEND(2)) then
 if(K>=FSTART(3).and. K<=FEND(3)) then

  RHS_UFOU(I,J,K,TN) = RHS_UFOU(I,J,K,TN) + FORCING_UFOU(NF)
  RHS_VFOU(I,J,K,TN) = RHS_VFOU(I,J,K,TN) + FORCING_VFOU(NF)
  RHS_WFOU(I,J,K,TN) = RHS_WFOU(I,J,K,TN) + FORCING_WFOU(NF)

 end if
 end if
 end if

end do







!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_FLUID(4) = CPU_FLUID(4) + TIME_END - TIME_START
end if


!!- Synchonize all task
call MPI_BARRIER(MPI_COMM_WORLD,IERR)


end subroutine ADD_FORCING

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function Gauss(ID)

  implicit none

  integer, intent(inout) :: ID

  real(kind=8) :: Gauss
  real(kind=8) :: FAC,RSQ,V1,V2
  real(kind=8) :: XRAND
  real(kind=8), external :: RAN1

  intrinsic LOG, SQRT

!!  common /GENE1/ ISET,GSET
  save ISET,GSET

  integer :: ISET = 0
  real(kind=8) :: GSET

  if (ISET /= 0) then
    Gauss = GSET
    ISET = 0
    return
  end if
  do
    ! Generateur aleatoire tire de Numerical receipes (routine interne)
!    call random_number(XRAND)
!    V1 = 2.*XRAND - 1.
!    call random_number(XRAND)
!    V2 = 2.*XRAND - 1.

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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
function ran1(IDUM)

  implicit none

  integer, parameter :: IA = 16807
  integer, parameter :: IM = 2147483647
  real(kind=8), parameter :: AM = 1./IM
  integer, parameter :: IQ = 127773
  integer, parameter :: IR = 2836
  integer, parameter :: NTAB = 32
  integer, parameter :: NDIV = 1+(IM-1)/NTAB
  real(kind=8), parameter :: EPS = 1.2e-7
  real(kind=8), parameter :: RNMX = 1.-EPS

  integer, intent(inout) :: IDUM

  real(kind=8) :: ran1

  integer :: J,K

  intrinsic MAX, MIN

!!  common /GENE2/ IV,IY
  save IV,IY

  integer, dimension(32) :: IV
  integer :: IY = 0

  data IV/NTAB*0/

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
end function ran1
