!!====================================================================
!!
!!  Initiation of Forcing scheme for statistically steady 
!!  turbulent flows
!!
!!====================================================================

subroutine INITIATION_FORCING

!!====================================================================
!! This routine works with ADD_FORCING
!!--------------------------------------------------------------------
!! The forcing scheme, proposed by Eswaran & Pope, J. Comp. & Fluids
!! (1988) is based on a stochastic force added at low-wavenumbers.
!!
!! The stochastic force is given by a Ornstein-Uhlenbeck process,
!! parameterized by a timescale (TIME_FORCE) and a variancce
!! (SIGMA_FORCE). The range of modified wavenumber is controled by
!! KFORCE_MIN and KFORCE_MAX.
!!
!! 
!!   KFORCE_MIN: minimum forced wavenumber
!!   KFORCE_MAX: maximum forced wavenumber
!!   TIME_FORCE: timescale of forcing
!!   SIGMA_FORCE: variance of forcing
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE
use PARAM_PHYS
use FORCING
use RANDOM_VAR

use P3DFFT

implicit none



!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
!- Wavenumber modulus
real(kind=8) :: KF

!- 
integer :: NF, NFM

!- Record size
integer :: RECSIZE

!- File name 
character(len=40) :: FILENAME


!- Statistics of the droped field
real(kind=8) :: STATX, STATY, STATZ
real(kind=8) :: STATXM, STATYM, STATZM

!- index
integer :: I, J, K
!---------------------------------------------------------------------


NFORCE_CPU = 0

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

   KF = (KX(I)**2 + KY(J)**2 + KZ(K)**2)**0.5

   if ((KF>=KFORCE_MIN).and.(KF<=KFORCE_MAX)) then
    NFORCE_CPU = NFORCE_CPU + 1
   end if

  end do
 end do
end do



!!- Sum over the domain
call ISUMCPU(NFORCE_CPU,NFORCE_FULL)




!!- allocation arrays containing stochastic acceleration
allocate(FORCING_UFOU(NFORCE_FULL))
allocate(FORCING_VFOU(NFORCE_FULL))
allocate(FORCING_WFOU(NFORCE_FULL))



if(MYID==0) write(*,10800) 'Forcing initiation --> Nwave= ',NFORCE_FULL



!!====================================================================
!!
!!====================================================================
if(INIT_FLUID_VELOCITY <= 2 .or.(INIT_FLUID_VELOCITY == 4)) then

 FORCING_UFOU(:) = cmplx(ZERO,ZERO)
 FORCING_VFOU(:) = cmplx(ZERO,ZERO)
 FORCING_WFOU(:) = cmplx(ZERO,ZERO)

 if(MYID==0) write(*,*)'Forcing coefficients initiation --> F=0'


!!====================================================================
!! 2. Read from file
!!====================================================================
elseif(INIT_FLUID_VELOCITY == 3) then


!!--------------------------------------------------------------------
!! 2.1. Multiple binary
!!--------------------------------------------------------------------
if (ISAVEFORCE==1) then
 !- Define file name
 FILENAME='forcing.ini'

 write(FILENAME,10101)'forcing',trim(FILE_EXT),'.ini'


 !- Open file containing the forcing coefficients
 open(unit =120, file=trim(FILENAME), form='unformatted')

 !- Number of forced waves
 read(120) NFORCE_FULL

 !- Random seed
 read(120) IDFORCE
 read(120) ISET, GSET
 read(120) XRANDF, YRANDF, ZRANDF


 !- Forcing coefficients
 read(120)(FORCING_UFOU(I),I=1,NFORCE_FULL)
 read(120)(FORCING_VFOU(I),I=1,NFORCE_FULL)
 read(120)(FORCING_WFOU(I),I=1,NFORCE_FULL)


 !- Close file
 close(120)

!!--------------------------------------------------------------------
!! 2.2. Direct access file
!!--------------------------------------------------------------------
elseif(ISAVEFORCE==2) then

!! RECSIZE = 8*2*NFORCE_FULL
 RECSIZE = (4 + 4 + 8 + 4 + 32*4) + 8*2*NFORCE_FULL

 FILENAME= 'forcing.ini'

!!- Open files
 open(unit=120,file=trim(FILENAME),access='direct',recl=RECSIZE,form='unformatted'   )

!!- Print data
 read(unit=120,rec=3*MYID+1)IDFORCE,ISET,GSET,IY,IV,FORCING_UFOU
 read(unit=120,rec=3*MYID+2)IDFORCE,ISET,GSET,IY,IV,FORCING_VFOU
 read(unit=120,rec=3*MYID+3)IDFORCE,ISET,GSET,IY,IV,FORCING_WFOU



 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

 close(120)


end if


 if(MYID==0) write(*,*) 'Forcing coefficients initiation --> Read from file'

end if


STATX = ZERO
STATY = ZERO
STATZ = ZERO

do I=1,NFORCE_FULL
 STATX = STATX + real(FORCING_UFOU(I)*conjg(FORCING_UFOU(I)))
 STATY = STATY + real(FORCING_VFOU(I)*conjg(FORCING_VFOU(I)))
 STATZ = STATZ + real(FORCING_WFOU(I)*conjg(FORCING_WFOU(I)))
end do

call RSUMCPU(STATX,STATXM)
call RSUMCPU(STATY,STATYM)
call RSUMCPU(STATZ,STATZM)



if(MYID==0) write(*,*) '       <fx> = ',STATXM / NFORCE_FULL
if(MYID==0) write(*,*) '       <fy> = ',STATYM / NFORCE_FULL
if(MYID==0) write(*,*) '       <fz> = ',STATZM / NFORCE_FULL



!!--------------------------------------------------------------------
10800 format (1x,A,I3,A)
10101 format (A,A,A)


end subroutine INITIATION_FORCING

