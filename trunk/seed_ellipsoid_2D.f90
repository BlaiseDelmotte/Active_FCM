 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating the particle positions and orientation 
!!
!! Date :  10/07/2013
!!
!!
!!> @author 
!!> Blaise Delmotte, Pascal Fede
!!====================================================================

subroutine SEED_ELLIPSOID_2D( NP,&
                             SIGMA,&
                             MAX_RAD, &
                             XP,&
                             YP,&
                             ZP,&
                             QUAT1,&
                             QUAT2,&
                             QUAT3,&
                             QUAT4,&
                             A )


!!====================================================================
!! Random values are obtained from "init_random_seed" defined at the bottom
!!====================================================================
!!
!!====================================================================
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_FORCING_VARIABLE
use FCM_PART_VARIABLE

implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,                                intent(in) :: NP
real(kind=8), dimension(3,3),           intent(in) :: SIGMA
real(kind=8),                           intent(in) :: MAX_RAD

!!----------------------------------------------------------------------
real(kind=8),                           intent(out) :: XP
real(kind=8),                           intent(out) :: YP
real(kind=8),                           intent(out) :: ZP
real(kind=8),                           intent(out) :: QUAT1
real(kind=8),                           intent(out) :: QUAT2
real(kind=8),                           intent(out) :: QUAT3
real(kind=8),                           intent(out) :: QUAT4
real(kind=8), dimension(3,3),           intent(out) :: A


!!======================================================================

!- Local variables

real(kind=8), dimension(3) :: AXIS
real(kind=8) ::  ANGLE

real(kind=8), dimension(3,3) :: U
real(kind=8), dimension(3,3) :: UT
real(kind=8), dimension(3,3) :: TMP


real(kind=8),dimension(3) :: RANDNUM_POS
real(kind=8),dimension(3) :: RANDNUM_AXIS
real(kind=8)              :: RANDNUM_ANGLE

real(kind=8)              :: NORM_QUAT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call init_random_seed()


!!-----------------------------------------------------
!! Center of mass position
!!----------------------------------------------------- 
 
call random_number(RANDNUM_POS) 

!!- x-direction
if (FCM_BC==2) then
 !- mirror plane at LXMAX/2
 !- Particles cannot be closer than 1.1*RAD_MAX to the planes
 XP = (LXMAX/2.0-2.2*MAX_RAD)*RANDNUM_POS(1) + 1.1*MAX_RAD
else if (FCM_BC==1) then
 XP = LXMAX*RANDNUM_POS(1)
end if
!!- y-direction
YP = LYMAX /2.0
!!- z-direction
ZP = LZMAX*RANDNUM_POS(3)
 


!!-----------------------------------------------------
!! Orientation
!!-----------------------------------------------------

! All particles orientation are set in the x-direction
if (FCM_INIT_PART_ORIENT==0) then

 ANGLE = PPI/2.0

 AXIS = 1.0

 QUAT1 = cos(ANGLE/2.0)
 QUAT2 = AXIS(1) * sin(ANGLE/2.0)
 QUAT3 = AXIS(2) * sin(ANGLE/2.0)
 QUAT4 = AXIS(3) * sin(ANGLE/2.0)

 NORM_QUAT = dsqrt( QUAT1**2 &
		  + QUAT2**2 &
		  + QUAT3**2 &
		  + QUAT4**2 )
		  
 QUAT1 = QUAT1/NORM_QUAT
 QUAT2 = QUAT2/NORM_QUAT
 QUAT3 = QUAT3/NORM_QUAT
 QUAT4 = QUAT4/NORM_QUAT


! All particles orientation are randomly set 
elseif (FCM_INIT_PART_ORIENT==1) then

 call random_number(RANDNUM_ANGLE)
 ANGLE = 2.D0*PPI*RANDNUM_ANGLE
 
 AXIS = 0.0
 AXIS(2) = 1.0

 AXIS = AXIS/sqrt( AXIS(1)**2 + AXIS(2)**2 + AXIS(3)**2 )

 QUAT1 = cos(ANGLE/2.0)
 QUAT2 = AXIS(1) * sin(ANGLE/2.0)
 QUAT3 = AXIS(2) * sin(ANGLE/2.0)
 QUAT4 = AXIS(3) * sin(ANGLE/2.0)
 
end if



U(1,1) = 2.0*( QUAT1**2 + QUAT2**2 - 0.5 )
U(1,2) = 2.0*( QUAT2*QUAT3 - QUAT1*QUAT4 )
U(1,3) = 2.0*( QUAT2*QUAT4 + QUAT1*QUAT3 )
                                  
U(2,1) = 2.0*( QUAT2 * QUAT3 + QUAT1 * QUAT4 )
U(2,2) = 2.0*( QUAT1**2 + QUAT3**2 - 0.5 )
U(2,3) = 2.0*( QUAT3 * QUAT4 - QUAT1 * QUAT2 )
					  
U(3,1) = 2.0*( QUAT2 * QUAT4 - QUAT1 * QUAT3 )
U(3,2) = 2.0*( QUAT3 * QUAT4 + QUAT1 * QUAT2 )
U(3,3) = 2.0*( QUAT1**2 + QUAT4**2 - 0.5 )       


UT = transpose( U )

TMP = matmul( SIGMA, UT )

A = matmul( U, TMP )



end subroutine SEED_ELLIPSOID_2D
!
!~ 
!~ subroutine init_random_seed()
!~ implicit none
!~ integer, allocatable :: seed(:)
!~ integer :: i, n, un, istat, dt(8), pid, t(2), s
!~ integer(8) :: count, tms
!~ 
!~ call random_seed(size = n)
!~ allocate(seed(n))
!~ ! First try if the OS provides a random number generator
!~ open(newunit=un, file="/dev/urandom", access="stream", &
!~ 	 form="unformatted", action="read", status="old", iostat=istat)
!~ if (istat == 0) then
!~    read(un) seed
!~    close(un)
!~ else
!~    ! Fallback to XOR:ing the current time and pid. The PID is
!~    ! useful in case one launches multiple instances of the same
!~    ! program in parallel.
!~    call system_clock(count)
!~    if (count /= 0) then
!~ 	  t = transfer(count, t)
!~    else
!~ 	  call date_and_time(values=dt)
!~ 	  tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
!~ 		   + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
!~ 		   + dt(3) * 24 * 60 * 60 * 60 * 1000 &
!~ 		   + dt(5) * 60 * 60 * 1000 &
!~ 		   + dt(6) * 60 * 1000 + dt(7) * 1000 &
!~ 		   + dt(8)
!~ 	  t = transfer(tms, t)
!~    end if
!~    s = ieor(t(1), t(2))
!~    pid = 1 + 1099279 ! Add a prime
!~    s = ieor(s, pid)
!~    if (n >= 3) then
!~ 	  seed(1) = t(1) + 36269
!~ 	  seed(2) = t(2) + 72551
!~ 	  seed(3) = pid
!~ 	  if (n > 3) then
!~ 		 seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
!~ 	  end if
!~    else
!~ 	  seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
!~    end if
!~ end if
!~ call random_seed(put=seed)
!~ end subroutine init_random_seed
