 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating swimming particle velocities
!!
!! Date :  26/02/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PARTICLE_VELOCITY


use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE



implicit none

character(len=40) :: FILENAME
integer :: I

!- Vsw is defined in terms of radii/s, it thus remains constant whichever the radius of swimmers

if (FCM_NSPHERE>0) then

FCM_VSW(1:FCM_NSWIM(1)) = FCM_VSW_REF * FCM_SPHERE_RADP(1:FCM_NSWIM(1))
 
! write(FILENAME,10201) 'Velocities_Nsw_',FCM_NSWIM(1),'.txt'
! open(unit=666, file=FILENAME, status='old')
!   do I = 1,FCM_NSWIM(1)
!    read(666,*) FCM_VSW(I)
!    FCM_VSW(I) = FCM_VSW(I)/5.0*FCM_SPHERE_RADP(I)/dsqrt(3.0)
!    !print*,MYID, FCM_VSW(I)/FCM_SPHERE_RADP(I)
!   end do
!  close(666)
else

 FCM_VSW_REF = FCM_VSW_REF * maxval(FCM_ELLIPSOID_RADP)
 
 FCM_VSW(1:FCM_NSWIM(1)) = FCM_VSW_REF * maxval(FCM_ELLIPSOID_RADP)

end if


10201 format (A,I4.4,A)

end subroutine FCM_INITIATION_PARTICLE_VELOCITY
