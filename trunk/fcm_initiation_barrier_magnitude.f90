 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating FBLEVEL
!!
!! Date :  26/02/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_BARRIER_MAGNITUDE


use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE



implicit none


! Scale rep force with drag : FB*6*pi*mu*Vsw 
! (the multiplication by the radius is done in Barrier Force computation)
if (FCM_NSPHERE>0) then

 if (FCM_NSWIM(1)>0) then
  FCM_FBLEVEL = FCM_FBLEVEL * 6.0 * PPI * FCM_VSW_REF * maxval(FCM_SPHERE_RADP)
 else if (FCM_NSWIM(1).eq.0) then
  !FCM_FBLEVEL = FCM_FBLEVEL * 6.0 * PPI * maxval(FCM_SPHERE_RADP)
  
  ! TESTS FOR BIDISPERSITY
  if (KBT.gt.0.0) then
   FCM_FBLEVEL = FCM_FBLEVEL * KBT/(minval(FCM_SPHERE_RADP)**2)
  else 
   FCM_FBLEVEL = FCM_FBLEVEL * 6.0 * PPI 
  end if
 end if


!~  print*,'FCM_FBLEVEL = ', FCM_FBLEVEL
!~  read(*,*)

!!!!--- TO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
else
 
 if (FCM_NSWIM(1)>0) then
  !- Define GB parameters
  POW_RHO = 4.0
  THRES_RHO = 2.0
  POW_REP = 3.0
  SIGMA_MIN = 2.0*minval(FCM_ELLIPSOID_RADP)

  ! Scale Repulsive potential amplitude 
  FCM_FBLEVEL = FCM_FBLEVEL &
               * PPI * maxval(FCM_ELLIPSOID_RADP)**2 * FCM_VSW_REF &
              / ((POW_REP - 1.0)*4.0) 
 !!!!--- TO CHECK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 else if (FCM_NSWIM(1).eq.0) then
  !- Define GB parameters
  POW_RHO = 4.0
  THRES_RHO = 2.0
  POW_REP = 3.0
  SIGMA_MIN = 2.0*minval(FCM_ELLIPSOID_RADP)

  ! Scale Repulsive potential amplitude 
  FCM_FBLEVEL = FCM_FBLEVEL &
               * PPI * maxval(FCM_ELLIPSOID_RADP)**2 &
              / ((POW_REP - 1.0)*4.0) 
 end if
end if


end subroutine FCM_INITIATION_BARRIER_MAGNITUDE
