!!====================================================================
!!
!! This routine print informations about the simulation in file
!!                           'run.info'
!!
!!====================================================================

subroutine FCM_WRITE_INFO

!!====================================================================
!! 
!!====================================================================

use DNS_DIM
use PARAM_PHYS 
use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- index
integer :: I, J, N
!---------------------------------------------------------------------




!!=====================================================================
!! RUN INFORMATIONS
!!=====================================================================

!- Opening file
open(unit=157, file='fcm_run.info', action='write', status='replace')

write(157,*), LXMAX
write(157,*), LYMAX
write(157,*), LZMAX
write(157,*), NX
write(157,*), NY
write(157,*), NZ
write(157,*), DTIME
write(157,*), NCYCLEMAX
write(157,*), FOUT1
write(157,*), FOUT3
write(157,*), NPART_FULL
write(157,*), FCM_NSPHERE, FCM_NSPHERE_1
write(157,*), FCM_NELLIPSOID
write(157,*), FCM_VSW_REF
if (FCM_NSPHERE>0) then
 write(157,*), maxval(FCM_SPHERE_RADP)
else
 write(157,*), maxval(FCM_ELLIPSOID_RADP)
end if
write(157,*), FCM_USE_QUAT
if (FCM_NSPHERE>0) then
 write(157,*), FCM_SPHERE_SIZE
else
 write(157,*), maxval(FCM_ELLIPSOID_SIZE)
end if
write(157,*), FCM_FLOW_TYPE
write(157,*), FCM_BC
write(157,*), FCM_WALL_RANGE
write(157,*), FCM_NSWIM
write(157,*), KBT
write(157,*), FCM_RUNTUMBLE
write(157,*), FCM_TAU_RUN

write(157,*), ' ! LX '
write(157,*), ' ! LY '
write(157,*), ' ! LZ '
write(157,*), ' ! NX '
write(157,*),  ' ! NY '
write(157,*),  ' ! NZ '
write(157,*), ' ! DTIME'
write(157,*), ' ! NCYCLEMAX'
write(157,*), ' ! FLUID DUMP FREQ'
write(157,*), ' ! PARTICLE DUMP FREQ'
write(157,*),  ' ! NPART_FULL'
write(157,*),  ' ! FCM_NSPHERE, FCM_NSPHERE_1'
write(157,*),  ' ! FCM_NELLIPSOID'
write(157,*),  ' ! FCM_VSW'
write(157,*),  ' ! RADMAX'
write(157,*),  ' ! FCM_USE_QUAT'
write(157,*),  ' ! MAX SIZE'
write(157,*),  ' ! FLOW TYPE'
write(157,*),  ' ! BC'
write(157,*),  ' ! WALL RANGE'
write(157,*),  ' ! NSWIM'
write(157,*),  ' ! KBT'
write(157,*),  ' ! FCM_RUNTUMBLE'
write(157,*),  ' ! FCM_TAU_RUN'
close(157)

end subroutine FCM_WRITE_INFO
