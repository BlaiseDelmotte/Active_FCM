
!!====================================================================
!!
!! 
!!> @brief
!!> Save particle characteristics for post-processing
!!
!! Date :  10/07/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_SAVE_PARTICLE_CHARACTERISTICS

!!=====================================================================
!! Save particle :
!!			- radii
!!
!!
!!------------------------------
!!
!!=====================================================================

use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use DNS_DIM
use MPI_structures

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Temporary arrays

real(kind=8), dimension(NPART_FULL,3) :: FCM_RADII


!- File names 
character(len=40) :: FCM_FILENAME_RADII


!- Index
integer :: IP, IND
!---------------------------------------------------------------------

IND = 0

do IP = 1,NPART_FULL
 
 
 if (IP.le.FCM_NSPHERE) then
  FCM_RADII(IP,1) = FCM_SPHERE_RADP(IP)
  FCM_RADII(IP,2) = FCM_SPHERE_RADP(IP)
  FCM_RADII(IP,3) = FCM_SPHERE_RADP(IP)
 else
  IND = IND + 1
  FCM_RADII(IP,1:3) = FCM_ELLIPSOID_RADP(IND,1:3)
 end if

end do

!!====================================================================
!!  MPI I/O 
!!====================================================================

 FCM_FILENAME_RADII = 'FCM_PART_RADII.end'
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_RADII,FCM_FILENAME_RADII)
 
 if(MYID==0) write(*,*) 'FCM : save particle characteristics --> OK'
 if(MYID==0) write(*,*) '    + MPI I/O'
 
 

!!====================================================================
!!  ASCII
!!====================================================================
!~ 
!~ !- Define file name
!~  write(FCM_FILENAME_POS,10301)'FCM_PART_POS.end'
!~  write(FCM_FILENAME_VEL,10301)'FCM_PART_VEL.end'
!~  write(FCM_FILENAME_ROT,10301)'FCM_PART_ROT.end'
!~  write(FCM_FILENAME_ORIENT,10301)'FCM_PART_ORIENT.end'
!~ 
!~  if (FCM_SWIMMING==1) then
!~   write(FCM_FILENAME_SWIM,10301)'FCM_PART_SWIM.end'
!~  end if
!~  
!~ !- Open file containing the last particle positions
!~  open(unit = 150, file = trim(FCM_FILENAME_POS),form='formatted')
!~  !- Open file containing the last particle translational velocities
!~  open(unit = 151, file = trim(FCM_FILENAME_VEL),form='formatted')
!~  !- Open file containing the last particle rotation rates
!~  open(unit = 152, file = trim(FCM_FILENAME_ROT),form='formatted')
!~  !- Open file containing the last particle orientations
!~  open(unit = 153, file = trim(FCM_FILENAME_ORIENT),form='formatted')
!~  !- Open file containing the last particle swimming directions
!~  if (FCM_SWIMMING==1) then 
!~   open(unit = 154, file = trim(FCM_FILENAME_SWIM),form='formatted')
!~  end if
!~  
!~ 
!~  do IP =1, NPART_FULL
!~ 
!~ !!====================================================================
!~ !! 1. POSITIONS
!~ !!====================================================================
!~   write(150,3) FCM_XP(IP), FCM_YP(IP), FCM_ZP(IP)
!~  
!~ !!====================================================================
!~ !! 2. TRANSLATIONAL VELOCITIES
!~ !!====================================================================
!~   write(151,3) FCM_UP(IP,1), FCM_VP(IP,1), FCM_WP(IP,1)
!~  
!~ !!====================================================================
!~ !! 3. ROTATIONAL VELOCITIES
!~ !!====================================================================
!~   write(152,3) FCM_OMPX(IP), FCM_OMPY(IP), FCM_OMPZ(IP)
!~  
!~ !!====================================================================
!~ !! 4. ORIENTATIONS, I.E. LOCAL FRAME (QUATERNION)
!~ !!====================================================================
!~   write(153,4) FCM_QUAT(IP,1), FCM_QUAT(IP,2), FCM_QUAT(IP,3), FCM_QUAT(IP,4)
!~  
!~ !!====================================================================
!~ !! 5. SWIMMING DIRECTIONS
!~ !!====================================================================
!~  if (FCM_SWIMMING==1) then 
!~   write(154,3) FCM_PSWIM(IP,1), FCM_PSWIM(IP,2), FCM_PSWIM(IP,3)
!~  end if
!~ 
!~  end do !!- Loop: 1, NPART_FULL
!~ 
!~ !- Close files
!~  close(150)
!~  close(151)
!~  close(152)
!~  close(153)
!~  if (FCM_SWIMMING==1) then 
!~   close(154)
!~  end if
 
 
 






!~ write(*,*) 'FCM : Final particle position, velocities and orientation --> Saved'
!~ write(*,*) '    + 1 text file per kinematic quantity'
!~ write(*,*) '    + all saved by PROC #', MYID


!!--------------------------------------------------------------------
10301 format (A)
3 format(3(e12.5, 3x))
4 format(4(e12.5, 3x))


end subroutine FCM_SAVE_PARTICLE_CHARACTERISTICS
