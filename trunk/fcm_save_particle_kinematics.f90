
!!====================================================================
!!
!! 
!!> @brief
!!> Save particle kinematics for restart
!!
!! Date :  26/03/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_SAVE_PARTICLE_KINEMATICS(NCYCLE)

!!=====================================================================
!! Save particle :
!!			- positions
!!			- velocities (translationnal & rotationnal)
!!			- orientation (pswim)
!!
!!
!!
!! Notes for further development:  quaternions !!
!!------------------------------
!!
!!=====================================================================

use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use DNS_DIM
use FCM_FORCING_VARIABLE, only:FCM_SWIMMING
use MPI_structures

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!- Cycle number
integer, intent(in) :: NCYCLE

!- Temporary arrays
real(kind=8), dimension(NPART_FULL,3) :: FCM_POS
real(kind=8), dimension(NPART_FULL,3) :: FCM_POS_NOPER
real(kind=8), dimension(NPART_FULL,3) :: FCM_VEL
real(kind=8), dimension(NPART_FULL,3) :: FCM_ROT



!- File names 
character(len=40) :: FCM_FILENAME_POS
character(len=40) :: FCM_FILENAME_POS_NOPER
character(len=40) :: FCM_FILENAME_VEL
character(len=40) :: FCM_FILENAME_ROT
character(len=40) :: FCM_FILENAME_ORIENT
character(len=40) :: FCM_FILENAME_SWIM
character(len=40) :: FCM_FILENAME_P2
character(len=40) :: FCM_FILENAME_P3



!- Index
integer :: IP
!---------------------------------------------------------------------


do IP = 1,NPART_FULL

 FCM_POS(IP,1) = FCM_XP(IP)
 FCM_POS(IP,2) = FCM_YP(IP)
 FCM_POS(IP,3) = FCM_ZP(IP)
 
 FCM_POS_NOPER(IP,1) = FCM_XP_NOPER(IP)
 FCM_POS_NOPER(IP,2) = FCM_YP_NOPER(IP)
 FCM_POS_NOPER(IP,3) = FCM_ZP_NOPER(IP)
 
 FCM_VEL(IP,1) = FCM_UP(IP,1)
 FCM_VEL(IP,2) = FCM_VP(IP,1)
 FCM_VEL(IP,3) = FCM_WP(IP,1)
 
 FCM_ROT(IP,1) = FCM_OMPX(IP)
 FCM_ROT(IP,2) = FCM_OMPY(IP)
 FCM_ROT(IP,3) = FCM_OMPZ(IP) 


end do




!!====================================================================
!!  MPI I/O 
!!====================================================================

 
 !!- Print file name
 if(NCYCLE>0) then
  write(FCM_FILENAME_POS,10201)'FCM_PART_POS_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_POS = 'FCM_PART_POS.end'
 end if  
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_POS,FCM_FILENAME_POS)
 
 if(NCYCLE>0) then
  write(FCM_FILENAME_POS_NOPER,10201)'FCM_PART_POS_NOPER_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_POS_NOPER = 'FCM_PART_POS_NOPER.end'
 end if 
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_POS_NOPER,FCM_FILENAME_POS_NOPER)
   
 if(NCYCLE>0) then
  write(FCM_FILENAME_VEL,10201)'FCM_PART_VEL_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_VEL= 'FCM_PART_VEL.end'
 end if 
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_VEL,FCM_FILENAME_VEL)

 if(NCYCLE>0) then
  write(FCM_FILENAME_ROT,10201)'FCM_PART_ROT_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_ROT= 'FCM_PART_ROT.end'
 end if 
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_ROT,FCM_FILENAME_ROT)

 if(NCYCLE>0) then
  write(FCM_FILENAME_SWIM,10201)'FCM_PART_SWIM_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_SWIM= 'FCM_PART_SWIM.end'
 end if
 call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_PSWIM,FCM_FILENAME_SWIM)

!
! if (FCM_USE_QUAT == 1) then 
! 
!  if(NCYCLE>0) then
!   write(FCM_FILENAME_ORIENT,10201)'FCM_PART_ORIENT_t',NCYCLE,'.bin'
!  else
!   FCM_FILENAME_ORIENT= 'FCM_PART_ORIENT.end'
!  end if 
!  call FCM_SAVE_VAR_MPIIO(NPART_FULL,4,1,FCM_QUAT,FCM_FILENAME_ORIENT)
!  
! else
! 
!  if(NCYCLE>0) then
!   write(FCM_FILENAME_P2,10201)'FCM_PART_P2_t',NCYCLE,'.bin'
!   write(FCM_FILENAME_P3,10201)'FCM_PART_P3_t',NCYCLE,'.bin'
!  else
!   FCM_FILENAME_P2= 'FCM_PART_P2.end'
!   FCM_FILENAME_P3= 'FCM_PART_P3.end'
!  end if 
!  call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_P2,FCM_FILENAME_P2)
!  call FCM_SAVE_VAR_MPIIO(NPART_FULL,3,1,FCM_P3,FCM_FILENAME_P3) 
!   
! end if


 
!~  if(MYID==0) write(*,*) 'FCM : save particle kinematics --> OK'
!~  if(MYID==0) write(*,*) '    + MPI I/O'
 
 

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
10201 format (A,I8.8,A)
10301 format (A)
3 format(3(e12.5, 3x))
4 format(4(e12.5, 3x))


end subroutine FCM_SAVE_PARTICLE_KINEMATICS
