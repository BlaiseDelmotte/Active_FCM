
!!====================================================================
!!
!! 
!!> @brief
!!> Save particle dynamics for restart
!!
!! Date :  21/06/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_SAVE_PARTICLE_DYNAMICS(NCYCLE)

!!=====================================================================
!! Save particle :
!!			- Stresslet coefficents
!!			- ROS
!!
!!
!!=====================================================================

use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use DNS_DIM
use MPI_structures

implicit none

!- Cycle number
integer, intent(in) :: NCYCLE
!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- File names 
character(len=40) :: FCM_FILENAME_STRESSLET
!~ character(len=40) :: FCM_FILENAME_ROS



!- Index
integer :: IP, IC
!---------------------------------------------------------------------

!!====================================================================
!!  MPI I/O 
!!====================================================================

!!- Print file name
 if(NCYCLE>0) then
  write(FCM_FILENAME_STRESSLET,10201)'FCM_PART_STRESSLET_t',NCYCLE,'.bin'
 else
  FCM_FILENAME_STRESSLET = 'FCM_PART_STRESSLET.end'
 end if
 if (FCM_NSWIM(1)>0) then
  call FCM_SAVE_VAR_MPIIO(FCM_NSWIM(1),5,1,FCM_SIJ+FCM_SPIJ,FCM_FILENAME_STRESSLET)
 else 
  call FCM_SAVE_VAR_MPIIO(FCM_ACTIVATE_STRESSLET,5,1,FCM_SIJ,FCM_FILENAME_STRESSLET)
 end if
 
!~  FCM_FILENAME_ROS = 'FCM_PART_ROS.end'
!~  call FCM_SAVE_VAR_MPIIO(NPART_FULL,5,1,FCM_EIJ,FCM_FILENAME_ROS)


!~   if(MYID==0) write(*,*) 'FCM : save final particle stresslets --> OK'
!~   if(MYID==0) write(*,*) '    + MPI I/O'


!!====================================================================
!!  ASCII
!!====================================================================
!~ !- Define file name
!~ write(FCM_FILENAME_STRESSLET,10301)'FCM_PART_STRESSLET.end'
!~ write(FCM_FILENAME_ROS,10301)'FCM_PART_ROS.end'
!~ 
!~ !- Open file containing the last particle stresslet coeffs
!~ open(unit = 155, file = trim(FCM_FILENAME_STRESSLET),form='formatted')
!- Open file containing the last particle ROS coeffs
!~ open(unit = 156, file = trim(FCM_FILENAME_ROS),form='formatted')
!~ 
!~ do IP =1, NPART_FULL
!~ 
!~ 
!~ !!====================================================================
!~ !! 1. PARTICLE STRESSLETS
!~ !!====================================================================
!~  write(155,5) (FCM_SIJ(IP,IC), IC = 1,5)
!~  
!!====================================================================
!! 2. PARTICLE ROS
!!====================================================================
!~  write(156,5) (FCM_EIJ(IP,IC), IC = 1,5)
!~ 
!~ 
!~ end do !!- Loop: 1, NPART_FULL
!~ 
!~ !- Close files
!~ close(155)
!~ close(156)
!~ 
!~ write(*,*) 'FCM : Final particle stresslets --> Saved'
!~ write(*,*) '    + 1 text file per dynamic quantity'
!~ write(*,*) '    + all saved by PROC #', MYID



!!--------------------------------------------------------------------
10201 format (A,I8.8,A)
10301 format (A)
3 format(3(e12.5, 3x))
4 format(4(e12.5, 3x))
5 format(5(e12.5, 3x))


end subroutine FCM_SAVE_PARTICLE_DYNAMICS
