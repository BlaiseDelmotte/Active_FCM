 !!====================================================================
!!
!! 
!!> @brief
!!> Create symmetric forcing field with respect to x-direction 
!!
!! Date :  28/06/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_MIRROR_FORCES_X

!!====================================================================
!! Create symmetric forcing field with respect to x = LXMAX/2
!! This is for the Stress-free condition
!!====================================================================
!!
!!====================================================================

use DNS_DIM	        !- Dimension
use PARAM_PHYS          !- Physical & numerical parameters
use FLUID_VARIABLE      !- Fluid velocity
use FCM_FORCING_VARIABLE
use FCM_PART_VARIABLE
use GEOMETRIC_VARIABLE 
use RHS_VARIABLES       !- Variable for cpu time checking
use WORK_ARRAYS
use CHECK_CPU	       

use MPI_STRUCTURES


implicit none

!!====================================================================
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------

!- Time control variable
real(kind=8) :: TIME_START, TIME_END
real(kind=8) :: TIME_START2

real(kind=8) :: MEASURE_START, MEASURE_END

integer :: I, J, K
!!
!!====================================================================

!- CPU check
if(MYID == 0) then
 TIME_START = MPI_WTIME()
 TIME_START2 = TIME_START
end if



!- Mirror boundary : x = LXMAX/2
do I = NX/2+2, NX
 FCM_FORCING_X(I,ISTART(2):IEND(2),ISTART(3):IEND(3)) =  &
     - FCM_FORCING_X(NX+2-I,ISTART(2):IEND(2),ISTART(3):IEND(3))
     
 FCM_FORCING_Y(I,ISTART(2):IEND(2),ISTART(3):IEND(3)) =  &
       FCM_FORCING_Y(NX+2-I,ISTART(2):IEND(2),ISTART(3):IEND(3))
 
 FCM_FORCING_Z(I,ISTART(2):IEND(2),ISTART(3):IEND(3)) =  &
       FCM_FORCING_Z(NX+2-I,ISTART(2):IEND(2),ISTART(3):IEND(3))
    
end do

!- No stress throuh lower plane
FCM_FORCING_X(1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = 0.0
!- No stress throuh upper plane
FCM_FORCING_X(NX/2+1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = 0.0


!- At the boundary, tangential stresses are twice bigger 
!- as we superimpose image systems
FCM_FORCING_Y(1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = &
                2.0*FCM_FORCING_Y(1,ISTART(2):IEND(2),ISTART(3):IEND(3))
FCM_FORCING_Z(1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = &
                2.0*FCM_FORCING_Z(1,ISTART(2):IEND(2),ISTART(3):IEND(3))
FCM_FORCING_Y(NX/2+1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = &
                2.0*FCM_FORCING_Y(NX/2+1,ISTART(2):IEND(2),ISTART(3):IEND(3))
FCM_FORCING_Z(NX/2+1,ISTART(2):IEND(2),ISTART(3):IEND(3)) = &
                2.0*FCM_FORCING_Z(NX/2+1,ISTART(2):IEND(2),ISTART(3):IEND(3))
                               
!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()

 CPU_FLUID(1) = CPU_FLUID(1) + TIME_END - TIME_START2
end if
   
end subroutine FCM_MIRROR_FORCES_X
