 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating the particle positions and orientation from elli
!!> psoids algorithm
!!
!! Date :  10/07/2013
!!
!!
!!> @author 
!!> Blaise Delmotte, Pascal Fede
!!====================================================================

subroutine INIT_ELLIPSOID_2D( NPART_START,& !Index of 1st particle to initiate
                           NPART_END,&   !Index of last particle to initiate
                           RAD1,&        !Radius in 1st principal direction
                           RAD2,&        !Radius in 2nd principal direction
                           RAD3,&        !Radius in 3rd principal direction
                           XP,&          !Particle x-positions
                           YP,&          !Particle y-positions
                           ZP,&          !Particle z-positions
                           QUAT1,&       !Particle orientations
                           QUAT2,&       !Particle orientations
                           QUAT3,&       !Particle orientations
                           QUAT4 )       !Particle orientations

!!====================================================================
!! Particles are initiated 
!! according to the prescribed size of the particle
!!====================================================================
!! Particle initiation: 
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE



implicit none



!!-------Arguments ----------------------------------------------------
integer,                             intent(in) :: NPART_START
integer,                             intent(in) :: NPART_END
real(kind=8), dimension(NPART_FULL), intent(in) :: RAD1
real(kind=8), dimension(NPART_FULL), intent(in) :: RAD2
real(kind=8), dimension(NPART_FULL), intent(in) :: RAD3
real(kind=8), dimension(NPART_FULL), intent(inout) :: XP
real(kind=8), dimension(NPART_FULL), intent(inout) :: YP
real(kind=8), dimension(NPART_FULL), intent(inout) :: ZP
real(kind=8), dimension(NPART_FULL), intent(inout) :: QUAT1
real(kind=8), dimension(NPART_FULL), intent(inout) :: QUAT2
real(kind=8), dimension(NPART_FULL), intent(inout) :: QUAT3
real(kind=8), dimension(NPART_FULL), intent(inout) :: QUAT4

!!-------Local variables-----------------------------------------------
real(kind=8), dimension(NPART_FULL,3,3) :: MAT_SIGMA
real(kind=8), dimension(NPART_FULL,3,3) :: MAT_U
real(kind=8), dimension(NPART_FULL,3,3) :: MAT_A

!- Temporary variables
real(kind=8), dimension(3,3) :: U, UT, TMP, A, SIGMA


!- Maximum radius per particle, to define the intersection region
real(kind=8), dimension(NPART_FULL)  :: MAX_RAD

!- Tests for while loops
logical :: CONT1
logical :: CONT2


!- Arrays for intersection detection
!- NL is the size of the stored symmetric matrix
integer, parameter :: NL = (3*(3+1))/2

!- FLAG determine the initiation type:
!       = 1 : random positions + orientations
integer :: FLAG

!- Max number of iterations to initiate a particle
integer :: NCYCLEMAX

!- Indices for loops
integer :: I,J, NCYCLE, NP
integer :: NBINTERSEC


!!============== START ROUTINE =========================================

FLAG = 1

NCYCLEMAX = 500000


!!===================================================================
!! 1. Initiation of particle sigma matrix
!!===================================================================
MAT_SIGMA = 0.0

do I = 1,NPART_FULL
 MAT_SIGMA(I,1,1) = 1.0/RAD1(I)**2.0
 MAT_SIGMA(I,2,2) = 1.0/RAD2(I)**2.0
 MAT_SIGMA(I,3,3) = 1.0/RAD3(I)**2.0
  
 MAX_RAD(I) = max( RAD1(I),RAD2(I),RAD3(I) )
end do

!!===================================================================
!! 2. Compute ellipsoid matrix A for already initated particle
!!===================================================================
if (NPART_START.gt.1) then
 do I= 1 , NPART_START 
 
  U(1,1) = 2.0*( QUAT1(I)**2 + QUAT2(I)**2 - 0.5 )
  U(1,2) = 2.0*( QUAT2(I)*QUAT3(I) - QUAT1(I)*QUAT4(I) )
  U(1,3) = 2.0*( QUAT2(I)*QUAT4(I) + QUAT1(I)*QUAT3(I) )
                                  
  U(2,1) = 2.0*( QUAT2(I) * QUAT3(I) + QUAT1(I) * QUAT4(I) )
  U(2,2) = 2.0*( QUAT1(I)**2 + QUAT3(I)**2 - 0.5 )
  U(2,3) = 2.0*( QUAT3(I) * QUAT4(I) - QUAT1(I) * QUAT2(I) )
					  
  U(3,1) = 2.0*( QUAT2(I) * QUAT4(I) - QUAT1(I) * QUAT3(I) )
  U(3,2) = 2.0*( QUAT3(I) * QUAT4(I) + QUAT1(I) * QUAT2(I) )
  U(3,3) = 2.0*( QUAT1(I)**2 + QUAT4(I)**2 - 0.5 )       


  SIGMA = MAT_SIGMA(I,1:3,1:3)
  
  UT = transpose( U )

  TMP = matmul( SIGMA, UT )

  A = matmul( U, TMP )
  
  MAT_A(I,1:3,1:3) = A
  
 end do
end if





!!===================================================================
!! 3. Spatial  distribution of ellipsoids
!!===================================================================
if (NPART_START.gt.1) then
 NP = NPART_START
else
 NP = 0
end if

CONT1 = .true.
CONT2 = .true.




do while(CONT1.and.CONT2) 

 NBINTERSEC = 0
 NP = NP +1
 !print*,'NP = ',NP
 

!!-------------------------------------------------------------------
!! 2.1 First step: initiate position and orientation
!!-------------------------------------------------------------------

 if (FLAG == 1) then
!!- Random ellipsoid
  call SEED_ELLIPSOID_2D( NP,&
                         MAT_SIGMA(NP,1:3,1:3),&
                         XP(NP),&
                         YP(NP),&
                         ZP(NP),&
                         QUAT1(NP),&
                         QUAT2(NP),&
                         QUAT3(NP),&
                         QUAT4(NP),&
                         MAT_A(NP,1:3,1:3)      )
 end if 





!!-------------------------------------------------------------------
!! 2.2 Second step: Check for intersection
!!-------------------------------------------------------------------
 if(NP>NPART_START) call INTERSEC_ELLIPSOID( NL,&
                                   NP,&
                                   NP,&
                                   XP,&
                                   YP,&
                                   ZP,&
                                   MAT_A,&
                                   MAX_RAD,&
                                   NBINTERSEC      )


!!-------------------------------------------------------------------
!! 2.3 Third step: in case of intersection rondly find another ellipsoid
!!-------------------------------------------------------------------
 NCYCLE = 0
 

 do while ((NBINTERSEC>0).and.(NCYCLE<=NCYCLEMAX))

  NCYCLE = NCYCLE + 1
  


!!- Random ellipsoid
  call SEED_ELLIPSOID_2D( NP,&
                         MAT_SIGMA(NP,1:3,1:3),&
                         XP(NP),&
                         YP(NP),&
                         ZP(NP),&
                         QUAT1(NP),&
                         QUAT2(NP),&
                         QUAT3(NP),&
                         QUAT4(NP),&
                         MAT_A(NP,1:3,1:3)  )
!!- Check intersection  
  call INTERSEC_ELLIPSOID( NL,&
                           NP,&
                           NP,&
                           XP,&
                           YP,&
                           ZP,&
                           MAT_A,&
                           MAX_RAD,&
                           NBINTERSEC      )

 end do
 

 

 if(NCYCLE>NCYCLEMAX/10) print*, 'ITERATIONS FOR ELLIPSOID INITIATIONS = ', NCYCLE

 if(NCYCLE>=NCYCLEMAX)     CONT1 =.false.

 if(NP==NPART_END)         CONT2 =.false.


end do

!!- The last ellipsoid does not satisfies non intersection criteria
if(.not.CONT1) NP = NP-1


if(.not.CONT1)then
 write(*,*)'Not converged'
 write(*,*)'Npart=',NP
 write(*,*)'NCYCLE=',NCYCLE
else
 write(*,*)'Converged'
 write(*,*)'Npart=',NP
 write(*,*)'NCYCLE=',NCYCLE
end if

end subroutine INIT_ELLIPSOID_2D
