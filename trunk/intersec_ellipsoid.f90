 !!====================================================================
!!
!! 
!!> @brief
!!> Routine which detects if ellipsoids are overlapping
!!
!! Date :  10/07/2013
!!
!!
!!> @author 
!!> Blaise Delmotte, Pascal Fede
!!====================================================================
subroutine INTERSEC_ELLIPSOID( NL,&
                               NCHECK_START,&
                               NCHECK_END,&
                               XP,&
                               YP,&
                               ZP,&
                               MAT_A,&
                               MAX_RAD,&
                               NBINT     )

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE

implicit none



!!====================================================================
!! Intersection is detected with the Pope Routine "ell_pair_separate"
!! now called "ell_pair_intersec"
!!====================================================================
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,                                   intent(in) :: NL
integer,                                   intent(in) :: NCHECK_START
integer,                                   intent(in) :: NCHECK_END
real(kind=8), dimension(NPART_FULL),       intent(inout) :: XP
real(kind=8), dimension(NPART_FULL),       intent(inout) :: YP
real(kind=8), dimension(NPART_FULL),       intent(inout) :: ZP
real(kind=8), dimension(NPART_FULL,3,3),   intent(in) :: MAT_A
real(kind=8), dimension(NPART_FULL),       intent(in) :: MAX_RAD

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer,                                   intent(out) :: NBINT

!!======================================================================

!- Local variables

real(kind=8), dimension(NL) :: L1, L2
real(kind=8), dimension(NL) :: G1, G2

real(kind=8), dimension(3,3) :: A1, A2
real(kind=8), dimension(3)      :: C1, C2

!!- Periodic images of particles
real(kind=8), dimension(3)      :: C1_IM

real(kind=8) :: DIST_CENTERS

logical      :: INTERSECT 


integer :: I, J, K
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
NBINT = 0

do I = NCHECK_START, NCHECK_END


!!- ith particle

 
 A1 = MAT_A(I,1:3,1:3)

 C1(1) = XP(I)
 C1(2) = YP(I)
 C1(3) = ZP(I)
 
!  Cholesky representation
 call ELL_FULL2LOW( 3, A1, L1 )  ! g contains lower triangle of aa in packed format
 call ELL_LOW2CHOL( 3, L1 ,G1 )  ! g contains lower Cholesky triangle in packed format

  
 do J = 1, I-1 

!!- jth particle
  A2 = MAT_A(J,1:3,1:3)

  C2(1) = XP(J)
  C2(2) = YP(J)
  C2(3) = ZP(J)
  
!  Cholesky representation
  call ELL_FULL2LOW( 3, A2, L2 )  ! g contains lower triangle of aa in packed format
  call ELL_LOW2CHOL( 3, L2 ,G2 )  ! g contains lower Cholesky triangle in packed format
 

 !  Compute the periodic images in each direction
  if ((C1(1)-C2(1)).ge.LXMAX/2.0) then
   C1_IM(1) = C1(1) - LXMAX
  elseif ((C2(1)-C1(1)).ge.LXMAX/2.0) then
   C1_IM(1) = C1(1) + LXMAX
  else
   C1_IM(1) = C1(1)
  end if
 
  if ((C1(2)-C2(2)).ge.LYMAX/2.0) then
   C1_IM(2) = C1(2) - LYMAX
  elseif ((C2(2)-C1(2)).ge.LYMAX/2.0) then
   C1_IM(2) = C1(2) + LYMAX
  else
   C1_IM(2) = C1(2)
  end if
 
  if ((C1(3)-C2(3)).ge.LZMAX/2.0) then
   C1_IM(3) = C1(3) - LZMAX
  elseif ((C2(3)-C1(3)).ge.LZMAX/2.0) then
   C1_IM(3) = C1(3) + LZMAX
  else
   C1_IM(3) = C1(3)
  end if

  DIST_CENTERS = sqrt( (C1_IM(1)-C2(1))**2 + (C1_IM(2)-C2(2))**2 + (C1_IM(3)-C2(3))**2 )
      

  if (DIST_CENTERS.le.(MAX_RAD(I)+MAX_RAD(J))) then
   
   call ELL_PAIR_INTERSEC( 3, C1_IM, G1, C2, G2, INTERSECT )
   
   if(INTERSECT) then
    NBINT = NBINT + 1    
    !! - If one image of the i-th particle intersects the j-th, stop loop on images
    exit  
   end if !if(INTERSECT) then
    
  end if !if (DIST_CENTERS.le.(2.0*(MAX_RAD(I)+MAX_RAD(J))) then
   
 end do !do J = I+1, NCHECK

end do !do I=NCHECK_START, NCHECK_END


end subroutine INTERSEC_ELLIPSOID
