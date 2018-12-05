subroutine ell_pair_intersec( n, c1, gg1, c2, gg2, intersect )

!  Determine if two non-concentric ellipsoids, E1 and E2, intersect.

!  E1 is defined by:  (x-c1)^T * G1 * G1^T * (x-c1) <= 1.  The array
!  gg1 contains the lower triangular n x n matrix G1 in packed format.
!  Similarly for E2.  


!  Method:
!
!Phase I
!  1/ transform to y-space in which E1 is the unit ball: y = G1^T * ( x - c1 )
!  2/ find the point y2 in E2 that is closest to the origin
!  3/ if |y2| <= 1, then E1 and E2 intersect

!  S.B. Pope  6/13/04, 12/28/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
logical, intent(out)    :: intersect


real(k_dp) :: v(n)
integer    :: i, j, k, iter
real(k_dp) :: g1(n,n), g2(n,n), c2y(n), y2(n), y2norm, g2y(n,n)
            

! unpack gg1 and gg2  -----------------------------------------------------------
k  = 0
g1 = 0.d0
g2 = 0.d0
do j = 1, n
   do i = j, n
      k = k + 1
	  g1(i,j) = gg1(k)
	  g2(i,j) = gg2(k)	  
   end do
end do

! transform to y-space:  y = G1^T * (x-c1);  G2y = G1^{-1} * G2
g2y = g2
call dtrsm( 'L', 'L', 'N', 'N', n, n, 1.d0, g1, n, g2y, n )

c2y = c2 - c1  !  c2y = G1^T * ( c2 - c1 )
call dtrmv( 'L', 'T', 'N', n, g1, n, c2y, 1 )

!  find point y2 in E2 closest to the origin
v = 0.d0  ! origin
call ellu_pt_near_far( 1, n, c2y, g2y, v, y2 ) 

y2norm = sum( y2*y2 )  !  if y2 in unit ball then E1 and E2 intersect

if( y2norm <= 1.d0 ) then
   intersect = .true.           ! E1 and E2 intersect: all done
   return         
endif

intersect = .false.

return

end subroutine ell_pair_intersec
