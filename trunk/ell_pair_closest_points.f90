subroutine ell_pair_closest_points( n, c1, gg1, c2, gg2, x1, x2, intersect )

!  Determine if two non-concentric ellipsoids, E1 and E2, intersect.
!  If they do not intersect, determine a separating hyperplane.
!  If they do intersect, determine the hyperplane which is the 
!     perpendicular bisector of the line of centers.

!  E1 is defined by:  (x-c1)^T * G1 * G1^T * (x-c1) <= 1.  The array
!  gg1 contains the lower triangular n x n matrix G1 in packed format.
!  Similarly for E2.  

!  The hyperplane is defined by v^T * ( x - xh  ) = 0, where v is
!  a unit vector.  The quantity s(x) = v^T * ( x - xh  ) is the signed
!  distance of the point x from the hyperplane: s(c1) is negative, and 
!  s(c2) is positive.  If E1 and E2 do not intersect, then s(x) is negative
!  for all points x in E1, and positive for all points x in E2.  

!  Method:
!
!Phase I
!  1/ transform to y-space in which E1 is the unit ball: y = G1^T * ( x - c1 )
!  2/ find the point y2 in E2 that is closest to the origin
!  3/ if |y2| <= 1, then E1 and E2 intersect; otherwise...
!  4/ define y1 = y2/|y2|;  yh = (y1+y2)/2;  vy = y2-y1; then vy^T * ( y - yh )
!     is a separating hyperplane
!  5/ transform back to x-space to obtain xh and v (and x1 and x2)
!  6/ all done if qual <=0 or max_it <=0
!
!Phase II
!  The objective of phase II is to improve the quality of the separating hyperplane.
!  The quality q ( q <= 1 ) is defined by: 
!     q= (distance between supporting hyperplanes)/|x1-x2|.
!  Starting from the value of x2 obtained in phase I, iterations are performed to:
!  7/ determine the point x1 in E1 closest to x2
!  8/ determine the point x2 in E2 closest to x1
!  9/ xh = (x1+x2)/2;  v = (x2-x1)/|x2-x1|
!  These operations (7-9) are performed iteratively until q >= qual, or the number
!  of iterations exceeds max_it.  The values of xh and v returned correspond to the
!  hyperplane of greatest separation encountered during (or before) the iterations
!  (which usually occurs on the final iteration). 
!
!  Notes:
!  1/ For 7 and 8, the same algorithm is used as in ell_pt_near_far, but there are 
!     efficiency gains in re-coding it here.
!  2/ A quadratic programming routine could be used to solve the whole problem.
!  3/ Convergence is slow if E1 and E2 nearly intersect; but then the separating
!     hyperplane obtained from phase I may be adequate.
!  4/ Aitken extrapolation could be used to speed up convergence.

!  S.B. Pope  6/13/04, 12/28/04

implicit none

integer, parameter      :: k_dp = kind(1.d0)
integer, intent(in)     :: n
real(k_dp), intent(in)  :: c1(n), gg1((n*(n+1))/2), c2(n), gg2((n*(n+1))/2)
real(k_dp), intent(out) :: x1(n),x2(n)
logical, intent(out)    :: intersect


integer    :: itmax = 1000  !  max. iterations in dgqt
integer    :: lu_diag = -1 !  logical unit for diagnostics (set < 0 to suppress)
real(k_dp) :: atol = 1.d-5, rtol = 1.d-4  !  tolerances for dgqt

integer    :: info, i, j, k, iter
real(k_dp) :: g1(n,n), g2(n,n), c1y(n), c2y(n), y2(n), y2norm, g1i(n,n), g2i(n,n), &
              dist, smin1, smax1, smin2, smax2, q, qual_tol, sep, sep_max, &
			  v_best(n), xh_best(n), gi(n,n), delta, par, par2, f, &
              a1(n,n), a2(n,n), b(n), z(2*n), wa1(2*n), wa2(2*n)
              
real(k_dp) :: afull(2*n,2*n), bfull(2*n), a11(n,n), a12(n,n), a21(n,n), a22(n,n)  
real(k_dp) :: ysol(2*n)
              

call ell_pair_intersec( n, c1, gg1, c2, gg2, intersect )

print*, 'intersect = ', intersect

if (intersect) then
 return
end if

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

!  G1^{-1} 
g1i = g1
call dtrtri( 'L', 'N', n, g1i, n, info )  

!  G2^{-1} 
g2i = g2
call dtrtri( 'L', 'N', n, g2i, n, info ) 



a11 = 0.d0
a12 = 0.d0
a21 = 0.d0
a22 = 0.d0



do j = 1, n
   a11(j:n,j) = g1i(j:n,j)
   a12(j:n,j) = -g1i(j:n,j)
   a21(j:n,j) = -g2i(j:n,j)
   a22(j:n,j) = g2i(j:n,j)
end do



! a11 = G1^{-1} * G1^{-T}
call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, g1i, n, a11, n )

! a12 = -G1^{-1} * G2^{-T}
call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, g2i, n, a12, n )

! a21 = -G2^{-1} * G1^{-T}
call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, g1i, n, a21, n )

! a22 = G2^{-1} * G2^{-T}
call dtrmm ( 'R', 'L', 'T', 'N', n, n, 1.d0, g2i, n, a22, n )


print*,' a12 = ', a12
read(*,*)
print*,' a21 = ', a21
read(*,*)


afull = 0.0d0
afull(1:n,1:n) = a11
afull(1:n,n+1:2*n) = a12
afull(n+1:2*n,1:n) = a21
afull(n+1:2*n,n+1:2*n) = a22

c1y = c1 - c2  !  c1y =  ( c1 - c2 ) * G1^-T 
call dtrmv( 'L', 'N', 'N', n, g1i, n, c1y, 1 )

c2y = c2 - c1  !  c2y =  ( c2 - c1 ) * G2^-T 
call dtrmv( 'L', 'N', 'N', n, g2i, n, c2y, 1 )


bfull(1:n) = c1y
bfull(n+1:2*n) = c2y



delta = 1.d0
par   = 0.d0

call dgqt(2*n,afull,2*n,bfull,delta,rtol,atol,itmax,par,f,ysol,info,z,wa1,wa2)

x1 = ysol(1:n)
x2 = ysol(n+1:2*n)


! x1 = c1 + G1^{-T}*y1 
call dtrmv( 'L', 'T', 'N', n, g1i, n, x1, 1 )
call dtrmv( 'L', 'T', 'N', n, g2i, n, x2, 1 )

x1 = x1 + c1
x2 = x2 + c2

return

end subroutine ell_pair_closest_points
