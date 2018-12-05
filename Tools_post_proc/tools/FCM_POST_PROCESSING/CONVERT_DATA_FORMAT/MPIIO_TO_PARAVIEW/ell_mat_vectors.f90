!!====================================================================
!!
!!
!!====================================================================

subroutine ELL_MAT_VECTORS(NPART, P1, P2, P3, RADII, ELL_MAT_VEC)

!!====================================================================
!!
!!
!!====================================================================


implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of swimmers
integer, intent(in) :: NPART
! Particle orientatiion
real(kind=8), dimension(NPART,3), intent(in) :: P1
real(kind=8), dimension(NPART,3), intent(in) :: P2
real(kind=8), dimension(NPART,3), intent(in) :: P3
real(kind=8), dimension(NPART,3), intent(in) :: RADII
real(kind=8), dimension(NPART,9), intent(inout) :: ELL_MAT_VEC

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
real(kind=8), dimension(3,3) :: RAD_MAT
real(kind=8), dimension(3,3) :: ROT_MAT
real(kind=8), dimension(3,3) :: ELL_MAT

real(kind=8) :: RADMAX
!- Index
integer :: IB, IND



!~ print*,' '
!~ print*,'-------------------COMPUTE ELLIPSOID MATRICES WITH VECTORS----------------------- '
!~ print*,' '


!=====================================================================
! 1.COMPUTE ELIIPSOIDS MATRICES
!=====================================================================
RAD_MAT = 0.0
RADMAX = maxval(RADII)
do IB = 1, NPART

!~  ROT_MAT(1,1:3) = P1(IB,1:3) 
!~  ROT_MAT(2,1:3) = P2(IB,1:3) 
!~  ROT_MAT(3,1:3) = P3(IB,1:3)
 
 ROT_MAT(1:3,1) = P1(IB,1:3) 
 ROT_MAT(1:3,2) = P2(IB,1:3) 
 ROT_MAT(1:3,3) = P3(IB,1:3)
 
 RAD_MAT(1,1) = RADII(IB,1)/RADMAX
 RAD_MAT(2,2) = RADII(IB,2)/RADMAX
 RAD_MAT(3,3) = RADII(IB,3)/RADMAX 

 
 ELL_MAT = matmul(ROT_MAT, RAD_MAT)
 ELL_MAT = matmul(ROT_MAT,transpose(ELL_MAT))
 
 ELL_MAT_VEC(IB,1:3) = ELL_MAT(1,1:3)
 ELL_MAT_VEC(IB,4:6) = ELL_MAT(2,1:3)
 ELL_MAT_VEC(IB,7:9) = ELL_MAT(3,1:3)
 
  
!~  print*,'P1(IB,1:3) = ', P1(IB,1:3) 
!~  print*,'P2(IB,1:3) = ', P2(IB,1:3) 
!~  print*,'P3(IB,1:3) = ', P3(IB,1:3) 
!~  print*,'RAD_MAT = ', RAD_MAT
!~  print*,'ROT_MAT = ', ROT_MAT
!~  print*,'ELL_MAT = ', ELL_MAT
!~  print*,'ELL_MAT_VEC(IB,:) = ', ELL_MAT_VEC(IB,:)
!~  read(*,*)
 
end do




!~ print*,' '
!~ print*,'-------------------END COMPUTE ROTATION MATRICES WITH VECTORS --------------------------- '
!~ print*,' '

 

!!====================================================================
1999 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "n", "p"')
2001 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',e12.5)

1998 format ('VARIABLES = "xp noper", "yp noper", "zp noper"')
2002 format ('VARIABLES = "xp", "yp", "zp"')
2003 format ('VARIABLES = "up", "vp", "wp"')
2004 format ('VARIABLES = "ompx", "ompy", "ompz"')
2005 format ('VARIABLES = "quat1", "quat2", "quat3", "quat4"')
2006 format ('VARIABLES = "pswimx", "pswimy", "pswimz"')
2007 format ('VARIABLES = "Sxx", "Sxy", "Sxz", "Syy", "Syz"')
2008 format ('VARIABLES = "a1", "a2", "a3"')

2010 format ('NPART_FULL = ', i4)
2011 format ('NPART_FULL = ')

10200 format (A)
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10205 format (I8.8)
10101 format (A,A,A)

end subroutine ELL_MAT_VECTORS
