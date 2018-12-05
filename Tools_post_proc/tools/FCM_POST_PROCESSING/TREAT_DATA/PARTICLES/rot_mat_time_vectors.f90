!!====================================================================
!!
!!
!!====================================================================

subroutine ROT_MAT_TIME_VECTORS(NSAVES, NPART, P1, P2, P3, ROT_MAT)

!!====================================================================
!!
!!
!!====================================================================


implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer, intent(in) :: NSAVES
! Number of swimmers
integer, intent(in) :: NPART
! Particle orientatiion
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: P1
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: P2
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: P3
real(kind=8), dimension(NSAVES,NPART,3,3), intent(inout) :: ROT_MAT

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Index
integer :: IB, IND



print*,' '
print*,'-------------------COMPUTE ROTATION MATRICES WITH VECTORS----------------------- '
print*,' '


!=====================================================================
! 1.COMPUTE ROTATION MATRICES
!=====================================================================
! p1 = ey
! p2 = ez
! p3 = ex
do IND = 1, NSAVES 
 do IB = 1, NPART
  ROT_MAT(IND,IB,3,1:3) = P1(IND,IB,1:3)
  
  ROT_MAT(IND,IB,1,1:3) = P2(IND,IB,1:3)
  
  ROT_MAT(IND,IB,2,1:3) = P3(IND,IB,1:3)

!~   print*,'NORM(P1(IND,IB,1:3)) = ', dsqrt(P1(IND,IB,1)**2 + P1(IND,IB,2)**2 + P1(IND,IB,3)**2 )
!~   print*,'NORM(P2(IND,IB,1:3)) = ', dsqrt(P2(IND,IB,1)**2 + P2(IND,IB,2)**2 + P2(IND,IB,3)**2 )
!~   print*,'NORM(P3(IND,IB,1:3)) = ', dsqrt(P3(IND,IB,1)**2 + P3(IND,IB,2)**2 + P3(IND,IB,3)**2 )
!~ 
!~   print*,'(P1.P2) = ', P1(IND,IB,1)*P2(IND,IB,1)+ P1(IND,IB,2)*P2(IND,IB,2)+P1(IND,IB,3)*P2(IND,IB,3)
!~ 
!~   read(*,*)



!~   ROT_MAT(IND,IB,1,1:3) = P1(IND,IB,1:3)
!~   
!~   ROT_MAT(IND,IB,2,1:3) = P2(IND,IB,1:3)
!~   
!~   ROT_MAT(IND,IB,3,1:3) = P3(IND,IB,1:3)


 end do
end do



print*,' '
print*,'-------------------END COMPUTE ROTATION MATRICES WITH VECTORS --------------------------- '
print*,' '

 

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

end subroutine ROT_MAT_TIME_VECTORS
