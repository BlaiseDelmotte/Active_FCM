!!====================================================================
!!
!!
!!====================================================================

subroutine MEAN_VEL(NSAVES, NPART, TIME_VEC, VEL)

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
! Physical Time
real(kind=8), dimension(NSAVES), intent(in) :: TIME_VEC
! Particle velocities
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: VEL

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

real(kind=8),  dimension(NSAVES) :: MEAN_V
real(kind=8),  dimension(NSAVES) :: MEAN_VX
real(kind=8),  dimension(NSAVES) :: MEAN_VY
real(kind=8),  dimension(NSAVES) :: MEAN_VZ

!- Norm of velocities
real(kind=8),  dimension(NSAVES,NPART) :: VELNORM


!- File name 
character(len=40) :: FILENAME

!- Index
integer :: I, J, K, IND


print*,' '
print*,'-------------------START MEAN VEL------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE NORM OF VELOCITIES
!=====================================================================
print*, 'COMPUTE NORM OF VELOCITIES'
do IND = 1, NSAVES
 do I = 1, NPART
 
  VELNORM(IND,I) = dsqrt( VEL(IND,I,1)**2 &
                        + VEL(IND,I,2)**2 &
                        + VEL(IND,I,3)**2 )
 
 end do
end do
 
print*, 'COMPUTE NORM OF VELOCITIES---> OK'



!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE MEAN VEL  '



do IND = 1, NSAVES 

 MEAN_V(IND) = 0.0
 MEAN_VX(IND) = 0.0
 MEAN_VY(IND) = 0.0
 MEAN_VZ(IND) = 0.0
 
 do I = 1, NPART
  MEAN_V(IND) = MEAN_V(IND) + VELNORM(IND,I)
  MEAN_VX(IND) = MEAN_VX(IND) + VEL(IND,I,1)
  MEAN_VY(IND) = MEAN_VY(IND) + VEL(IND,I,2)
  MEAN_VZ(IND) = MEAN_VZ(IND) + VEL(IND,I,3)  
 end do
 
 MEAN_V(IND) = MEAN_V(IND)/NPART
 MEAN_VX(IND) = MEAN_VX(IND)/NPART
 MEAN_VY(IND) = MEAN_VY(IND)/NPART
 MEAN_VZ(IND) = MEAN_VZ(IND)/NPART
 
end do

print*,'COMPUTE MEAN VEL --->  OK '

!=====================================================================
! 3. WRITE STAT
!=====================================================================

!---------- norm-----------
!!-Print filename
write(FILENAME,10200)'MEAN_VELNORM.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, NSAVES
 write(301,'(2(e17.7))') TIME_VEC(J), MEAN_V(J)
end do

!- close file
close(301)

!---------- x-component-----------
!!-Print filename
write(FILENAME,10200)'MEAN_VELX.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, NSAVES
 write(302,'(2(e17.7))') TIME_VEC(J), MEAN_VX(J)
end do

!- close file
close(302)


!---------- y-component-----------
!!-Print filename
write(FILENAME,10200)'MEAN_VELY.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, NSAVES
 write(303,'(2(e17.7))') TIME_VEC(J), MEAN_VY(J)
end do

!- close file
close(303)

!---------- z-component-----------
!!-Print filename
write(FILENAME,10200)'MEAN_VELZ.dat'

!!- ASCII
open(unit=304,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, NSAVES
 write(304,'(2(e17.7))') TIME_VEC(J), MEAN_VZ(J)
end do

!- close file
close(304)

print*,'SAVE MEAN VEL--->  OK '


print*,' '
print*,'-------------------END STAT VELOCITIES--------------------------- '
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

end subroutine MEAN_VEL
