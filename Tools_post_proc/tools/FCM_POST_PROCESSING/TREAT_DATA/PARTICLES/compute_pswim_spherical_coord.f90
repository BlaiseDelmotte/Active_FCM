!!====================================================================
!!
!!
!!====================================================================

subroutine COMPUTE_PSWIM_SPHERICAL_COORD(NSAVES, &
                       PART_START, &
                       PART_END, &
                       PSWIM, &
                       PSWIM_SPH)
                       

!!====================================================================
!!
!!
!!====================================================================

use MPI

implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer, intent(in) :: NSAVES
! Number of swimmers
integer, intent(in) :: PART_START, PART_END
! Swimmers orientation
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: PSWIM
! Swimmers orientation along time in spherical coordinates (theta,phi)
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,2), intent(out) :: PSWIM_SPH


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Pi
real(kind=8) :: PPI

!- File name 
character(len=40) :: FILENAME
character(len=20) :: MYFMT

integer :: I, J, K, IND


PPI = 4.0*datan(1.0d0)

!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE ORIENTATION VECTOR IN SPHERICAL COORDINATES
!=====================================================================
print*,' '
print*,'-------------------START PSWIM_SPH--------------------------- '
print*,' '


PSWIM_SPH = 0.0

do IND = 1, NSAVES
 do I = 1, PART_END-PART_START+1 
 
  ! In spherical coordinates (theta,phi)
  PSWIM_SPH(IND,I,1) = dacos( PSWIM(IND,I,3) ) 
  
  PSWIM_SPH(IND,I,2) = datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) )
  if (PSWIM_SPH(IND,I,2)<0) then
   PSWIM_SPH(IND,I,2) = PSWIM_SPH(IND,I,2) + 2.0*PPI
  end if
 end do
 
end do

print*,' '
print*,'-------------------END PSWIM_SPH------------------------------ '
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
10201 format (A,I1,A)
10202 format (A,I2,A)
10203 format (A,I3,A)
10204 format (A,I4.4,A)
10205 format (I8.8)
10101 format (A,A,A)

end subroutine COMPUTE_PSWIM_SPHERICAL_COORD
