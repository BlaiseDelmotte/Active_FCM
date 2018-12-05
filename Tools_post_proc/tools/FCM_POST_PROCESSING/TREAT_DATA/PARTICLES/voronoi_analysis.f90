!!====================================================================
!!
!!
!!====================================================================

subroutine VORONOI_ANALYSIS(NSAVES, &
                            NPART, &
                            TIME_VEC, &
                            PSWIM, &
                            MEAN_PSWIM, &
                            PSWIM_SPH, &
                            MEAN_PSWIM_SPH )

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
integer, intent(in) :: NPART
! Physical Time
real(kind=8), dimension(NSAVES), intent(in) :: TIME_VEC
! Swimmers orientation
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: PSWIM
! Mean Swimmers orientation along time
real(kind=8), dimension(NSAVES,3), intent(out) :: MEAN_PSWIM
! Swimmers orientation along time in spherical coordinates (theta,phi)
real(kind=8), dimension(NSAVES,NPART,2), intent(out) :: PSWIM_SPH
! Mean Swimmers orientation along time in spherical coordinates (theta,phi)
real(kind=8), dimension(NSAVES,2), intent(out) :: MEAN_PSWIM_SPH

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

!- Polar order
real(kind=8), dimension(NSAVES) :: POLAR_ORDER_TIME

!- Module of Mean Pswim for one time
real(kind=8) :: MOD_MEAN_PSWIM

!- Pi
real(kind=8) :: PPI

!- File name 
character(len=40) :: FILENAME


!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_POLAR



PPI = 4.0*datan(1.0)

!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE POLAR ORDER ALONG TIME
!=====================================================================
print*,' '
print*,'-------------------START POLAR ORDER--------------------------- '
print*,' '
PSWIM_SPH = 0.0
do IND = 1, NSAVES

 MEAN_PSWIM(IND,1:3) = 0.0
 MEAN_PSWIM_SPH(IND,1:2) = 0.0
 
 do I = 1, NPART 
  MEAN_PSWIM(IND,1:3) = MEAN_PSWIM(IND,1:3) + PSWIM(IND,I,1:3)
!~   print*,'I = ', I
!~   print*,'PSWIM(IND,I,1:3) = ', PSWIM(IND,I,1:3) 
!~   print*,'MEAN_PSWIM(IND,1:3) = ', MEAN_PSWIM(IND,1:3) 
  ! In spherical coordinates (theta,phi)
  PSWIM_SPH(IND,I,1) = dacos( PSWIM(IND,I,3) )  
  
!~   PSWIM_SPH(IND,I,2) = datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) )
  PSWIM_SPH(IND,I,2) = datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) )
  if (PSWIM_SPH(IND,I,2)<0) then
   PSWIM_SPH(IND,I,2) = PSWIM_SPH(IND,I,2) + 2.0*PPI
  end if
!~   print*,'PSWIM(IND,I,1), PSWIM(IND,I,2) = ', PSWIM(IND,I,1), PSWIM(IND,I,2)
!~   print*,'datan2( PSWIM(IND,I,1), PSWIM(IND,I,2) ) = ', datan2( PSWIM(IND,I,1), PSWIM(IND,I,2) )
!~   print*,'datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) ) = ', datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) )
  
!~   print*,'PSWIM_SPH(IND,I,1:2) = ', PSWIM_SPH(IND,I,1:2) 
!~   read(*,*)
  
  MEAN_PSWIM_SPH(IND,1:2) = MEAN_PSWIM_SPH(IND,1:2) + PSWIM_SPH(IND,I,1:2)
  
  
 end do
 
 MOD_MEAN_PSWIM = dsqrt( MEAN_PSWIM(IND,1)**2 &
                       + MEAN_PSWIM(IND,2)**2 &
                       + MEAN_PSWIM(IND,3)**2 )
                       
 POLAR_ORDER_TIME(IND) = MOD_MEAN_PSWIM/ real(NPART)
 
 MEAN_PSWIM(IND,1:3) = MEAN_PSWIM(IND,1:3)/ real(NPART)
 
!~  read(*,*)
 MEAN_PSWIM_SPH(IND,1:2) = MEAN_PSWIM_SPH(IND,1:2)/ real(NPART)
 
!~  print*,'IND = ', IND
!~  print*,'MEAN_PSWIM_SPH(IND,1:2) = ', MEAN_PSWIM_SPH(IND,1:2)
!~  print*,'MOD_MEAN_PSWIM = ', MOD_MEAN_PSWIM
!~  print*,'NPART = ', NPART
!~  print*,'MOD_MEAN_PSWIM/ real(NPART) = ', MOD_MEAN_PSWIM/ real(NPART)
!~  read(*,*)
 
end do

print*,'COMPUTE POLAR_ORDER_TIME--->  OK '


!!-Print filename
write(FILENAME,10200)'POLAR_ORDER.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(301,'(2(e17.7))') TIME_VEC(IND), POLAR_ORDER_TIME(IND)
end do

!- close file
close(301)


!!-Print filename
write(FILENAME,10200)'MEAN_PSWIM.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(302,'(4(e17.7))') TIME_VEC(IND), MEAN_PSWIM(IND,1), MEAN_PSWIM(IND,2), MEAN_PSWIM(IND,3)
end do

!- close file
close(302)

!!-Print filename
write(FILENAME,10200)'MEAN_PSWIM_SPH.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(303,'(3(e17.7))') TIME_VEC(IND), MEAN_PSWIM_SPH(IND,1), MEAN_PSWIM_SPH(IND,2)
end do

!- close file
close(303)

print*,'SAVE POLAR_ORDER_TIME--->  OK '

print*,' '
print*,'-------------------END POLAR ORDER------------------------------ '
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

end subroutine VORONOI_ANALYSIS
