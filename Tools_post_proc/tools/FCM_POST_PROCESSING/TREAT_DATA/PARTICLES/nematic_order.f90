!!====================================================================
!!
!!
!!====================================================================

subroutine NEMATIC_ORDER(NSAVES, NPART, TIME_VEC, PSWIM)

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

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- NEMATIC order
real(kind=8), dimension(NSAVES,6) :: NEMATIC_ORDER_TIME

!- Mean Pswim for one time
real(kind=8), dimension(6) :: MEAN_NEMATIC

!- File name 
character(len=40) :: FILENAME


!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_NEMATIC





!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE NEMATIC ORDER ALONG TIME
!=====================================================================
print*,' '
print*,'-------------------START NEMATIC ORDER--------------------------- '
print*,' '

do IND = 1, NSAVES

 MEAN_NEMATIC = 0.0
 
 do I = 1, NPART
  MEAN_NEMATIC(1)= MEAN_NEMATIC(1) + PSWIM(IND,I,1)**2 -1.0/3.0
  MEAN_NEMATIC(2)= MEAN_NEMATIC(2) + PSWIM(IND,I,1)*PSWIM(IND,I,2)
  MEAN_NEMATIC(3)= MEAN_NEMATIC(3) + PSWIM(IND,I,1)*PSWIM(IND,I,3)
  MEAN_NEMATIC(4)= MEAN_NEMATIC(4) + PSWIM(IND,I,2)**2 -1.0/3.0
  MEAN_NEMATIC(5)= MEAN_NEMATIC(5) + PSWIM(IND,I,2)*PSWIM(IND,I,3)
  MEAN_NEMATIC(6)= MEAN_NEMATIC(6) + PSWIM(IND,I,3)**3 -1.0/3.0
 end do
 

 NEMATIC_ORDER_TIME(IND,1:6) = MEAN_NEMATIC/ real(NPART)
 
!~  print*,'IND = ', IND
!~  print*,'MEAN_PSWIM = ', MEAN_PSWIM
!~  print*,'MOD_MEAN_PSWIM = ', MOD_MEAN_PSWIM
!~  print*,'NPART = ', NPART
!~  print*,'MOD_MEAN_PSWIM/ real(NPART) = ', MOD_MEAN_PSWIM/ real(NPART)
!~  read(*,*)
 
end do

print*,'COMPUTE NEMATIC_ORDER_TIME--->  OK '


!!-Print filename
write(FILENAME,10200)'NEMATIC_ORDER.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(301,'(7(e17.7))') TIME_VEC(IND), (NEMATIC_ORDER_TIME(IND,I) , I=1,6)
end do

!- close file
close(301)



print*,'SAVE NEMATIC_ORDER_TIME--->  OK '

print*,' '
print*,'-------------------END NEMATIC ORDER------------------------------ '
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

end subroutine NEMATIC_ORDER
