!!====================================================================
!!
!!
!!====================================================================

subroutine VEL_AUTOCORREL(NSAVES, SAVE_START, NPART, DTIME, VEL)

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
! Where to start
integer, intent(in) :: SAVE_START
! Number of swimmers
integer, intent(in) :: NPART
! Time step
real(kind=8), intent(in) :: DTIME
! Particle velocities
real(kind=8), dimension(NSAVES-SAVE_START,NPART,3), intent(in) :: VEL
! Particle Rotations

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Vel autocorrel
real(kind=8), dimension(NPART,3) :: MEAN_AUTOCORREL_T

real(kind=8) :: MEAN_DIFF_T
real(kind=8) :: MEAN_DIFF_R

!- File name 
character(len=40) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2

!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_DIFF

print*,' '
print*,'-------------------START AUTOCORREL------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE DIFF_T AND DIFF_R
!=====================================================================



MEAN_AUTOCORREL_T = 0.0

print*, 'NSAVES = ', NSAVES
print*, 'NPART = ', NPART
do K = 1,NSAVES-SAVE_START
 do I = 1, NPART


  MEAN_AUTOCORREL_T(I,1) = MEAN_AUTOCORREL_T(I,1) &
                         + VEL(K,I,1)**2
  MEAN_AUTOCORREL_T(I,2) = MEAN_AUTOCORREL_T(I,2) &
                         + VEL(K,I,2)**2
  MEAN_AUTOCORREL_T(I,3) = MEAN_AUTOCORREL_T(I,3) &
                         + VEL(K,I,3)**2
 end do

end do

print*, 'real(NSAVES-SAVE_START-1)= ', real(NSAVES-SAVE_START-1)


MEAN_AUTOCORREL_T = MEAN_AUTOCORREL_T/real(NSAVES-SAVE_START-1)

print*, 'COMPUTE AUTOCORREL---> OK'


!=====================================================================
! 3. WRITE AUTOCORREL
!=====================================================================

!---------- MEAN_AUTOCORREL_T -----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILENAME,10102)'MEAN_AUTOCORREL_T_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII


do I = 1,NPART
 write(301,'(4(e17.7))') real(I), &
                         MEAN_AUTOCORREL_T(I,1), &
                         MEAN_AUTOCORREL_T(I,2), &
                         MEAN_AUTOCORREL_T(I,3) 
end do

!- close file
close(301)


print*,'SAVE DIFF_T AND DIFF_R--->  OK '


print*,' '
print*,'-------------------END DIFF_T AND DIFF_R--------------------------- '
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
10102 format (A,A,A,A,A)
end subroutine VEL_AUTOCORREL
