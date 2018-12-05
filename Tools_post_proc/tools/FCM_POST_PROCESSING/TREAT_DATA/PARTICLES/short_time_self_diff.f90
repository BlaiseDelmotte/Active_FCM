!!====================================================================
!!
!!
!!====================================================================

subroutine SHORT_TIME_SELF_DIFF(NSAVES, &
                                SAVE_START, &
                                PART_START, &
                                PART_END, &
                                DTIME, &
                                VEL, &
                                ROT)

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
integer, intent(in) :: PART_START, PART_END
! Time step
real(kind=8), intent(in) :: DTIME
! Particle velocities
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: VEL
! Particle Rotations
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: ROT

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Short Time Self_Diff
real(kind=8), dimension(NSAVES-SAVE_START+1) :: DIFF_T
real(kind=8), dimension(PART_END-PART_START+1,3) :: MEAN_AUTOCORREL_T
real(kind=8), dimension(NSAVES-SAVE_START+1) :: DIFF_R
real(kind=8) :: MEAN_DIFF_T
real(kind=8) :: MEAN_DIFF_R

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4
!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_DIFF

print*,' '
print*,'-------------------START SHORT TIME SELF DIFF------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE DIFF_T AND DIFF_R
!=====================================================================
print*, 'COMPUTE DIFF_T AND DIFF_R'

DIFF_T = 0.0
DIFF_R = 0.0

MEAN_DIFF_T = 0.0
MEAN_DIFF_R = 0.0

MEAN_AUTOCORREL_T = 0.0

print*, 'NSAVES = ', NSAVES
print*, 'PART_START, PART_END = ', PART_START, PART_END

K = 0
do IND = SAVE_START, NSAVES
 K = K+1
 do I = 1, PART_END-PART_START+1
  DIFF_T(K) = DIFF_T(K) &
            + VEL(IND,I,1)**2 &
            + VEL(IND,I,2)**2 &
            + VEL(IND,I,3)**2
            
  DIFF_R(K) = DIFF_R(K) &
            + ROT(IND,I,1)**2 &
            + ROT(IND,I,2)**2 &
            + ROT(IND,I,3)**2

  MEAN_AUTOCORREL_T(I,1) = MEAN_AUTOCORREL_T(I,1) &
                         + VEL(IND,I,1)**2
  MEAN_AUTOCORREL_T(I,2) = MEAN_AUTOCORREL_T(I,2) &
                         + VEL(IND,I,2)**2
  MEAN_AUTOCORREL_T(I,3) = MEAN_AUTOCORREL_T(I,3) &
                         + VEL(IND,I,3)**2
 end do
 DIFF_T(K) = DIFF_T(K)/(6.0*real(PART_END-PART_START+1))*DTIME
 DIFF_R(K) = DIFF_R(K)/(6.0*real(PART_END-PART_START+1))*DTIME
 
 MEAN_DIFF_T = MEAN_DIFF_T + DIFF_T(K)
 MEAN_DIFF_R = MEAN_DIFF_R + DIFF_R(K)
 
end do

print*, 'real(K-1)= ', real(K-1)

MEAN_DIFF_T = MEAN_DIFF_T/real(K-1)
MEAN_DIFF_R = MEAN_DIFF_R/real(K-1)

MEAN_AUTOCORREL_T = MEAN_AUTOCORREL_T/real(K-1)

print*, 'COMPUTE DIFF_T AND DIFF_R---> OK'


!=====================================================================
! 3. WRITE DIFF_T AND DIFF_R
!=====================================================================

!---------- DIFF_T-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

write(FILENAME,10103)'SHORT_DIFF_T_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

write(301,'(2(e17.7))') real(NSAVES-SAVE_START+1), MEAN_DIFF_T
K = 0
do IND = SAVE_START, NSAVES
 K=K+1
 write(301,'(2(e17.7))') real(IND), DIFF_T(K)
end do

!- close file
close(301)

!---------- DIFF_R-----------
!!-Print filename
write(FILENAME,10103)'SHORT_DIFF_R_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

write(302,'(2(e17.7))') 9999999.0, MEAN_DIFF_R
K = 0
do IND = SAVE_START, NSAVES
 K=K+1
 write(302,'(2(e17.7))') real(IND), DIFF_R(K)
end do

!- close file
close(302)


!---------- MEAN_AUTOCORREL_T -----------
!!-Print filename
write(FILENAME,10103)'MEAN_AUTOCORREL_T_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII


do I = 1,PART_END-PART_START+1
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
10103 format (A,A,A,A,A,A,A,A,A)

end subroutine SHORT_TIME_SELF_DIFF
