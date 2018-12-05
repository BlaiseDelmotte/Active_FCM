!!====================================================================
!!
!!
!!====================================================================

subroutine MSD_ONLY(SAVE_START, &
               SAVE_END, &               
               PART_START, &
               PART_END, &
               NSAVES_MSD, &
               POS_NOPER)

!!====================================================================
!!
!!
!!====================================================================



implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! First save to use
integer, intent(in) :: SAVE_START
! Last save to use
integer, intent(in) :: SAVE_END
! Number of swimmers
integer, intent(in) :: PART_START, PART_END
! Time window over which to compute MSD
integer, intent(in) :: NSAVES_MSD
! Particle positions with no periodicity
real(kind=8), dimension(SAVE_END-SAVE_START+1,PART_END-PART_START+1,3), intent(in) :: POS_NOPER


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------



!- MSD and Displacement
real(kind=8), dimension(NSAVES_MSD) :: MEAN_MSD
real(kind=8), dimension(SAVE_END-SAVE_START+1,3) :: POSCOM ! Center of mass position
real(kind=8), dimension(NSAVES_MSD,PART_END-PART_START+1,3) :: DISP
real(kind=8), dimension(NSAVES_MSD) :: NDT

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4

!- Index
integer :: JNSAVES

!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_DIFF

print*,' '
print*,'-------------------START MEANS SQUARED DISPLACEMENT------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE MSD
!=====================================================================
print*, 'COMPUTE MSD'


MEAN_MSD = 0.0
DISP = 0.0
NDT = 0.0
POSCOM = 0.0


print*, 'NSAVES_MSD = ', NSAVES_MSD

print*, 'SAVE_START, SAVE_END = ', SAVE_START, SAVE_END
print*, 'PART_START, PART_END = ', PART_START, PART_END


do K = 1, SAVE_END-SAVE_START+1 


 do I = 1, PART_END-PART_START+1
  POSCOM(K,1) = POSCOM(K,1) + POS_NOPER(K,I,1)  
  POSCOM(K,2) = POSCOM(K,2) + POS_NOPER(K,I,2) 
  POSCOM(K,3) = POSCOM(K,3) + POS_NOPER(K,I,3) 
 end do

 POSCOM(K,1:3) = POSCOM(K,1:3)/real(PART_END-PART_START+1)
end do

do K = 1, SAVE_END-SAVE_START
 JNSAVES = min(K+NSAVES_MSD,SAVE_END-SAVE_START+1)
 
 do J = K+1, JNSAVES
 
!  print*,'K,J = ', K,J 
 
  do I = 1, PART_END-PART_START+1
  

			DISP(J-K,I,1:3) = POS_NOPER(J,I,1:3)-POSCOM(J,1:3) &
			              - ( POS_NOPER(K,I,1:3)-POSCOM(K,1:3) )
			       
!   print*, 'DISP(J-K,I,1:3) = ', DISP(J-K,I,1:3)
!   read(*,*)
   MEAN_MSD(J-K) = MEAN_MSD(J-K) &
               + DISP(J-K,I,1)**2 + DISP(J-K,I,2)**2 + DISP(J-K,I,3)**2   
  end do
  
  NDT(J-K) = NDT(J-K) + real(PART_END-PART_START+1)  
  
  
!   print*,'MEAN_MSD(J-K) = ', MEAN_MSD(J-K)
!   
!   print*,'NDT(J-K) = ', NDT(J-K)
!   read(*,*)
 end do
 
end do



MEAN_MSD = MEAN_MSD/NDT

print*, 'COMPUTE MSD---> OK'

!=====================================================================
! 3. WRITE DIFF_T AND DIFF_R
!=====================================================================

!---------- DIFF_T-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) SAVE_END
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

write(FILENAME,10103)'MSD_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1,  NSAVES_MSD
 write(301,'(2(e17.7))') real(I), MEAN_MSD(I)
end do

!- close file
close(301)



print*,'SAVE MSD --> OK '

print*,' '
print*,'-------------------END MSD AND DSIP--------------------------- '
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

end subroutine MSD_ONLY
