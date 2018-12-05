!!====================================================================
!!
!!
!!====================================================================

subroutine RDF_ONLY(SAVE_START, &
               SAVE_END, &               
               PART_START, &
               PART_END, &
               L, &
               RAD, &
               POSI)

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
! Domain size
real(kind=8), intent(in):: L
! MAximal particle radius
real(kind=8), intent(in):: RAD
! Particle positions with no periodicity
real(kind=8), dimension(SAVE_END-SAVE_START+1,PART_END-PART_START+1,3), intent(in) :: POSI


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- r-component
real(kind=8), allocatable, dimension(:) :: RANGE_R, VOL_SHELL
!- g(r)
real(kind=8), allocatable, dimension(:,:) :: RDF_R

real(kind=8), allocatable, dimension(:) :: RDF_R_MEAN

!-Min/Maximal distance accepted for RDF computation
real(kind=8) ::  MIN_DIST, MAX_DIST, RADMAX, TWORADMAX

real(kind=8),dimension(3) :: POS_I, POS_J, POS_IJ

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4

!- Range_r discretization
integer :: NR
real(kind=8) :: DR
real(kind=8) :: RIJ

real(kind=8) :: PPI
!- Index
integer :: I, J, K, M, IND, IND_I, IND_J, IND_R

!- Integers to signalize problems when writing
integer :: ERR_FILE_DIFF

print*,' '
print*,'-------------------START RDF------------------------ '
print*,' '

PPI =  4.0*atan(1.0)

RADMAX = RAD
TWORADMAX = 2.0*RADMAX
MIN_DIST = TWORADMAX*0.5
MAX_DIST = L/2.0
DR = 0.1*RADMAX
NR = ceiling((MAX_DIST - MIN_DIST)/DR)

print*,'NR = '
print*, NR

allocate(RANGE_R(NR))
allocate(VOL_SHELL(NR))
allocate(RDF_R(SAVE_END-SAVE_START+1,NR))
allocate(RDF_R_MEAN(NR))
!---------------------------------------------------------------------
!=====================================================================
! 1. COMUTE RDF
!=====================================================================
print*, 'COMPUTE RDF'


print*, 'SAVE_START, SAVE_END = ', SAVE_START, SAVE_END
print*, 'PART_START, PART_END = ', PART_START, PART_END

RANGE_R(1) = MIN_DIST
do J = 2, NR
 RANGE_R(J) = RANGE_R(J-1) + DR
end do

do I = 1, NR
 VOL_SHELL(I) = &
     4.0/3.0*PPI*( (RANGE_R(I)+DR)**3 - RANGE_R(I)**3 )
end do

RDF_R = 0.0
RDF_R_MEAN = 0.0



do M = 1, SAVE_END-SAVE_START+1

 do IND_I = 1, PART_END-PART_START
 
  POS_I = POSI(M,IND_I,1:3)
  do IND_J = IND_I+1, PART_END-PART_START+1
			POS_J = POSI(M,IND_J,1:3)
			POS_IJ = POS_I - POS_J   
			POS_IJ(1) = POS_IJ(1) - L* real(int(POS_IJ(1)/(L/2.0)))   
			POS_IJ(2) = POS_IJ(2) - L* real(int(POS_IJ(2)/(L/2.0)))  
			POS_IJ(3) = POS_IJ(3) - L* real(int(POS_IJ(3)/(L/2.0)))
			
			RIJ = dsqrt(POS_IJ(1)**2 + POS_IJ(2)**2 + POS_IJ(3)**2)
! 			print*,'M, IND_I, IND_J = ', M, IND_I, IND_J
! 			print*,'POS_IJ = ', POS_IJ
! 			print*,'RIJ = ', RIJ
! 			read(*,*)
			
			if ((RIJ>=MIN_DIST).and.(RIJ<MAX_DIST)) then
		 	IND_R = ceiling( ( RIJ-RANGE_R(1) )/DR ) 
		 	RDF_R(M,IND_R) = RDF_R(M,IND_R) + 1.0 
			end if
		end do			
 end do
 RDF_R(M,:) = RDF_R(M,:)/(real(PART_END-PART_START+1)*real(PART_END-PART_START)/2.0*VOL_SHELL/L**(3))

 RDF_R_MEAN = RDF_R_MEAN + RDF_R(M,:)
end do

RDF_R_MEAN =  RDF_R_MEAN/real(SAVE_END-SAVE_START + 1)
print*,'maxval(abs(RDF_R_MEAN) ) = ',maxval(abs(RDF_R_MEAN) )



!=====================================================================
! 3. WRITE RDF
!=====================================================================

!---------- RDF_MEAN-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) SAVE_END
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

write(FILENAME,10103)'RDF_MEAN_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1,  NR
 write(301,'(2(e17.7))') RANGE_R(I), RDF_R_MEAN(I)
end do

!- close file
close(301)

!---------- RDF TIME-----------
write(FILENAME,10103)'RDF_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII
do J = 1, SAVE_END - SAVE_START +1
	do I = 1,  NR
		write(301,'(3(e17.7))'), real(J), RANGE_R(I), RDF_R(J,I)
	end do
end do
!- close file
close(301)



print*,'SAVE RDF --> OK '

print*,' '
print*,'-------------------END RDF--------------------------- '
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

end subroutine RDF_ONLY
