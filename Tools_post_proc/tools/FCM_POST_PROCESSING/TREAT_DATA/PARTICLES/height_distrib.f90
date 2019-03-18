!!====================================================================
!!
!!
!!====================================================================

subroutine HEIGHT_DISTRIB(NSAVES, &
                          SAVE_START, &
                          PART_START, &
                          PART_END, &
                          LX, &
                          RADMAX, &
                          POS )

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
! Height of the domain
real(kind=8), intent(in) :: LX
! Max radius
real(kind=8), intent(in) :: RADMAX
! Particle velocities
real(kind=8), dimension(NSAVES-SAVE_START,PART_END-PART_START+1,3), intent(in) :: POS

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Height distrib and its standars deviation for error bars
real(kind=8), allocatable, dimension(:,:) :: H_DISTRIB_TIME
real(kind=8), allocatable, dimension(:) :: H_DISTRIB
real(kind=8), allocatable, dimension(:) :: STD_H_DISTRIB
real(kind=8), allocatable, dimension(:) :: RANGE_H


!- Step size to discretize the interval
real(kind=8) :: DH

!- Number of steps to discretize the interval 
integer :: NH

!- Indices for distrib
integer :: IND_H

!- File name 
character(len=100) :: FILENAME
!- String for saves
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4


!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_HEIGHT

print*,' '
print*,'-------------------START HEIGHT DISTRIB------------------------ '
print*,' '
!---------------------------------------------------------------------

 
!=====================================================================
! 2. DEFINE DISCRETIZATION OF H_DISTRIB
!=====================================================================
print*, 'DISCRETIZE RANGE OF VELOCITIES FOR PDF'
DH = 0.02*RADMAX
NH = (LX/2.0 - 0.0)/DH

print*,'LX = ', LX
print*,'DH = ', DH
print*,'RADMAX = ', RADMAX
print*,'NH = '
print*, NH
        

allocate(H_DISTRIB_TIME(NSAVES-SAVE_START,NH))
allocate(H_DISTRIB(NH))
allocate(STD_H_DISTRIB(NH))
allocate(RANGE_H(NH))


RANGE_H(1) = 0.0
do J = 2, NH
 RANGE_H(J) = RANGE_H(J-1) + DH
end do


print*, 'DISCRETIZE RANGE OF HEIGHT FOR DISTRIB---> OK'
print*, 'NSAVES-SAVE_START = ', NSAVES-SAVE_START


!=====================================================================
! 2. COMPUTE DSITRIB
!=====================================================================
print*,'COMPUTE DISTRIB HEIGHT  '

H_DISTRIB = 0.0
STD_H_DISTRIB = 0.0
H_DISTRIB_TIME = 0.0



!~ do IND = SAVE_START, NSAVES 
do K = 1,NSAVES-SAVE_START
 do I = 1, PART_END-PART_START+1
  IND_H = floor( ( POS(K,I,1)-RANGE_H(1) )/DH ) + 1
 
  
  if ((IND_H<1).or.(IND_H>NH)) then
   print*, 'POS(K,I,1) = ', POS(K,I,1)
   print*, 'RANGE_H(1) = ', RANGE_H(1)
   print*, '( POS(K,I,1)-RANGE_H(1) )/DH = ', ( POS(K,I,1)-RANGE_H(1) )/DH
   print*, 'IND_H = ', IND_H   
   print*, 'I = ', I
   print*, 'SAVE # = ', K + SAVE_START
   read(*,*)
  else 
   H_DISTRIB_TIME(K,IND_H) = H_DISTRIB_TIME(K,IND_H) + 1.0
  end if
 end do
 
 H_DISTRIB_TIME(K,:) = H_DISTRIB_TIME(K,:)/(PART_END-PART_START+1)
 H_DISTRIB = H_DISTRIB + H_DISTRIB_TIME(K,:)
end do

H_DISTRIB = H_DISTRIB/(NSAVES-SAVE_START+1)

do K = 1,NSAVES-SAVE_START
 STD_H_DISTRIB = STD_H_DISTRIB + ( H_DISTRIB_TIME(K,:) - H_DISTRIB )**2
end do

STD_H_DISTRIB = STD_H_DISTRIB/(NSAVES-SAVE_START)
STD_H_DISTRIB = dsqrt( STD_H_DISTRIB/(NSAVES-SAVE_START+1) )
 
print*,'COMPUTE H_DISTRIB  --->  OK '

!=====================================================================
! 3. WRITE H_DISTRIB
!=====================================================================

!---------- Height-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

write(FILENAME,10103)'H_DISTRIB_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII


do J = 1, NH
 write(301,'(3(e17.7))') RANGE_H(J)/RADMAX, H_DISTRIB(J), STD_H_DISTRIB(J)
end do


!- close file
close(301)

print*,'SAVE H_DISTRIB--->  OK '


deallocate(H_DISTRIB)
deallocate(RANGE_H)


print*,' '
print*,'-------------------END PDF VELOCITIES--------------------------- '
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
end subroutine HEIGHT_DISTRIB
