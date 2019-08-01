!!====================================================================
!!
!!
!!====================================================================

subroutine PDF_PSWIM_SPHERICAL_COORD(NSAVES, &
                            SAVE_START, &
                            PART_START, &
                            PART_END, &
                            P1_SPH)

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
! Particle Orientations in spherical coord (theta,phi)
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,2), intent(in) :: P1_SPH
!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

!- PDF of theta-component
real(kind=8), allocatable, dimension(:) :: PDF_P1TH
real(kind=8), allocatable, dimension(:) :: RANGE_P1TH

!- PDF of phi-component
real(kind=8), allocatable, dimension(:) :: PDF_P1PHI
real(kind=8), allocatable, dimension(:) :: RANGE_P1PHI

!- PDF of theta-phi-component
real(kind=8), allocatable, dimension(:,:) ::  PDF_P1TH_PHI

!- PDF of theta-phi-component along time
real(kind=8), allocatable, dimension(:,:,:) ::  PDF_P1TH_PHI_TIME


!- Step size to discretize the interval of orientations
real(kind=8) ::  DSCAL, DP1TH, DP1PHI

real(kind=8) :: PPI

!- Number of steps to discretize the interval of orientations
integer ::   N_STEP_P1TH, N_STEP_P1PHI

!- Indices for pdf
integer :: IND_TH, IND_PHI

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4


!- Index
integer :: I, J, K, IND

!- Jump in time to save pdf along time
integer :: JUMP_SAVES

!- Integers to signalize problems when writing
integer :: ERR_FILE_PDF

print*,' '
print*,'-------------------START PDF Orientations------------------------ '
print*,' '

PPI = 4.0*datan(1.0d0)
!---------------------------------------------------------------------
!=====================================================================
! 1. DEFINE DISCRETIZATION OF PDF
!=====================================================================
print*, 'DISCRETIZE RANGE OF Orientations FOR PDF'


N_STEP_P1TH = 60
N_STEP_P1PHI = 60

DP1TH = PPI/(real(N_STEP_P1TH))
DP1PHI = 2.0*PPI/(real(N_STEP_P1PHI))    


print*,' N_STEP_P1TH ,N_STEP_P1PHI = '
print*,  N_STEP_P1TH ,N_STEP_P1PHI
        

allocate(PDF_P1TH(N_STEP_P1TH))
allocate(PDF_P1PHI(N_STEP_P1PHI))
allocate(PDF_P1TH_PHI(N_STEP_P1TH,N_STEP_P1PHI))
allocate(PDF_P1TH_PHI_TIME(NSAVES,N_STEP_P1TH,N_STEP_P1PHI))

allocate(RANGE_P1TH(N_STEP_P1TH))
allocate(RANGE_P1PHI(N_STEP_P1PHI))

RANGE_P1TH(1) = 0.0 
do J = 2, N_STEP_P1TH
 RANGE_P1TH(J) = RANGE_P1TH(J-1) + DP1TH
end do

RANGE_P1PHI(1) = 0.0
do J = 2, N_STEP_P1PHI
 RANGE_P1PHI(J) = RANGE_P1PHI(J-1) + DP1PHI
end do

print*, 'DISCRETIZE RANGE OF Orientations FOR PDF---> OK'


!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE Orientations PDF  '

PDF_P1TH = 0.0
PDF_P1PHI = 0.0
PDF_P1TH_PHI = 0.0
PDF_P1TH_PHI_TIME = 0.0
K = 0


do IND = 1, NSAVES 
!~ do IND = 1, 1

 if (IND>SAVE_START-1) then
  K = K+1
 end if
 
 do I = 1, PART_END-PART_START+1

  
  
  if (P1_SPH(IND,I,1)==RANGE_P1TH(1)) then
   IND_TH = 1
  else
   IND_TH = ceiling( ( P1_SPH(IND,I,1)-RANGE_P1TH(1) )/DP1TH ) 
  end if
  if (P1_SPH(IND,I,2)==RANGE_P1PHI(1)) then
   IND_PHI = 1
  else
   IND_PHI = ceiling( ( P1_SPH(IND,I,2)-RANGE_P1PHI(1) )/DP1PHI ) 
  end if
  
   
 
  
  if ((IND_TH<1).or.(IND_TH>N_STEP_P1TH)) then
   print*, 'P1_SPH(IND,I,1) = ', P1_SPH(IND,I,1)
   print*, 'RANGE_P1TH(1) = ', RANGE_P1TH(1)
   print*, '( P1_SPH(IND,I,1)-RANGE_P1TH(1) )/DP1TH = ', ( P1_SPH(IND,I,1)-RANGE_P1TH(1) )/DP1TH
   print*, 'IND_TH = ', IND_TH
!~    read(*,*)
  end if
  
  if ((IND_PHI<1).or.(IND_PHI>N_STEP_P1PHI)) then
   print*, 'P1_SPH(IND,I,2) = ', P1_SPH(IND,I,2)
   print*, 'RANGE_P1PHI(1) = ', RANGE_P1PHI(1)
   print*, '( P1_SPH(IND,I,2)-RANGE_P1PHI(1) )/DP1PHI = ', ( P1_SPH(IND,I,2)-RANGE_P1PHI(1) )/DP1PHI
   print*, 'IND_PHI = ', IND_PHI
!~    read(*,*)
  end if
  
  if (IND>SAVE_START-1) then
  
   PDF_P1TH(IND_TH) = PDF_P1TH(IND_TH) + 1.0
   PDF_P1PHI(IND_PHI) = PDF_P1PHI(IND_PHI) + 1.0
   PDF_P1TH_PHI(IND_TH,IND_PHI) =  PDF_P1TH_PHI(IND_TH,IND_PHI) + 1.0 

   
  end if
  
  PDF_P1TH_PHI_TIME(IND,IND_TH,IND_PHI) =  PDF_P1TH_PHI_TIME(IND,IND_TH,IND_PHI) + 1.0
  
 end do
   
end do

print*,'COMPUTE Orientations PDF --->  OK '

!=====================================================================
! 3. WRITE PDF
!=====================================================================
!---------- Scalar product with mean vector-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

!---------- theta-component-----------
!!-Print filename
write(FILENAME,10103) 'PDF_P1THETA_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'
!!- ASCII
open(unit=305,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_P1TH
 write(305,'(2(e17.7))') RANGE_P1TH(J), PDF_P1TH(J)
end do

!- close file
close(305)

!---------- phi-component-----------
!!-Print filename
write(FILENAME,10103) 'PDF_P1PHI_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=305,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_P1PHI
 write(305,'(2(e17.7))') RANGE_P1PHI(J), PDF_P1PHI(J)
end do

!- close file
close(305)

!---------- theta-phi-component-----------
!!-Print filename
write(FILENAME,10103) 'PDF_P1TH_PHI_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'
!!- ASCII
open(unit=305,file=trim(FILENAME))
!!- Ecriture ASCII
do I = 1, N_STEP_P1TH
 do J = 1, N_STEP_P1PHI
  write(305,'(3(e17.7))') RANGE_P1TH(I), RANGE_P1PHI(J), PDF_P1TH_PHI(I,J)
 end do
end do
!- close file

!~ print*,' PDF_P1TH_PHI(:,1) =', PDF_P1TH_PHI(:,1)

close(305)

 !---------- theta-phi-component along time-----------
!!-Print filename
write(FILENAME,10103) 'PDF_P1TH_PHI_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'
JUMP_SAVES = 10
!!- ASCII
open(unit=305,file=trim(FILENAME))
!!- Ecriture ASCII
write(305,'(4(I8))') JUMP_SAVES, NSAVES, N_STEP_P1TH, N_STEP_P1PHI

do IND =  1, NSAVES, JUMP_SAVES 
!~  print*,' PDF_P1TH_PHI_TIME(IND,:,1) =', PDF_P1TH_PHI_TIME(IND,:,1)
!~  read(*,*)
 if (IND<NSAVES) then
  !- Average over jumps
  do K = 1, JUMP_SAVES-1
   PDF_P1TH_PHI_TIME(IND,:,:) = PDF_P1TH_PHI_TIME(IND,:,:) + PDF_P1TH_PHI_TIME(IND + K,:,:)
  end do
 end if
 
 do I = 1, N_STEP_P1TH
  do J = 1, N_STEP_P1PHI
   write(305,'(4(e17.7))') real(IND), RANGE_P1TH(I), RANGE_P1PHI(J), PDF_P1TH_PHI_TIME(IND,I,J)
  end do
 end do
end do
!- close file
close(305) 

print*,'SAVE PDF Orientations--->  OK '


deallocate(PDF_P1TH)
deallocate(PDF_P1PHI)
deallocate(RANGE_P1TH)
deallocate(RANGE_P1PHI)

print*,' '
print*,'-------------------END PDF Orientations--------------------------- '
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
10103 format (A,A,A,A,A,A,A,A,A)

end subroutine PDF_ORIENTATIONS
