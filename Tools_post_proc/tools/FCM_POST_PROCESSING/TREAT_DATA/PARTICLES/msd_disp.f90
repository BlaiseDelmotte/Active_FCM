!!====================================================================
!!
!!
!!====================================================================

subroutine MSD_DISP(NSAVES, &
               SAVE_START, &
               PART_START, &
               PART_END, &
               NSAVES_MSD, &
               RAD, &
               LXMAX, &
               LYMAX, &
               LZMAX, &
               POS_NOPER)

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
! Time interval over compute disp and msd
integer, intent(in) :: NSAVES_MSD
! Radius of particles considered
real(kind=8), intent(in) :: RAD
! Domain size
real(kind=8), intent(in) :: LXMAX
real(kind=8), intent(in) :: LYMAX
real(kind=8), intent(in) :: LZMAX
! Particle positions with no periodicity
real(kind=8), dimension(NSAVES-SAVE_START+1,PART_END-PART_START+1,3), intent(inout) :: POS_NOPER


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

!- PDF of norm
real(kind=8), allocatable, dimension(:,:) :: PDF_DISPNORM
real(kind=8), allocatable, dimension(:) :: RANGE_DISPNORM
!- PDF of x-component
real(kind=8), allocatable, dimension(:,:) :: PDF_DISPX
real(kind=8), allocatable, dimension(:) :: RANGE_DISPX
!- PDF of y-component
real(kind=8), allocatable, dimension(:,:) :: PDF_DISPY
real(kind=8), allocatable, dimension(:) :: RANGE_DISPY
!- PDF of z-component
real(kind=8), allocatable, dimension(:,:) :: PDF_DISPZ
real(kind=8), allocatable, dimension(:) :: RANGE_DISPZ

!- MSD and Displacement
real(kind=8), dimension(NSAVES_MSD,3) :: MEAN_MSD_DIR
real(kind=8), dimension(NSAVES_MSD) :: MEAN_MSD
real(kind=8), dimension(NSAVES_MSD) :: NDT
real(kind=8), dimension(3) :: DISP_TEMP
real(kind=8) :: DISPNORM_TEMP

!- Step size to discretize the interval of displacements
real(kind=8) :: DDISPNORM, DDISPX, DDISPY, DDISPZ
real(kind=8) :: MAX_DISPNORM, MAX_DISPX, MAX_DISPY, MAX_DISPZ
!- Number of steps to discretize the interval of displacements
integer :: N_STEP_DISPNORM, N_STEP_DISPX, N_STEP_DISPY, N_STEP_DISPZ

!- Indices for pdf
integer :: IND_NORM, IND_X, IND_Y, IND_Z

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4
!- Index
integer :: I, J, K, IND, JNSAVES, NMISS

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


MEAN_MSD_DIR = 0.0
DISP_TEMP = 0.0
DISPNORM_TEMP = 0.0
MEAN_MSD = 0.0
NDT = 0



print*, 'NSAVES_MSD = ', NSAVES_MSD

print*, 'SAVE_START, SAVE_END = ', SAVE_START, NSAVES
print*, 'PART_START, PART_END = ', PART_START, PART_END



do K = 1, NSAVES-SAVE_START
 JNSAVES = min(K+NSAVES_MSD,NSAVES-SAVE_START+1)
 
 do J = K+1, JNSAVES
 
!  print*,'K,J = ', K,J 
 
  do I = 1, PART_END-PART_START+1
   
   if (K>1) then
    if(abs(POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1))>LXMAX/2.0) then
!!~      print*,'K = ', K
!!~      print*,'I = ', I
!!~      print*,'POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1) = ', POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1)
!!~      print*,'POS_NOPER(K,I,1) = ', POS_NOPER(K,I,1)
!!~      print*,'POS_NOPER(K-1,I,1) = ', POS_NOPER(K-1,I,1)
     POS_NOPER(K,I,1) = POS_NOPER(K,I,1) - sign(LXMAX,POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1))
!!~      print*,'POS_NOPER(K,I,1) = ', POS_NOPER(K,I,1)
!!~      print*,'POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1) = ', POS_NOPER(K,I,1)-POS_NOPER(K-1,I,1)
!!~      read(*,*)
    end if
    if(abs(POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2))>LYMAX/2.0) then
!!~      print*,'K = ', K
!!~      print*,'I = ', I
!!~      print*,'POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2) = ', POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2)
!!~      print*,'POS_NOPER(K,I,2) = ', POS_NOPER(K,I,2)
!!~      print*,'POS_NOPER(K-1,I,2) = ', POS_NOPER(K-1,I,2)
     POS_NOPER(K,I,2) = POS_NOPER(K,I,2) - sign(LYMAX,POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2))
!!~      print*,'POS_NOPER(K,I,2) = ', POS_NOPER(K,I,2)
!!~      print*,'POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2) = ', POS_NOPER(K,I,2)-POS_NOPER(K-1,I,2)
!!~      read(*,*)
    end if
    if(abs(POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3))>LZMAX/2.0) then
!!~      print*,'K = ', K
!!~      print*,'I = ', I
!!~      print*,'POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3) = ', POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3)
!!~      print*,'POS_NOPER(K,I,3) = ', POS_NOPER(K,I,3)
!!~      print*,'POS_NOPER(K-1,I,3) = ', POS_NOPER(K-1,I,3)
     POS_NOPER(K,I,3) = POS_NOPER(K,I,3) - sign(LZMAX,POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3))
!!~      print*,'POS_NOPER(K,I,3) = ', POS_NOPER(K,I,3)
!!~      print*,'POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3) = ', POS_NOPER(K,I,3)-POS_NOPER(K-1,I,3)
!!~      read(*,*)
    end if
   end if 
   
    if(abs(POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1))>LXMAX/2.0) then
!~      print*,'J = ', J
!~      print*,'I = ', I
!~      print*,'POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1) = ', POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1)
!~      print*,'POS_NOPER(J,I,1) = ', POS_NOPER(J,I,1)
!~      print*,'POS_NOPER(J-1,I,1) = ', POS_NOPER(J-1,I,1)
     POS_NOPER(J,I,1) = POS_NOPER(J,I,1) - sign(LXMAX,POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1))
!~      print*,'POS_NOPER(J,I,1) = ', POS_NOPER(J,I,1)
!~      print*,'POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1) = ', POS_NOPER(J,I,1)-POS_NOPER(J-1,I,1)
!~      read(*,*)
    end if
    if(abs(POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2))>LYMAX/2.0) then
!~      print*,'J = ', J
!~      print*,'I = ', I
!~      print*,'POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2) = ', POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2)
!~      print*,'POS_NOPER(J,I,2) = ', POS_NOPER(J,I,2)
!~      print*,'POS_NOPER(J-1,I,2) = ', POS_NOPER(J-1,I,2)
     POS_NOPER(J,I,2) = POS_NOPER(J,I,2) - sign(LYMAX,POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2))
!~      print*,'POS_NOPER(J,I,2) = ', POS_NOPER(J,I,2)
!~      print*,'POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2) = ', POS_NOPER(J,I,2)-POS_NOPER(J-1,I,2)
!~      read(*,*)
    end if
    if(abs(POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3))>LZMAX/2.0) then
!~      print*,'J = ', J
!~      print*,'I = ', I
!~      print*,'POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3) = ', POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3)
!~      print*,'POS_NOPER(J,I,3) = ', POS_NOPER(J,I,3)
!~      print*,'POS_NOPER(J-1,I,3) = ', POS_NOPER(J-1,I,3)
     POS_NOPER(J,I,3) = POS_NOPER(J,I,3) - sign(LZMAX,POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3))
!~      print*,'POS_NOPER(J,I,3) = ', POS_NOPER(J,I,3)
!~      print*,'POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3) = ', POS_NOPER(J,I,3)-POS_NOPER(J-1,I,3)
!~      read(*,*)
   end if 

   DISP_TEMP(1:3) = POS_NOPER(J,I,1:3)- POS_NOPER(K,I,1:3)
   
   if (isnan(DISP_TEMP(1))) then
     print*,'DISP_TEMP(1) = ', DISP_TEMP(1)
     print*,'POS_NOPER(K,I,1:3) = ', POS_NOPER(K,I,1:3)
     print*,'POS_NOPER(J,I,1:3) = ', POS_NOPER(J,I,1:3)
     print*,'I = ', I
     print*,'K = ', K
     print*,'J = ', J     
     read(*,*)
   end if
   if (isnan(DISP_TEMP(2))) then
     print*,'DISP_TEMP(2) = ', DISP_TEMP(2)
     print*,'POS_NOPER(K,I,1:3) = ', POS_NOPER(K,I,1:3)
     print*,'POS_NOPER(J,I,1:3) = ', POS_NOPER(J,I,1:3)
     print*,'I = ', I
     print*,'K = ', K
     print*,'J = ', J     
     read(*,*)
   end if
   if (isnan(DISP_TEMP(3))) then
     print*,'DISP_TEMP(3) = ', DISP_TEMP(3)
     print*,'POS_NOPER(K,I,1:3) = ', POS_NOPER(K,I,1:3)
     print*,'POS_NOPER(J,I,1:3) = ', POS_NOPER(J,I,1:3)
     print*,'I = ', I
     print*,'K = ', K
     print*,'J = ', J     
     read(*,*)
   end if
			       
   MEAN_MSD_DIR(J-K,1) = MEAN_MSD_DIR(J-K,1) + DISP_TEMP(1)**2
   MEAN_MSD_DIR(J-K,2) = MEAN_MSD_DIR(J-K,2) + DISP_TEMP(2)**2
   MEAN_MSD_DIR(J-K,3) = MEAN_MSD_DIR(J-K,3) + DISP_TEMP(3)**2
   MEAN_MSD(J-K) = MEAN_MSD(J-K) &
                 + DISP_TEMP(1)**2 + DISP_TEMP(2)**2 + DISP_TEMP(3)**2   
  end do
  
  NDT(J-K) = NDT(J-K) + real(PART_END-PART_START+1)  
  
  
!   print*,'MEAN_MSD(J-K) = ', MEAN_MSD(J-K)
!   
!   print*,'NDT(J-K) = ', NDT(J-K)
!   read(*,*)
 end do
 
end do

MEAN_MSD = MEAN_MSD/NDT
MEAN_MSD_DIR(:,1) = MEAN_MSD_DIR(:,1)/NDT
MEAN_MSD_DIR(:,2) = MEAN_MSD_DIR(:,2)/NDT
MEAN_MSD_DIR(:,3) = MEAN_MSD_DIR(:,3)/NDT


print*, 'COMPUTE MSD---> OK'

!=====================================================================
! 2. COMPUTE PDF DISPLACEMENTS
!=====================================================================
print*, 'COMPUTE PDF DISP'

print*, 'DISCRETIZE RANGE OF DISPLACEMENT FOR PDF'
if (PART_START==1) then
 DDISPNORM = sqrt(3.0)*0.5*RAD
 DDISPX = 0.5*RAD
 DDISPY = DDISPX
 DDISPZ = DDISPX
!~ !  DDISPNORM = sqrt(3.0)*0.05*RAD
!~ ! DDISPX = 0.05*RAD
!~ ! DDISPY = DDISPX
!~ ! DDISPZ = DDISPX
else
!~  DDISPNORM = sqrt(3.0)*0.1*RAD
!~  DDISPX = 0.1*RAD
  DDISPNORM = sqrt(3.0)*0.2*RAD
 DDISPX = 0.2*RAD
 DDISPY = DDISPX
 DDISPZ = DDISPX
end if

if (PART_START==1) then
 MAX_DISPNORM = dsqrt(3.0d0)*100.0*RAD *1.5
 MAX_DISPX = 100.0*RAD *1.5
!~ !  MAX_DISPNORM = dsqrt(3.0)*10.0*RAD *1.5
!~ ! MAX_DISPX = 10.0*RAD *1.5
 MAX_DISPY = MAX_DISPX !200.0*maxval(dsqrt(MEAN_MSD))
 MAX_DISPZ = MAX_DISPX !200.0*maxval(dsqrt(MEAN_MSD))
else 
!~  MAX_DISPNORM = dsqrt(3.0)*20.0*RAD *1.5
!~  MAX_DISPX = 20.0*RAD *1.5
 MAX_DISPNORM = dsqrt(3.0d0)*100.0*RAD *1.5
 MAX_DISPX = 100.0*RAD *1.5
 MAX_DISPY = MAX_DISPX
 MAX_DISPZ = MAX_DISPX
end if

N_STEP_DISPNORM = ceiling(( MAX_DISPNORM &
                         - (-MAX_DISPNORM) )/DDISPNORM)
N_STEP_DISPX = ceiling(( MAX_DISPX &
                         - (-MAX_DISPX)  )/DDISPX)
N_STEP_DISPY = ceiling(( MAX_DISPY &
                         - (-MAX_DISPY) )/DDISPY)
N_STEP_DISPZ = ceiling(( MAX_DISPZ &
                         - (-MAX_DISPZ) )/DDISPZ)

print*,'N_STEP_DISPNORM ,N_STEP_DISPX ,N_STEP_DISPY ,N_STEP_DISPZ = '
print*, N_STEP_DISPNORM ,N_STEP_DISPX ,N_STEP_DISPY ,N_STEP_DISPZ


allocate(PDF_DISPNORM(NSAVES_MSD,N_STEP_DISPNORM))
allocate(PDF_DISPX(NSAVES_MSD,N_STEP_DISPX))
allocate(PDF_DISPY(NSAVES_MSD,N_STEP_DISPY))
allocate(PDF_DISPZ(NSAVES_MSD,N_STEP_DISPZ))
allocate(RANGE_DISPNORM(N_STEP_DISPNORM))
allocate(RANGE_DISPX(N_STEP_DISPX))
allocate(RANGE_DISPY(N_STEP_DISPY))
allocate(RANGE_DISPZ(N_STEP_DISPZ))

RANGE_DISPNORM(1) = -MAX_DISPNORM
do J = 2, N_STEP_DISPNORM
 RANGE_DISPNORM(J) = RANGE_DISPNORM(J-1) + DDISPNORM
end do

RANGE_DISPX(1) = -MAX_DISPX 
do J = 2, N_STEP_DISPX
 RANGE_DISPX(J) = RANGE_DISPX(J-1) + DDISPX
end do

RANGE_DISPY(1) = -MAX_DISPY
do J = 2, N_STEP_DISPY
 RANGE_DISPY(J) = RANGE_DISPY(J-1) + DDISPY
end do

RANGE_DISPZ(1) = -MAX_DISPZ
do J = 2, N_STEP_DISPZ
 RANGE_DISPZ(J) = RANGE_DISPZ(J-1) + DDISPZ
end do


print*, 'DISCRETIZE RANGE OF DISPLACEMENT FOR PDF---> OK'

PDF_DISPNORM = 0.0
PDF_DISPX = 0.0
PDF_DISPY = 0.0
PDF_DISPZ = 0.0

NMISS = 0

do K = 1, NSAVES-SAVE_START+1
 JNSAVES = min(K+NSAVES_MSD,NSAVES-SAVE_START+1)

 do J = K+1, JNSAVES
  do I = 1, PART_END-PART_START+1

   DISP_TEMP(1:3) = POS_NOPER(J,I,1:3)- POS_NOPER(K,I,1:3)
   DISPNORM_TEMP = dsqrt(DISP_TEMP(1)**2 + DISP_TEMP(2)**2 + DISP_TEMP(3)**2) 

   IND_NORM = floor( ( DISPNORM_TEMP-RANGE_DISPNORM(1) )/DDISPNORM ) + 1
   IND_X = floor( ( DISP_TEMP(1)-RANGE_DISPX(1) )/DDISPX ) + 1
   IND_Y = floor( ( DISP_TEMP(2)-RANGE_DISPY(1) )/DDISPY ) + 1
   IND_Z = floor( ( DISP_TEMP(3)-RANGE_DISPZ(1) )/DDISPZ ) + 1
   
   if ((IND_NORM<1).or.(IND_NORM>N_STEP_DISPNORM)) then
!~     print*, 'DISPNORM_TEMP = ', DISPNORM_TEMP
!~     print*, 'RANGE_DISPNORM(1), RANGE_DISPNORM(N_STEP_DISPNORM) = ', RANGE_DISPNORM(1), RANGE_DISPNORM(N_STEP_DISPNORM) 
!~     print*, 'N_STEP_DISPNORM = ', N_STEP_DISPNORM
!~     print*, 'IND_NORM = ', IND_NORM
    NMISS = NMISS + 1
  else 
   PDF_DISPNORM(J-K,IND_NORM) = PDF_DISPNORM(J-K,IND_NORM) + 1.0
  end if
   
   if ((IND_X<1).or.(IND_X>N_STEP_DISPX)) then
!~     print*, 'DISP_TEMP(1) = ', DISP_TEMP(1)
!~     print*, 'RANGE_DISPX(1), RANGE_DISPX(N_STEP_DISPX) = ', RANGE_DISPX(1), RANGE_DISPX(N_STEP_DISPX)
!~     print*, 'N_STEP_DISPX = ', N_STEP_DISPX
!~     print*, 'IND_X = ', IND_X
    NMISS = NMISS + 1
   else
    PDF_DISPX(J-K,IND_X) = PDF_DISPX(J-K,IND_X) + 1.0
   end if
  
   if ((IND_Y<1).or.(IND_Y>N_STEP_DISPX)) then
!~     print*, 'DISP_TEMP(2) = ', DISP_TEMP(2)
!~     print*, 'RANGE_DISPY(1), RANGE_DISPY(N_STEP_DISPX) = ', RANGE_DISPY(1), RANGE_DISPY(N_STEP_DISPY)
!~     print*, 'N_STEP_DISPY = ', N_STEP_DISPY
!~     print*, 'IND_Y = ', IND_Y
    NMISS = NMISS + 1
   else
    PDF_DISPY(J-K,IND_Y) = PDF_DISPY(J-K,IND_Y) + 1.0
   end if
  
   if ((IND_Z<1).or.(IND_Z>N_STEP_DISPZ)) then
!~     print*, 'DISP_TEMP(3) = ', DISP_TEMP(3)
!~     print*, 'RANGE_DISPZ(1), RANGE_DISPZ(N_STEP_DISPZ) = ', RANGE_DISPZ(1), RANGE_DISPZ(N_STEP_DISPZ)
!~     print*, 'N_STEP_DISPZ = ', N_STEP_DISPZ
!~     print*, 'IND_Z = ', IND_Z
    NMISS = NMISS + 1
   else
    PDF_DISPZ(J-K,IND_Z) = PDF_DISPZ(J-K,IND_Z) + 1.0
   end if

  end do
 end do
end do

print*,'NMISS = ', NMISS 
print*,'COMPUTE DISPLACEMENTS PDF --->  OK '
!=====================================================================
! 3. WRITE DIFF_T AND DIFF_R
!=====================================================================
print*,'SAVE MSD  '
!---------- MSD_DIR-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END

write(FILENAME,10103)'MSD_DIR_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1,  NSAVES_MSD
 write(301,'(4(e17.7))') real(I), (MEAN_MSD_DIR(I,J), J=1,3)
end do

!- close file
close(301)


!---------- MSD-----------
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
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


print*,'SAVE MSD --->  OK '
!=====================================================================
! 4. WRITE PDF
!=====================================================================

!---------- norm-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END


write(FILENAME,10103) 'PDF_DISPNORM_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII
write(301,'(4(e17.7))') real(NSAVES_MSD), real(N_STEP_DISPNORM), 0.0, 0.0

do K = 1, NSAVES_MSD
 do J = 1, N_STEP_DISPNORM
  write(301,'(4(e17.7))') real(K),NDT(K), RANGE_DISPNORM(J), PDF_DISPNORM(K,J)
 end do
end do

!- close file
close(301)

!---------- x-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_DISPX_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(4(e17.7))') real(NSAVES_MSD), real(N_STEP_DISPX), 0.0, 0.0

do K = 1, NSAVES_MSD
 do J = 1, N_STEP_DISPX
  write(302,'(4(e17.7))') real(K), NDT(K), RANGE_DISPX(J), PDF_DISPX(K,J)
 end do
end do
!- close file
close(302)


!---------- y-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_DISPY_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII
write(303,'(4(e17.7))') real(NSAVES_MSD), real(N_STEP_DISPY), 0.0, 0.0

do K = 1, NSAVES_MSD
 do J = 1, N_STEP_DISPY
  write(303,'(4(e17.7))') real(K), NDT(K), RANGE_DISPY(J), PDF_DISPY(K,J)
 end do
end do
!- close file
close(303)

!---------- z-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_DISPZ_TIME_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=304,file=trim(FILENAME))
!!- Ecriture ASCII
write(304,'(4(e17.7))') real(NSAVES_MSD), real(N_STEP_DISPZ), 0.0, 0.0

do K = 1, NSAVES_MSD
 do J = 1, N_STEP_DISPZ
  write(304,'(4(e17.7))') real(K), NDT(K), RANGE_DISPZ(J), PDF_DISPZ(K,J)
 end do
end do
!- close file
close(304)

print*,'SAVE PDF DISPLACEMENTS--->  OK '


deallocate(PDF_DISPNORM)
deallocate(PDF_DISPX)
deallocate(PDF_DISPY)
deallocate(PDF_DISPZ)
deallocate(RANGE_DISPNORM)
deallocate(RANGE_DISPX)
deallocate(RANGE_DISPY)
deallocate(RANGE_DISPZ)



print*,'SAVE MSD AND PDF DISP--->  OK '


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

end subroutine MSD_DISP
