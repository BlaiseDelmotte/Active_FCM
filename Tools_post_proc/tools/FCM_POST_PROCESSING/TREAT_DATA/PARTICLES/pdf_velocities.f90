!!====================================================================
!!
!!
!!====================================================================

subroutine PDF_VELOCITIES(NSAVES, &
                          SAVE_START, &
                          PART_START, &
                          PART_END, &
                          FCM_RADIUS, &
                          VEL)

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
! Particle radius
real(kind=8), intent(in) :: FCM_RADIUS
! Particle velocities
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: VEL

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- PDF of norm
real(kind=8), allocatable, dimension(:) :: PDF_VELNORM
real(kind=8), allocatable, dimension(:) :: RANGE_VELNORM
!- PDF of x-component
real(kind=8), allocatable, dimension(:) :: PDF_VELX
real(kind=8), allocatable, dimension(:) :: RANGE_VELX
!- PDF of y-component
real(kind=8), allocatable, dimension(:) :: PDF_VELY
real(kind=8), allocatable, dimension(:) :: RANGE_VELY
!- PDF of z-component
real(kind=8), allocatable, dimension(:) :: PDF_VELZ
real(kind=8), allocatable, dimension(:) :: RANGE_VELZ

! Particle velocities normalized by their radius
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3) :: VEL_RAD
!- Norm of velocities
real(kind=8),  dimension(NSAVES - SAVE_START + 1,PART_END-PART_START+1) :: VELNORM

!- Step size to discretize the interval of velocities
real(kind=8) :: DVELNORM, DVELX, DVELY, DVELZ

!- Number of steps to discretize the interval of velocities
integer :: N_STEP_VELNORM, N_STEP_VELX, N_STEP_VELY, N_STEP_VELZ

!- Indices for pdf
integer :: IND_NORM, IND_X, IND_Y, IND_Z

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT1, FILE_EXT2, FILE_EXT3, FILE_EXT4


!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_PDF

print*,' '
print*,'-------------------START PDF VELOCITIES------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE NORM OF VELOCITIES
!=====================================================================

! Normalize velocity by particle radius
VEL_RAD = VEL / FCM_RADIUS

print*, 'COMPUTE NORM OF VELOCITIES'
K = 0
do IND = SAVE_START, NSAVES
 K = K+1
 do I = 1, PART_END-PART_START+1
 
  VELNORM(K,I) = dsqrt( VEL_RAD(IND,I,1)**2 &
                        + VEL_RAD(IND,I,2)**2 &
                        + VEL_RAD(IND,I,3)**2 )
 
 end do
end do
 
print*, 'COMPUTE NORM OF VELOCITIES---> OK'
!=====================================================================
! 2. DEFINE DISCRETIZATION OF PDF
!=====================================================================
print*, 'DISCRETIZE RANGE OF VELOCITIES FOR PDF'
DVELNORM = 0.05
DVELX = 0.05
DVELY = 0.05
DVELZ = 0.05

N_STEP_VELNORM = ceiling(( maxval(VELNORM) &
                         - minval(VELNORM) )/DVELNORM)
N_STEP_VELX = ceiling(( maxval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,1)) &
                      - minval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,1)) )/DVELX)
N_STEP_VELY = ceiling(( maxval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,2)) &
                      - minval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,2)) )/DVELY)
N_STEP_VELZ = ceiling(( maxval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,3)) &
                      - minval(VEL_RAD(SAVE_START:NSAVES,1:PART_END-PART_START+1,3)) )/DVELZ)

print*,'N_STEP_VELNORM ,N_STEP_VELX ,N_STEP_VELY ,N_STEP_VELZ = '
print*, N_STEP_VELNORM ,N_STEP_VELX ,N_STEP_VELY ,N_STEP_VELZ
        

allocate(PDF_VELNORM(N_STEP_VELNORM))
allocate(PDF_VELX(N_STEP_VELX))
allocate(PDF_VELY(N_STEP_VELY))
allocate(PDF_VELZ(N_STEP_VELZ))
allocate(RANGE_VELNORM(N_STEP_VELNORM))
allocate(RANGE_VELX(N_STEP_VELX))
allocate(RANGE_VELY(N_STEP_VELY))
allocate(RANGE_VELZ(N_STEP_VELZ))

RANGE_VELNORM(1) = minval(VELNORM)
do J = 2, N_STEP_VELNORM
 RANGE_VELNORM(J) = RANGE_VELNORM(J-1) + DVELNORM
end do

RANGE_VELX(1) = minval(VEL_RAD(SAVE_START:NSAVES,:,1))
do J = 2, N_STEP_VELX
 RANGE_VELX(J) = RANGE_VELX(J-1) + DVELX
end do

RANGE_VELY(1) = minval(VEL_RAD(SAVE_START:NSAVES,:,2))
do J = 2, N_STEP_VELY
 RANGE_VELY(J) = RANGE_VELY(J-1) + DVELY
end do

RANGE_VELZ(1) = minval(VEL_RAD(SAVE_START:NSAVES,:,3))
do J = 2, N_STEP_VELZ
 RANGE_VELZ(J) = RANGE_VELZ(J-1) + DVELZ
end do


print*, 'DISCRETIZE RANGE OF VELOCITIES FOR PDF---> OK'


!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE VELOCITIES PDF  '

PDF_VELNORM = 0.0
PDF_VELX = 0.0
PDF_VELY = 0.0
PDF_VELZ = 0.0

K = 0

do IND = SAVE_START, NSAVES 
 K = K+1
 do I = 1, PART_END-PART_START+1
  IND_NORM = floor( ( VELNORM(K,I)-RANGE_VELNORM(1) )/DVELNORM ) + 1
  IND_X = floor( ( VEL_RAD(IND,I,1)-RANGE_VELX(1) )/DVELX ) + 1
  IND_Y = floor( ( VEL_RAD(IND,I,2)-RANGE_VELY(1) )/DVELY ) + 1
  IND_Z = floor( ( VEL_RAD(IND,I,3)-RANGE_VELZ(1) )/DVELZ ) + 1
  
  if (isnan(VELNORM(K,I))) then
     print*,'VELNORM(K,I) = ', VELNORM(K,I)
     print*,'I = ', I
     print*,'K = ', K
     read(*,*)
   end if
  
  
  if ((IND_NORM<1).or.(IND_NORM>N_STEP_VELNORM)) then
   print*, 'VELNORM(IND,I) = ', VELNORM(K,I)
   print*, 'RANGE_VELNORM(1) = ', RANGE_VELNORM(1)
   print*, '( VELNORM(IND,I)-RANGE_VELNORM(1) )/DVELNORM = ', ( VELNORM(K,I)-RANGE_VELNORM(1) )/DVELNORM
   print*, 'IND_NORM = ', IND_NORM
   read(*,*)
  end if
  
  if ((IND_X<1).or.(IND_X>N_STEP_VELX)) then
   print*, 'VEL(IND,I,1) = ', VEL_RAD(IND,I,1)
   print*, 'RANGE_VELX(1) = ', RANGE_VELX(1)
   print*, '( VEL(IND,I,1)-RANGE_VELX(1) )/DVELX = ', ( VEL_RAD(IND,I,1)-RANGE_VELX(1) )/DVELX
   print*, 'IND_X = ', IND_X
   read(*,*)
  end if
  
  if ((IND_Y<1).or.(IND_Y>N_STEP_VELY)) then
   print*, 'VEL(IND,I,2) = ', VEL_RAD(IND,I,2)
   print*, 'RANGE_VELY(1) = ', RANGE_VELY(1)
   print*, '( VEL(IND,I,2)-RANGE_VELY(1) )/DVELY = ', ( VEL_RAD(IND,I,2)-RANGE_VELY(1) )/DVELY
   print*, 'IND_Y = ', IND_Y
   read(*,*)
  end if
  
  if ((IND_Z<1).or.(IND_Z>N_STEP_VELZ)) then
   print*, 'VEL(IND,I,3) = ', VEL_RAD(IND,I,3)
   print*, 'RANGE_VELZ(1) = ', RANGE_VELZ(1)
   print*, '( VEL(IND,I,3)-RANGE_VELZ(1) )/DVELZ = ', ( VEL_RAD(IND,I,3)-RANGE_VELZ(1) )/DVELZ
   print*, 'IND_Z = ', IND_Z
   read(*,*)
  end if
  
  PDF_VELNORM(IND_NORM) = PDF_VELNORM(IND_NORM) + 1.0
  PDF_VELX(IND_X) = PDF_VELX(IND_X) + 1.0
  PDF_VELY(IND_Y) = PDF_VELY(IND_Y) + 1.0
  PDF_VELZ(IND_Z) = PDF_VELZ(IND_Z) + 1.0
 end do
end do

print*,'COMPUTE VELOCITIES PDF --->  OK '

!=====================================================================
! 3. WRITE PDF
!=====================================================================

!---------- norm-----------
!!-Print filename
write(FILE_EXT1,10205) SAVE_START
write(FILE_EXT2,10205) NSAVES
write(FILE_EXT3,10205) PART_START
write(FILE_EXT4,10205) PART_END


write(FILENAME,10103) 'PDF_VELNORM_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                      '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_VELNORM
 write(301,'(2(e17.7))') RANGE_VELNORM(J), PDF_VELNORM(J)
end do

!- close file
close(301)

!---------- x-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_VELX_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_VELX
 write(302,'(2(e17.7))') RANGE_VELX(J), PDF_VELX(J)
end do

!- close file
close(302)


!---------- y-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_VELY_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_VELY
 write(303,'(2(e17.7))') RANGE_VELY(J), PDF_VELY(J)
end do

!- close file
close(303)

!---------- z-component-----------
!!-Print filename
write(FILENAME,10103)'PDF_VELZ_PART_',trim(FILE_EXT3),'_',trim(FILE_EXT4),&
                     '_SAVE_',trim(FILE_EXT1),'_',trim(FILE_EXT2),'.dat'

!!- ASCII
open(unit=304,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_VELZ
 write(304,'(2(e17.7))') RANGE_VELZ(J), PDF_VELZ(J)
end do

!- close file
close(304)

print*,'SAVE PDF VELOCITIES--->  OK '


deallocate(PDF_VELNORM)
deallocate(PDF_VELX)
deallocate(PDF_VELY)
deallocate(PDF_VELZ)
deallocate(RANGE_VELNORM)
deallocate(RANGE_VELX)
deallocate(RANGE_VELY)
deallocate(RANGE_VELZ)

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

end subroutine PDF_VELOCITIES
