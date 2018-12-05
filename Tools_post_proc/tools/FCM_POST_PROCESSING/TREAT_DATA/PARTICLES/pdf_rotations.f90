!!====================================================================
!!
!!
!!====================================================================

subroutine PDF_ROTATIONS(NSAVES, SAVE_START, NPART, ROT)

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
! Particle ROTocities
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: ROT

!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- PDF of norm
real(kind=8), allocatable, dimension(:) :: PDF_ROTNORM
real(kind=8), allocatable, dimension(:) :: RANGE_ROTNORM
!- PDF of x-component
real(kind=8), allocatable, dimension(:) :: PDF_ROTX
real(kind=8), allocatable, dimension(:) :: RANGE_ROTX
!- PDF of y-component
real(kind=8), allocatable, dimension(:) :: PDF_ROTY
real(kind=8), allocatable, dimension(:) :: RANGE_ROTY
!- PDF of z-component
real(kind=8), allocatable, dimension(:) :: PDF_ROTZ
real(kind=8), allocatable, dimension(:) :: RANGE_ROTZ

!- Norm of ROTocities
real(kind=8),  dimension(NSAVES - SAVE_START + 1,NPART) :: ROTNORM

!- Step size to discretize the interval of ROTocities
real(kind=8) :: DROTNORM, DROTX, DROTY, DROTZ

!- Number of steps to discretize the interval of ROTocities
integer :: N_STEP_ROTNORM, N_STEP_ROTX, N_STEP_ROTY, N_STEP_ROTZ

!- Indices for pdf
integer :: IND_NORM, IND_X, IND_Y, IND_Z

!- File name 
character(len=40) :: FILENAME


!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_PDF

print*,' '
print*,'-------------------START PDF ROTATIONS------------------------ '
print*,' '
!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE NORM OF ROTOCITIES
!=====================================================================
print*, 'COMPUTE NORM OF ROTATIONS'
K = 0
do IND = SAVE_START, NSAVES
 K = K+1
 do I = 1, NPART
 
  ROTNORM(K,I) = dsqrt( ROT(IND,I,1)**2 &
                        + ROT(IND,I,2)**2 &
                        + ROT(IND,I,3)**2 )
 
 end do
end do
 
print*, 'COMPUTE NORM OF ROTATIONS---> OK'
!=====================================================================
! 2. DEFINE DISCRETIZATION OF PDF
!=====================================================================
print*, 'DISCRETIZE RANGE OF ROTATIONS FOR PDF'
DROTNORM = 0.01
DROTX = 0.01
DROTY = 0.01
DROTZ = 0.01

N_STEP_ROTNORM = ceiling(( maxval(ROTNORM) &
                         - minval(ROTNORM) )/DROTNORM)
N_STEP_ROTX = ceiling(( maxval(ROT(SAVE_START:NSAVES,:,1)) &
                      - minval(ROT(SAVE_START:NSAVES,:,1)) )/DROTX)
N_STEP_ROTY = ceiling(( maxval(ROT(SAVE_START:NSAVES,:,2)) &
                      - minval(ROT(SAVE_START:NSAVES,:,2)) )/DROTY)
N_STEP_ROTZ = ceiling(( maxval(ROT(SAVE_START:NSAVES,:,3)) &
                      - minval(ROT(SAVE_START:NSAVES,:,3)) )/DROTZ)

print*,'N_STEP_ROTNORM ,N_STEP_ROTX ,N_STEP_ROTY ,N_STEP_ROTZ = '
print*, N_STEP_ROTNORM ,N_STEP_ROTX ,N_STEP_ROTY ,N_STEP_ROTZ
        

allocate(PDF_ROTNORM(N_STEP_ROTNORM))
allocate(PDF_ROTX(N_STEP_ROTX))
allocate(PDF_ROTY(N_STEP_ROTY))
allocate(PDF_ROTZ(N_STEP_ROTZ))
allocate(RANGE_ROTNORM(N_STEP_ROTNORM))
allocate(RANGE_ROTX(N_STEP_ROTX))
allocate(RANGE_ROTY(N_STEP_ROTY))
allocate(RANGE_ROTZ(N_STEP_ROTZ))

RANGE_ROTNORM(1) = minval(ROTNORM)
do J = 2, N_STEP_ROTNORM
 RANGE_ROTNORM(J) = RANGE_ROTNORM(J-1) + DROTNORM
end do

RANGE_ROTX(1) = minval(ROT(SAVE_START:NSAVES,:,1)) 
do J = 2, N_STEP_ROTX
 RANGE_ROTX(J) = RANGE_ROTX(J-1) + DROTX
end do

RANGE_ROTY(1) = minval(ROT(SAVE_START:NSAVES,:,2)) 
do J = 2, N_STEP_ROTY
 RANGE_ROTY(J) = RANGE_ROTY(J-1) + DROTY
end do

RANGE_ROTZ(1) = minval(ROT(SAVE_START:NSAVES,:,3)) 
do J = 2, N_STEP_ROTZ
 RANGE_ROTZ(J) = RANGE_ROTZ(J-1) + DROTZ
end do


print*, 'DISCRETIZE RANGE OF ROTATIONS FOR PDF---> OK'


!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE ROTATIONS PDF  '

PDF_ROTNORM = 0.0
PDF_ROTX = 0.0
PDF_ROTY = 0.0
PDF_ROTZ = 0.0

K = 0

do IND = SAVE_START, NSAVES 
 K = K+1
 do I = 1, NPART
  IND_NORM = floor( ( ROTNORM(K,I)-RANGE_ROTNORM(1) )/DROTNORM ) + 1
  IND_X = floor( ( ROT(IND,I,1)-RANGE_ROTX(1) )/DROTX ) + 1
  IND_Y = floor( ( ROT(IND,I,2)-RANGE_ROTY(1) )/DROTY ) + 1
  IND_Z = floor( ( ROT(IND,I,3)-RANGE_ROTZ(1) )/DROTZ ) + 1
  
  if ((IND_NORM<1).or.(IND_NORM>N_STEP_ROTNORM)) then
   print*, 'ROTNORM(IND,I) = ', ROTNORM(K,I)
   print*, 'RANGE_ROTNORM(1) = ', RANGE_ROTNORM(1)
   print*, '( ROTNORM(IND,I)-RANGE_ROTNORM(1) )/DROTNORM = ', ( ROTNORM(K,I)-RANGE_ROTNORM(1) )/DROTNORM
   print*, 'IND_NORM = ', IND_NORM
   read(*,*)
  end if
  
  if ((IND_X<1).or.(IND_X>N_STEP_ROTX)) then
   print*, 'ROT(IND,I,1) = ', ROT(IND,I,1)
   print*, 'RANGE_ROTX(1) = ', RANGE_ROTX(1)
   print*, '( ROT(IND,I,1)-RANGE_ROTX(1) )/DROTX = ', ( ROT(IND,I,1)-RANGE_ROTX(1) )/DROTX
   print*, 'IND_X = ', IND_X
   read(*,*)
  end if
  
  if ((IND_Y<1).or.(IND_Y>N_STEP_ROTY)) then
   print*, 'ROT(IND,I,2) = ', ROT(IND,I,2)
   print*, 'RANGE_ROTY(1) = ', RANGE_ROTY(1)
   print*, '( ROT(IND,I,2)-RANGE_ROTY(1) )/DROTY = ', ( ROT(IND,I,2)-RANGE_ROTY(1) )/DROTY
   print*, 'IND_Y = ', IND_Y
   read(*,*)
  end if
  
  if ((IND_Z<1).or.(IND_Z>N_STEP_ROTZ)) then
   print*, 'ROT(IND,I,3) = ', ROT(IND,I,3)
   print*, 'RANGE_ROTZ(1) = ', RANGE_ROTZ(1)
   print*, '( ROT(IND,I,3)-RANGE_ROTZ(1) )/DROTZ = ', ( ROT(IND,I,3)-RANGE_ROTZ(1) )/DROTZ
   print*, 'IND_Z = ', IND_Z
   read(*,*)
  end if
  
  PDF_ROTNORM(IND_NORM) = PDF_ROTNORM(IND_NORM) + 1.0
  PDF_ROTX(IND_X) = PDF_ROTX(IND_X) + 1.0
  PDF_ROTY(IND_Y) = PDF_ROTY(IND_Y) + 1.0
  PDF_ROTZ(IND_Z) = PDF_ROTZ(IND_Z) + 1.0
 end do
end do

print*,'COMPUTE ROTATIONS PDF --->  OK '

!=====================================================================
! 3. WRITE PDF
!=====================================================================

!---------- norm-----------
!!-Print filename
write(FILENAME,10200)'PDF_ROTNORM.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_ROTNORM
 write(301,'(2(e17.7))') RANGE_ROTNORM(J), PDF_ROTNORM(J)
end do

!- close file
close(301)

!---------- x-component-----------
!!-Print filename
write(FILENAME,10200)'PDF_ROTX.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_ROTX
 write(302,'(2(e17.7))') RANGE_ROTX(J), PDF_ROTX(J)
end do

!- close file
close(302)


!---------- y-component-----------
!!-Print filename
write(FILENAME,10200)'PDF_ROTY.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_ROTY
 write(303,'(2(e17.7))') RANGE_ROTY(J), PDF_ROTY(J)
end do

!- close file
close(303)

!---------- z-component-----------
!!-Print filename
write(FILENAME,10200)'PDF_ROTZ.dat'

!!- ASCII
open(unit=304,file=trim(FILENAME))
!!- Ecriture ASCII

do J = 1, N_STEP_ROTZ
 write(304,'(2(e17.7))') RANGE_ROTZ(J), PDF_ROTZ(J)
end do

!- close file
close(304)

print*,'SAVE PDF ROTOCITIES--->  OK '


deallocate(PDF_ROTNORM)
deallocate(PDF_ROTX)
deallocate(PDF_ROTY)
deallocate(PDF_ROTZ)
deallocate(RANGE_ROTNORM)
deallocate(RANGE_ROTX)
deallocate(RANGE_ROTY)
deallocate(RANGE_ROTZ)

print*,' '
print*,'-------------------END PDF ROTATIONS-------------------------- '
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

end subroutine PDF_ROTATIONS
