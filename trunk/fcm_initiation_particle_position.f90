 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating particle positions with different ways
!!
!! Date :  15/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PARTICLE_POSITION

!!====================================================================
!! Here the particle positions are initiated according
!! to the variable INIT_PART_POSITION specified by user in the
!! parameter file: 'param.in'.
!!====================================================================
!! Particle position initiation: 
!!------------------------------
!! FCM_INIT_PART_POS=
!!
!!  0: The particle positions are specified in the Fortran file.
!!     The default position correspond to TO WRITE
!!
!!  1: Random particle distribution in the box. Take care that the
!!      minimal distance between particle must ensure that the 
!!     particles are not overlapping.
!!
!!  2: Particle positions are read from the binary file 'POSPART.ini'
!!
!!  3: Particle aligned according to routine 'place_spheres_aligned_x'.
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use MPI_STRUCTURES


implicit none

! Temporary array to read particle positions
real(kind=8), dimension(NPART_FULL,3) :: FCM_POS

! Default distance between particles
real(kind=8) :: DIST_INIT

!- File name 
character(len=40) :: FILENAME

! Process information variables
integer :: IERROR

!- Integers to signalize a non-existing file
integer :: ERR_FILE_POS

!- Indices 
integer :: IP

!- To delete
real(kind=8) ::  theta_0

! FCM_INIT_PART_POS  
!     =0 : positions are prescribed (NAPRT_FULL<=4)

if ((FCM_INIT_PART_POS ==0).or.(FCM_INIT_PART_POS ==4)) then
! For test_cases mostly

 if (NPART_FULL ==1) then
  if (FCM_BC==1) then
   FCM_XP(1) = LXMAX    /2.0
   FCM_YP(1) = LYMAX    /2.0
   FCM_ZP(1) = LZMAX    /2.0
  else if (FCM_BC ==2) then
   FCM_XP(1) = LXMAX    /4.0
   FCM_YP(1) = LYMAX    /2.0
   FCM_ZP(1) = LZMAX    /2.0
  end if
  
 else if (NPART_FULL ==2) then 
  FCM_XP(1) = LXMAX/2.0 - 1.1*FCM_SPHERE_RADP(1)
  FCM_YP(1) = LYMAX    /2.0
  FCM_ZP(1) = LZMAX    /2.0
  
  if (FCM_BC==2) then
   FCM_XP(2) = LXMAX - FCM_XP(1)
   FCM_YP(2) = LYMAX    /2.0
   FCM_ZP(2) = LZMAX    /2.0
   
  else if (FCM_BC==1) then
 
   if (FCM_NELLIPSOID==0) then

 !~    FCM_XP(2) = 3.0*LXMAX /4.0
 !~    FCM_YP(2) = LYMAX    /2.0
 !~    FCM_ZP(2) = LZMAX    /2.0 + 1.0*FCM_SPHERE_RADP(1)
 !~    
    theta_0 = 10.0*PPI/20.0
    FCM_XP(2) = LXMAX/2.0+1.1*FCM_SPHERE_RADP(2)
!~     FCM_XP(2) = FCM_XP(1) + 3.0*FCM_SPHERE_RADP(1)*(1.0 - dcos(PPI - theta_0)) !+ 2.2*FCM_SPHERE_RADP(1)
    FCM_YP(2) = FCM_YP(1) !+ 1.1*FCM_SPHERE_RADP(1)
    FCM_ZP(2) = FCM_ZP(1) !+ 3.0*FCM_SPHERE_RADP(1)*dsin(PPI - theta_0)
   else
    FCM_XP(2) = FCM_XP(1) &
 !~              + 0.1*maxval(FCM_ELLIPSOID_RADP) 
              + 1.1*( FCM_ELLIPSOID_RADP(1,1) + FCM_ELLIPSOID_RADP(2,1) )
    FCM_YP(2) = FCM_YP(1) !+ 1.1*FCM_SPHERE_RADP(1)
    FCM_ZP(2) = FCM_ZP(1) &
 !~              + 2.3*maxval(FCM_ELLIPSOID_RADP)
              + 1.0*FCM_ELLIPSOID_RADP(1,2) 
   end if
  end if !else if (FCM_BC==1) then
 
    if (MYID==0) then
   print*, 'FCM_XP/FCM_SPHERE_RADP(1)  = ', FCM_XP/FCM_SPHERE_RADP(1)
   print*, 'FCM_YP/FCM_SPHERE_RADP(1) = ', FCM_YP/FCM_SPHERE_RADP(1)
   print*, 'FCM_ZP/FCM_SPHERE_RADP(1) = ', FCM_ZP/FCM_SPHERE_RADP(1)
  end if

 else if (NPART_FULL ==4) then
 
!~   FCM_XP(1) = LXMAX    /2.0
!~   FCM_YP(1) = LYMAX    /4.0
!~   FCM_ZP(1) = LZMAX    /4.0
!~   
!~   FCM_XP(2) = LXMAX    /2.0
!~   FCM_YP(2) = LYMAX    /4.0 * 3.0
!~   FCM_ZP(2) = LZMAX    /4.0
!~   
!~   FCM_XP(3) = LXMAX    /2.0
!~   FCM_YP(3) = LYMAX    /4.0
!~   FCM_ZP(3) = LZMAX    /4.0 * 3.0
!~   
!~   FCM_XP(4) = LXMAX    /2.0
!~   FCM_YP(4) = LYMAX    /4.0 * 3.0
!~   FCM_ZP(4) = LZMAX    /4.0 * 3.0
  
  if (FCM_NELLIPSOID==0) then
   FCM_XP(1) = LXMAX    /2.0
   FCM_YP(1) = LYMAX    /2.0
   FCM_ZP(1) = LZMAX    /2.0
   
   FCM_XP(2) = FCM_XP(1) + 4.0*FCM_SPHERE_RADP(1) 
   FCM_YP(2) = FCM_YP(1) 
   FCM_ZP(2) = FCM_ZP(1) 
   
   FCM_XP(3) = FCM_XP(1) 
   FCM_YP(3) = FCM_YP(1) 
   FCM_ZP(3) = FCM_ZP(1) + 4.0*FCM_SPHERE_RADP(1) 
   
   FCM_XP(4) = FCM_XP(1) 
   FCM_YP(4) = FCM_YP(1) + 4.0*FCM_SPHERE_RADP(1)
   FCM_ZP(4) = FCM_ZP(1)  
   
  else
  
   DIST_INIT = 2.5*FCM_ELLIPSOID_RADP(1,1)
   FCM_XP(1) = LXMAX    /2.0
   FCM_YP(1) = LYMAX    /2.0
   FCM_ZP(1) = LZMAX    /2.0
   
   FCM_XP(2) = FCM_XP(1) + DIST_INIT/dsqrt(2.0D0)
   FCM_YP(2) = FCM_YP(1) 
   FCM_ZP(2) = FCM_ZP(1) + 0.8*DIST_INIT/dsqrt(2.0D0)
   
   FCM_XP(3) = FCM_XP(2) - DIST_INIT/dsqrt(2.0D0)
   FCM_YP(3) = FCM_YP(2) 
   FCM_ZP(3) = FCM_ZP(2) + 1.2*DIST_INIT/dsqrt(2.0D0)
   
   FCM_XP(4) = FCM_XP(3) - DIST_INIT/dsqrt(2.0D0)
   FCM_YP(4) = FCM_YP(3) 
   FCM_ZP(4) = FCM_ZP(3) - 1.3*DIST_INIT/dsqrt(2.0D0)
  
  end if
  if (MYID==0) then
   print*, 'FCM_XP/FCM_SPHERE_RADP(1)  = ', FCM_XP/FCM_SPHERE_RADP(1)
   print*, 'FCM_YP/FCM_SPHERE_RADP(1) = ', FCM_YP/FCM_SPHERE_RADP(1)
   print*, 'FCM_ZP/FCM_SPHERE_RADP(1) = ', FCM_ZP/FCM_SPHERE_RADP(1)
  end if
 end if
 
 
! FCM_INIT_PART_POS  
!     =2 : positions from MPIIO file
elseif (FCM_INIT_PART_POS ==2) then
  
  FILENAME='FCM_PART_POS.ini' 
  call FCM_READ_VAR_MPIIO(NPART_FULL,3,1,FCM_POS,FILENAME,ERR_FILE_POS)
 
  if (ERR_FILE_POS.eq.1) then
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!                     ERROR                    !!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!'
   write(*,*)'!! file     : FCM_PART_POS.ini'
   write(*,*)'!!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop
  end if
  
  do IP=1,NPART_FULL
   FCM_XP(IP) = FCM_POS(IP,1) - LXMAX*real(int(FCM_POS(IP,1)/LXMAX))
   FCM_YP(IP) = FCM_POS(IP,2) - LYMAX*real(int(FCM_POS(IP,2)/LYMAX))
   FCM_ZP(IP) = FCM_POS(IP,3) - LZMAX*real(int(FCM_POS(IP,3)/LZMAX))
  end do
  
  
 
 
 call MPI_BCAST(FCM_XP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_YP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_ZP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)


! FCM_INIT_PART_POS  
!     =3 : positions from ASCII file
elseif (FCM_INIT_PART_POS ==3) then
  
 if (MYID==0) then
  FILENAME='FCM_PART_POS.ini.ascii' 
  
  open(unit=200,file=FILENAME,status='old', iostat=ERR_FILE_POS)
 
  if (ERR_FILE_POS.eq.1) then
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!                     ERROR                    !!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!'
   write(*,*)'!! file     : FCM_PART_POS.ini.ascii'
   write(*,*)'!!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop
   
   
  else
   do IP = 1, NPART_FULL
    read(200,*) FCM_POS(IP,1:3)
   end do
  end if
  
  close(200)
  
  do IP=1,NPART_FULL
   FCM_XP(IP) = FCM_POS(IP,1) - LXMAX*real(int(FCM_POS(IP,1)/LXMAX))
   FCM_YP(IP) = FCM_POS(IP,2) - LYMAX*real(int(FCM_POS(IP,2)/LYMAX))
   FCM_ZP(IP) = FCM_POS(IP,3) - LZMAX*real(int(FCM_POS(IP,3)/LZMAX))
  end do
  
 end if
  
 call MPI_BCAST(FCM_XP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_YP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(FCM_ZP,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)

 
else

 print *, "CPU -- ", MYID, ":: ERROR: Wrong initialization mode!"
 call MPI_FINALIZE(IERROR)	




endif ! if FCM_INIT_PART_POS

! Initial non-periodic posistions = initial positions
FCM_XP_NOPER(1:NPART_FULL) = FCM_XP(1:NPART_FULL)
FCM_YP_NOPER(1:NPART_FULL) = FCM_YP(1:NPART_FULL)
FCM_ZP_NOPER(1:NPART_FULL) = FCM_ZP(1:NPART_FULL)
   
end subroutine FCM_INITIATION_PARTICLE_POSITION
