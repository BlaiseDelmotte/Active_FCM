!!====================================================================
!!
!!
!!====================================================================

Program MPIIO_TO_ASCII_POS_ONLY

!!====================================================================
!!
!!
!!====================================================================

use MPI

implicit none


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------

!- File name 
character(len=100) :: FILENAME
character(len=10) :: FILE_EXT, FILE_EXT2, FILE_EXT3



!!- Vraiables from direct access files
!-Mesh
real(kind=8), allocatable, dimension(:) :: XMESH
real(kind=8), allocatable, dimension(:) :: YMESH
real(kind=8), allocatable, dimension(:) :: ZMESH
!-Fluid field
real(kind=8), allocatable, dimension(:,:,:) :: UFLU
real(kind=8), allocatable, dimension(:,:,:) :: VFLU
real(kind=8), allocatable, dimension(:,:,:) :: WFLU

real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS
real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS_NOPER
real(kind=8), allocatable, dimension(:,:,:) :: FCM_VEL
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ROT
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ORIENT
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM
real(kind=8), allocatable, dimension(:,:,:) :: FCM_P2
real(kind=8), allocatable, dimension(:,:,:) :: FCM_P3
real(kind=8), allocatable, dimension(:,:,:) :: FCM_SIJ
real(kind=8), allocatable, dimension(:,:,:) :: FCM_RADII


real(kind=8) :: PPI
real(kind=8) :: RADMIN
!!- Dicrestization points in each direction
integer :: NX, NY, NZ
real(kind=8) :: LXMAX, LYMAX, LZMAX

!!- Particle number temporary
integer :: NPART, NSPHERE, NSPHERE_1, NELLIPSOID

!!- Particles to read
integer :: PART_START, PART_END
!!- Temporary dimension integers
integer :: DIM2, DIM3

!!- Number of time iterations
integer :: NCYCLEMAX

!!- Time-step
real(kind=8) :: DTIME

!!- Dump for fluid solution
integer :: FOUT1

!!- Dump for particle solution
integer :: FOUT3

!!- Flag if quaternions are used or orientation vectors
integer :: USE_QUAT

!!- Number of dumps to read
integer :: NB_DUMP_FILES_FLUID
integer :: NB_DUMP_FILES_PART

!- Time 
integer :: TIME

!- Save_start, Save_end, STEP_SAVE
integer :: SAVE_START, SAVE_END, STEP_SAVE

!- Index
integer :: I, J, K, IJK, IND

!- Integers to signalize a non-existing file
integer :: ERR_INFO, &
	   ERR_FILE_UFLU, &
           ERR_FILE_VFLU, &
           ERR_FILE_WFLU, &
           ERR_FILE_POS, &
	   ERR_FILE_POS_NOPER, &
           ERR_FILE_VEL, &
           ERR_FILE_ROT, &
           ERR_FILE_ORIENT, &
           ERR_FILE_SWIM, &
	   ERR_FILE_P2, &
           ERR_FILE_P3, &
           ERR_FILE_STRESS, &
           ERR_FILE_RADII

!- Integer to say if we convert fluid velocity into ASCII
integer ::  CONVERT_FLUID
!- if VORONOI = 1, only positions are saved  and in a way to be read by voro++
integer ::  VORONOI

!-  MPI Variables
integer :: IERR,NPROC,MYID



!!--------------------------------------------------------------------
!! 0.0 MPI World initiation
!!--------------------------------------------------------------------
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)


PPI = acos(-1.0)


!---------------------------------------------------------------------
!=====================================================================
! 0. READ PARAMETERS AND ALLOCATION
!=====================================================================


!- Open file
open(unit=200,file='fcm_run.info',status='old', iostat=ERR_INFO)

if(ERR_INFO /= 0) then

 print *, "CPU -- ", MYID, ":: ERROR: Input data file fcm_run.info open error!"
 call MPI_FINALIZE(ERR_INFO)
 call abort()
else

 read(200,*), LXMAX
 read(200,*), LYMAX
 read(200,*), LZMAX
 read(200,*), NX
 read(200,*), NY
 read(200,*), NZ
 read(200,*), DTIME
 read(200,*), NCYCLEMAX
 read(200,*), FOUT1
 read(200,*), FOUT3
 read(200,*), NPART
 read(200,*), NSPHERE, NSPHERE_1
 read(200,*), NELLIPSOID
 read(200,*)  ! VSW
 read(200,*)  ! RADMAX, but useless here
 read(200,*), USE_QUAT

end if
close(200)

!- Open file
open(unit=200,file='convert_ascii_choice.input',status='old', iostat=ERR_INFO)

if(ERR_INFO /= 0) then

 print *, "CPU -- ", MYID, ":: ERROR: Input data file convert_ascii_choice.input open error!"
 call MPI_FINALIZE(ERR_INFO)
 call abort()
else

 read(200,*), CONVERT_FLUID
 read(200,*), VORONOI
 read(200,*), SAVE_START
 read(200,*), SAVE_END
 read(200,*), STEP_SAVE
 read(200,*), PART_START
 read(200,*), PART_END
end if
close(200)

print*,' CONVERT_FLUID = ', CONVERT_FLUID
print*,' VORONOI = ', VORONOI
print*,' SAVE_START = ', SAVE_START
print*,' SAVE_END = ', SAVE_END
print*,' SAVE_STEP = ', STEP_SAVE
print*,' PART_START = ', PART_START
print*,' PART_END = ', PART_END
!~ LXMAX = 2.0*PPI
!~ LYMAX = 2.0*PPI
!~ LZMAX = 2.0*PPI

!~ ! Read dimensions from ending files
!~ FILENAME='uf.end'
!~ call READ_MPIIO_DIMS(NX,NY,NZ,FILENAME,ERR_FILE_UFLU)
!~ 
!~ ! Read dimensions from ending files
!~ FILENAME='FCM_PART_POS.end'
!~ call READ_MPIIO_DIMS(NPART,DIM2,DIM3,FILENAME,ERR_FILE_POS)

!~ NCYCLEMAX = 500
!~ FOUT1 = 1001
!~ FOUT3 = 2
if (CONVERT_FLUID==1) then
 NB_DUMP_FILES_FLUID = (NCYCLEMAX)/FOUT1
 print*,'NB_DUMP_FILES_FLUID = ' ,  NB_DUMP_FILES_FLUID
end if

NB_DUMP_FILES_PART= (NCYCLEMAX)/FOUT3
print*,'NB_DUMP_FILES_PART = ' ,  NB_DUMP_FILES_PART

if (CONVERT_FLUID==1) then
 allocate( XMESH(NX), YMESH(NY), ZMESH(NZ) )
 allocate( UFLU(NX,NY,NZ), VFLU(NX,NY,NZ), WFLU(NX,NY,NZ) )
end if
allocate( FCM_POS(NPART,3,1) )
allocate( FCM_RADII(NPART,3,1) )


if (CONVERT_FLUID==1) then
 !=====================================================================
 ! 1. READ ENDING FLUID VELOCITY
 !=====================================================================

 !!--------------------------------------------------------------------
 !! 1.1 Read Direct access file written with MPI/IO
 !!--------------------------------------------------------------------
 !!--------------------------------------------------------------------
 !!- x-component
 !!--------------------------------------------------------------------
 !!- Print file name
 FILENAME='uf.end'
 call READ_MPIIO(NX,NY,NZ,UFLU,FILENAME,ERR_FILE_UFLU)

 !!--------------------------------------------------------------------
 !! y-component
 !!--------------------------------------------------------------------
 !!- Print file name
 FILENAME='vf.end'
 call READ_MPIIO(NX,NY,NZ,VFLU,FILENAME,ERR_FILE_VFLU)

 !!--------------------------------------------------------------------
 !! z-component
 !!--------------------------------------------------------------------
 !!- Print file name
  FILENAME='wf.end'
 call READ_MPIIO(NX,NY,NZ,WFLU,FILENAME,ERR_FILE_WFLU)


 !=====================================================================
 ! 2. MESH GENERATION
 !=====================================================================
 do I = 1, NX
  XMESH(I) = (I-1)*LXMAX/NX
 end do

 do I = 1, NY
  YMESH(I) = (I-1)*LYMAX/NY
 end do

 do I = 1, NZ
  ZMESH(I) = (I-1)*LZMAX/NZ
 end do

 
 !=====================================================================
 ! 3. SAVE  FIELD + MESH + ENDING TIME IN ASCII 
 !=====================================================================
 TIME = NCYCLEMAX + 1

 !! - Write ASCII file only if the three files exist
 if ( (ERR_FILE_UFLU.eq.0).and.(ERR_FILE_VFLU.eq.0).and.(ERR_FILE_WFLU.eq.0) ) then

	 write(FILE_EXT,10205) TIME

	 !!-Print filename
	 write(FILENAME,10101)'uf_t',trim(FILE_EXT),'.dat'

	 !!- ASCII
	 open(unit=300,file=trim(FILENAME))
	 !!- Ecriture ASCII
	 write(300,1999)
	 write(300,2001)NX,NY,NZ,TIME
	 do K = 1, NZ
	  do J = 1, NY
	   do I = 1, NX
	    write(300,'(6(e17.7))')XMESH(I), YMESH(J), ZMESH(K), &
				   UFLU(I,J,K), &
				   VFLU(I,J,K), &
				   WFLU(I,J,K)
									    
	   end do
	  end do
	 end do
	 !- close file
	 close(300)

 end if

 !=====================================================================
 ! 4. READ INTERMEDIATE FLUID VELOCITY
 !=====================================================================

 do IND=1,NB_DUMP_FILES_FLUID

  TIME = (IND-1)*FOUT1 + 1

  write(FILE_EXT,10205) TIME
  
  write(FILENAME,10101)'uf_t',trim(FILE_EXT),'.bin'
  call READ_MPIIO(NX,NY,NZ,UFLU,FILENAME,ERR_FILE_UFLU) 

  write(FILENAME,10101)'vf_t',trim(FILE_EXT),'.bin'
  call READ_MPIIO(NX,NY,NZ,VFLU,FILENAME,ERR_FILE_VFLU) 

  write(FILENAME,10101)'wf_t',trim(FILE_EXT),'.bin'
  call READ_MPIIO(NX,NY,NZ,WFLU,FILENAME,ERR_FILE_WFLU)
  
  !=====================================================================
  ! 5. SAVE MESH +  INTERMEDIATE FLUID VELOCITY
  !=====================================================================

  !! - Write Tecplot file only if the three files exist
  if ( (ERR_FILE_UFLU.eq.0).and.(ERR_FILE_VFLU.eq.0).and.(ERR_FILE_WFLU.eq.0) ) then
  
	 write(FILE_EXT,10205) TIME

	 !!-Print filename
	 write(FILENAME,10101)'uf_t',trim(FILE_EXT),'.dat'

	 !!- ASCII
	 open(unit=300,file=trim(FILENAME))
	 !!- Ecriture ASCII
	 write(300,1999)
	 write(300,2001)NX,NY,NZ,TIME
	 do K = 1, NZ
	  do J = 1, NY
	   do I = 1, NX
	    write(300,'(6(e17.7))')XMESH(I), YMESH(J), ZMESH(K), &
				   UFLU(I,J,K), &
				   VFLU(I,J,K), &
				   WFLU(I,J,K)
									    
	   end do
	  end do
	 end do
	 !- close file
	 close(300)
  end if

 end do
  
end if
!=====================================================================
! 6. READ LAGRANGIAN DATA AT ENDING TIME
!=====================================================================
 
print*, '-------------------------------------------------------------'
print*, '                    START SAVING FILES                            '
print*, '-------------------------------------------------------------'

FILENAME='FCM_PART_RADII.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_RADII,FILENAME,ERR_FILE_RADII)

! Periodicity
if (ERR_FILE_RADII.eq.0) then

	!!-Print filename
	write(FILENAME,10200)'FCM_PART_RADII.dat'

	!!- ASCII
	open(unit=307,file=trim(FILENAME))
	!!- Ecriture ASCII

	do I = 1, NPART
	 write(307,'(3(e17.7))') (FCM_RADII(I,J,1), J=1,3)
	end do

	!- close file
	close(307)
end if

RADMIN = minval(FCM_RADII)

if (SAVE_END == NCYCLEMAX + 1) then
 FILENAME='FCM_PART_POS.end'
 call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS)


 !=====================================================================
 ! 7. SAVE LAGRANGIAN DATA AT ENDING TIME IN ASCII 
 !=====================================================================

 TIME = NCYCLEMAX + 1

 !!--------------------------------------------------------------------
 !! 7.1 Positions
 !!--------------------------------------------------------------------
 if (VORONOI == 1) then
  ! Periodicity
  if (ERR_FILE_POS.eq.0) then
	  !!-Print time
	  write(FILE_EXT,10205) TIME

	  !!-Print filename
	  write(FILENAME,10101)'FCM_PART_POS_VORO_t',trim(FILE_EXT),'.dat'

	  !!- ASCII
	  open(unit=301,file=trim(FILENAME))
	  !!- Ecriture ASCII
	  do I = PART_START, PART_END
	   write(301,2012) I, (FCM_POS(I,J,1)/RADMIN, J=1,3)
	  end do

	  !- close file
	  close(301)
  end if
 else
  ! Periodicity
  if (ERR_FILE_POS.eq.0) then
	  !!-Print time
	  write(FILE_EXT,10205) TIME

	  !!-Print filename
	  write(FILENAME,10101)'FCM_PART_POS_t',trim(FILE_EXT),'.dat'

	  !!- ASCII
	  open(unit=301,file=trim(FILENAME))
	  !!- Ecriture ASCII
	  write(301,2002)
	  write(301,2011)
	  write(301,'(I4.4)') PART_END - PART_START +1

	  do I =  PART_START, PART_END
	   write(301,'(3(e17.7))') (FCM_POS(I,J,1), J=1,3)
	  end do

	  !- close file
	  close(301)
  end if

  end if
end if

!=====================================================================
! 8. READ INTERMEDIATE LAGRANGIAN DATA
!=====================================================================

do IND=1,NB_DUMP_FILES_PART,STEP_SAVE


 TIME = (IND-1)*FOUT3 + 1
 
 if ( (TIME>=SAVE_START).and.(TIME<=SAVE_END) ) then
 
  write(FILE_EXT,10205) TIME
  
  write(FILENAME,10101)'FCM_PART_POS_t',trim(FILE_EXT),'.bin'
  call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS) 

  
  !=====================================================================
  ! 9. SAVE INTERMEDIATE LAG DATA IN ASCII
  !=====================================================================

  !!--------------------------------------------------------------------
  !! 9.1 Positions
  !!--------------------------------------------------------------------
  if (VORONOI == 1) then
  ! Periodicity
   if (ERR_FILE_POS.eq.0) then
	   !!-Print time
	   write(FILE_EXT,10205) TIME

	   !!-Print filename
	   write(FILENAME,10101)'FCM_PART_POS_VORO_t',trim(FILE_EXT),'.dat'

	   !!- ASCII
	   open(unit=301,file=trim(FILENAME))
	   !!- Ecriture ASCII
	   do I = PART_START, PART_END
	    write(301,2012) I, (FCM_POS(I,J,1)/RADMIN, J=1,3)
	   end do

	   !- close file
	   close(301)
   end if
  else
   ! Periodicity
   if (ERR_FILE_POS.eq.0) then
	   !!-Print time
	   write(FILE_EXT,10205) TIME

	   !!-Print filename
	   write(FILENAME,10101)'FCM_PART_POS_t',trim(FILE_EXT),'.dat'

	   !!- ASCII
	   open(unit=301,file=trim(FILENAME))
	   !!- Ecriture ASCII
	   write(301,2002)
	   write(301,2011)
	   write(301,'(I4.4)') PART_END - PART_START +1

	   do I = PART_START, PART_END
	    write(301,'(3(e17.7))') (FCM_POS(I,J,1), J=1,3)
	   end do

	   !- close file
	   close(301)
   end if

  end if ! IF VORONOI
 end if !SAVE
end do



print*, '-------------------------------------------------------------'
print*, '                    END SAVING FILES                            '
print*, '-------------------------------------------------------------'

 
!- Free MPI environement
call MPI_FINALIZE (ierr)

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

2012 format (I5.5,3(e17.7))

10200 format (A)
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I5)
10205 format (I8.8)
10101 format (A,A,A)

end program MPIIO_TO_ASCII_POS_ONLY
