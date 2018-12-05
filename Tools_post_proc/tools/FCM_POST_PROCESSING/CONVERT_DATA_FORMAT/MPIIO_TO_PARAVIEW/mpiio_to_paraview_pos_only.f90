!!====================================================================
!!
!!
!!====================================================================

Program MPIIO_TO_PARAVIEW_POS_ONLY

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
character(len=40) :: FILENAME



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
real(kind=8), allocatable, dimension(:,:,:) :: FCM_VEL
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ROT
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ORIENT
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM
real(kind=8), allocatable, dimension(:,:,:) :: FCM_P2
real(kind=8), allocatable, dimension(:,:,:) :: FCM_P3
real(kind=8), allocatable, dimension(:,:,:) :: FCM_SIJ
real(kind=8), allocatable, dimension(:,:,:) :: FCM_RADII

real(kind=8), allocatable, dimension(:,:) :: ELL_MAT_VEC


real(kind=8), allocatable, dimension(:,:) :: FCM_PMEAN
real(kind=8), allocatable, dimension(:,:) :: P_SCAL_PMEAN

!!-  Vswim
real(kind=8) :: FCM_VSW

real(kind=8) :: PPI

!!- Dicrestization points in each direction
integer :: NX, NY, NZ
real(kind=8) :: LXMAX, LYMAX, LZMAX, RADMAX

!!- Particle number temporary
integer :: NPART, NSPHERE, NELLIPSOID

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

!- Time 
integer :: TIME

!!- Number of dumps to read
integer :: NB_DUMP_FILES_FLUID
integer :: NB_DUMP_FILES_PART

!!- Flag if quaternions are used or orientation vectors
integer :: USE_QUAT
!- Save_start, Save_end, STEP_SAVE
integer :: SAVE_START, SAVE_END, STEP_SAVE
! Flag to decide to use mean_pswim or not
integer :: USE_MEAN_PSWIM
!- Index
integer :: I, J, K, IND

!-  MPI Variables
integer :: IERR,NPROC,MYID

!- Integers to signalize a non-existing file
integer :: ERR_INFO, &
           ERR_FILE_UFLU, &
           ERR_FILE_VFLU, &
           ERR_FILE_WFLU, &
           ERR_FILE_POS, &
           ERR_FILE_VEL, &
           ERR_FILE_ROT, &
           ERR_FILE_ORIENT, &
           ERR_FILE_P2, &
           ERR_FILE_P3, &
           ERR_FILE_SWIM, &
           ERR_FILE_STRESS, &
           ERR_FILE_RADII
!- To decide whether particles are written in TCPLT or ASCII
logical :: BINFLAG

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
 read(200,*), NSPHERE
 read(200,*), NELLIPSOID
 read(200,*), FCM_VSW
 read(200,*), RADMAX
 read(200,*), USE_QUAT
end if

!- Open file
open(unit=200,file='convert_paraview_choice.input',status='old', iostat=ERR_INFO)

if(ERR_INFO /= 0) then

 print *, "CPU -- ", MYID, ":: ERROR: Input data file convert_ascii_choice.input open error!"
 call MPI_FINALIZE(ERR_INFO)
 call abort()
else
 read(200,*), SAVE_START
 read(200,*), SAVE_END
 read(200,*), STEP_SAVE
 read(200,*), USE_MEAN_PSWIM
end if
close(200)

NB_DUMP_FILES_PART= (NCYCLEMAX)/FOUT3
print*,'NB_DUMP_FILES_PART = ' ,  NB_DUMP_FILES_PART

!~ allocate( XMESH(NX), YMESH(NY), ZMESH(NZ) )
!~ allocate( UFLU(NX,NY,NZ), VFLU(NX,NY,NZ), WFLU(NX,NY,NZ) )
allocate( FCM_POS(NPART,3,1) )




!=====================================================================
! 4. READ LAGRANGIAN DATA 
!=====================================================================


BINFLAG = .false.
 
 
FILENAME='FCM_PART_POS.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS)





TIME = NCYCLEMAX + 1





 if (ERR_FILE_POS.eq.0) then

  
   call PRINT_PARAVIEW_PART_POS_ONLY(TIME, &
                             NPART, &
                             FCM_POS(:,:,1), &
                             RADMAX )
 end if
 

!=====================================================================
! 5. READ INTERMEDIATE LAGRANGIAN DATA
!=====================================================================

do IND=1,NB_DUMP_FILES_PART,STEP_SAVE

 TIME = (IND-1)*FOUT3 + 1
 
 if ( (TIME>=SAVE_START).and.(TIME<=SAVE_END) ) then
 
 write(FILENAME,10205)'FCM_PART_POS_t',TIME,'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS) 

  
 !=====================================================================
 ! 6. SAVE INTERMEDIATE LAG DATA IN PARAVIEW
 !=====================================================================
  if (ERR_FILE_POS.eq.0) then

  
   call PRINT_PARAVIEW_PART_POS_ONLY(TIME, &
                             NPART, &
                             FCM_POS(:,:,1), &
                             RADMAX )
  end if
 end if !SAVE
end do




 
 !- Free MPI environement
 call MPI_FINALIZE (ierr)

!!====================================================================
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10205 format (A,I8.8,A)
10101 format (A,A,A)
10001 format (I8.8)


end program MPIIO_TO_PARAVIEW_POS_ONLY
