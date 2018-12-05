!!====================================================================
!!
!!
!!====================================================================

Program MPIIO_TO_TECPLOT

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
real(kind=8), allocatable, dimension(:,:,:) :: FCM_SIJ
real(kind=8), allocatable, dimension(:,:,:) :: FCM_RADII


real(kind=8) :: PPI

!!- Dicrestization points in each direction
integer :: NX, NY, NZ
real(kind=8) :: LXMAX, LYMAX, LZMAX

!!- Particle number temporary
integer :: NPART, NSPHERE, NELLIPSOID

!!- Maximal radius
real(kind=8) :: RADMAX

!!- Vsw
real(kind=8) :: VSW

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

 print *, "CPU -- ", MYID, ":: ERROR: Input data file fcm_param.in open error!"
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
 read(200,*), VSW
 read(200,*), RADMAX
 read(200,*), !USE_QUAT, useless here

end if

NB_DUMP_FILES_FLUID = (NCYCLEMAX)/FOUT1

print*,'NB_DUMP_FILES_FLUID = ' ,  NB_DUMP_FILES_FLUID

NB_DUMP_FILES_PART= (NCYCLEMAX)/FOUT3
print*,'NB_DUMP_FILES_PART = ' ,  NB_DUMP_FILES_PART

allocate( XMESH(NX), YMESH(NY), ZMESH(NZ) )
allocate( UFLU(NX,NY,NZ), VFLU(NX,NY,NZ), WFLU(NX,NY,NZ) )
allocate( FCM_POS(NPART,3,1) )
allocate( FCM_VEL(NPART,3,1) )
!~ allocate( FCM_ROT(NPART,3,1) )
!~ allocate( FCM_ORIENT(NPART,4,1) )
allocate( FCM_PSWIM(NPART,3,1) )
!~ allocate( FCM_SIJ(NPART,5,1) )
allocate( FCM_RADII(NPART,3,1) )




!=====================================================================
! 1. MESH GENERATION
!=====================================================================
do I = 1, NX
 XMESH(I) = (I-1)*LXMAX/NX/RADMAX
end do

do I = 1, NY
 YMESH(I) = (I-1)*LYMAX/NY/RADMAX
end do

do I = 1, NZ
 ZMESH(I) = (I-1)*LZMAX/NZ/RADMAX
end do

!=====================================================================
! 2. FLUID VELOCITY
!=====================================================================

!!--------------------------------------------------------------------
!! 2.1 Read Direct access file written with MPI/IO at ending time
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

!!--------------------------------------------------------------------
! 2.2 SAVE MESH + ENDING FIELDS IN TECPLOT FORMAT
!!--------------------------------------------------------------------
TIME = NCYCLEMAX + 1

!! - Write Tecplot file only if the three files exist
if ( (ERR_FILE_UFLU.eq.0).and.(ERR_FILE_VFLU.eq.0).and.(ERR_FILE_WFLU.eq.0) ) then
 call PRINT_TECPLOT(TIME,NX,NY,NZ,XMESH,YMESH,ZMESH,UFLU/VSW,VFLU/VSW,WFLU/VSW)
end if


!!--------------------------------------------------------------------
!! 2.3 Read Direct access file written with MPI/IO at intermediate times
!!--------------------------------------------------------------------

do I=1,NB_DUMP_FILES_FLUID

 TIME = (I-1)*FOUT1 + 1



 write(FILENAME,10205)'uf_t',TIME,'.bin'
 call READ_MPIIO(NX,NY,NZ,UFLU,FILENAME,ERR_FILE_UFLU)
 
 write(FILENAME,10205)'vf_t',TIME,'.bin'
 call READ_MPIIO(NX,NY,NZ,VFLU,FILENAME,ERR_FILE_VFLU)
 
 write(FILENAME,10205)'wf_t',TIME,'.bin'
 call READ_MPIIO(NX,NY,NZ,WFLU,FILENAME,ERR_FILE_WFLU)
 
 !!--------------------------------------------------------------------
 ! 2.4 SAVE MESH + INTERMEDIATE FIELDS IN TECPLOT FORMAT
 !!--------------------------------------------------------------------


 !! - Write Tecplot file only if the three files exist
 if ( (ERR_FILE_UFLU.eq.0).and.(ERR_FILE_VFLU.eq.0).and.(ERR_FILE_WFLU.eq.0) ) then
  call PRINT_TECPLOT(TIME,NX,NY,NZ,XMESH,YMESH,ZMESH,UFLU/VSW,VFLU/VSW,WFLU/VSW)
 end if

end do




!=====================================================================
! 4. READ LAGRANGIAN DATA 
!=====================================================================

BINFLAG = .false.
 
 
FILENAME='FCM_PART_POS.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS)

FILENAME='FCM_PART_VEL.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_VEL,FILENAME,ERR_FILE_VEL)
!~ 
!~ FILENAME='FCM_PART_ROT.end'
!~ call READ_VAR_MPIIO(NPART,3,1,FCM_ROT,FILENAME,ERR_FILE_ROT)

!~ FILENAME='FCM_PART_ORIENT.end'
!~ call READ_VAR_MPIIO(NPART,4,1,FCM_ORIENT,FILENAME,ERR_FILE_ORIENT)

FILENAME='FCM_PART_SWIM.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_PSWIM,FILENAME,ERR_FILE_SWIM)

!~ FILENAME='FCM_PART_STRESSLET.end'
!~ call READ_VAR_MPIIO(NPART,5,1,FCM_SIJ,FILENAME,ERR_FILE_STRESS)

FILENAME='FCM_PART_RADII.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_RADII,FILENAME,ERR_FILE_RADII)



TIME = NCYCLEMAX + 1
if ( (ERR_FILE_POS.eq.0).and.(ERR_FILE_SWIM.eq.0) ) then
 call PRINT_TECPLOT_PART_POS_PSWIM(TIME,NPART, &
      FCM_POS(:,1,1)/RADMAX,FCM_POS(:,2,1)/RADMAX,FCM_POS(:,3,1)/RADMAX, &
      FCM_PSWIM(:,1,1),FCM_PSWIM(:,2,1),FCM_PSWIM(:,3,1),BINFLAG)
end if

if ( (ERR_FILE_POS.eq.0).and.(ERR_FILE_VEL.eq.0) ) then
 call PRINT_TECPLOT_PART_POS_VEL(TIME,NPART, &
      FCM_POS(:,1,1)/RADMAX,FCM_POS(:,2,1)/RADMAX,FCM_POS(:,3,1)/RADMAX, &
      FCM_VEL(:,1,1)/VSW,FCM_VEL(:,2,1)/VSW,FCM_VEL(:,3,1)/VSW,BINFLAG)
end if



!=====================================================================
! 5. READ INTERMEDIATE LAGRANGIAN DATA
!=====================================================================

do IND=1,NB_DUMP_FILES_PART

 TIME = (IND-1)*FOUT3 + 1
 

 
 write(FILENAME,10205)'FCM_PART_POS_t',TIME,'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS) 

 write(FILENAME,10205)'FCM_PART_VEL_t',TIME,'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_VEL,FILENAME,ERR_FILE_VEL)
 
!~  write(FILENAME,10101)'FCM_PART_ROT_t',TIME,'.bin'
!~  call READ_VAR_MPIIO(NPART,3,1,FCM_ROT,FILENAME,ERR_FILE_ROT)
 
 
 write(FILENAME,10205)'FCM_PART_SWIM_t',TIME,'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_PSWIM,FILENAME,ERR_FILE_SWIM)

  
 !=====================================================================
 ! 6. SAVE INTERMEDIATE LAG DATA IN ASCII
 !=====================================================================

 !!--------------------------------------------------------------------
 !! 6.1 Positions + swimming direction
 !!--------------------------------------------------------------------

  if ( (ERR_FILE_POS.eq.0).and.(ERR_FILE_SWIM.eq.0) ) then
   call PRINT_TECPLOT_PART_POS_PSWIM(TIME,NPART, &
        FCM_POS(:,1,1)/RADMAX,FCM_POS(:,2,1)/RADMAX,FCM_POS(:,3,1)/RADMAX, &
        FCM_PSWIM(:,1,1),FCM_PSWIM(:,2,1),FCM_PSWIM(:,3,1),BINFLAG)
  end if
  
  
  !!--------------------------------------------------------------------
 !! 6.2 Positions + velocity
 !!--------------------------------------------------------------------
  if ( (ERR_FILE_POS.eq.0).and.(ERR_FILE_VEL.eq.0) ) then
   call PRINT_TECPLOT_PART_POS_VEL(TIME,NPART, &
        FCM_POS(:,1,1)/RADMAX,FCM_POS(:,2,1)/RADMAX,FCM_POS(:,3,1)/RADMAX, &
        FCM_VEL(:,1,1)/VSW,FCM_VEL(:,2,1)/VSW,FCM_VEL(:,3,1)/VSW,BINFLAG)
  end if


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


end program MPIIO_TO_TECPLOT
