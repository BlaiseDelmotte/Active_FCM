!!====================================================================
!!
!!
!!====================================================================

Program TREAT_PARTICLES_BIDISPERSE

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
!- String for time
character(len=10) :: FILE_EXT

! Physical Time
real(kind=8), allocatable, dimension(:) :: TIME_VEC

!!- Vraiables from direct access files
! Read Lagrangian Data
real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS
real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS_NOPER
real(kind=8), allocatable, dimension(:,:,:) :: FCM_VEL
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ROT
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM
real(kind=8), allocatable, dimension(:,:,:) :: FCM_SIJ
real(kind=8), allocatable, dimension(:,:,:) :: FCM_RADII

! Lagrangian Data along time
real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS_TIME
real(kind=8), allocatable, dimension(:,:,:) :: FCM_POS_NOPER_TIME
real(kind=8), allocatable, dimension(:,:,:) :: FCM_VEL_TIME
real(kind=8), allocatable, dimension(:,:,:) :: FCM_ROT_TIME
real(kind=8), allocatable, dimension(:,:,:,:) :: FCM_ROT_MAT_TIME
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM_TIME
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM_TIME_SPH
real(kind=8), allocatable, dimension(:,:,:) :: FCM_PSWIM_TIME_SPH_PMEAN
real(kind=8), allocatable, dimension(:,:,:) :: FCM_SIJ_TIME

real(kind=8), allocatable, dimension(:,:) :: MEAN_PSWIM_TIME
real(kind=8) :: PPI

!!- Dicrestization points in each direction
integer :: NX, NY, NZ
real(kind=8) :: LXMAX, LYMAX, LZMAX

!!- Swilmming velocity (radii/s)
real(kind=8) :: VSW

!!- Minimal/Maximal radius, to compute dimensionless time
real(kind=8) :: RADMIN
real(kind=8) :: RADMAX

real(kind=8), dimension(2) :: SPH_SIZE

!!- Volume fraction
real(kind=8) :: VOL_FRAC

!!- Asymptotic Value Mean Orientation vector
real(kind=8), dimension(3) :: PINF

!! Flags for Brownian + Wall
integer :: FLOW_TYPE, BC

!!- Number of boxes in ones direction to compute local statistics
integer :: NBOX_DIR

!!- Particle number temporary
integer :: NPART, NPSTAT, NSPHERE, NSPHERE_1, NELLIPSOID
integer :: PART_START_1, PART_END_1, PART_START_2, PART_END_2
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

!!- Number of dumps to read
integer :: NB_DUMP_FILES_FLUID
integer :: NB_DUMP_FILES_PART
!! Starting dump to treat correlations
integer :: SAVE_START
integer:: NSAVES_MSD
!!- Flag if quaternions are used or orientation vectors
integer :: USE_QUAT

!- Time 
integer :: TIME

!- Index
integer :: I, J, K, IJK, IND

!- Integers to signalize a non-existing file
integer :: ERR_STAT, &
           ERR_INFO, &
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
	   
!!- Integers to chose which treatment to perform
integer :: TREAT_VEL
integer :: TREAT_STRESS
integer :: TREAT_CORREL_POS
integer :: TREAT_CORREL_ORIENT
integer :: TREAT_DIFF

!-  MPI Variables
integer :: IERR,NPROC,MYID



!!--------------------------------------------------------------------
!! 0.0 MPI World initiation
!!--------------------------------------------------------------------
call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)


PPI = 4.0*atan(1.0)

!---------------------------------------------------------------------
!=====================================================================
! 0. READ INPUT, PARAMETERS AND ALLOCATION
!=====================================================================

!- Open file
open(unit=100,file='stat_choice.in',status='old', iostat=ERR_STAT)

if(ERR_STAT /= 0) then

 print *, "CPU -- ", MYID, ":: ERROR: Input data file stat_choice.in open error!"
 call MPI_FINALIZE(ERR_STAT)
 call abort()
else

 read(100,*), TREAT_DIFF
 read(100,*), TREAT_VEL
 read(100,*), TREAT_STRESS
 read(100,*), TREAT_CORREL_POS
 read(100,*), TREAT_CORREL_ORIENT
 read(100,*), NPSTAT
 read(100,*), SAVE_START
 read(100,*), PINF(1), PINF(2), PINF(3)
 read(100,*), NBOX_DIR
if (TREAT_DIFF==1) then
 read(100,*), NSAVES_MSD
end if

end if
close(100)



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
 read(200,*), VSW
 read(200,*), RADMAX
 read(200,*), USE_QUAT
 read(200,*), SPH_SIZE
 read(200,*), FLOW_TYPE
 read(200,*), BC

end if

close(200)




NB_DUMP_FILES_PART= (NCYCLEMAX)/FOUT3
print*,'NB_DUMP_FILES_PART + 1 = ' ,  NB_DUMP_FILES_PART + 1


allocate( TIME_VEC(NB_DUMP_FILES_PART+1) )

allocate( FCM_RADII(NPART,3,1) )

allocate( FCM_POS(NPART,3,1) )
allocate( FCM_POS_TIME(NB_DUMP_FILES_PART+1,NPART,3) )
allocate( FCM_PSWIM(NPART,3,1) )

allocate( FCM_PSWIM_TIME(NB_DUMP_FILES_PART+1,NPART,3) )
allocate( MEAN_PSWIM_TIME(NB_DUMP_FILES_PART+1,3) )

allocate( FCM_PSWIM_TIME_SPH(NB_DUMP_FILES_PART+1,NPART,2) )
allocate( FCM_PSWIM_TIME_SPH_PMEAN(NB_DUMP_FILES_PART+1,NPART,2) )

if (TREAT_DIFF==1) then
 allocate( FCM_POS_NOPER(NPART,3,1) )
 allocate( FCM_POS_NOPER_TIME(NB_DUMP_FILES_PART+1,NPART,3) )
end if


if (TREAT_VEL==1) then
 allocate( FCM_VEL(NPART,3,1) )
 allocate( FCM_ROT(NPART,3,1) )
 allocate( FCM_VEL_TIME(NB_DUMP_FILES_PART+1,NPART,3) )
 allocate( FCM_ROT_TIME(NB_DUMP_FILES_PART+1,NPART,3) )
end if



if (TREAT_STRESS==1) then
 allocate( FCM_SIJ(NPART,5,1) )
 allocate( FCM_SIJ_TIME(NB_DUMP_FILES_PART+1,NPART,5) )
end if






print*,' '
print*,'-------------------------START READING---------------------- '
print*,' '
!=====================================================================
! 1. READ LAGRANGIAN DATA AT ENDING TIME
!=====================================================================

FILENAME='FCM_PART_RADII.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_RADII,FILENAME,ERR_FILE_RADII)

if (ERR_FILE_RADII==0) then
 RADMIN = minval(FCM_RADII)
 RADMAX = maxval(FCM_RADII)

 VOL_FRAC = 0.0
 do I = 1, NPART 
  VOL_FRAC = VOL_FRAC &
           + 4.0/3.0*PPI*FCM_RADII(I,1,1)*FCM_RADII(I,2,1)*FCM_RADII(I,3,1) &
           / (LXMAX*LYMAX*LZMAX)
 end do

 print*,'VOLUME FRACTION = ', VOL_FRAC
else
 VOL_FRAC = 0.0
 RADMIN = VSW/2.0
 RADMAX = VSW
end if

FILENAME='FCM_PART_POS.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS)
FILENAME='FCM_PART_SWIM.end'
call READ_VAR_MPIIO(NPART,3,1,FCM_PSWIM,FILENAME,ERR_FILE_SWIM)

FCM_POS_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_POS(:,:,1)
FCM_PSWIM_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_PSWIM(:,:,1)

if (TREAT_DIFF==1) then
 FILENAME='FCM_PART_POS_NOPER.end'
 call READ_VAR_MPIIO(NPART,3,1,FCM_POS_NOPER,FILENAME,ERR_FILE_POS_NOPER)
  
 FCM_POS_NOPER_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_POS_NOPER(:,:,1)
end if

if (TREAT_VEL==1) then
 FILENAME='FCM_PART_VEL.end'
 call READ_VAR_MPIIO(NPART,3,1,FCM_VEL,FILENAME,ERR_FILE_VEL)
 FILENAME='FCM_PART_ROT.end'
 call READ_VAR_MPIIO(NPART,3,1,FCM_ROT,FILENAME,ERR_FILE_ROT)
 
 FCM_VEL_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_VEL(:,:,1)
 FCM_ROT_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_ROT(:,:,1)
end if


if (TREAT_STRESS==1) then
 FILENAME='FCM_PART_STRESSLET.end'
 call READ_VAR_MPIIO(NPART,5,1,FCM_SIJ,FILENAME,ERR_FILE_STRESS)
 
 FCM_SIJ_TIME(NB_DUMP_FILES_PART+1,:,:) = FCM_SIJ(:,:,1)
end if


 

TIME_VEC(NB_DUMP_FILES_PART+1) = (real(NCYCLEMAX)-1.0)*DTIME*VSW


!=====================================================================
! 2. READ INTERMEDIATE LAGRANGIAN DATA
!=====================================================================

do IND=1,NB_DUMP_FILES_PART


 if (NB_DUMP_FILES_PART.ge.5) then
  if ((mod(IND,NB_DUMP_FILES_PART/5)==0).or.(IND==NB_DUMP_FILES_PART)) then
   print*, IND, '/', NB_DUMP_FILES_PART+1, 'FILES READ'
  end if
 end if

 TIME = (IND-1)*FOUT3 + 1
 TIME_VEC(IND) = real(TIME)*DTIME*VSW
 
 write(FILE_EXT,10205) TIME
 
 write(FILENAME,10101)'FCM_PART_POS_t',trim(FILE_EXT),'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_POS,FILENAME,ERR_FILE_POS) 
 write(FILENAME,10101)'FCM_PART_SWIM_t',trim(FILE_EXT),'.bin'
 call READ_VAR_MPIIO(NPART,3,1,FCM_PSWIM,FILENAME,ERR_FILE_SWIM)

 FCM_POS_TIME(IND,:,:) = FCM_POS(:,:,1)
 FCM_PSWIM_TIME(IND,:,:) = FCM_PSWIM(:,:,1)

 if (TREAT_DIFF==1) then
  write(FILENAME,10101)'FCM_PART_POS_NOPER_t',trim(FILE_EXT),'.bin'
  call READ_VAR_MPIIO(NPART,3,1,FCM_POS_NOPER,FILENAME,ERR_FILE_POS_NOPER) 
  
!~   if (IND == 1024) then
!~    print*, ' FCM_POS_NOPER(1181+12,:,1) = ',   FCM_POS_NOPER(1181+12,:,1)
!~    read(*,*)
!~   end if 
  FCM_POS_NOPER_TIME(IND,:,:) = FCM_POS_NOPER(:,:,1)
 end if

 if (TREAT_VEL==1) then
  write(FILENAME,10101)'FCM_PART_VEL_t',trim(FILE_EXT),'.bin'
  call READ_VAR_MPIIO(NPART,3,1,FCM_VEL,FILENAME,ERR_FILE_VEL)  
  write(FILENAME,10101)'FCM_PART_ROT_t',trim(FILE_EXT),'.bin'
  call READ_VAR_MPIIO(NPART,3,1,FCM_ROT,FILENAME,ERR_FILE_ROT)

  FCM_VEL_TIME(IND,:,:) = FCM_VEL(:,:,1)
  FCM_ROT_TIME(IND,:,:) = FCM_ROT(:,:,1)
 end if
 
 
 
 if (TREAT_STRESS==1) then 
  write(FILENAME,10101)'FCM_PART_STRESSLET_t',trim(FILE_EXT),'.bin'
  call READ_VAR_MPIIO(NPART,5,1,FCM_SIJ,FILENAME,ERR_FILE_STRESS)
  
  FCM_SIJ_TIME(IND,:,:) = FCM_SIJ(:,:,1)
 end if
 

 
end do


print*,' '
print*,'-------------------------END READING---------------------- '
print*,' '


!=====================================================================
! 3. TREAT DATA
!=====================================================================

print*,'SAVE_START = ', SAVE_START
print*,'NPART = ', NPART
print*,'SPH_SIZE = ', SPH_SIZE
print*,'LXMAX = ',LXMAX

if ( (SPH_SIZE(2).eq.SPH_SIZE(1)).or.(NSPHERE_1.eq.NPART) ) then
 print*,'Monodisperse stat'
 PART_START_1 = 1
 PART_END_1 = NSPHERE_1
 PART_START_2 = -1
 PART_END_2 = -1
else if (SPH_SIZE(1).ne.SPH_SIZE(2)) then
 print*,'Bidisperse stat' 
 PART_START_1 = 1
 PART_END_1 = NSPHERE_1
 PART_START_2 = NSPHERE_1+1
 PART_END_2 = NPART
end if !SPH_SIZE

print*, 'PART_START_1 = ', PART_START_1
print*, 'PART_END_1 = ', PART_END_1
print*, 'PART_START_2 = ', PART_START_2
print*, 'PART_END_2 = ', PART_END_2

if (TREAT_DIFF==1) then
 print*, 'NSAVES_MSD = ', NSAVES_MSD
end if

 if (BC == 2) then
  call HEIGHT_DISTRIB(NB_DUMP_FILES_PART+1, &
                      SAVE_START, & 
                      PART_START_1, &
                      PART_END_1, &
                      LXMAX, &
                      RADMAX, &
                      FCM_POS_TIME(SAVE_START:NB_DUMP_FILES_PART,PART_START_1:PART_END_1,:) )
  if (PART_START_2>0) then                     
  call HEIGHT_DISTRIB(NB_DUMP_FILES_PART+1, &
                      SAVE_START, &                     
                      PART_START_2, &
                      PART_END_2, &
                      LXMAX, &
                      RADMAX, &
                      FCM_POS_TIME(SAVE_START:NB_DUMP_FILES_PART,PART_START_2:PART_END_2,:) )
  end if                                            
 end if
 
  if (TREAT_VEL==1) then
  call PDF_VELOCITIES(NB_DUMP_FILES_PART+1, &
                      SAVE_START, &
                      PART_START_1, &
                      PART_END_1, &
                      FCM_RADII(PART_START_1,1,1), &
                      FCM_VEL_TIME(:,PART_START_1:PART_END_1,:) )
  if (PART_START_2>0) then                    
   call PDF_VELOCITIES(NB_DUMP_FILES_PART+1, &
                      SAVE_START, &
                      PART_START_2, &
                      PART_END_2, &
                      FCM_RADII(PART_START_2,1,1), &
                      FCM_VEL_TIME(:,PART_START_2:PART_END_2,:) )   
  end if
  if (TREAT_DIFF==1) then
   call SHORT_TIME_SELF_DIFF(NB_DUMP_FILES_PART+1, &
                      SAVE_START, &
                      PART_START_1, &
                      PART_END_1, &
                      DTIME, &
                      FCM_VEL_TIME(:,PART_START_1:PART_END_1,:), &
                      FCM_ROT_TIME(:,PART_START_1:PART_END_1,:)  )
   if (PART_START_2>0) then                    
    call SHORT_TIME_SELF_DIFF(NB_DUMP_FILES_PART+1, &
                      SAVE_START, &
                      PART_START_2, &
                      PART_END_2, &
                      DTIME, &
                      FCM_VEL_TIME(:,PART_START_2:PART_END_2,:), &
                      FCM_ROT_TIME(:,PART_START_2:PART_END_2,:)  ) 
   end if
  end if !(TREAT_DIFF==1)

 if (TREAT_DIFF==1) then
 
   call MSD_DISP(NB_DUMP_FILES_PART+1,&
             SAVE_START, &
             PART_START_1, &
             PART_END_1, &
             NSAVES_MSD, &
             FCM_RADII(PART_START_1,1,1), &
             LXMAX, &
             LYMAX, &
             LZMAX, &
              FCM_POS_NOPER_TIME(SAVE_START:NB_DUMP_FILES_PART+1,PART_START_1:PART_END_1,:) )  
!~   !!           FCM_POS_TIME(SAVE_START:NB_DUMP_FILES_PART+1,PART_START_1:PART_END_1,:) )  


   if (PART_START_2>0) then    
                      
    call MSD_DISP(NB_DUMP_FILES_PART+1,&
             SAVE_START, &
             PART_START_2, &
             PART_END_2, &
             NSAVES_MSD, &
             FCM_RADII(PART_START_2,1,1), &
             LXMAX, &
             LYMAX, &
             LZMAX, &
             FCM_POS_NOPER_TIME(SAVE_START:NB_DUMP_FILES_PART+1,PART_START_2:PART_END_2,:) )      
!~    !!          FCM_POS_TIME(SAVE_START:NB_DUMP_FILES_PART+1,PART_START_2:PART_END_2,:) ) 
   end if 
                                                     
 end if  


!~                                            
!~   call PDF_ROTATIONS(NB_DUMP_FILES_PART+1, &
!~                     SAVE_START, &
!~                     PART_START_1, &
!~                     PART_END_1, &
!~                     FCM_ROT_TIME(:,PART_START_1:PART_END_1,:) )
!~   if (PART_START_2>0) then                     
!~    call PDF_ROTATIONS(NB_DUMP_FILES_PART+1, &
!~                     SAVE_START, &
!~                     PART_START_2, &
!~                     PART_END_2, &
!~                     FCM_ROT_TIME(:,PART_START_2:PART_END_2,:) ) 
!~   end if
!~                                        
!~   call MEAN_VEL(NB_DUMP_FILES_PART+1, &
!~                 PART_START_1, &
!~                 PART_END_1, &
!~                 TIME_VEC, &
!~                 FCM_VEL_TIME(:,PART_START_1:PART_END_1,:) )
!~   call MEAN_VEL(NB_DUMP_FILES_PART+1, &
!~                 PART_START_2, &
!~                 PART_END_2, &
!~                 TIME_VEC, &
!~                 FCM_VEL_TIME(:,PART_START_2:PART_END_2,:) )                
 end if
!~ 
 call POLAR_ORDER_ONLY(NB_DUMP_FILES_PART+1, &
                  PART_START_1, &
                  PART_END_1, &
                  NBOX_DIR, &
                  LXMAX, &
                  LYMAX, &
                  LZMAX, &
                  TIME_VEC, &
                  FCM_PSWIM_TIME(:,PART_START_1:PART_END_1,:), &
                   FCM_POS_TIME(:,PART_START_1:PART_END_1,:), &             
                  MEAN_PSWIM_TIME )
                  
                
                  
!~  call PDF_ORIENTATIONS(NB_DUMP_FILES_PART+1, &
!~ 																								SAVE_START, &
!~ 																								PART_START_1, &
!~ 																								PART_END_1, &
!~ 																								FCM_PSWIM_TIME(:,PART_START_1:PART_END_1,:), &
!~ 																								MEAN_PSWIM_TIME, &
!~ 																								FCM_PSWIM_TIME_SPH(:,PART_START_1:PART_END_1,:), &
!~ 																								FCM_PSWIM_TIME_SPH_PMEAN(:,PART_START_1:PART_END_1,:), &
!~ 																								PINF)

 
if ((TREAT_CORREL_POS==1).and.(TREAT_VEL==1)&
                         .and.(TREAT_CORREL_ORIENT==1)) then 
                         
 if (PART_START_2>0) then                         
                                    
 call POSITION_VELOCITY_ORIENT_CORREL_BIDISPERSE(NB_DUMP_FILES_PART+1, &
                                                 SAVE_START, &
                                                 PART_START_1, &
                                                 PART_END_1, &
                                                 PART_START_2, &
                                                 PART_END_2,&
                                                 NPART, &
                                                 NPSTAT, &
                                                 LXMAX, &
                                                 FCM_RADII(:,1,1), &
                                                 FCM_PSWIM_TIME, &
                                                 FCM_POS_TIME, &
                                                 FCM_VEL_TIME )  
 else
  call POSITION_VELOCITY_ORIENT_CORREL(NB_DUMP_FILES_PART+1, &
                                                 SAVE_START, &
                                                 NPART, &
                                                 NPSTAT, &
                                                 LXMAX, &
                                                 maxval(FCM_RADII(:,1,1)), &
                                                 FCM_ROT_MAT_TIME, &
                                                 FCM_POS_TIME, &
                                                 FCM_VEL_TIME ) 
 end if
end if 
!~ 
!~  if (NELLIPSOID>0) then
!~ 
!~   if ((TREAT_CORREL_POS==1).and.(TREAT_VEL==1)&
!~                            .and.(TREAT_CORREL_ORIENT==1)) then  
!~                                       
!~    call POSITION_VELOCITY_ORIENT_CORREL_PROLATE(NB_DUMP_FILES_PART+1, &
!~                                                 SAVE_START, &
!~                                                 NPART, &
!~                                                 NPSTAT, &
!~                                                 LXMAX, &
!~                                                 RADMIN, &
!~                                                 RADMAX, &
!~                                                 FCM_ROT_MAT_TIME, &
!~                                                 FCM_POS_TIME, &
!~                                                 FCM_VEL_TIME ) 
!~                                              
!~ 
!~   end if
!~   
!~  else
!~  
!~  end if



print*,' '
print*,'PARTICLE POST PROCESSING--->DONE'
print*,' '
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

10200 format (A)
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10205 format (I8.8)
10101 format (A,A,A)

end program TREAT_PARTICLES_BIDISPERSE
