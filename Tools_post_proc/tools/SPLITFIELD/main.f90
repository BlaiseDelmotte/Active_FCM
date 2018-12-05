!!====================================================================
!! 
!!  
!!====================================================================

program INTERP_SOLUTION

!!====================================================================
!!
!!
!!====================================================================

use DNS_DIM	       !- Dimension
use MPI_STRUCTURES
use FLUID_VARIABLE
use GEOMETRIC


implicit none


integer :: DXCPU, DYCPU


!!- File name 
character(len=40) :: FILENAME
character(len=40), dimension(3,2) :: FILELIST


!!- 
real(kind=8) :: MAX_VAR, MAX_VARI
real(kind=8) :: MEAN_VAR, MEAN_VARI
real(kind=8) :: MEAN_VAR2, MEAN_VARI2

integer :: IORDER

integer :: NBFILE



integer :: I, J ,K, N
!---------------------------------------------------------------------


NBFILE = 3

FILELIST(1,1) = 'uf.end'
FILELIST(1,2) = 'uf.int'
FILELIST(2,1) = 'vf.end'
FILELIST(2,2) = 'vf.int'
FILELIST(3,1) = 'wf.end'
FILELIST(3,2) = 'wf.int'


IORDER = 4

!- Size of the file that is read
NX = 128
NY = NX
NZ = NY


NXI = 1024
NYI = NXI
NZI = NYI


LXMAX = 0.128
LYMAX = LXMAX
LZMAX = LXMAX



call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD,NPROC,IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD,MYID,IERR)


!!--------------------------------------------------------------------
!! Domain splitting
!!--------------------------------------------------------------------

!- Dimentionality of the cpu decomposition (=1: slab, =2: squared)
NDIM = 2

if(NDIM == 1) then
 DIMS(1) = 1
 DIMS(2) = NPROC
else if(NDIM == 2) then
! if (MYID==0) print *, 'Creating proc. grid with mpi_dims_create'
 DIMS(1) = 0
 DIMS(2) = 0
 call MPI_DIMS_CREATE(NPROC,2,DIMS,IERR)
 if(DIMS(1) > DIMS(2)) then
  DIMS(1) = DIMS(2)
  DIMS(2) = NPROC / DIMS(1)
 endif
endif

IPROC = DIMS(1)
JPROC = DIMS(2)

!if(MYID == 0)write(*,*)'Using processor grid ',iproc,' x ',jproc




!write(*,1608)'NDIM = ',NDIM 
!write(*,1608)'IPROC = ',IPROC
!write(*,1608)'JPROC = ',JPROC
!write(*,*)'----------------------------'
!write(*,1608)'Nx = ',NX
!write(*,1608)'Ny = ',NY
!write(*,1608)'Nz = ',NZ



call GET_DIMS_NEW(MYID,NX ,NY ,NZ ,IPROC,JPROC,ISTART ,IEND ,ISIZE )
call GET_DIMS_NEW(MYID,NXI,NYI,NZI,IPROC,JPROC,ISTARTI,IENDI,ISIZEI)


!!- FFt Initiation
!!call P3DFFT_SETUP(DIMS,NX,NY,NZ,.TRUE.)


!!- Split the geometry
!!call P3DFFT_GET_DIMS(ISTARTI,IENDI,ISIZEI,1)
!!call P3DFFT_GET_DIMS(FSTARTI,FENDI,FSIZEI,2)
!!call GET_DIMS(ISTART,IEND,ISIZE,1)
!!call GET_DIMS(FSTARTI,FENDI,FSIZEI,2)


!write(*,1609)'My Id=',MYID,' IS = ',ISTART(1),' ISI = ',ISTARTI(1),' IE = ',IEND(1),' IEI = ',IENDI(1)
!write(*,1609)'My Id=',MYID,' JS = ',ISTART(2),' JSI = ',ISTARTI(2),' JE = ',IEND(2),' JEI = ',IENDI(2)
!write(*,1609)'My Id=',MYID,' KS = ',ISTART(3),' KSI = ',ISTARTI(3),' KE = ',IEND(3),' KEI = ',IENDI(3)


!- Clean P3DFFT
!!call P3DFFT_CLEAN





!!--------------------------------------------------------------------
!! Allocate array
!!--------------------------------------------------------------------

allocate(TEMP(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )



allocate(TEMPI(ISTARTI(1)         :IENDI(1)             &
              ,ISTARTI(2)-NGHTCELL:IENDI(2)+NGHTCELL    &
              ,ISTARTI(3)-NGHTCELL:IENDI(3)+NGHTCELL )  )


!- Allocate arrays for mesh
allocate(XMESH(ISTART(1):IEND(1)))
allocate(YMESH(ISTART(2):IEND(2)))
allocate(ZMESH(ISTART(3):IEND(3)))

!- Allocate arrays for mesh
allocate(XMESHI(ISTARTI(1):IENDI(1)))
allocate(YMESHI(ISTARTI(2):IENDI(2)))
allocate(ZMESHI(ISTARTI(3):IENDI(3)))



!!- Create the two meshes
call MESHING




!write(*,1610)'My Id=',MYID,' xs = ', XMESH(ISTART(1))    &
!                          ,' xsi = ',XMESHI(ISTARTI(1)) &
!                          ,' xe = ', XMESH(IEND(1))      &
!                          ,' xei = ',XMESHI(IENDI(1))

!write(*,1610)'My Id=',MYID,' ys = ', YMESH(ISTART(2))    &
!                          ,' ysi = ',YMESHI(ISTARTI(2)) &
!                          ,' ye = ', YMESH(IEND(2))      &
!                          ,' yei = ',YMESHI(IENDI(2))

!write(*,1610)'My Id=',MYID,' zs = ',ZMESH(ISTART(3))    &
!                          ,' zsi = ',ZMESHI(ISTARTI(3)) &
!                          ,' ze = ',ZMESH(IEND(3))      &
!                          ,' zei = ',ZMESHI(IENDI(3))






! Create MPI structures for ghost cells MPI-exchange
call CREATE_MPI_VECTOR

! Find neighbouring for each processor
call NEIGHBOURING




!!- Loop on the file number
do N = 1, NBFILE

if(MYID==0) write(*,*)
if(MYID==0) write(*,*)
if(MYID==0) write(*,*)'============================================='
if(MYID==0) write(*,*)'== Processing file=',N
if(MYID==0) write(*,*)'============================================='
if(MYID==0) write(*,*)



FILENAME=trim(FILELIST(N,1))

call READ_MPIIO(TEMP,FILENAME)
if(MYID==0) write(*,*)'Read initial solution --> OK'
if(MYID==0) write(*,*)'             FILENAME --> ',trim(FILENAME)



!!- Halo updating
call FLUIDCOMM(TEMP)
if(MYID==0) write(*,*)'Halo updating --> OK'



!!=====================================
!! 0th order Lagrangian interpolation
!!=====================================
if(IORDER==0) then

 call INTERP_LAG0( XMESH,YMESH,ZMESH,   &
                   TEMP,                 &
                   XMESHI,YMESHI,ZMESHI, &
                   TEMPI                )


!!=====================================
!! 1st order Lagrangian interpolation
!!=====================================
elseif(IORDER==1) then

 call INTERP_LAG1( XMESH,YMESH,ZMESH,   &
                    TEMP,                 &
                    XMESHI,YMESHI,ZMESHI, &
                    TEMPI                )


!!=====================================
!! 2nd order Lagrangian interpolation
!!=====================================
elseif(IORDER==2) then

 call INTERP_LAG2( XMESH, YMESH, ZMESH, &
                    TEMP,               &
                  XMESHI,YMESHI,ZMESHI, &
                   TEMPI                )


!!=====================================
!! 3rd order Lagrangian interpolation
!!=====================================
elseif(IORDER==3) then

 call INTERP_LAG3( XMESH, YMESH, ZMESH, &
                    TEMP,               &
                  XMESHI,YMESHI,ZMESHI, &
                   TEMPI                )

!!=====================================
!! 4th order Lagrangian interpolation
!!=====================================
elseif(IORDER==4) then

 call INTERP_LAG4(XMESH,YMESH,ZMESH,   &
                     TEMP,               &
                   XMESHI,YMESHI,ZMESHI, &
                    TEMPI                )

!!=====================================
!! Security
!!=====================================
else
 IORDER=0
end if


call RSUMCPU(SUM( TEMP( ISTART(1): IEND(1), ISTART(2): IEND(2), ISTART(3): IEND(3)))   ,MEAN_VAR)
call RSUMCPU(SUM( TEMP( ISTART(1): IEND(1), ISTART(2): IEND(2), ISTART(3): IEND(3))**2),MEAN_VAR2)
call RMAXCPU(abs( TEMP( ISTART(1): IEND(1), ISTART(2): IEND(2), ISTART(3): IEND(3)))   ,MAX_VAR)

call RSUMCPU(SUM(TEMPI(ISTARTI(1):IENDI(1),ISTARTI(2):IENDI(2),ISTARTI(3):IENDI(3)))   ,MEAN_VARI)
call RSUMCPU(SUM(TEMPI(ISTARTI(1):IENDI(1),ISTARTI(2):IENDI(2),ISTARTI(3):IENDI(3))**2),MEAN_VARI2)
call RMAXCPU(abs(TEMPI(ISTARTI(1):IENDI(1),ISTARTI(2):IENDI(2),ISTARTI(3):IENDI(3)))   ,MAX_VARI)

MEAN_VAR = MEAN_VAR /NX/NY/NZ
MEAN_VAR2 = MEAN_VAR2 /NX/NY/NZ

MEAN_VARI = MEAN_VARI /NXI/NYI/NZI
MEAN_VARI2 = MEAN_VARI2 /NXI/NYI/NZI



if (MYID==0) then
write(*,*)'Before interpolation'
write(*,*)'max(var) = ',MAX_VAR
write(*,*)'   <var> = ',MEAN_VAR
write(*,*)'  <var2> = ',MEAN_VAR2
write(*,*)
write(*,*)'After interpolation'
write(*,*)'max(var_inter) = ',MAX_VARI
write(*,*)'   <var_inter> = ',MEAN_VARI
write(*,*)'  <var_inter2> = ',MEAN_VARI2
write(*,*)
write(*,*)'max(var)/max(var_inter) = ',MAX_VAR/MAX_VARI
!!write(*,*)'      <var>/<var_inter> = ',MEAN_VAR/MEAN_VARI
write(*,*)'    <var2>/<var_inter2> = ',MEAN_VAR2/MEAN_VARI2
end if




!TEMPI(:,:,:) = MYID
!FILENAME='split.int'
!call SAVE_MPIIO(TEMPI,FILENAME)

FILENAME=trim(FILELIST(N,2))
call SAVE_MPIIO(TEMPI,FILENAME)
if(MYID==0) write(*,*)'Printing interpolated solution --> OK'
if(MYID==0) write(*,*)'                      FILENAME --> ',trim(FILENAME)



end do !!- do N=1,NBFILE



!- Free MPI environment
 call MPI_FINALIZE (ierr)


!!----------------------------------------------------------------------
1600 format (1x,A,I3,A,I3,A,I3)
1608 format (1x,A,I5)
1609 format (1x,A,I5,A,I5,A,I5,A,I5,A,I5,A,I5)
1610 format (1x,A,I5,A,E13.5,A,E13.5,A,E13.5,A,E13.5,A,E13.5)

end program INTERP_SOLUTION

