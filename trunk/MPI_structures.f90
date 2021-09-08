   module MPI_structures

   use DNS_DIM 
   use GEOMETRIC_VARIABLE
   use COLLISION_VARIABLE

  implicit none

  integer :: XZ_2PLANES,XY_2PLANES,CORNER
  INTEGER, PARAMETER               :: NB_NEIGHBORS = 8
  INTEGER, DIMENSION(NB_NEIGHBORS) :: NEIGHBOR
  INTEGER, PARAMETER               :: FJ=1,BJ=2,FK=3,BK=4 ! Neighbors : Forward_J, Backward_J, ... 
  INTEGER, PARAMETER               :: FJBK=5,BJFK=6,BJBK=7,FJFK=8 ! Common edge neighbor

  INTEGER :: COUNT_STAY
  integer, dimension(:,:), allocatable :: COUNTER ! nb of particles leaving in each direction  

  integer, dimension(:,:), allocatable :: IND_LEAV
  integer, dimension(:), allocatable :: IND_STAY



contains

  SUBROUTINE NEIGHBOURING

    implicit none

    !*******************************************************************
    ! Neighbouring table Initiation
    NEIGHBOR(:) = MPI_PROC_NULL

    NEIGHBOR(FJ) = MYID+1
    NEIGHBOR(BJ) = MYID-1
    if(mod(MYID+1,JPROC)==0) NEIGHBOR(FJ) = MYID+1-JPROC
    if(mod(MYID  ,JPROC)==0) NEIGHBOR(BJ) = MYID-1+JPROC

    NEIGHBOR(FK) = mod(MYID+JPROC,NPROC)
    NEIGHBOR(BK) = mod(MYID-JPROC+NPROC,NPROC)
  
    NEIGHBOR(FJFK) = NEIGHBOR(FK) + 1
    if(mod(NEIGHBOR(FJFK),JPROC)==0) NEIGHBOR(FJFK) = NEIGHBOR(FK)+1-JPROC
    NEIGHBOR(FJBK) = NEIGHBOR(BK) + 1
    if(mod(NEIGHBOR(FJBK),JPROC)==0) NEIGHBOR(FJBK) = NEIGHBOR(BK)+1-JPROC

    NEIGHBOR(BJFK) = NEIGHBOR(FK) - 1
    if(mod(NEIGHBOR(FK),JPROC)==0) NEIGHBOR(BJFK) = NEIGHBOR(FK)-1+JPROC
    NEIGHBOR(BJBK) = NEIGHBOR(BK) - 1
    if(mod(NEIGHBOR(BK),JPROC)==0) NEIGHBOR(BJBK) = NEIGHBOR(BK)-1+JPROC

  end subroutine NEIGHBOURING


  subroutine CREATE_MPI_VECTOR

   implicit none 

   ! xOz 2 planes 
   CALL MPI_TYPE_VECTOR(IEND(3)-ISTART(3)+1, & ! blocks number
                        (IEND(1)-ISTART(1)+1)*NGHTCELL, & ! block size
                        (IEND(1)-ISTART(1)+1)*(IEND(2)-ISTART(2)+NGHTCELL*2+1),& ! step between 2 blocks
                           MPI_DOUBLE_PRECISION,XZ_2PLANES,IERR) ! name of MPI_TYPE
   CALL MPI_TYPE_COMMIT(XZ_2PLANES,IERR)

   ! xOy 2 planes
   CALL MPI_TYPE_VECTOR(2, & ! blocks number
                        (IEND(1)-ISTART(1)+1)*(IEND(2)-ISTART(2)+1), & ! block size
                        (IEND(1)-ISTART(1)+1)*(IEND(2)-ISTART(2)+NGHTCELL*2+1), & ! step between 2 blocks 
                           MPI_DOUBLE_PRECISION,XY_2PLANES,IERR) ! name of MPI_TYPE
   CALL MPI_TYPE_COMMIT(XY_2PLANES,IERR)

   ! Corner 
   CALL MPI_TYPE_VECTOR(2, & ! blocks number
                        (IEND(1)-ISTART(1)+1)*NGHTCELL, & ! block size
                        (IEND(1)-ISTART(1)+1)*(IEND(2)-ISTART(2)+NGHTCELL*2+1), & ! step between 2 blocks 
                           MPI_DOUBLE_PRECISION,CORNER,IERR) ! name of MPI_TYPE
   CALL MPI_TYPE_COMMIT(CORNER,IERR)


  end subroutine CREATE_MPI_VECTOR

subroutine FLUIDCOMM(VAR)

    implicit none

    real (kind=8), dimension(ISTART(1):IEND(1),ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL,&
                             ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL) :: VAR


    ! FJ => BJ
    CALL MPI_SENDRECV(VAR(ISTART(1),IEND(2)-NGHTCELL+1,ISTART(3)),1,XZ_2PLANES,NEIGHBOR(FJ),&
                      101,VAR(ISTART(1),ISTART(2)-NGHTCELL,ISTART(3)),1,XZ_2PLANES,NEIGHBOR(BJ),&
                      101,MPI_COMM_WORLD,STATUT,IERR)

    ! BJ => FJ
    CALL MPI_SENDRECV(VAR(ISTART(1),ISTART(2),ISTART(3)),1,XZ_2PLANES,NEIGHBOR(BJ),&
                      102,VAR(ISTART(1),IEND(2)+1,ISTART(3)),1,XZ_2PLANES,NEIGHBOR(FJ),&
                      102,MPI_COMM_WORLD,STATUT,IERR)

    ! FK => BK 
    CALL MPI_SENDRECV(VAR(ISTART(1),ISTART(2),IEND(3)-NGHTCELL+1),1,XY_2PLANES,NEIGHBOR(FK),&
                      104,VAR(ISTART(1),ISTART(2),ISTART(3)-NGHTCELL),1,XY_2PLANES,NEIGHBOR(BK),&
                      104,MPI_COMM_WORLD,STATUT,IERR)

    ! BK => FK
    CALL MPI_SENDRECV(VAR(ISTART(1),ISTART(2),ISTART(3)),1,XY_2PLANES,NEIGHBOR(BK),&
                      103,VAR(ISTART(1),ISTART(2),IEND(3)+1),1,XY_2PLANES,NEIGHBOR(FK),&
                      103,MPI_COMM_WORLD,STATUT,IERR)

    ! FJFK => BJBK
    CALL MPI_SENDRECV(VAR(ISTART(1),IEND(2)-NGHTCELL+1,IEND(3)-NGHTCELL+1),1,CORNER,NEIGHBOR(FJFK),&
                      105,VAR(ISTART(1),ISTART(2)-NGHTCELL,ISTART(3)-NGHTCELL),1,CORNER,NEIGHBOR(BJBK),&
                      105,MPI_COMM_WORLD,STATUT,IERR)

    ! BJBK => FJFK
    CALL MPI_SENDRECV(VAR(ISTART(1),ISTART(2),ISTART(3)),1,CORNER,NEIGHBOR(BJBK),&
                      106,VAR(ISTART(1),IEND(2)+1,IEND(3)+1),1,CORNER,NEIGHBOR(FJFK),&
                      106,MPI_COMM_WORLD,STATUT,IERR)
   
    ! BJFK => FJBK
    CALL MPI_SENDRECV(VAR(ISTART(1),ISTART(2),IEND(3)-NGHTCELL+1),1,CORNER,NEIGHBOR(BJFK),&
                      107,VAR(ISTART(1),IEND(2)+1,ISTART(3)-NGHTCELL),1,CORNER,NEIGHBOR(FJBK),&
                      107,MPI_COMM_WORLD,STATUT,IERR)

    ! FJBK => BJFK
    CALL MPI_SENDRECV(VAR(ISTART(1),IEND(2)-NGHTCELL+1,ISTART(3)),1,CORNER,NEIGHBOR(FJBK),&
                      108,VAR(ISTART(1),ISTART(2)-NGHTCELL,IEND(3)+1),1,CORNER,NEIGHBOR(BJFK),&
                      108,MPI_COMM_WORLD,STATUT,IERR)

  end subroutine FLUIDCOMM
  

subroutine SAVE_MPIIO(VAR,FILENAME)
  

  real (kind=8),dimension(ISTART(1)         :IEND(1)             &
                         ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                         ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL) :: VAR
  character(len=40) :: FILENAME

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  integer :: FILETYPE

  ! Temporary array
  real (kind=8), dimension(ISTART(1):IEND(1)  &
                          ,ISTART(2):IEND(2)  &
                          ,ISTART(3):IEND(3)) :: VAR_TEMP
  
  ! Dimension of array and subarray
  integer, dimension(3) :: VARCOUNT,VARSTART,DIMSUIDS


  ! File opening with writing permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )

  ! Total dimensions of the array to write
  DIMSUIDS(1) = NX ! NX
  DIMSUIDS(2) = NY ! NY
  DIMSUIDS(3) = NZ ! NZ

  ! Position in the file (in bytes)
  POS_FILE = 0

  ! Write the box dimensions (nx,ny,nz)
  if(MYID == 0) then
    call MPI_File_write(DESCRIPTEUR, DIMSUIDS, 3, MPI_INTEGER, MPI_STATUS_IGNORE, IERR)
  end if
 
  ! Position offset after 3 integers (3 * 4 bytes)
  POS_FILE = sizeof(DIMSUIDS) 

  ! Starting index for each proc -1 to start from 0 (MPI-IO C-language)
  VARSTART(1) = 0
  VARSTART(2) = ISTART(2)-1
  VARSTART(3) = ISTART(3)-1

  ! Number of cells in each direction for each proc 
  VARCOUNT(1) = IEND(1) - ISTART(1) + 1
  VARCOUNT(2) = IEND(2) - ISTART(2) + 1
  VARCOUNT(3) = IEND(3) - ISTART(3) + 1

  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(3, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_REAL8, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_REAL8, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)

  ! Temporary copy of VAR to avoid ghostcells treatment 
  VAR_TEMP(:,:,:) = VAR(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))

  ! Collective data writing 
  call MPI_File_write_all(DESCRIPTEUR, VAR_TEMP, product(VARCOUNT), MPI_REAL8, MPI_STATUS_IGNORE, IERR)

  ! Close file
  call MPI_File_close(DESCRIPTEUR, IERR)

  ! Derivated type suppression 
  call MPI_Type_free(FILETYPE, IERR)

end subroutine SAVE_MPIIO

subroutine SAVE_MPIIO_RHS(VAR,FILENAME)
  
  double complex,dimension(FSTART(1):FEND(1)    &
                          ,FSTART(2):FEND(2)    &
                          ,FSTART(3):FEND(3),3) :: VAR
  character(len=40) :: FILENAME

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  integer :: FILETYPE

  ! Dimension of array and subarray
  integer, dimension(4) :: VARCOUNT,VARSTART,DIMSUIDS
  integer, dimension(3) :: TNTAB


  ! File opening with writing permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )

  ! Total dimensions of the array to write
  DIMSUIDS(1) = (IEND(1) - ISTART(1) + 1) / 2 + 1 ! NX / 2 + 1
  DIMSUIDS(2) = IEND(1) - ISTART(1) + 1 ! NX
  DIMSUIDS(3) = IEND(1) - ISTART(1) + 1 ! NX
  DIMSUIDS(4) = 3 
  TNTAB(1) = TN ; TNTAB(2) = TNM1 ; TNTAB(3) = TNM2
  
  ! Position in the file (in bytes)
  POS_FILE = 0

  ! Write the box dimensions (nx,ny,nz)
  if(MYID == 0) then
    call MPI_File_write(DESCRIPTEUR, DIMSUIDS, 4, MPI_INTEGER, MPI_STATUS_IGNORE, IERR)
    call MPI_File_write(DESCRIPTEUR, TNTAB, 3,MPI_INTEGER, MPI_STATUS_IGNORE, IERR) 
  end if
 
  ! Position offset after 6 integers (3 * 4 bytes)
  POS_FILE = sizeof(DIMSUIDS)+sizeof(TNTAB) 

  ! Starting index for each proc -1 to start from 0 (MPI-IO C-language)
  VARSTART(1) = FSTART(1)-1
  VARSTART(2) = FSTART(2)-1
  VARSTART(3) = FSTART(3)-1
  VARSTART(4) = 0

  ! Number of cells in each direction for each proc 
  VARCOUNT(1) = FEND(1) - FSTART(1) + 1
  VARCOUNT(2) = FEND(2) - FSTART(2) + 1
  VARCOUNT(3) = FEND(3) - FSTART(3) + 1
  VARCOUNT(4) = 3 
 
  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(4, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_COMPLEX16, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_COMPLEX16, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)

  ! Collective data writing 
  call MPI_File_write_all(DESCRIPTEUR, VAR, product(VARCOUNT), MPI_COMPLEX16, MPI_STATUS_IGNORE, IERR)

  ! Close file
  call MPI_File_close(DESCRIPTEUR, IERR)

  ! Derivated type suppression 
  call MPI_Type_free(FILETYPE, IERR)

end subroutine SAVE_MPIIO_RHS

subroutine READ_MPIIO(VAR,FILENAME)


  real (kind=8),dimension(ISTART(1)         :IEND(1)             &
                         ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
                         ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL) :: VAR
  character(len=40) :: FILENAME

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  real (kind=8), dimension(ISTART(1):IEND(1)  &
                          ,ISTART(2):IEND(2)  &
                          ,ISTART(3):IEND(3)) :: VAR_TEMP

  integer :: FILETYPE
  integer, dimension(3) :: VARCOUNT,VARSTART,DIMSUIDS


  ! File opening with reading permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
                           MPI_MODE_RDONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )


  ! Position in the file (in bytes)
  POS_FILE = 0

  ! Write the box dimensions
  call MPI_File_read(DESCRIPTEUR, DIMSUIDS, 3, MPI_INTEGER, MPI_STATUS_IGNORE, IERR)

  POS_FILE = sizeof(DIMSUIDS) !offset because we just wrote 3 integers
  VARSTART(1) = 0
  VARSTART(2) = ISTART(2)-1
  VARSTART(3) = ISTART(3)-1

  VARCOUNT(1) = IEND(1) - ISTART(1) + 1
  VARCOUNT(2) = IEND(2) - ISTART(2) + 1
  VARCOUNT(3) = IEND(3) - ISTART(3) + 1
 
  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(3, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_REAL8, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_REAL8, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)


  ! Collective data reading 
  call MPI_File_read_all(DESCRIPTEUR, VAR_TEMP, product(VARCOUNT), MPI_REAL8, MPI_STATUS_IGNORE, IERR)

  VAR(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) =  VAR_TEMP(:,:,:)

  call MPI_File_close(DESCRIPTEUR, IERR)

  call MPI_Type_free(FILETYPE, IERR)

end subroutine READ_MPIIO
!***************
subroutine READ_MPIIO_RHS(VAR,FILENAME)
  
  double complex, dimension(FSTART(1):FEND(1)    &
                           ,FSTART(2):FEND(2)    &
                           ,FSTART(3):FEND(3),3) :: VAR
  character(len=40) :: FILENAME

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  integer :: FILETYPE

  ! Dimension of array and subarray
  integer, dimension(4) :: VARCOUNT,VARSTART,DIMSUIDS
  integer, dimension(3) :: TNTAB


  ! File opening with writing permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
                           MPI_MODE_RDONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )

  ! Total dimensions of the array to write
  DIMSUIDS(1) = (IEND(1) - ISTART(1) + 1) / 2 + 1 ! NX / 2 + 1
  DIMSUIDS(2) = IEND(1) - ISTART(1) + 1 ! NX 
  DIMSUIDS(3) = IEND(1) - ISTART(1) + 1 ! NX
  DIMSUIDS(4) = 3 
  
  TNTAB(1) = TN ; TNTAB(2) = TNM1 ; TNTAB(3) = TNM2
  
  ! Position in the file (in bytes)
  POS_FILE = 0

  ! Write the box dimensions (nx,ny,nz)
  call MPI_File_read(DESCRIPTEUR, DIMSUIDS, 4, MPI_INTEGER, MPI_STATUS_IGNORE, IERR)
  call MPI_File_read(DESCRIPTEUR, TNTAB, 3,MPI_INTEGER, MPI_STATUS_IGNORE, IERR) 
 
  ! Position offset after 6 integers (3 * 4 bytes)
  POS_FILE = sizeof(DIMSUIDS)+sizeof(TNTAB) 

  ! Starting index for each proc -1 to start from 0 (MPI-IO C-language)
  VARSTART(1) = FSTART(1)-1
  VARSTART(2) = FSTART(2)-1
  VARSTART(3) = FSTART(3)-1
  VARSTART(4) = 0

  ! Number of cells in each direction for each proc 
  VARCOUNT(1) = FEND(1) - FSTART(1) + 1
  VARCOUNT(2) = FEND(2) - FSTART(2) + 1
  VARCOUNT(3) = FEND(3) - FSTART(3) + 1
  VARCOUNT(4) = 3 
 
  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(4, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_COMPLEX16, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_COMPLEX16, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)

  ! Collective data writing 
  call MPI_File_read_all(DESCRIPTEUR, VAR, product(VARCOUNT), MPI_COMPLEX16, MPI_STATUS_IGNORE, IERR)

  ! Close file
  call MPI_File_close(DESCRIPTEUR, IERR)

  ! Derivated type suppression 
  call MPI_Type_free(FILETYPE, IERR)

end subroutine READ_MPIIO_RHS

subroutine FCM_SAVE_VAR_MPIIO(DIM1,DIM2,DIM3,VAR,FILENAME)
  
  integer,                                 intent(in) :: DIM1,DIM2,DIM3
  real (kind=8), dimension(DIM1,DIM2,DIM3),intent(in) :: VAR
  character(len=40),                       intent(in) :: FILENAME
  
  real (kind=8), dimension(:,:,:), allocatable :: TEMP_VAR
  

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE, POS_OFFSET
  integer :: FILETYPE
  
  integer :: NOI, NOD, NOL
  integer :: DIM1_START, DIM1_END
  integer :: LOAD_SIZE, LOAD_LOC_1, RES_LOAD
  integer :: K
  
  ! File opening with writing permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
           MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )
                                    


  ! Position in the file (in bytes)
  POS_FILE = 0
  
  
  ! Write the dimensions 
  if(MYID == 0) then
    call MPI_File_write(DESCRIPTEUR, DIM1, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
    call MPI_File_write(DESCRIPTEUR, DIM2, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
    call MPI_File_write(DESCRIPTEUR, DIM3, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
  end if

  call mpi_type_size(MPI_INTEGER,NOI,IERR)
  call mpi_type_size(MPI_LOGICAL,NOL,IERR)
  call mpi_type_size(MPI_DOUBLE_PRECISION,NOD,IERR)

  POS_OFFSET = 3 * NOI 
  
  
  ! Load balancing between procs
  LOAD_LOC_1 = int(DIM1/NPROC)  
  RES_LOAD = mod(DIM1,NPROC)   
  DIM1_START = (MYID)*LOAD_LOC_1 + min(MYID,RES_LOAD) +1

  
   if ( (MYID+1).le.RES_LOAD ) then
   LOAD_LOC_1 = LOAD_LOC_1 + 1
  end if
  
  
  allocate(TEMP_VAR(LOAD_LOC_1, DIM2, DIM3))
  
  ! Define DIM1 ending boudary
  DIM1_END = DIM1_START + LOAD_LOC_1 - 1
  
  if ((DIM1_END >  DIM1) .or. ( MYID ==  NPROC-1)) DIM1_END=DIM1
  
  LOAD_SIZE = LOAD_LOC_1 * DIM2 * DIM3
 
!~   print*, MYID, ' DIM1_START = ',DIM1_START
!~   print*, MYID, ' DIM1_END = ',DIM1_END

  POS_FILE = POS_OFFSET + ((DIM1_START-1)*DIM2*DIM3) * ( NOD )
!~   print*, MYID, ' LOAD_LOC_1 = ', LOAD_LOC_1
!~   print*, MYID, ' VAR(DIM1_START:DIM1_END,DIM2,DIM3) =  ', VAR(DIM1_START:DIM1_END,DIM2,DIM3)
  
  TEMP_VAR = VAR(DIM1_START:DIM1_END,1:DIM2,1:DIM3)  
  
    
!~      print*, MYID, ' Start writing'
  do K=1,LOAD_LOC_1  
   call mpi_file_write_at( DESCRIPTEUR,POS_FILE,TEMP_VAR(K,1:DIM2,1:DIM3),&
                         DIM2*DIM3,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR )
   POS_FILE = POS_FILE + (DIM2*DIM3 ) * ( NOD )   
  end do
!~    print*, MYID, ' End writing'

  call MPI_File_close(DESCRIPTEUR, IERR)
  
  deallocate(TEMP_VAR)
 
end subroutine FCM_SAVE_VAR_MPIIO





subroutine FCM_READ_VAR_MPIIO(DIM1,DIM2,DIM3,VAR,FILENAME,ERR_FILE)
  
  implicit none
  
  integer,                                 intent(in) :: DIM1,DIM2,DIM3
  integer,                                 intent(out) :: ERR_FILE
  real (kind=8), dimension(DIM1,DIM2,DIM3),intent(out) :: VAR
  character(len=40),                       intent(in) :: FILENAME
  
  real (kind=8), dimension(DIM1,DIM2,DIM3) :: TEMP_VAR
  

  !- File descriptor
  integer :: DESCRIPTEUR, IERR, ERR_CODE
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE, POS_OFFSET
  integer :: FILETYPE
  
  integer :: NOI, NOD, NOL
  integer :: K
  
  integer :: DIM1_READ, DIM2_READ, DIM3_READ
  
  

  ! File opening with reading permission
  call MPI_FILE_OPEN(       MPI_COMM_WORLD, &
                            trim(FILENAME), &
                           MPI_MODE_RDONLY, &
                             MPI_INFO_NULL, &
                               DESCRIPTEUR, &
                                    IERR  )
                                    
                                  
      
  ERR_FILE = IERR 
                                 
  if (IERR.ne.0) then
   print*,'"',trim(FILENAME),'"', 'DOES NOT EXIST--> SKIP FILE'  
   return                           
  end if                              
  

  ! Position in the file (in bytes)
  POS_FILE = 0

  
  ! Write the dimensions 
  if (MYID==0) then
   call MPI_File_read(DESCRIPTEUR, DIM1_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
   call MPI_File_read(DESCRIPTEUR, DIM2_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
   call MPI_File_read(DESCRIPTEUR, DIM3_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)

  
   if ( (DIM1.ne.DIM1_READ).or.(DIM2.ne.DIM2_READ).or.(DIM3.ne.DIM3_READ) ) then
    print*, 'PROBLEM WITH ', '"',trim(FILENAME),'"'
    print*, 'ERROR : WRITTEN DIMENSIONS ARE DIFFERENT FROM INIDICATED IN THE MAIN PROG!!'
    print*, 'EXIT PROGRAM'
    ERR_CODE = -1
    call mpi_abort(MPI_COMM_WORLD,ERR_CODE,IERR)
   end if
  
  end if


  call mpi_type_size(MPI_INTEGER,NOI,IERR)
  call mpi_type_size(MPI_LOGICAL,NOL,IERR)
  call mpi_type_size(MPI_DOUBLE_PRECISION,NOD,IERR)

  POS_OFFSET = 3 * NOI 
  

  POS_FILE = POS_OFFSET 
  
  if (MYID==0) then
   do K=1,DIM1 
    call mpi_file_read_at( DESCRIPTEUR,POS_FILE,TEMP_VAR(K,1:DIM2,1:DIM3),&
                          DIM2*DIM3,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR )
    POS_FILE = POS_FILE + (DIM2*DIM3 ) * ( NOD )   
   end do
  end if
  
  VAR = TEMP_VAR

  call MPI_File_close(DESCRIPTEUR, IERR)
  
  
  
 
end subroutine FCM_READ_VAR_MPIIO


end module MPI_structures

