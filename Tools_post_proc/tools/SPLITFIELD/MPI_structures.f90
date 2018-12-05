module MPI_structures

use DNS_DIM 

integer :: XZ_2PLANES,XY_2PLANES,CORNER
INTEGER, PARAMETER		 :: NB_NEIGHBORS = 8
INTEGER, DIMENSION(NB_NEIGHBORS) :: NEIGHBOR
INTEGER, PARAMETER		 :: FJ=1,BJ=2,FK=3,BK=4 ! Neighbors : Forward_J, Backward_J, ... 
INTEGER, PARAMETER		 :: FJBK=5,BJFK=6,BJBK=7,FJFK=8 ! Common edge neighbor


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!     Print 3D-field with MPIIO
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SAVE_MPIIO(VAR,FILENAME)

  real (kind=8),dimension(ISTARTI(1)         :IENDI(1)             &
                         ,ISTARTI(2)-NGHTCELL:IENDI(2)+NGHTCELL    &
                         ,ISTARTI(3)-NGHTCELL:IENDI(3)+NGHTCELL) :: VAR
  character(len=40) :: FILENAME

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  integer :: FILETYPE

  ! Temporary array
  real (kind=8), dimension(ISTARTI(1):IENDI(1)  &
                          ,ISTARTI(2):IENDI(2)  &
                          ,ISTARTI(3):IENDI(3)) :: VAR_TEMP
  
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
  DIMSUIDS(1:3) = IENDI(1) - ISTARTI(1) + 1 ! NX

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
  VARSTART(2) = ISTARTI(2)-1
  VARSTART(3) = ISTARTI(3)-1

  ! Number of cells in each direction for each proc 
  VARCOUNT(1) = IENDI(1) - ISTARTI(1) + 1
  VARCOUNT(2) = IENDI(2) - ISTARTI(2) + 1
  VARCOUNT(3) = IENDI(3) - ISTARTI(3) + 1
 
  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(3, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_REAL8, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_REAL8, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)

  ! Temporary copy of VAR to avoid ghostcells treatment 
  VAR_TEMP(:,:,:) = VAR(ISTARTI(1):IENDI(1),ISTARTI(2):IENDI(2),ISTARTI(3):IENDI(3))

  ! Collective data writing 
  call MPI_File_write_all(DESCRIPTEUR, VAR_TEMP, product(VARCOUNT), MPI_REAL8, MPI_STATUS_IGNORE, IERR)

  ! Close file
  call MPI_File_close(DESCRIPTEUR, IERR)

  ! Derivated type suppression 
  call MPI_Type_free(FILETYPE, IERR)

end subroutine SAVE_MPIIO

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


!!write(*,*)'MYID=',MYID,' FILENAME=',FILENAME


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

end module MPI_structures
