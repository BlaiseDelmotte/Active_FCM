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
  DIMSUIDS(1:3) = IEND(1) - ISTART(1) + 1 ! NX

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
