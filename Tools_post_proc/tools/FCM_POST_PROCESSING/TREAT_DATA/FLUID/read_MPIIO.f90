subroutine READ_MPIIO(NX,NY,NZ,VAR,FILENAME)


  use MPI

  implicit none 
  
  integer, intent(in):: NX, NY, NZ  ! Var size
  
  real (kind=8),dimension(NX    &
                         ,NY    &
                         ,NZ ):: VAR
                         
  character(len=40), intent(in) :: FILENAME
  
   

  !- File descriptor
  integer :: DESCRIPTEUR
  integer(kind=MPI_OFFSET_KIND) :: POS_FILE
  real (kind=8), dimension(NX  &
                          ,NY  &
                          ,NZ) :: VAR_TEMP

  integer :: FILETYPE
  integer :: IERR
  
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
  VARSTART(2) = 0
  VARSTART(3) = 0

  VARCOUNT(1) = NX
  VARCOUNT(2) = NY
  VARCOUNT(3) = NZ
 
  ! Derivated type creation : a subarray different for each proc 
  call MPI_Type_create_subarray(3, DIMSUIDS, VARCOUNT, VARSTART, &
                                MPI_ORDER_FORTRAN, MPI_REAL8, FILETYPE, IERR)
  call MPI_Type_commit(FILETYPE, IERR)

  ! Set view in the file 
  call MPI_File_set_view(DESCRIPTEUR, POS_FILE, MPI_REAL8, &
                         FILETYPE, "native", MPI_INFO_NULL, IERR)


  ! Collective data reading 
  call MPI_File_read_all(DESCRIPTEUR, VAR_TEMP, product(VARCOUNT), MPI_REAL8, MPI_STATUS_IGNORE, IERR)

  VAR(1:NX,1:NY,1:NZ) =  VAR_TEMP(:,:,:)

  call MPI_File_close(DESCRIPTEUR, IERR)

  call MPI_Type_free(FILETYPE, IERR)

end subroutine READ_MPIIO
