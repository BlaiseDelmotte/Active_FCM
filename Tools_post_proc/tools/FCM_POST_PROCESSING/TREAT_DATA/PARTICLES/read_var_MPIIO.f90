subroutine READ_VAR_MPIIO(DIM1,DIM2,DIM3,VAR,FILENAME,ERR_FILE)

  use MPI
  
  
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
  
  call MPI_File_read(DESCRIPTEUR, DIM1_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
  call MPI_File_read(DESCRIPTEUR, DIM2_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)
  call MPI_File_read(DESCRIPTEUR, DIM3_READ, 1, MPI_INTEGER,MPI_STATUS_IGNORE,IERR)

  
  if ( (DIM1.ne.DIM1_READ).or.(DIM2.ne.DIM2_READ).or.(DIM3.ne.DIM3_READ) ) then
   print*, 'PROBLEM WITH ', '"',trim(FILENAME),'"'
   print*, 'ERROR : WRITTEN DIMENSIONS ARE DIFFERENT FROM INIDICATED IN THE MAIN PROG!!'
   print*, 'DIM1 = ', DIM1
   print*, 'DIM1_READ = ', DIM1_READ
   print*, 'DIM2 = ', DIM2
   print*, 'DIM2_READ = ', DIM2_READ
   print*, 'DIM3 = ', DIM3
   print*, 'DIM3_READ = ', DIM3_READ
   
!~    print*, 'EXIT PROGRAM'
!~    ERR_CODE = -1
!~    call mpi_abort(MPI_COMM_WORLD,ERR_CODE,IERR)
  end if


  call mpi_type_size(MPI_INTEGER,NOI,IERR)
  call mpi_type_size(MPI_LOGICAL,NOL,IERR)
  call mpi_type_size(MPI_DOUBLE_PRECISION,NOD,IERR)

  POS_OFFSET = 3 * NOI 
  

  POS_FILE = POS_OFFSET 
  
  do K=1,DIM1 
   call mpi_file_read_at( DESCRIPTEUR,POS_FILE,TEMP_VAR(K,1:DIM2,1:DIM3),&
                         DIM2*DIM3,MPI_DOUBLE_PRECISION,MPI_STATUS_IGNORE,IERR )
   POS_FILE = POS_FILE + (DIM2*DIM3 ) * ( NOD )   
  end do
  
  VAR = TEMP_VAR

  call MPI_File_close(DESCRIPTEUR, IERR)
  
  if ( (DIM1.ne.DIM1_READ).or.(DIM2.ne.DIM2_READ).or.(DIM3.ne.DIM3_READ) ) then
   print*,'VAR = ', VAR(1:10,1,1)
   read(*,*)
  end if
 
end subroutine READ_VAR_MPIIO
