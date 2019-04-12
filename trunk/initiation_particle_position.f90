!!====================================================================
!!
!!          Particle positions initiation
!!
!!====================================================================

subroutine INITIATION_PARTICLE_POSITION

!!====================================================================
!! Here the particle position and velocity are initiated according
!! to the variable INIT_PART_POSITION specified by user in the
!! parameter file: 'param.in'.
!! Note that the fluid velocity at the particle position is computed
!! at the end of the subroutine because it is needed by time-advancing
!! numerical scheme.
!!====================================================================
!! Particle position initiation: 
!!------------------------------
!! INIT_PART_POSITION=
!!
!!  0: The particle positions are spefied in the Fortran file.
!!     The default position correspond of one particle in a cell of
!!     fluid mesh with a small shift from the cell center. This position
!!     distribution allows to check the accuracy of the interpolation
!!     scheme of fluid velocity at the particle position.
!!
!!  1: Random particle distribution in the box. Take care that a
!!     factor 0.99 im used in order to ensure that the particles are
!!     inside the box.
!!
!!  2: Particle positions are read from the binary file 'POSPART.ini'
!!
!!  3: Particle injection at one edge. This is like an injector.
!!
!!--------------------
!!
!! Warning: The particle number in each CPU is defined with the
!!          particle position
!!
!!--------------------
!!
!! Warning: Particle velocities and positions are read simultaneously
!!
!!====================================================================

use DNS_DIM

use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use FLUID_VARIABLE
use MPI_STRUCTURES


implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!---------------------------------------------------------------
!- Read data variables
!---------------------------------------------------------------
!- File name 
character(len=40) :: FILENAME

integer :: IDUMMY

integer :: RECSIZE, SIZE_INT, SIZE_REAL, NVARIABLE

!- 
integer :: NP_READ, ND_READ

!---------------------------------------------------------------
!- Variable for random and uniform initiation
!---------------------------------------------------------------
real(kind=8) :: X0, Y0, Z0
real(kind=8) :: U0, V0, W0
real(kind=8) :: XMAX, YMAX, ZMAX

!- 
real(kind=8) :: XRAND, YRAND, ZRAND, GAUSS
integer :: ID

real(kind=8) :: DD

!---------------------------------------------------------------
!- Injection parameter
!---------------------------------------------------------------
!- Position of injection center
real(kind=8) :: XINJ, YINJ 

!- Radius and width of injection
real(kind=8) :: RINJ, DRINJ

!- Angle and velocity of injection
real(kind=8) :: PHINJ, UINJ

!- 
real(kind=8) :: RADIUS, THETA


!- Shift
integer :: SHIFT

!- Index
integer :: I, J, K, L, M, N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!- Number of particle in each CPU
do J = 1, NIG
 NPART_LOC(J) = int(NPART_FULL/NPROC)
end do



!!====================================================================
!! 1. Specific position given by user
!!====================================================================
if(INIT_PART_POSITION == 0 .or. INIT_PART_POSITION == 1) then

!- Shift 
! SHIFT = 1: 1 cell  between two particles
!       = 2: 2 cells between two particles 
SHIFT = 2

X0 = XMESH(ISTART(1)) + DX*0.5D0
Y0 = YMESH(ISTART(2)) + DY*0.5D0
Z0 = ZMESH(ISTART(3)) + DZ*0.5D0

!- positions are shifted
! X0 = XMESH(ISTART(1)) + DX/6.1D0
! Y0 = YMESH(ISTART(2)) - DY/3.5D0
! Z0 = ZMESH(ISTART(3)) + DZ/2.85D0

L = 0
M = 0
N = 0

do J = 1, NIG
 do I = 1, NPART_LOC(J)

  PART(I,J)%XP = X0 + L*DX
  PART(I,J)%YP = Y0 + M*DY
  PART(I,J)%ZP = Z0 + N*DZ

  L = L + SHIFT

  if((PART(I,J)%XP + DX) > XMESH(IEND(1))) then
   L = 0
   M = M + SHIFT

   if ((PART(I,J)%YP + DY) > YMESH(IEND(2))) then
    M = 0
    N = N + SHIFT

    if ((PART(I,J)%ZP + DZ) > ZMESH(IEND(3))) then
      N = 0
    end if
   end if
  end if

  end do
end do



if(MYID==0) write(*,*)'Particle position initiation --> Uniform'



!!====================================================================
!! 2. Random position
!!====================================================================
elseif(INIT_PART_POSITION == 2) then
 

do J = 1, NIG
 do I = 1, NPART_LOC(J)

  call random_number(XRAND)
  PART(I,J)%XP = XMESH(ISTART(1)) &
	       + XRAND*(XMESH(IEND(1))+DX-XMESH(ISTART(1)))

  call random_number(XRAND)
  PART(I,J)%YP = YMESH(ISTART(2)) &
	       + XRAND*(YMESH(IEND(2))+DY-YMESH(ISTART(2)))

  call random_number(XRAND)
  PART(I,J)%ZP = ZMESH(ISTART(3)) &
	       + XRAND*(ZMESH(IEND(3))+DZ-ZMESH(ISTART(3)))

  PART(I,J)%COLOR = MYID
  PART(I,J)%COLOR = 1.0

 end do
end do



if(MYID==0) write(*,*)'Particle position initiation: Uniformly randomized --> OK'


!!====================================================================
!! 3. Stored particle position 
!!====================================================================
elseif(INIT_PART_POSITION == 3) then


!!--------------------------------------------------------------------
!! 3.1 Multiple binary files
!!--------------------------------------------------------------------
if(ISAVEPART == 1) then

!do J = 1, NIG
!
! !- Define file name
! write(FILENAME,10101)'PART_c',J,trim(FILE_EXT),'.ini'
!
! !- Open file containing the last particle position and velocity
! open(unit = 150,          &
!       file = trim(FILENAME),&
!       status='old',        &
!       form='unformatted')
!
! read(150)NP_READ
!
! if(NP_READ>NPMAX_LOC) then
!  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!  write(*,*)'!!                     ERROR                    !!'
!  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!  write(*,*)'!!'
!  write(*,*)'!! file     : init_part.f90'
!  write(*,*)'!!'
!  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!  stop
! end if  
!
!!write(*,*)'ID:',MYID,' Class:',J,'  NPREAD=',NP_READ
!
! NPART_LOC(J) = NP_READ
!
! !- Read block data 
! read(150)(PART(I,J)%XP,I=1,NP_READ)
! read(150)(PART(I,J)%YP,I=1,NP_READ)
! read(150)(PART(I,J)%ZP,I=1,NP_READ)
!
! read(150)(PART(I,J)%UP,I=1,NP_READ)
! read(150)(PART(I,J)%VP,I=1,NP_READ)
! read(150)(PART(I,J)%WP,I=1,NP_READ)
! 
! if(SOLVE_SCALAR) read(150)(PART(I,J)%TP,I=1,NP_READ)
!
! !- Close file
! close(150)
!
!end do

!- Define file name
write(FILENAME,10303)'PART',trim(FILE_EXT),'.ini'

!- Open file containing the last particle position and velocity
open(unit = 150,            &
       file = trim(FILENAME), &
       status ='old',         &
       form ='unformatted')

do J = 1, NIG

 read(150)NP_READ

 if(NP_READ>NPMAX_LOC) then
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!                     ERROR                    !!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!'
  write(*,*)'!! file     : initiation_particle_position.f90'
  write(*,*)'!!'
  write(*,*)'!!   NP_READ=',NP_READ
  write(*,*)'!! NPMAX_LOC=', NPMAX_LOC
  write(*,*)'!!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  stop
 end if  

!!write(*,*)'ID:',MYID,' Class:',J,'  NPREAD=',NP_READ

 NPART_LOC(J) = NP_READ

 !- Read block data 
 read(150)(PART(I,J)%XP,I=1,NP_READ)
 read(150)(PART(I,J)%YP,I=1,NP_READ)
 read(150)(PART(I,J)%ZP,I=1,NP_READ)

 read(150)(PART(I,J)%UP,I=1,NP_READ)
 read(150)(PART(I,J)%VP,I=1,NP_READ)
 read(150)(PART(I,J)%WP,I=1,NP_READ)
 
 read(150)(PART(I,J)%UFAP,I=1,NP_READ)
 read(150)(PART(I,J)%VFAP,I=1,NP_READ)
 read(150)(PART(I,J)%WFAP,I=1,NP_READ)
 
 if(SOLVE_SCALAR) then
   read(150)(PART(I,J)%TP,I=1,NP_READ)
   read(150)(PART(I,J)%TFAP,I=1,NP_READ)
 end if

end do

!- Close file
close(150)



!- Define file name
write(FILENAME,10303)'rhs_part',trim(FILE_EXT),'.ini'

!- Open file containing the last particle position and velocity
open(unit = 150,            &
       file = trim(FILENAME), &
       status ='old',         &
       form ='unformatted')

do J = 1, NIG

 read(150)NP_READ
 read(150)(PART(I,J)%UP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%UP_NM2,I=1,NP_READ)
 read(150)(PART(I,J)%VP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%VP_NM2,I=1,NP_READ)
 read(150)(PART(I,J)%WP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%WP_NM2,I=1,NP_READ)
 
 read(150)(PART(I,J)%UFAP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%UFAP_NM2,I=1,NP_READ)
 read(150)(PART(I,J)%VFAP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%VFAP_NM2,I=1,NP_READ)
 read(150)(PART(I,J)%WFAP_NM1,I=1,NP_READ)
 read(150)(PART(I,J)%WFAP_NM2,I=1,NP_READ)

 if(SOLVE_SCALAR) then
   read(150)(PART(I,J)%TP_NM1,I=1,NP_READ)
   read(150)(PART(I,J)%TP_NM2,I=1,NP_READ)
   read(150)(PART(I,J)%TFAP_NM1,I=1,NP_READ)
   read(150)(PART(I,J)%TFAP_NM2,I=1,NP_READ)
 end if


end do

!- Close file
close(150)















if(MYID==0) write(*,*)'Particle position initiation: Read from file --> OK'
if(MYID==0) write(*,*)'     + Multiple Binary Files'



!!--------------------------------------------------------------------
!! 3.2 Direct access file
!!--------------------------------------------------------------------
elseif(ISAVEPART == 2) then

!! 
call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,SIZE_REAL,IERR)
call MPI_TYPE_SIZE(MPI_INTEGER,SIZE_INT,IERR)


!!-------------------------------------
!! Read binary file
!!-------------------------------------
!!- Number of variables
!!  Basically 6 variables are read: xp, yp, zp, up, vp, wp
NVARIABLE = 6

!!- in case of add temperature
if(SOLVE_SCALAR) NVARIABLE = NVARIABLE + 1

 
!!- Compute the length of the record 
RECSIZE = (SIZE_INT +  NPMAX_LOC*SIZE_REAL*NVARIABLE)*NIG

!- Define file name
FILENAME = 'traj.ini'


open(unit=150,          &
      file=trim(FILENAME),&
      status='old',       &
      access='direct',    &
      action='read',      &
      recl=RECSIZE,       &
      form='unformatted') ! Had to change 'binary' to 'unformatted' to make it compile with mpifort 

if(SOLVE_SCALAR) then
read(unit=150,rec=1+MYID)  &
           (NPART_LOC(J), &
	    (PART(I,J)%XP, &
              PART(I,J)%YP, &
              PART(I,J)%ZP, &
              PART(I,J)%UP, &
              PART(I,J)%VP, &
              PART(I,J)%WP, &
              PART(I,J)%TP, &
	      I=1,NPMAX_LOC), J=1,NIG)

else
read(unit=150,rec=1+MYID)  &
           (NPART_LOC(J), &
	    (PART(I,J)%XP, &
              PART(I,J)%YP, &
              PART(I,J)%ZP, &
              PART(I,J)%UP, &
              PART(I,J)%VP, &
              PART(I,J)%WP,I=1,NPMAX_LOC), J=1,NIG)
end if



call MPI_BARRIER(MPI_COMM_WORLD,IERR)
close(150)






if(DEBUG) then
! Total particles number checking
do J=1, NIG
call ISUMCPU(NPART_LOC(J),IDUMMY)
if (MYID==0) write(*,*) 'Read Part -> Class:',J,' Full number of particles ',IDUMMY
end do
!do J=1, NIG
!write(650+MYID,*)MYID,'IG=',J,'  NPART_LOC=',NPART_LOC(J)
!end do
end if


if(MYID==0) write(*,*)'Particle position initiation: Read from file --> OK'
if(DEBUG) then
 if(MYID==0)write(*,*) 'Read Part -> NVARIABLE = ',NVARIABLE
 if(MYID==0)write(*,*) 'Read Part ->   RECSIZE = ',RECSIZE
 if(MYID==0)write(*,*) 'Read Part -> NPMAX_LOC = ',NPMAX_LOC
 if(MYID==0)write(*,*) 'Read Part -> SIZE_REAL = ',SIZE_REAL
 if(MYID==0)write(*,*) 'Read Part ->  SIZE_INT = ',SIZE_INT
end if
if(MYID==0)write(*,*) '    + Direct access file'
if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)



!!--------------------------------------------------------------------
!! 3.2 MPI I/O
!!--------------------------------------------------------------------
elseif(ISAVEPART == 3) then


FILENAME = 'traj.ini'


call READ_PART_MPIIO(PART,FILENAME)


if(MYID==0) write(*,*)'Particle position initiation: Read from file --> OK'
if(MYID==0)write(*,*) '    + MPI I/O'
if(MYID==0)write(*,10300) '     + File name = ',trim(FILENAME)



end if !!- If: ISAVEPART

!!- Switch particle velocity initiation flag
INIT_PART_VELOCITY = 3



!!====================================================================
!! 4. Injection at the edge (x,y)
!!====================================================================
elseif(INIT_PART_POSITION == 4) then


 if(MYID==0)write(*,*)'Particle position initiation: Injection --> OK'

 !- For this case the particle position and velocity initiation are
 !  performed simultaneously. So the velocity initiation variable
 !  is switched to 100 in order to avoid a another particle velocity
 !  initiation.
!! INIT_PART_VELOCITY = 100

 !- Position of injector
!! XINJ = XMESH(NX/2)
!! YINJ = YMESH(NY/2)

 !- Radius of injector
!! RINJ = 1.E-3 
!! DRINJ = 1.E-3

 !- Angle of injection
!! PHINJ = 2*PPI/180*20

 !- Injection velocity 
!! UINJ = 5. !- m/s

!! do IG = 1, NIG
!!  do NP = 1, NPMAX

   !- Random number for radial position
!!   call random_number(XRAND)

   !- Random number for angular position
!!   call random_number(YRAND)

   !- Random number for angular position
!!   call random_number(ZRAND)


!- Random radial distribution
!   RADIUS = RINJ + DRINJ*(1.-XRAND)
!!   RADIUS = RINJ*(1.-XRAND)

!- Random azimuthal position
!!   THETA = 2*PPI*YRAND

!- Random angula distribution 
!!   PHINJ = (1.-ZRAND)*PPI/6.


!- Particle position
!!   XPART(NP,IG) = XINJ + RADIUS*cos(THETA)
!!   YPART(NP,IG) = YINJ + RADIUS*sin(THETA)
!!   ZPART(NP,IG) = 0.

!- Particle velocity
!!   UPART(NP,IG) = UINJ*sin(PHINJ)*cos(THETA)
!!   VPART(NP,IG) = UINJ*sin(PHINJ)*sin(THETA)
!!   WPART(NP,IG) = UINJ*cos(PHINJ)
    

!!  end do
!! end do


end if


!!- Proc ID managing the particles
!PART(:,:)%PROC_ID = MYID




!!- Boundary conditions
call BOUNDARY_PARTICLE






if(MYID==0) write(*,*) 'Particle position initiation --> OK'
! Total particles number checking
do J=1, NIG
call ISUMCPU(NPART_LOC(J),IDUMMY)
if (MYID==0) write(*,10102) '   + C',J,' Nbr of Init. Part. = ',IDUMMY
end do



!!--------------------------------------------------------------------
10100 format (A,I2.2,A)
10101 format (A,I2.2,A,A)
10102 format (A,I2.2,A,I8)
10300 format (A,A)
10303 format (A,A,A)


end subroutine INITIATION_PARTICLE_POSITION
