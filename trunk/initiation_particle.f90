!!====================================================================
!!
!!          Particle position and velocity initiation
!!
!!====================================================================

subroutine INITIATION_PARTICLE

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
!!--------------------------------------------------------------------
!! Particle velocity initiation: 
!!------------------------------
!! INIT_PART_VELOCITY=
!!
!!  0: The particle velocity is equal to zero.
!!
!!  1: The particle velocity is randomly chosen in a Gaussian
!!     distribution.
!!
!!  2: The particle velocity is equal to the fluid velocity.
!!
!!  3: Particle velocities are read from the binary file 'VELPART.ini'
!!====================================================================

use dns_dim

use geometric_variable
use param_phys
use particle_parallel
use fluid_variable

implicit none


!!====================================================================
! ARRAYS STATEMENT
!---------------------------------------------------------------------
real(kind=8) :: X0, Y0, Z0
real(kind=8) :: U0, V0, W0
real(kind=8) :: XMAX, YMAX, ZMAX

!- 
real(kind=8) :: XRAND, YRAND, ZRAND, GAUSS
integer :: ID

real(kind=8) :: DD

integer :: IDUMMY

!- 
integer :: NP_READ, ND_READ

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

!- File name 
character(len=40) :: FILENAME

!- Shift
integer :: SHIFT

!- Index
integer :: I, J, K, L, M, N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!====================================================================
!! 1. Particle position
!!====================================================================

!!- Number of particle in each CPU
do J = 1, NIG
 NPART_LOC(J) = int(NPART_FULL/NPROC)
end do



!!-------------------------------------------------------------------
!! 1.0. Specific position given by user
!!-------------------------------------------------------------------
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

   if ((PART(I,J)%YP + DY) > XMESH(IEND(2))) then
    M = 0
    N = N + SHIFT

    if ((PART(I,J)%ZP + DZ) > XMESH(IEND(3))) then
      N = 0
    end if
   end if
  end if

  end do
end do


if(MYID==0) write(*,*)'Particle position initiation --> Uniform'



!!-------------------------------------------------------------------
!! 1.2. Random position
!!-------------------------------------------------------------------
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

 end do
end do



if(MYID==0) write(*,*)'Particle position initiation: Uniformly randomized --> OK'


!!-------------------------------------------------------------------
!! 1.3. Stored particle position 
!!-------------------------------------------------------------------
!!
elseif(INIT_PART_POSITION == 3) then


do J = 1, NIG

 !- Define file name
 write(FILENAME,10101)'PART_c',J,'_',trim(FILE_EXT),'.ini'

 !- Open file containing the last particle position and velocity
 open(unit = 150, file = trim(FILENAME),form='unformatted')

 read(150)NP_READ

 if(NP_READ>NPMAX_LOC) then
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!                     ERROR                    !!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  write(*,*)'!!'
  write(*,*)'!! file     : init_part.f90'
  write(*,*)'!!'
  write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  stop
 end if  

 NPART_LOC(J) = NP_READ

 !- Read block data 
 read(150)(PART(I,J)%XP,I=1,NP_READ)
 read(150)(PART(I,J)%YP,I=1,NP_READ)
 read(150)(PART(I,J)%ZP,I=1,NP_READ)

 read(150)(PART(I,J)%UP,I=1,NP_READ)
 read(150)(PART(I,J)%VP,I=1,NP_READ)
 read(150)(PART(I,J)%WP,I=1,NP_READ)

 !- Close file
 close(150)

end do


if(MYID==0) write(*,*)'Particle position initiation: Read from file --> OK'





!!-------------------------------------------------------------------
!! 1.3. Injection at the edge (x,y)
!!-------------------------------------------------------------------
!!
elseif(INIT_PART_POSITION == 3) then

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



!!====================================================================
!! 2. Particle velocity
!!====================================================================


!!-------------------------------------------------------------------
!! 2.0. Equal to zero Up = 0 
!!-------------------------------------------------------------------
if(INIT_PART_VELOCITY == 0) then

PART(:,:)%UP = ZERO
PART(:,:)%VP = ZERO
PART(:,:)%WP = ZERO

if(MYID==0) write(*,*)'Particle velocity initiation: Up(t=0)=0 --> OK'



!!-------------------------------------------------------------------
!! 2.1. Equal to fluid velocity 
!!-------------------------------------------------------------------
elseif(INIT_PART_VELOCITY == 1) then


do J = 1, NIG

 !- U-velocity Interpolation
 call INTERPH(INTERP_SCHEME,                              &
                      XMESH, YMESH, ZMESH,                &
                       UFLU,                              &
               NPART_LOC(J),                              &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%UFAP			          )

  !- V-velocity Interpolation
 call INTERPH(INTERP_SCHEME,				   &
		      XMESH, YMESH, ZMESH,		   &
		       VFLU,				   &
	       NPART_LOC(J),				   &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%VFAP 			   )


 !- W-velocity Interpolation
 call INTERPH(INTERP_SCHEME,				   &
		      XMESH, YMESH, ZMESH,		   &
		       WFLU,				   &
	       NPART_LOC(J),				   &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%WFAP 			   )


  PART(:,J)%UP = PART(:,J)%UFAP
  PART(:,J)%VP = PART(:,J)%VFAP
  PART(:,J)%WP = PART(:,J)%WFAP

 end do



if(MYID==0) write(*,*)'Particle velocity initiation: Up(t=0)=Uf@p(t=0) --> Bof'

 

!!-------------------------------------------------------------------
!! 2.1. Random (from Gaussian distribution 
!!-------------------------------------------------------------------
elseif(INIT_PART_VELOCITY == 2) then

 do J = 1, NIG
  do I = 1, NPART_LOC(J)

   call random_number(XRAND)
   PART(I,J)%XP = XRAND - 0.5

   call random_number(XRAND)
   PART(I,J)%YP = XRAND - 0.5

   call random_number(XRAND)
   PART(I,J)%ZP = XRAND - 0.5

  end do
 end do

 write(*,*)'Particle velocity initiation: Random --> OK'



!!-------------------------------------------------------------------
!! 2.2. Read from stored file 
!!-------------------------------------------------------------------
elseif(INIT_PART_VELOCITY == 3) then

if(MYID==0) write(*,*)'Particle velocity initiation: Read from file --> OK'

!!- already done previously



end if



!!===================================================================
!! 4. Interpolate fluid velocity at the particle position
!!===================================================================
do J = 1, NIG

 !- U-velocity Interpolation
 call INTERPH(INTERP_SCHEME,                              &
                      XMESH, YMESH, ZMESH,                &
                       UFLU,                              &
               NPART_LOC(J),                              &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%UFAP			          )

  !- V-velocity Interpolation
 call INTERPH(INTERP_SCHEME,				   &
		      XMESH, YMESH, ZMESH,		   &
		       VFLU,				   &
	       NPART_LOC(J),				   &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%VFAP 			   )


 !- W-velocity Interpolation
 call INTERPH(INTERP_SCHEME,				   &
		      XMESH, YMESH, ZMESH,		   &
		       WFLU,				   &
	       NPART_LOC(J),				   &
 	        PART(1:NPART_LOC(J),J)%XP, PART(1:NPART_LOC(J),J)%YP, PART(1:NPART_LOC(J),J)%ZP, &
 	        PART(1:NPART_LOC(J),J)%WFAP 			   )

end do

!!- Proc ID managing the particles
PART(:,:)%PROC_ID = MYID


! Total particles number checking
do J=1, NIG
call ISUMCPU(NPART_LOC(J),IDUMMY)
if (MYID==0) then
write(*,*) 'Class:',J,' Number of particle in each CPU ',NPART_LOC(J)
write(*,*) '            Full number of particles ',IDUMMY
end if
end do



do I=1, NPART_LOC(1)
 if(    (PART(I,1)%XP<ZERO ) &
    .or.(PART(I,1)%XP>LXMAX) &
    .or.(PART(I,1)%YP<ZERO ) &
    .or.(PART(I,1)%YP>LYMAX) &
    .or.(PART(I,1)%ZP<ZERO ) &
    .or.(PART(I,1)%ZP>LXMAX) ) then

PART(I,1)%PROC_ID = 10

write(*,*)'Oupss'
end if
end do



if(MYID==0) write(*,*) 'Particle initiation --> OK'

10101 format (A,I2.2,A,A,A)


end subroutine INITIATION_PARTICLE
