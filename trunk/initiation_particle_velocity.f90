!!====================================================================
!!
!!          Particle velocities initiation
!!
!!====================================================================

subroutine INITIATION_PARTICLE_VELOCITY

!!====================================================================
!! Here the particle position and velocity are initiated according
!! to the variable INIT_PART_VELOCITY specified by user in the
!! parameter file: 'param.in'.
!! Note that the fluid velocity at the particle position is computed
!! at the end of the subroutine because it is needed by time-advancing
!! numerical scheme.
!!====================================================================
!! Particle velocity initiation: 
!!------------------------------
!!
!! INIT_PART_VELOCITY=
!!
!!  0: The particle velocity is equal to zero.
!!
!!  1: The particle velocity is randomly chosen in a Gaussian
!!     distribution.
!!
!!  2: The particle velocity is equal to the fluid velocity.
!!
!!  3: Particle velocities are read from the binary file
!!      (done in subroutine INIT_PARTICLE_POSITION)
!!====================================================================

use DNS_DIM

use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use FLUID_VARIABLE

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

real(kind=8) :: XTMP, YTMP,ZTMP,UTMP,VTMP,WTMP,COLORTMP

!- File name 
character(len=40) :: FILENAME

!- Shift
integer :: SHIFT

!- Index
integer :: I, J, K, L, M, N, NP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!====================================================================
!! 1. Interpolate fluid velocity at particle position
!!====================================================================
if(SOLVE_FLUID>0) then 

do J = 1, NIG

 !- U-velocity Interpolation
 call INTERPH(INTERP_SCHEME,                &
                      XMESH, YMESH, ZMESH,  &
                       UFLU,                &
               NPART_LOC(J),                &
 	        PART(1:NPART_LOC(J),J)%XP,  &
		PART(1:NPART_LOC(J),J)%YP,  &
		PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%UFAP )

  !- V-velocity Interpolation
 call INTERPH(INTERP_SCHEME,                &
		      XMESH, YMESH, ZMESH,  &
		       VFLU,		    &
	       NPART_LOC(J),		    &
 	        PART(1:NPART_LOC(J),J)%XP,  &
		PART(1:NPART_LOC(J),J)%YP,  &
		PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%VFAP )


 !- W-velocity Interpolation
 call INTERPH(INTERP_SCHEME,                &
		      XMESH, YMESH, ZMESH,  &
		       WFLU,                &
	       NPART_LOC(J),                &
 	        PART(1:NPART_LOC(J),J)%XP,  &
		PART(1:NPART_LOC(J),J)%YP,  &
		PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%WFAP )
end do

end if !!- if(SOLVE_FLUID>0)

!!====================================================================
!! 2. Particle velocity equal to zero Up = 0 
!!====================================================================
if(INIT_PART_VELOCITY == 0) then

PART(:,:)%UP = ZERO
PART(:,:)%VP = ZERO
PART(:,:)%WP = ZERO




!!==================================================================
!! Part #1
!!==================================================================
!NP = 1
!XTMP = 1.D0/5.D0*LXMAX
!YTMP = YMESH(ISTART(2)) + (YMESH(IEND(2))+DY-YMESH(ISTART(2)))*0.5
!ZTMP = ZMESH(ISTART(3)) + (ZMESH(IEND(3))+DZ-ZMESH(ISTART(3)))*0.5
!if(ZTMP<LZMAX*0.5) then
! YTMP = YTMP+3e-3
! UTMP = ZERO
! VTMP = ZERO
! WTMP = 0.10D0
! COLORTMP = 10.0
!else
! UTMP = ZERO
! VTMP = ZERO
! WTMP = -0.1D0
! COLORTMP = -10.0
!end if
!PART(NP,:)%XP = XTMP
!PART(NP,:)%YP = YTMP
!PART(NP,:)%ZP = ZTMP
!PART(NP,:)%UP = UTMP
!PART(NP,:)%VP = VTMP
!PART(NP,:)%WP = WTMP
!
!PART(NP,:)%COLOR = COLORTMP


!!==================================================================
!! Part #2
!!==================================================================
!NP = NP + 1
!
!XTMP = 4.D0/5.D0*LXMAX
!YTMP = YMESH(ISTART(2)) + (YMESH(IEND(2))+DY-YMESH(ISTART(2)))*0.5
!ZTMP = ZMESH(ISTART(3)) + (ZMESH(IEND(3))+DZ-ZMESH(ISTART(3)))*0.5
!
!if(ZTMP<LZMAX*0.5) then
! UTMP = ZERO
! VTMP = ZERO
! WTMP = 0.1D0
! COLORTMP = 10.0
!else
! UTMP = ZERO
! VTMP = ZERO
! WTMP = -0.1D0
! COLORTMP = -10.0
!end if
!
!PART(NP,:)%XP = XTMP
!PART(NP,:)%YP = YTMP
!PART(NP,:)%ZP = ZTMP
!PART(NP,:)%UP = UTMP
!PART(NP,:)%VP = VTMP
!PART(NP,:)%WP = WTMP

!PART(NP,:)%COLOR = COLORTMP

!NP=1
!write(*,20000)MYID,NP,PART(NP,1)%XP,PART(NP,1)%YP,PART(NP,1)%ZP,PART(NP,1)%UP,PART(NP,1)%VP,PART(NP,1)%WP,PART(NP,1)%COLOR
!NP=2
!write(*,20000)MYID,NP,PART(NP,1)%XP,PART(NP,1)%YP,PART(NP,1)%ZP,PART(NP,1)%UP,PART(NP,1)%VP,PART(NP,1)%WP,PART(NP,1)%COLOR


if(MYID==0) write(*,*)'Particle velocity initiation --> Up(t=0)=0'



!!====================================================================
!! 3. Equal to fluid velocity 
!!====================================================================
elseif(INIT_PART_VELOCITY == 1) then


do J = 1, NIG

  PART(:,J)%UP = PART(:,J)%UFAP
  PART(:,J)%VP = PART(:,J)%VFAP
  PART(:,J)%WP = PART(:,J)%WFAP

end do



if(MYID==0) write(*,*)'Particle velocity initiation --> Up(t=0)=Uf@p(t=0)'

 

!!====================================================================
!! 4. Random (from uniform distribution)
!!====================================================================
elseif(INIT_PART_VELOCITY == 2) then

 do J = 1, NIG
  do I = 1, NPART_LOC(J)

   call random_number(XRAND)
   PART(I,J)%UP = XRAND - 0.5

   call random_number(XRAND)
   PART(I,J)%VP = XRAND - 0.5

   call random_number(XRAND)
   PART(I,J)%WP = XRAND - 0.5

  end do
 end do

if(MYID==0)  write(*,*)'Particle velocity initiation: Random --> OK'



!!====================================================================
!! 5. Read from stored file 
!!====================================================================
elseif(INIT_PART_VELOCITY == 3) then

if(MYID==0) write(*,*)'Particle velocity initiation: Read from file --> OK'

!!- already done previously in INITIATION_PARTICLE_POSITION


end if




if(MYID==0) write(*,*) 'Particle velocity initiation --> OK'



!!--------------------------------------------------------------------
10101 format (A,I2.2,A,A)
20000 format (I3,2x,I3,10(e17.7))


end subroutine INITIATION_PARTICLE_VELOCITY
