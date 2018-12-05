!!====================================================================
!!
!!          Particle scalar intiation
!!
!!====================================================================

subroutine INITIATION_PARTICLE_SCALAR

!!====================================================================
!!
!!====================================================================

use DNS_DIM

use GEOMETRIC_VARIABLE
use PARAM_PHYS
use PARTICLE_PARALLEL
use FLUID_VARIABLE
use MPI_STRUCTURES
use SCALAR_VARIABLE


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


!- Index
integer :: I, J, K, L, M, N

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!====================================================================
!! Interpolation of scalar at particle position
!!====================================================================
do J = 1, NIG

 !- Interpolation
 call INTERPH(INTERP_SCHEME,                &
                      XMESH, YMESH, ZMESH,  &
                      THETA,                &
               NPART_LOC(J),                &
 	        PART(1:NPART_LOC(J),J)%XP,  &
		PART(1:NPART_LOC(J),J)%YP,  &
		PART(1:NPART_LOC(J),J)%ZP,  &
 	        PART(1:NPART_LOC(J),J)%TFAP )

end do !!- End loop:J = 1, NIG



!!====================================================================
!! 0. Specific position given by user
!!====================================================================
if(INIT_PART_SCALAR == 0) then


PART(:,:)%TP = ZERO

if(MYID==0) write(*,*)'Particle scalar initiation: Theta(t=0)=0 --> OK'


!!====================================================================
!! 1. Equal to fluid scalar 
!!====================================================================
elseif(INIT_PART_SCALAR == 1) then


do J = 1, NIG
 PART(:,J)%TP = PART(:,J)%TFAP
end do !!- End loop:J = 1, NIG



if(MYID==0) write(*,*)'Particle scalar initiation: Theta_p(t=0)=Theta_f@p(t=0) --> Ok'




!!====================================================================
!! 2. Random (from uniform distribution)
!!====================================================================
elseif(INIT_PART_SCALAR == 2) then


 write(*,*)'Particle scalar initiation: Random --> not done'



!!====================================================================
!! 3. Read from stored file 
!!====================================================================
elseif(INIT_PART_VELOCITY == 3) then


if(MYID==0) write(*,*)'Particle scalar initiation: Read from file --> OK'

!!- already done previously in INITIATION_PARTICLE_SCALAR



end if






if(MYID==0) write(*,*) 'Particle scalar initiation --> OK'




end subroutine INITIATION_PARTICLE_SCALAR
