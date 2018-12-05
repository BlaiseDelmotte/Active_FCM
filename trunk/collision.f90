!!====================================================================
!!
!!          Particle-particle interaction
!!
!!====================================================================

subroutine COLLISION(NCYCLE,TIME)

!!====================================================================
!!
!!
!!====================================================================

use STATISTICS 
use PARTICLE_PARALLEL
use MPI_STRUCTURES
use DNS_DIM               
use PARAM_PHYS
use CHECK_CPU
use COLLISION_VARIABLE

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!!- cycle number
integer, intent(in) :: NCYCLE

!- Curent time
real(kind=8), intent(in) :: TIME


!!====================================================================
integer, parameter :: NPCLOSEMAX = 150000
integer, dimension(NPCLOSEMAX) :: HOME
integer, dimension(NPCLOSEMAX) :: NEAR

integer :: NCLOSE

integer :: NP1, NP2
!!====================================================================



!!- Index
integer :: I, J ,K
integer :: NS, NF 

!!====================================================================



do J = 1, NIG

!!!! A EVALUER ET A DEPLACER AILLEURS !!!!!
!  DELTA_COLL = LXMAX/50.0 
  DELTA_COLL = 1.2*DPART(J)


!!====================================================================
!! 1. Gathering of particles located in neighboring CPU
!!====================================================================

if(NPROC>1) then

!! neighbor detection 
 call PARTICLES_COUNTING_COLL(J)

!! concatenation
 call EXCHANGE_P_COLL(J,NPNEIGH_COLL)


!! Y/Z periodicity shifts for neightbor particles
 NS = NPART_LOC(J)+1
 NF = NPART_LOC(J)+NPNEIGH_COLL


 do I = NS, NF


 if(YMESH(ISTART(2)) <= ZERO) then
  if (PART(I,J)%YP > LYMAX - DELTA_COLL .and. PART(I,J)%YP <= LYMAX) then
    PART(I,J)%YP = PART(I,J)%YP - LYMAX
  end if
 end if
 if(YMESH(IEND(2))+DY >= LYMAX) then
  if (PART(I,J)%YP >= ZERO .and. PART(I,J)%YP < DELTA_COLL) then
    PART(I,J)%YP = PART(I,J)%YP + LYMAX
  end if
 end if

 if(ZMESH(ISTART(3)) <= ZERO) then
  if (PART(I,J)%ZP > LZMAX - DELTA_COLL .and. PART(I,J)%ZP <= LZMAX) then
    PART(I,J)%ZP = PART(I,J)%ZP - LZMAX
  end if
 end if
 if(ZMESH(IEND(3))+DZ >= LZMAX) then
  if (PART(I,J)%ZP >= ZERO .and. PART(I,J)%ZP < DELTA_COLL) then
    PART(I,J)%ZP = PART(I,J)%ZP + LZMAX
  end if
 end if


 end do


else

 NPNEIGH_COLL = 0

 NS = NPART_LOC(J)+1
 NF = NPART_LOC(J)+NPNEIGH_COLL

end if !!- end if NPROC>1




!!====================================================================
!! 2. Detection of neighboring particles
!!====================================================================
call COLLISION_DETECTION(NF,                             &
                         XMESH(ISTART(1))   -DELTA_COLL, & !- XMIN
                         XMESH(  IEND(1))+DX+DELTA_COLL, & !- XMAX
                         YMESH(ISTART(2))   -DELTA_COLL, & !- YMIN
                         YMESH(  IEND(2))+DY+DELTA_COLL, & !- YMAX
                         ZMESH(ISTART(3))   -DELTA_COLL, & !- ZMIN
                         ZMESH(  IEND(3))+DZ+DELTA_COLL, & !- ZMAX
                         PART(1:NF,J)%XP, &
                         PART(1:NF,J)%YP, &
                         PART(1:NF,J)%ZP, &
                                  NCLOSE, & 
                              NPCLOSEMAX, &
                                    HOME, &
                                    NEAR  )




!!====================================================================
!! 3. Particle-particle interaction
!!====================================================================
call COLLISION_INTERACTION(TIME,NCYCLE,J,NCLOSE,HOME,NEAR)



end do !- end loop  J = 1, NIG



!!====================================================================
!! 3. Periodic boundary conditions
!!====================================================================
!!- Boundary conditions if particles move
! call BOUNDARY_PARTICLE


!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine COLLISION

