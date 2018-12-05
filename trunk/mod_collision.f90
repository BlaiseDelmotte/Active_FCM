module COLLISION_VARIABLE

 implicit none 

 real (kind=8) :: DELTA_COLL
 integer :: NPNEIGH_COLL

 integer, parameter :: NPDFMAX = 40
 real(kind=8), dimension(NPDFMAX,5) :: PDF
 real(kind=8)  :: NMOM

end module COLLISION_VARIABLE


