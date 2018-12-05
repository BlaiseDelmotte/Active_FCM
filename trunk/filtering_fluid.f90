!!====================================================================
!!
!!         Filtering of Fluid Velocity Field
!!
!!====================================================================

subroutine FILTERING_FLUID(IG)

!!====================================================================

use DNS_DIM	        !- Dimension
use PARAM_PHYS          !- Physical & numerical parameters
use FLUID_VARIABLE      !- Fluid velocity
use GEOMETRIC_VARIABLE
use WORK_ARRAYS
use PARTICLE_PARALLEL
use MPI_STRUCTURES


use P3DFFT

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!!------------------------------------------------------------------
!! ARGUMENTS
!!------------------------------------------------------------------
integer, intent(in) :: IG
!!
!!------------------------------------------------------------------
!! LOCAL ARRAYS STATEMENT
!!------------------------------------------------------------------
real(kind=8), dimension(ISTART(1)         :IEND(1)   ,  &
                        ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL, &
                        ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL   ) :: UFILTER

real(kind=8), dimension(ISTART(1)         :IEND(1)        ,  &
                             ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL, &
                             ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL   ) :: DUF

!!- Squared wavenumber
real(kind=8) :: K2

!!- Squared wavenumber for cutoff
real(kind=8) :: K2CUT

!!- Index
integer :: I, J, K
!!====================================================================


!!- Define squared wavenumber for cutoff
K2CUT = KCUT**2


!!====================================================================
!! 1. x-component
!!====================================================================

!!- Fill temporary array with the eulerian value
TMPFOU(:,:,:) =UFOU(:,:,:)

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

  K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2
    
  if(K2 >= K2CUT)TMPFOU(I,J,K)=cmplx(ZERO,ZERO) 
  
  end do
 end do
end do


!!- Tranform back in physical space
call P3DFFT_BTRAN_C2R(TMPFOU,UFILTER(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)

!!- Subgrid velocity
DUF(:,:,:) = UFLU(:,:,:) - UFILTER(:,:,:)

!!- Halo updating for interpolation
if(NGHTCELL>0) call FLUIDCOMM(UFILTER)
if(NGHTCELL>0) call FLUIDCOMM(DUF)

!!- Filtered fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                 &
		      XMESH,                  &
                      YMESH,                  &
                      ZMESH,                  &
                    UFILTER,                  &
                NPART_LOC(IG),                &
 	        PART(1:NPART_LOC(IG),IG)%XP,  &
                PART(1:NPART_LOC(IG),IG)%YP,  &
                PART(1:NPART_LOC(IG),IG)%ZP,  &
 	        PART(1:NPART_LOC(IG),IG)%UFAP )


!!- Subgrid fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                 &
		      XMESH,                  &
                      YMESH,                  &
                      ZMESH,                  &
                        DUF,                  &
                NPART_LOC(IG),                &
 	        PART(1:NPART_LOC(IG),IG)%XP,  &
                PART(1:NPART_LOC(IG),IG)%YP,  &
                PART(1:NPART_LOC(IG),IG)%ZP,  &
 	        PART(1:NPART_LOC(IG),IG)%DUFAP )


!!====================================================================
!! 2. y-component
!!====================================================================

!!- Fill temporary array with the eulerian value
TMPFOU(:,:,:) =VFOU(:,:,:)

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

  K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2
    
  if(K2 >= K2CUT)TMPFOU(I,J,K)=cmplx(ZERO,ZERO)
  
  end do
 end do
end do


!!- Tranform back in physical space
call P3DFFT_BTRAN_C2R(TMPFOU,UFILTER(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)


!!- Subgrid velocity
DUF(:,:,:) = VFLU(:,:,:) - UFILTER(:,:,:)


!!- Halo updating for interpolation
if(NGHTCELL>0) call FLUIDCOMM(UFILTER)
if(NGHTCELL>0) call FLUIDCOMM(DUF)


!!- Filtered fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                 &
		      XMESH,                &
                      YMESH,                &
                      ZMESH,                &
                       UFILTER,             &
                NPART_LOC(IG),               &
 	        PART(1:NPART_LOC(IG),IG)%XP,  &
                PART(1:NPART_LOC(IG),IG)%YP,  &
                PART(1:NPART_LOC(IG),IG)%ZP,  &
 	        PART(1:NPART_LOC(IG),IG)%VFAP )

!!- Subgrid fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                 &
		      XMESH,                  &
                      YMESH,                  &
                      ZMESH,                  &
                        DUF,                  &
                NPART_LOC(IG),                &
 	        PART(1:NPART_LOC(IG),IG)%XP,  &
                PART(1:NPART_LOC(IG),IG)%YP,  &
                PART(1:NPART_LOC(IG),IG)%ZP,  &
 	        PART(1:NPART_LOC(IG),IG)%DVFAP )

!!====================================================================
!! 3. z-component
!!====================================================================

!!- Fill temporary array with the eulerian value
TMPFOU(:,:,:) =WFOU(:,:,:)

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)

  K2 = KX(I)**2 + KY(J)**2 + KZ(K)**2
    
  if(K2 >= K2CUT)TMPFOU(I,J,K)=cmplx(ZERO,ZERO) 
  end do
 end do
end do


!!- Tranform back in physical space
call P3DFFT_BTRAN_C2R(TMPFOU,UFILTER(:,ISTART(2):IEND(2),ISTART(3):IEND(3)),FFTFLAG)

!!- Subgrid velocity
DUF(:,:,:) = WFLU(:,:,:) - UFILTER(:,:,:)


!!- Halo updating for interpolation
if(NGHTCELL>0) call FLUIDCOMM(UFILTER)
if(NGHTCELL>0) call FLUIDCOMM(DUF)


!!- Filtered fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                &
		     XMESH,                  &
                     YMESH,                  &
                     ZMESH,                  &
                   UFILTER,                  &
               NPART_LOC(IG),                &
 	       PART(1:NPART_LOC(IG),IG)%XP,  &
               PART(1:NPART_LOC(IG),IG)%YP,  &
               PART(1:NPART_LOC(IG),IG)%ZP,  &
 	       PART(1:NPART_LOC(IG),IG)%WFAP )


!!- Subgrid fluid velocity at particle position
call INTERPH(INTERP_SCHEME,                 &
		      XMESH,                  &
                      YMESH,                  &
                      ZMESH,                  &
                        DUF,                  &
                NPART_LOC(IG),                &
 	        PART(1:NPART_LOC(IG),IG)%XP,  &
                PART(1:NPART_LOC(IG),IG)%YP,  &
                PART(1:NPART_LOC(IG),IG)%ZP,  &
 	        PART(1:NPART_LOC(IG),IG)%DWFAP )


!!====================================================================

end subroutine FILTERING_FLUID
