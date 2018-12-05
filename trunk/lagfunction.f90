!!====================================================================
!!
!! Compute the Lagrangian correlation function for the filtered and
!! the subgrid the fluid velocity at the particle position
!!
!!====================================================================
!!
subroutine LAGFUNCTION(IFLAG,NCYCLE,                      &
                         UP_MEAN,    VP_MEAN,    WP_MEAN, &
                        UUP_MEAN,   VVP_MEAN,   WWP_MEAN, &
                       UFAP_MEAN,  VFAP_MEAN,  WFAP_MEAN, &
                      UUFAP_MEAN, VVFAP_MEAN, WWFAP_MEAN  )
!!
!!====================================================================
!! The normalized autocorrelation function is computed as
!!
!!                 <u(t0).u(t0+t)>
!! R(t) = ---------------------------------
!!        sqrt(<u(t0)^2>).sqrt(<u(t0+t)^2>)
!!
!!--------------------------------------------------------------------
!! IFLAG  <2: Computation
!!       ==2: Printings
!!====================================================================

use dns_dim
use param_phys
use particle_parallel
use statistics

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------------------------------------------------------------
!- Flag
integer, intent(in) :: IFLAG

!- Cycle number
integer, intent(in) :: NCYCLE

!- 
real(kind=8), dimension(NIG), intent(in) :: UP_MEAN
real(kind=8), dimension(NIG), intent(in) :: VP_MEAN
real(kind=8), dimension(NIG), intent(in) :: WP_MEAN
real(kind=8), dimension(NIG), intent(in) :: UUP_MEAN
real(kind=8), dimension(NIG), intent(in) :: VVP_MEAN
real(kind=8), dimension(NIG), intent(in) :: WWP_MEAN
real(kind=8), dimension(NIG), intent(in) :: UFAP_MEAN
real(kind=8), dimension(NIG), intent(in) :: VFAP_MEAN
real(kind=8), dimension(NIG), intent(in) :: WFAP_MEAN
real(kind=8), dimension(NIG), intent(in) :: UUFAP_MEAN
real(kind=8), dimension(NIG), intent(in) :: VVFAP_MEAN
real(kind=8), dimension(NIG), intent(in) :: WWFAP_MEAN



!--------------------------------------------------------------------
!- Local arrays
!--------------------------------------------------------------------
!- Lagrangian function
real(kind=8), dimension(DIMLGR) ::   RPX,   RPY,   RPZ
real(kind=8), dimension(DIMLGR) :: RFAPX, RFAPY, RFAPZ

!- True time step (accouting for frequency)
real(kind=8) :: DTIME0

!- 
integer :: NT0, NTLGR 

!- Index
integer :: I, J, LGR, K
!---------------------------------------------------------------------


DTIME0 = DTIME*FOUT2
NT0 = 1



if(IFLAG == 1) then


!!====================================================================
!! 1. Storage Lagrangian fluctuating velocity at t=t0
!!====================================================================

if(NCYCLE == NT0*FOUT2) then

if(MYID==0) write(*,*)'Lagrangian function computing starts'

do J = 1, NIG

!- Particle velocity 
 PART(:,J)%UPT0 = PART(:,J)%UP - UP_MEAN(J)
 PART(:,J)%VPT0 = PART(:,J)%VP - VP_MEAN(J)
 PART(:,J)%WPT0 = PART(:,J)%WP - WP_MEAN(J)

!- Fluid velocity at particle position
 PART(:,J)%UFAPT0 = PART(:,J)%UFAP - UFAP_MEAN(J)
 PART(:,J)%VFAPT0 = PART(:,J)%VFAP - VFAP_MEAN(J)
 PART(:,J)%WFAPT0 = PART(:,J)%WFAP - WFAP_MEAN(J)

!- Variance
 UUPT0_MEAN(J) = UUP_MEAN(J)
 VVPT0_MEAN(J) = VVP_MEAN(J)
 WWPT0_MEAN(J) = WWP_MEAN(J)

!- Variance
 UUFAPT0_MEAN(J) = UUFAP_MEAN(J)
 VVFAPT0_MEAN(J) = VVFAP_MEAN(J)
 WWFAPT0_MEAN(J) = WWFAP_MEAN(J)


end do

end if !- If: (NCYCLE == NT0)




!!====================================================================
!! 2. Correlation at time t
!!====================================================================
NTLGR = int( (NCYCLE - NT0*FOUT2)/FOUT2) + 1

if(NTLGR<DIMLGR) then

do J = 1, NIG

RPX_LOC(NTLGR,J) = 0.
RPY_LOC(NTLGR,J) = 0.
RPY_LOC(NTLGR,J) = 0.

do I = 1, NPART_LOC(J)

!!-------------------------------------------------------------------
!! 2.1 Particle-particle
!!-------------------------------------------------------------------
!- <up.up>
 RPX_LOC(NTLGR,J) = RPX_LOC(NTLGR,J)                             &
              + (PART(I,J)%UP-UP_MEAN(J))*PART(I,J)%UPT0 &
                 /sqrt(UUPT0_MEAN(J)*UUP_MEAN(J))

!- <vp.vp>
 RPY_LOC(NTLGR,J) = RPY_LOC(NTLGR,J)                             &
              + (PART(I,J)%VP-VP_MEAN(J))*PART(I,J)%VPT0 &
                 /sqrt(VVPT0_MEAN(J)*VVP_MEAN(J))

!- <wp.wp>
 RPZ_LOC(NTLGR,J) = RPZ_LOC(NTLGR,J)                             &
              + (PART(I,J)%WP-WP_MEAN(J))*PART(I,J)%WPT0 &
                 /sqrt(WWPT0_MEAN(J)*WWP_MEAN(J))


!!-------------------------------------------------------------------
!! 2.2 Fluid-particle
!!-------------------------------------------------------------------

!- <ufap.ufap>
 RFAPX_LOC(NTLGR,J) = RFAPX_LOC(NTLGR,J)                       &
              + (PART(I,J)%UFAP-UFAP_MEAN(J))*PART(I,J)%UFAPT0 &
                 /sqrt(UUFAPT0_MEAN(J)*UUFAP_MEAN(J))

!- <vfap.vfap>
 RFAPY_LOC(NTLGR,J) = RFAPY_LOC(NTLGR,J)                       &
              + (PART(I,J)%VFAP-VFAP_MEAN(J))*PART(I,J)%VFAPT0 &
                 /sqrt(VVFAPT0_MEAN(J)*VVFAP_MEAN(J))

!- <wfap.wfap>
 RFAPZ_LOC(NTLGR,J) = RFAPZ_LOC(NTLGR,J)                       &
              + (PART(I,J)%WFAP-WFAP_MEAN(J))*PART(I,J)%WFAPT0 &
                 /sqrt(WWFAPT0_MEAN(J)*WWFAP_MEAN(J))


end do !- Loop: I = 1, NPART_LOC(J)


!!-------------------------------------------------------------------
!! 2.3 Normalization
!!-------------------------------------------------------------------
 RPX_LOC(NTLGR,J) = RPX_LOC(NTLGR,J) / real(NPART_LOC(J))
 RPY_LOC(NTLGR,J) = RPY_LOC(NTLGR,J) / real(NPART_LOC(J))
 RPZ_LOC(NTLGR,J) = RPZ_LOC(NTLGR,J) / real(NPART_LOC(J))

 RFAPX_LOC(NTLGR,J) = RFAPX_LOC(NTLGR,J) / real(NPART_LOC(J))
 RFAPY_LOC(NTLGR,J) = RFAPY_LOC(NTLGR,J) / real(NPART_LOC(J))
 RFAPZ_LOC(NTLGR,J) = RFAPZ_LOC(NTLGR,J) / real(NPART_LOC(J))


end do !- Loop:  J = 1, NIG


end if !- If: (NTLGR<DIMLGR)

!!====================================================================
!! 3. Print and so on ...
!!====================================================================
else !- elseIf: (IFLAG == 1)


do J = 1, NIG
!!- Add all Lag function from each cpu
if (NPROC>1) then
 call MPI_ALLREDUCE(RPX_LOC(:,J),RPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPY_LOC(:,J),RPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPZ_LOC(:,J),RPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

 call MPI_ALLREDUCE(RFAPX_LOC(:,J),RPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RFAPY_LOC(:,J),RPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RFAPZ_LOC(:,J),RPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!- Normalization
 RPX = RPX / real(NPROC)
 RPY = RPY / real(NPROC)
 RPZ = RPZ / real(NPROC)

 RFAPX = RFAPX / real(NPROC)
 RFAPY = RFAPY / real(NPROC)
 RFAPZ = RFAPZ / real(NPROC)
end if



!!--------------------------------------------------------------------
!! 3.1 Print in file
!!--------------------------------------------------------------------
if(MYID==0) then

!- Print in file
do I = 1,DIMLGR 
  write(560+J,10000) (I-1)*DTIME0,                     &
                           RPX(I),   RPY(I), RPZ(I),   &
        (RPX(I)+RPY(I)+RPZ(I))/3.,                     &
                         RFAPX(I), RFAPY(I), RFAPZ(I), &
        (RFAPX(I)+RFAPY(I)+RFAPZ(I))/3.0                     
end do 

end if !- If: (MYID==0)

end do  !!- Loop: J = 1, NIG



end if !- If: (IFLAG == 1)


!!--------------------------------------------------------------------
10600 format (2x,A,I2)
10601 format (2x,A,e12.6,A)
10000 format (15(e17.7))
10001 format (50(e17.7))

end subroutine LAGFUNCTION 
