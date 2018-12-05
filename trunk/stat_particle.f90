!!====================================================================
!!
!! This fortran file compute all statistics on the paticles
!!
!!====================================================================

subroutine STAT_PARTICLE(NCYCLE,TIME)

!!====================================================================
!! We compute 4 levels of statistics, with graduate computational cost.
!!
!! LEVEL1_STPAR: + Mean velocity
!!               + Kinetic stress
!!
!! LEVEL2_STPAR: + Mean velocity
!!               + Kinetic stress
!!
!!--------------------------------------------------------------------
!! Time-averaged statistics are performed using a macro array
!! called MEAN_TIME_PART. The averaging is done as
!!   --> MEAN_TIME_PART = MEAN_TIME_PART + MEAN_PART
!!
!!--------------------------------------------------------------------
!!  MEAN_PART( 1): <up>
!!             2 : <vp>
!!             3 : <wp>
!!             4 : <up.up>
!!             5 : <vp.vp>
!!             6 : <wp.wp>
!!             7 : <up.vp>
!!             8 : <up.wp>
!!             9 : <vp.wp>
!!            10 :  qp = (<up.up>+<vp.vp>+<wp.wp>)/2
!!
!!            11 : <uf@p>
!!            12 : <vf@p>
!!            13 : <wf@p>
!!            14 : <uf@p.uf@p>
!!            15 : <vf@p.vf@p>
!!            16 : <wf@p.wf@p>
!!            17 : <uf@p.vf@p>
!!            18 : <uf@p.wf@p>
!!            19 : <vf@p.wf@p>
!!            20 :  qf@p = (<uf@p.uf@p>+<vf@p.vf@p>+<wf@p.wf@p>)/2
!!
!!            21 : <up.uf@p>
!!            22 : <vp.vf@p>
!!            23 : <wp.wf@p>
!!            24 : <up.vf@p>
!!            25 : <up.wf@p>
!!            26 : <vp.wf@p>
!!            27 :  qfp = <up.uf@p>+<vp.vf@p>+<wp.wf@p>
!!
!!            28 : <Vr,fp>
!!            29 : <Rep>
!!            30 : <Cd>
!!            31 : <1/tau_p>
!!-- Temperature
!!            35 : <Tp>
!!            36 : <Tf@p>
!!            37 : <Tp^2>
!!            38 : <Tf@p^2>
!!            39 : <Tp.up>
!!            40 : <Tp.vp>
!!            41 : <Tp.wp>
!!            42 : <Tp.uf@p>
!!            43 : <Tp.vf@p>
!!            44 : <Tp.wf@p>
!!            45 : <Tf@p.up>
!!            46 : <Tf@p.vp>
!!            47 : <Tf@p.wp>
!!            48 : <Tf@p.uf@p>
!!            49 : <Tf@p.vf@p>
!!            50 : <Tf@p.wf@p>
!!            51 : <Nup>
!!            52 : <1/tau_theta>
!!            53 : <Tp.Tf@p>
!!            54 : <Tp.dTp/dt>
!!            55 : <Tp.dvp/dt>
!!            56 : <vp.dTp/dt>
!!            57 : <Tf@p.dTp/dt>
!!            58 : <Tp.dTf@p/dt>
!!            59 : 
!!            60 : 
!!            61 : <np>
!!            62 : <np'^2>
!!            63 : <np^2>
!!            64 : 
!!            65 : 
!!            66 : 
!!            67 : 
!!            68 : 
!!            69 : 
!!            70 : <dufap>
!!            71 : <dvfap>
!!            72 : <dwfap>
!!            73 : <dufap.dufap>
!!            74 : <dvfap.dvfap>
!!            75 : <dwfap.dwfap>
!!            76 : <dufap.dvfap>
!!            77 : <dufap.dwfap>
!!            78 : <dvfap.dwfap>
!!            79 : dqfap
!!            80 : <dufap.ufap>
!!            81 : <dvfap.vfap>
!!            82 : <dwfap.wfap>
!!            83 : <dufap.vfap>
!!            84 : <dufap.wfap>
!!            85 : <dvfap.wfap>
!!            86 : <dufap.up>
!!            87 : <dufap.vp>
!!            88 : <dvfap.wp>
!!
!!====================================================================

use DNS_DIM	       !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use PARTICLE_PARALLEL
use GEOMETRIC_VARIABLE !- 
use STATISTICS         !- Statistics
use CHECK_CPU	       !- CPU time checks

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Global arrays
!---------------------------------------------------------------------
!- cycle number
integer, intent(in) :: NCYCLE

!- Curent time
real(kind=8), intent(in) :: TIME

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Mean velocities
real(kind=8) :: UPM, VPM, WPM
real(kind=8) :: UFAPM, VFAPM, WFAPM
real(kind=8) :: DUFAPM, DVFAPM, DWFAPM

!- Mean Temperature
real(kind=8) :: TPM, TFAPM

!- Fluctuating particle velocities
real(kind=8) :: UPFLC, VPFLC, WPFLC

!- Fluctuating fluid velocities
real(kind=8) :: UFAPFLC, VFAPFLC, WFAPFLC

!- Fluctuating if subgrid fluid velocities
real(kind=8) :: DUFAPFLC, DVFAPFLC, DWFAPFLC

!- Fluctuating temperatures
real(kind=8) :: TPFLC, TFAPFLC

!- Inverse of dynamic particle response time
real(kind=8) :: INVTAUP

!- Fluid-Particle relative velocity
real(kind=8) :: VRNRM

!- Particle Reynolds number
real(kind=8) :: REP

!- Drag coefficient
real(kind=8) :: CDRAG

!- Nusslet Number
real(kind=8) :: NUP

!- Inverse of heat particle response time
real(kind=8) :: INVTAUP_THETA

!- Number of particle
integer :: NPART




!- Particle concentration
real(kind=8) :: DXCP, DYCP,DZCP, CMEAN


!- 
integer :: NTLGR 

!- 
real(kind=8) :: RDUMMY

integer :: IFLAG


!- Time control variable
real(kind=8) :: TIME_START, TIME_END

integer :: I,J,K, IFLAG1, IPART, JPART, KPART, IPDF, NLGR
!---------------------------------------------------------------------

!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


MEAN_PART(:,:) = ZERO
MEAN_PART_LOC(:,:) = ZERO


!!====================================================================
!! 1. 1st level of  statistic
!!====================================================================
if(LEVEL1_STPAR) then


do J = 1, NIG


!!--------------------------------------------------------------------
!! 1.1. Mean variable
!!--------------------------------------------------------------------

!!- Full number of particle 
call ISUMCPU(NPART_LOC(J),NPART)


 do I = 1, NPART_LOC(J)
!!--------------------------------------------------------------------
!!- Particle velocities
!!--------------------------------------------------------------------
!!- <up>
  MEAN_PART_LOC(1,J) = MEAN_PART_LOC(1,J) + PART(I,J)%UP

!!- <vp>
  MEAN_PART_LOC(2,J) = MEAN_PART_LOC(2,J) + PART(I,J)%VP

!!- <wp>
  MEAN_PART_LOC(3,J) = MEAN_PART_LOC(3,J) + PART(I,J)%WP

!!--------------------------------------------------------------------
!!- Fluid velocities at particle position
!!--------------------------------------------------------------------
!!- <uf@p>
  MEAN_PART_LOC(11,J) = MEAN_PART_LOC(11,J) + PART(I,J)%UFAP

!!- <vf@p>
  MEAN_PART_LOC(12,J) = MEAN_PART_LOC(12,J) + PART(I,J)%VFAP

!!- <wf@p>
  MEAN_PART_LOC(13,J) = MEAN_PART_LOC(13,J) + PART(I,J)%WFAP
 end do



!!- Compute the mean velocities over the whole domain
 call RSUMCPU(MEAN_PART_LOC( 1,J),UPM)
 call RSUMCPU(MEAN_PART_LOC( 2,J),VPM)
 call RSUMCPU(MEAN_PART_LOC( 3,J),WPM)
 
 call RSUMCPU(MEAN_PART_LOC(11,J),UFAPM)
 call RSUMCPU(MEAN_PART_LOC(12,J),VFAPM)
 call RSUMCPU(MEAN_PART_LOC(13,J),WFAPM)


 UPM = UPM/ real(NPART)
 VPM = VPM/ real(NPART)
 WPM = WPM/ real(NPART)
 
 UFAPM = UFAPM/ real(NPART)
 VFAPM = VFAPM/ real(NPART)
 WFAPM = WFAPM/ real(NPART)



!!--------------------------------------------------------------------
!! 1.2. Particle kinetic stress
!!--------------------------------------------------------------------
 do I = 1, NPART_LOC(J)

!!- u'p = up-<up>
  UPFLC = PART(I,J)%UP - UPM

!!- v'p = vp-<vp>
  VPFLC = PART(I,J)%VP - VPM

!!- w'p = wp-<wp>
  WPFLC = PART(I,J)%WP - WPM


!!- <up*up>
  MEAN_PART_LOC(4,J) = MEAN_PART_LOC(4,J) + UPFLC*UPFLC

!!- <vp*vp>
  MEAN_PART_LOC(5,J) = MEAN_PART_LOC(5,J) + VPFLC*VPFLC

!!- <wp*wp>
  MEAN_PART_LOC(6,J) = MEAN_PART_LOC(6,J) + WPFLC*WPFLC

!!- <up*vp>
  MEAN_PART_LOC(7,J) = MEAN_PART_LOC(7,J) + UPFLC*VPFLC

!!- <up*wp>
  MEAN_PART_LOC(8,J) = MEAN_PART_LOC(8,J) + UPFLC*WPFLC

!!- <vp*wp>
  MEAN_PART_LOC(9,J) = MEAN_PART_LOC(9,J) + VPFLC*WPFLC

 end do

!!- qp
 MEAN_PART_LOC(10,J) = 0.5*(  MEAN_PART_LOC(4,J) &
                            + MEAN_PART_LOC(5,J) &
                            + MEAN_PART_LOC(6,J) )


!!--------------------------------------------------------------------
!! 1.3. Fluid at particle position kinetic stress
!!--------------------------------------------------------------------
 do I = 1, NPART_LOC(J)

!!- u'f@p = uf@p-<uf@p>
  UFAPFLC = PART(I,J)%UFAP - UFAPM

!!- v'f@p = vf@p-<vf@p>
  VFAPFLC = PART(I,J)%VFAP - VFAPM

!!- w'f@p = wf@p-<wf@p>
  WFAPFLC = PART(I,J)%WFAP - WFAPM


!!- <uf@p*uf@p>
  MEAN_PART_LOC(14,J) = MEAN_PART_LOC(14,J) + UFAPFLC*UFAPFLC

!!- <vf@p*vf@p>
  MEAN_PART_LOC(15,J) = MEAN_PART_LOC(15,J) + VFAPFLC*VFAPFLC

!!- <wf@p*wf@p>
  MEAN_PART_LOC(16,J) = MEAN_PART_LOC(16,J) + WFAPFLC*WFAPFLC

!!- <uf@p*vf@p>
  MEAN_PART_LOC(17,J) = MEAN_PART_LOC(17,J) + UFAPFLC*VFAPFLC

!!- <uf@p*wf@p>
  MEAN_PART_LOC(18,J) = MEAN_PART_LOC(18,J) + UFAPFLC*WFAPFLC

!!- <vf@p*wf@p>
  MEAN_PART_LOC(19,J) = MEAN_PART_LOC(19,J) + VFAPFLC*WFAPFLC

 end do

!!- qfap
 MEAN_PART_LOC(20,J) = 0.5*(MEAN_PART_LOC(14,J) &
                           +MEAN_PART_LOC(15,J) &
                           +MEAN_PART_LOC(16,J))

!!--------------------------------------------------------------------
!! 1.4. Fluid-particle correlation
!!--------------------------------------------------------------------
 do I = 1, NPART_LOC(J)

!!- u'p = up-<up>
  UPFLC = PART(I,J)%UP - UPM

!!- v'p = vp-<vp>
  VPFLC = PART(I,J)%VP - VPM

!!- w'p = wp-<wp>
  WPFLC = PART(I,J)%WP - WPM

!!- u'f@p = uf@p-<uf@p>
  UFAPFLC = PART(I,J)%UFAP - UFAPM

!!- v'f@p = vf@p-<vf@p>
  VFAPFLC = PART(I,J)%VFAP - VFAPM

!!- w'f@p = wf@p-<wf@p>
  WFAPFLC = PART(I,J)%WFAP - WFAPM



!!- <uf@p*up>
  MEAN_PART_LOC(21,J) = MEAN_PART_LOC(21,J) + UFAPFLC*UPFLC

!!- <vf@p*vp>
  MEAN_PART_LOC(22,J) = MEAN_PART_LOC(22,J) + VFAPFLC*VPFLC

!!- <wf@p*wp>
  MEAN_PART_LOC(23,J) = MEAN_PART_LOC(23,J) + WFAPFLC*WPFLC

!!- <uf@p*vp>
  MEAN_PART_LOC(24,J) = MEAN_PART_LOC(24,J) + UFAPFLC*VPFLC

!!- <uf@p*wp>
  MEAN_PART_LOC(25,J) = MEAN_PART_LOC(25,J) + UFAPFLC*WPFLC

!!- <vf@p*wp>
  MEAN_PART_LOC(26,J) = MEAN_PART_LOC(26,J) + VFAPFLC*WPFLC

 end do

!!- qfp
 MEAN_PART_LOC(27,J) = MEAN_PART_LOC(21,J) &
                     + MEAN_PART_LOC(22,J) &
                     + MEAN_PART_LOC(23,J)


!!--------------------------------------------------------------------
!! 1.5. Fluid-particle correlation
!!--------------------------------------------------------------------
if(PARTDEF(J) >1) then

 do I = 1, NPART_LOC(J)
 
 call GASPART_TRANSFER(         J,  & ! <- Particle class
                       PART(I,J)%UP,  & ! <- Part x-vel
		       PART(I,J)%VP,  & ! <- Part y-vel
		       PART(I,J)%WP,  & ! <- Part z-vel
		     PART(I,J)%UFAP,  & ! <- Fluid x-vel
		     PART(I,J)%VFAP,  & ! <- Fluid y-vel
		     PART(I,J)%WFAP,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

!!- <1/taup>
  MEAN_PART_LOC(28,J) = MEAN_PART_LOC(28,J) + INVTAUP

!!- <VR>
  MEAN_PART_LOC(29,J) = MEAN_PART_LOC(29,J) + VRNRM

!!- <Cd>
  MEAN_PART_LOC(30,J) = MEAN_PART_LOC(30,J) + CDRAG

!!- <Rep>
  MEAN_PART_LOC(31,J) = MEAN_PART_LOC(31,J) + REP


 end do

else

!!- <1/taup>
  MEAN_PART_LOC(28,J) = 1.0

!!- <VR>
  MEAN_PART_LOC(29,J) = ZERO

!!- <Cd>
  MEAN_PART_LOC(30,J) = ZERO

!!- <Rep>
  MEAN_PART_LOC(31,J) = ZERO



end if !!- If:(PARTDEF(J) >1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!- Temperature
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
if(SOLVE_SCALAR) then

!!--------------------------------------------------------------------
!! 1.1. Mean variable
!!--------------------------------------------------------------------
 do I = 1, NPART_LOC(J)
!!- <Tp>
  MEAN_PART_LOC(35,J) = MEAN_PART_LOC(35,J) + PART(I,J)%TP

!!- <Tf@p>
  MEAN_PART_LOC(36,J) = MEAN_PART_LOC(36,J) + PART(I,J)%TFAP

 end do
 

!- Compute the mean velocities over the whole domain
 call RSUMCPU(MEAN_PART_LOC(35,J),TPM)
 call RSUMCPU(MEAN_PART_LOC(36,J),TFAPM)


 TPM = TPM / real(NPART)
 TFAPM = TFAPM / real(NPART)
 


!!--------------------------------------------------------------------
!! 1.2. Temperature-velocity correlation
!!--------------------------------------------------------------------
 do I = 1, NPART_LOC(J)
 
 !!- T'p = Tp-<Tp>
  TPFLC = PART(I,J)%TP - TPM

!!- T'f@p = Tf@p-<Tf@p>
  TFAPFLC = PART(I,J)%TFAP - TFAPM

!!- u'p = up-<up>
  UPFLC = PART(I,J)%UP - UPM

!!- v'p = vp-<vp>
  VPFLC = PART(I,J)%VP - VPM

!!- w'p = wp-<wp>
  WPFLC = PART(I,J)%WP - WPM

!!- u'f@p = uf@p-<uf@p>
  UFAPFLC = PART(I,J)%UFAP - UFAPM

!!- v'f@p = vf@p-<vf@p>
  VFAPFLC = PART(I,J)%VFAP - VFAPM

!!- w'f@p = wf@p-<wf@p>
  WFAPFLC = PART(I,J)%WFAP - WFAPM



!!- <T'p^2>
  MEAN_PART_LOC(37,J) = MEAN_PART_LOC(37,J) + TPFLC*TPFLC

!!- <T'f@p^2>
  MEAN_PART_LOC(38,J) = MEAN_PART_LOC(38,J) + TFAPFLC*TFAPFLC


!!- <T'p.u'p>
  MEAN_PART_LOC(39,J) = MEAN_PART_LOC(39,J) + TPFLC*UPFLC
  
!!- <T'p.v'p>
  MEAN_PART_LOC(40,J) = MEAN_PART_LOC(40,J) + TPFLC*VPFLC
  
!!- <T'p.w'p>
  MEAN_PART_LOC(41,J) = MEAN_PART_LOC(41,J) + TPFLC*WPFLC



!!- <T'p.u'f@p>
  MEAN_PART_LOC(42,J) = MEAN_PART_LOC(42,J) + TPFLC*UFAPFLC
  
!!- <T'p.v'f@p>
  MEAN_PART_LOC(43,J) = MEAN_PART_LOC(43,J) + TPFLC*VFAPFLC
  
!!- <T'p.w'f@p>
  MEAN_PART_LOC(44,J) = MEAN_PART_LOC(44,J) + TPFLC*WFAPFLC
  
  
  
 !!- <T'f@p.u'p>
  MEAN_PART_LOC(45,J) = MEAN_PART_LOC(45,J) + TFAPFLC*UPFLC
  
!!- <T'f@p.v'p>
  MEAN_PART_LOC(46,J) = MEAN_PART_LOC(46,J) + TFAPFLC*VPFLC
  
!!- <T'f@p.w'p>
  MEAN_PART_LOC(47,J) = MEAN_PART_LOC(47,J) + TFAPFLC*WPFLC



!!- <T'f@p.u'f@p>
  MEAN_PART_LOC(48,J) = MEAN_PART_LOC(48,J) + TFAPFLC*UFAPFLC
  
!!- <T'f@p.v'f@p>
  MEAN_PART_LOC(49,J) = MEAN_PART_LOC(49,J) + TFAPFLC*VFAPFLC
  
!!- <T'f@p.w'f@p>
  MEAN_PART_LOC(50,J) = MEAN_PART_LOC(50,J) + TFAPFLC*WFAPFLC

!!- <T'p.T'p>
  MEAN_PART_LOC(53,J) = MEAN_PART_LOC(53,J) + TPFLC*TFAPFLC
  


!!- <Tp.dTp/dt>
  MEAN_PART_LOC(54,J) = MEAN_PART_LOC(54,J) &
                     + TPFLC*(PART(I,J)%TP-PART(I,J)%TP_NM1)/DTIME
!!- <Tp.dvp/dt>
  MEAN_PART_LOC(55,J) = MEAN_PART_LOC(55,J) &
                     + TPFLC*(PART(I,J)%VP-PART(I,J)%VP_NM1)/DTIME
!!- <vp.dTp/dt>
  MEAN_PART_LOC(56,J) = MEAN_PART_LOC(56,J) &
                     + VPFLC*(PART(I,J)%TP-PART(I,J)%TP_NM1)/DTIME
!!- <Tp.dTfap/dt>
  MEAN_PART_LOC(57,J) = MEAN_PART_LOC(57,J) &
                     + TPFLC*(PART(I,J)%TFAP-PART(I,J)%TFAP_NM1)/DTIME
 
!!- <Tfap.dTp/dt>
  MEAN_PART_LOC(58,J) = MEAN_PART_LOC(58,J) &
                     + TFAPFLC*(PART(I,J)%TP-PART(I,J)%TP_NM1)/DTIME
    

 end do


!!- Fluid-particle correlation
if(PARTDEF(J) >1) then

 do I = 1, NPART_LOC(J)
 
 call GASPART_TRANSFER(         J,  & ! <- Particle class
                       PART(I,J)%UP,  & ! <- Part x-vel
		       PART(I,J)%VP,  & ! <- Part y-vel
		       PART(I,J)%WP,  & ! <- Part z-vel
		     PART(I,J)%UFAP,  & ! <- Fluid x-vel
		     PART(I,J)%VFAP,  & ! <- Fluid y-vel
		     PART(I,J)%WFAP,  & ! <- Fluid z-vel
                            INVTAUP,  & ! -> Particle response time
			      VRNRM,  & ! -> |Up - Uf@p|
			      CDRAG,  & ! -> Drag coefficient
			        REP,  & ! -> Part Reynolds number
		      INVTAUP_THETA,  & ! -> Heat response time
		                NUP   ) ! -> Nusselt Number

!!- <1/taup_theta>
  MEAN_PART_LOC(51,J) = MEAN_PART_LOC(51,J) + INVTAUP_THETA

!!- <Nup>
  MEAN_PART_LOC(52,J) = MEAN_PART_LOC(52,J) + NUP

 end do

else

  MEAN_PART_LOC(51,J) = 1.0
  MEAN_PART_LOC(52,J) = ZERO
  
end if !!- If:(PARTDEF(J) >1)


 end if !!- End if SOLVE_SCALAR




!!--------------------------------------------------------------------
!! 1.3. Subgrid turbulence seen by solid particles
!!--------------------------------------------------------------------
if(FILTERING) then

 do I = 1, NPART_LOC(J)
!!- <dufap>
  MEAN_PART_LOC(70,J) = MEAN_PART_LOC(70,J) + PART(I,J)%DUFAP

!!- <dvfap>
  MEAN_PART_LOC(71,J) = MEAN_PART_LOC(71,J) + PART(I,J)%DVFAP

!!- <dwfap>
  MEAN_PART_LOC(72,J) = MEAN_PART_LOC(72,J) + PART(I,J)%DWFAP

 end do
 

!- Compute the mean velocities over the whole domain
  call RSUMCPU(MEAN_PART_LOC(70,J),DUFAPM)
  call RSUMCPU(MEAN_PART_LOC(71,J),DVFAPM)
  call RSUMCPU(MEAN_PART_LOC(71,J),DWFAPM)

!- Normalization
  DUFAPM = DUFAPM / real(NPART)
  DVFAPM = DVFAPM / real(NPART)
  DWFAPM = DWFAPM / real(NPART)


  do I = 1, NPART_LOC(J)

!!- du'f@p = duf@p-<duf@p>
  DUFAPFLC = PART(I,J)%DUFAP - DUFAPM

!!- dv'f@p = dvf@p-<dvf@p>
  DVFAPFLC = PART(I,J)%DVFAP - DVFAPM

!!- dw'f@p = dwf@p-<dwf@p>
  DWFAPFLC = PART(I,J)%DWFAP - DWFAPM

!!- u'p = up-<up>
  UPFLC = PART(I,J)%UP - UPM

!!- v'p = vp-<vp>
  VPFLC = PART(I,J)%VP - VPM

!!- w'p = wp-<wp>
  WPFLC = PART(I,J)%WP - WPM

!!- u'f@p = uf@p-<uf@p>
  UFAPFLC = PART(I,J)%UFAP - UFAPM

!!- v'f@p = vf@p-<vf@p>
  VFAPFLC = PART(I,J)%VFAP - VFAPM

!!- w'f@p = wf@p-<wf@p>
  WFAPFLC = PART(I,J)%WFAP - WFAPM



!!- SGS-SGS fluid-fluid correlation
  MEAN_PART_LOC(73,J) = MEAN_PART_LOC(73,J) + DUFAPFLC*DUFAPFLC !!- <duf@p.duf@p>
  MEAN_PART_LOC(74,J) = MEAN_PART_LOC(74,J) + DVFAPFLC*DVFAPFLC !!- <dvf@p.dvf@p>
  MEAN_PART_LOC(75,J) = MEAN_PART_LOC(75,J) + DWFAPFLC*DWFAPFLC !!- <dwf@p.dwf@p>
  MEAN_PART_LOC(76,J) = MEAN_PART_LOC(76,J) + DUFAPFLC*DVFAPFLC !!- <duf@p.dvf@p>
  MEAN_PART_LOC(77,J) = MEAN_PART_LOC(77,J) + DUFAPFLC*DWFAPFLC !!- <duf@p.dwf@p>
  MEAN_PART_LOC(78,J) = MEAN_PART_LOC(78,J) + DVFAPFLC*DWFAPFLC !!- <dvf@p.dwf@p>

  MEAN_PART_LOC(79,J) = MEAN_PART_LOC(79,J) &
                      + 0.5*(DUFAPFLC*DUFAPFLC+DVFAPFLC*DVFAPFLC+DWFAPFLC*DWFAPFLC)
  
!!- SGS-Filtered fluid-fluid correlation
  MEAN_PART_LOC(80,J) = MEAN_PART_LOC(80,J) +  DUFAPFLC*UFAPFLC!!- <duf@p.duf@p>
  MEAN_PART_LOC(81,J) = MEAN_PART_LOC(81,J) +  DVFAPFLC*VFAPFLC !!- <dvf@p.dvf@p>
  MEAN_PART_LOC(82,J) = MEAN_PART_LOC(82,J) +  DWFAPFLC*WFAPFLC !!- <dwf@p.dwf@p>

!!- SGS-velocity fluid-particle correlation
  MEAN_PART_LOC(83,J) = MEAN_PART_LOC(83,J) + DUFAPFLC*UPFLC !!- <duf@p.dvf@p>
  MEAN_PART_LOC(84,J) = MEAN_PART_LOC(84,J) + DVFAPFLC*VPFLC !!- <duf@p.dwf@p>
  MEAN_PART_LOC(85,J) = MEAN_PART_LOC(85,J) + DWFAPFLC*WPFLC !!- <dvf@p.dwf@p>

 end do


end if !!- end if (FILTERING)





end do!!- End loop over NIG



!!====================================================================
!! Summation overall domain and normalization
!!====================================================================
!if (NPROC>1) then
do J = 1, NIG

!!- Full number of particle 
 call ISUMCPU(NPART_LOC(J),NPART)

 call MPI_ALLREDUCE(MEAN_PART_LOC(:,J),MEAN_PART(:,J),NSTAT,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!- Normalization by the full number of particle
 MEAN_PART(:,J) = MEAN_PART(:,J) / real(NPART)

end do



do J = 1, NIG
!!--------------------------------------------------------------------
!! 1.4. Print in file
!!--------------------------------------------------------------------
 if(MYID==0) write(500+J,10000)   &
    TIME, MEAN_PART( 1,J),  & !- <up>
          MEAN_PART( 2,J),  & !- <vp>
          MEAN_PART( 3,J),  & !- <wp>
          MEAN_PART(11,J),  & !- <uf@p>
          MEAN_PART(12,J),  & !- <vf@p>
          MEAN_PART(13,J),  & !- <wf@p>
	  MEAN_PART( 4,J),  & !- <up.up>
          MEAN_PART( 5,J),  & !- <vp.vp>
          MEAN_PART( 6,J),  & !- <wp.wp>
          MEAN_PART( 7,J),  & !- <up.vp>
          MEAN_PART( 8,J),  & !- <up.wp>
          MEAN_PART( 9,J),  & !- <vp.wp>
          MEAN_PART(14,J),  & !- <uf@p.uf@p>
          MEAN_PART(15,J),  & !- <vf@p.vf@p>
          MEAN_PART(16,J),  & !- <wf@p.wf@p>
          MEAN_PART(17,J),  & !- <uf@p.vf@p>
          MEAN_PART(18,J),  & !- <uf@p.wf@p>
          MEAN_PART(19,J),  & !- <vf@p.wf@p>
          MEAN_PART(21,J),  & !- <uf@p.up>
          MEAN_PART(22,J),  & !- <vf@p.vp>
          MEAN_PART(23,J),  & !- <wf@p.wp>
          MEAN_PART(24,J),  & !- <uf@p.vp>
          MEAN_PART(25,J),  & !- <uf@p.wp>
          MEAN_PART(26,J)     !- <vf@p.wp>

  

 if(MYID==0) write(520+J,10000) &
    TIME, MEAN_PART(10,J),      & !- qp
          MEAN_PART(27,J),      & !- qfp
          MEAN_PART(20,J),      & !- qf@p
          MEAN_PART(28,J),      & !- <1/tp>
          MEAN_PART(29,J),      & !- <Rep>
          MEAN_PART(30,J),      & !- <Cd>
          MEAN_PART(31,J)         !- <Vr>

!!- Print scalar's particle statistics
 if((MYID==0).and.SOLVE_SCALAR)write(560+J,10000)TIME, MEAN_PART(35:58,J)
 

 end do !!- End loop over NIG


end if



!!======================================================================
!! 4. Lagrangian correlation function
!!======================================================================
!!
!! The normalized autocorrelation function is computed as
!!
!!                 <u(t0).u(t0+t)>
!! R(t) = ---------------------------------
!!        sqrt(<u(t0)^2>).sqrt(<u(t0+t)^2>)
!!
!!====================================================================
!!- The following procedure works only if 
!!  the end of computation is normal

if(LEVEL2_STPAR) then


do NLGR = 1, NBLGRMAX

 
!- Index
  NTLGR = int( (NCYCLE - NT0(NLGR)*FOUT2)/FOUT2) + 1
  
  
  
  
if(NCYCLE >= NT0(NLGR)*FOUT2) then

!!--------------------------------------------------------------------
!! 4.1 Storage Lagrangian fluctuating velocity at t=t0
!!--------------------------------------------------------------------
if(NCYCLE == NT0(NLGR)*FOUT2) then

if(MYID==0) write(*,*)'Lagrangian function computing starts at NCYCLE= ',NT0(NLGR)*FOUT2

do J = 1, NIG

!- Particle velocity 
 PART(1:NPART_LOC(J),J)%UPT0(NLGR) = PART(1:NPART_LOC(J),J)%UP - MEAN_PART(1,J)
 PART(1:NPART_LOC(J),J)%VPT0(NLGR) = PART(1:NPART_LOC(J),J)%VP - MEAN_PART(2,J)
 PART(1:NPART_LOC(J),J)%WPT0(NLGR) = PART(1:NPART_LOC(J),J)%WP - MEAN_PART(3,J)

!- Fluid velocity at particle position
 PART(1:NPART_LOC(J),J)%UFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%UFAP - MEAN_PART(11,J)
 PART(1:NPART_LOC(J),J)%VFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%VFAP - MEAN_PART(12,J)
 PART(1:NPART_LOC(J),J)%WFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%WFAP - MEAN_PART(13,J)

 if(SOLVE_SCALAR) then
  PART(1:NPART_LOC(J),J)%TPT0(NLGR)   = PART(1:NPART_LOC(J),J)%TP   - MEAN_PART(35,J)
  PART(1:NPART_LOC(J),J)%TFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%TFAP - MEAN_PART(36,J)
 end if !-  SOLVE_SCALAR
 
 if(FILTERING) then
  PART(1:NPART_LOC(J),J)%DUFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%DUFAP - MEAN_PART(70,J)
  PART(1:NPART_LOC(J),J)%DVFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%DVFAP - MEAN_PART(71,J)
  PART(1:NPART_LOC(J),J)%DWFAPT0(NLGR) = PART(1:NPART_LOC(J),J)%DWFAP - MEAN_PART(72,J)
 end if
 

end do !->  J = 1, NIG

end if !-> (NCYCLE == NT0*FOUT2)





!!--------------------------------------------------------------------
!! 4.2 correlation at time t
!!--------------------------------------------------------------------
!if(NTLGR<=DIMLGR) then

do J = 1, NIG

do I = 1, NPART_LOC(J)

!!-------------------------------------------------------------------
!! 4.1 Particle-particle
!!-------------------------------------------------------------------
!- <up(t0).up(t0+t)>
 RPX_LOC(NTLGR,J,NLGR) = RPX_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%UP - MEAN_PART(1,J))*PART(I,J)%UPT0(NLGR)/ real(NPART_LOC(J))

!- <vp(t0).vp(t0+t)>
 RPY_LOC(NTLGR,J,NLGR) = RPY_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%VP - MEAN_PART(2,J))*PART(I,J)%VPT0(NLGR)/ real(NPART_LOC(J))

!- <wp(t0).wp(t0+t)>
 RPZ_LOC(NTLGR,J,NLGR) = RPZ_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%WP - MEAN_PART(3,J))*PART(I,J)%WPT0(NLGR)/ real(NPART_LOC(J))


!!-------------------------------------------------------------------
!! 4.2 Fluid-particle
!!-------------------------------------------------------------------
!- <ufap(t0).ufap(t0+t)>
 RFAPX_LOC(NTLGR,J,NLGR) = RFAPX_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%UFAP-MEAN_PART(11,J))*PART(I,J)%UFAPT0(NLGR) / real(NPART_LOC(J))

!- <vfap(t0).vfap(t0+t)>
 RFAPY_LOC(NTLGR,J,NLGR) = RFAPY_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%VFAP-MEAN_PART(12,J))*PART(I,J)%VFAPT0(NLGR) / real(NPART_LOC(J))

!- <wfap(t0).wfap(t0+t)>
 RFAPZ_LOC(NTLGR,J,NLGR) = RFAPZ_LOC(NTLGR,J,NLGR) &
              + (PART(I,J)%WFAP-MEAN_PART(13,J))*PART(I,J)%WFAPT0(NLGR) / real(NPART_LOC(J))


!!-------------------------------------------------------------------
!! 4.3 Temperature
!!-------------------------------------------------------------------
 if(SOLVE_SCALAR) then
!- <Tp(t0).Tp(t0+t)>
  RTP_LOC(NTLGR,J,NLGR) = RTP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%TP-MEAN_PART(35,J))*PART(I,J)%TPT0(NLGR) / real(NPART_LOC(J))
 
!- <Tf@p(t0).Tf@p(t0+t)>
  RTFAP_LOC(NTLGR,J,NLGR) = RTFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%TFAP-MEAN_PART(36,J))*PART(I,J)%TFAPT0(NLGR) / real(NPART_LOC(J))

 
!- <vf@p(t0).Tf@p(t0+t)>
  RVFAPTFAP_LOC(NTLGR,J,NLGR) = RVFAPTFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%VFAP-MEAN_PART(12,J))*PART(I,J)%TFAPT0(NLGR) / real(NPART_LOC(J))

!- <Tf@p(t0).vf@p(t0+t)>
  RTFAPVFAP_LOC(NTLGR,J,NLGR) = RTFAPVFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%TFAP-MEAN_PART(36,J))*PART(I,J)%VFAPT0(NLGR) / real(NPART_LOC(J))

 end if !-> SOLVE_SCALAR


!!-------------------------------------------------------------------
!! 4.4 Subgrid turbulence
!!-------------------------------------------------------------------
 if(FILTERING) then
!- <dufap(t0).dufap(t0+t)>
  RDUFAPDUFAP_LOC(NTLGR,J,NLGR) = RDUFAPDUFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DUFAP-MEAN_PART(70,J))*PART(I,J)%DUFAPT0(NLGR) / real(NPART_LOC(J))

!- <dvfap(t0).dvfap(t0+t)>
  RDVFAPDVFAP_LOC(NTLGR,J,NLGR) = RDVFAPDVFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DVFAP-MEAN_PART(71,J))*PART(I,J)%DVFAPT0(NLGR) / real(NPART_LOC(J))
 
!- <dwfap(t0).dwfap(t0+t)>
  RDWFAPDWFAP_LOC(NTLGR,J,NLGR) = RDWFAPDWFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DWFAP-MEAN_PART(72,J))*PART(I,J)%DWFAPT0(NLGR) / real(NPART_LOC(J))


!- <dufap(t0).ufap(t0+t)>
  RDUFAPUFAP_LOC(NTLGR,J,NLGR) = RDUFAPUFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DUFAP-MEAN_PART(70,J))*PART(I,J)%UFAPT0(NLGR) / real(NPART_LOC(J))

!- <dvfap(t0).vfap(t0+t)>
  RDVFAPVFAP_LOC(NTLGR,J,NLGR) = RDVFAPVFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DVFAP-MEAN_PART(71,J))*PART(I,J)%VFAPT0(NLGR) / real(NPART_LOC(J))
 
!- <dwfap(t0).wfap(t0+t)>
  RDWFAPWFAP_LOC(NTLGR,J,NLGR) = RDWFAPWFAP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DWFAP-MEAN_PART(72,J))*PART(I,J)%WFAPT0(NLGR) / real(NPART_LOC(J))

!- <dufap(t0).up(t0+t)>
  RDUFAPUP_LOC(NTLGR,J,NLGR) = RDUFAPUP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DUFAP-MEAN_PART(70,J))*PART(I,J)%UPT0(NLGR) / real(NPART_LOC(J))

!- <dvfap(t0).vp(t0+t)>
  RDVFAPVP_LOC(NTLGR,J,NLGR) = RDVFAPVP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DVFAP-MEAN_PART(71,J))*PART(I,J)%VPT0(NLGR) / real(NPART_LOC(J))
 
!- <dwfap(t0).wp(t0+t)>
  RDWFAPWP_LOC(NTLGR,J,NLGR) = RDWFAPWP_LOC(NTLGR,J,NLGR) &
               + (PART(I,J)%DWFAP-MEAN_PART(72,J))*PART(I,J)%WPT0(NLGR) / real(NPART_LOC(J))

 end if !-> FILTERING




end do !-> I = 1, NPART_LOC(J)


end do !->  J = 1, NIG



end if !!- end if(NCYCLE >= NT0(NLGR)*FOUT2)


end do !!- end do NLGR=1, NBLGRMAX


end if !!-> (LEVEL2_STPAR)


!!======================================================================
!! 5. Spatial correlation function
!!======================================================================
if(LEVEL3_STPAR) then

DXCP = DX*RATIOCP
DYCP = DY*RATIOCP
DZCP = DZ*RATIOCP

!if(MYID==0) then
! write(*,*)'DXCP=',DXCP,' NXCP=',NXCP
! write(*,*)'DYCP=',DYCP,' NYCP=',NYCP
! write(*,*)'DZCP=',DZCP,' NZCP=',NZCP
!end if

!if(MYID==0) then
! write(*,*)'CMEAN=',CMEAN
!end if

do J = 1, NIG

!!----------------------------------------------------------------------
!! 5.1 Concentration of particle
!!----------------------------------------------------------------------
CONCP(:,:,:) = ZERO
do I = 1, NPART_LOC(J)
 IPART = int((PART(I,J)%XP-XMESH(ISTART(1)))/DXCP) + 1
 JPART = int((PART(I,J)%YP-YMESH(ISTART(2)))/DYCP) + 1
 KPART = int((PART(I,J)%ZP-ZMESH(ISTART(3)))/DZCP) + 1
 
 CONCP(IPART,JPART,KPART) = CONCP(IPART,JPART,KPART) + 1.0
end do

!!----------------------------------------------------------------------
!! 5.2 Statistic on concentration
!!----------------------------------------------------------------------
!- Mean concentration
 MEAN_PART_LOC(61,J) = sum(CONCP(:,:,:)) / NXCP / NYCP / NZCP
 
 
 call RSUMCPU(MEAN_PART_LOC(61,J),CMEAN)
 CMEAN = CMEAN / real(NPROC)

 
!- Mean variance
 MEAN_PART_LOC(62,J) = sum((CONCP(:,:,:)-CMEAN)**2)/ NXCP / NYCP / NZCP

!- Mean squared concentration
 MEAN_PART_LOC(63,J) = sum(CONCP(:,:,:)**2)/ NXCP / NYCP / NZCP


 call RSUMCPU(MEAN_PART_LOC(61,J),MEAN_PART(61,J))
 call RSUMCPU(MEAN_PART_LOC(62,J),MEAN_PART(62,J))
 call RSUMCPU(MEAN_PART_LOC(63,J),MEAN_PART(63,J))


 MEAN_PART(61,J) =  MEAN_PART(61,J) / real(NPROC)
 MEAN_PART(62,J) =  MEAN_PART(62,J) / real(NPROC)
 MEAN_PART(63,J) =  MEAN_PART(63,J) / real(NPROC)


!if(MYID==0) then
! write(*,*)'Class: ',J,'CMEAN=',CMEAN
!end if


!!----------------------------------------------------------------------
!! 5.3 PDF of particle concentration
!!----------------------------------------------------------------------
do KPART = 1, NZCP
 do JPART = 1, NYCP
  do IPART = 1, NXCP
   IPDF = int(CONCP(IPART,JPART,KPART))  + 1
   if(IPDF>NPDFCP) IPDF = NPDFCP
   PDFCP_LOC(IPDF,J) = PDFCP_LOC(IPDF,J) + 1.0
  end do
 end do
end do





end do

end if


!!======================================================================
!! 6. Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_PART = MEAN_TIME_PART + MEAN_PART
end if 



!!======================================================================
!! 7. Print statistics in "info" files
!!======================================================================
!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
 CPU_PART(6) = CPU_PART(6) + TIME_END - TIME_START
end if

!!----------------------------------------------------------------------
10000 format (30(e17.7))
10601 format (2x,A,E13.6)
10602 format (2x,A,E13.6,1x,A,E13.6)
10603 format (2x,A,E13.6,2x,A,E13.6,2x,A,E13.6)

end subroutine STAT_PARTICLE
