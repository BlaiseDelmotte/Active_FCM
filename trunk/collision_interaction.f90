!!====================================================================
!!
!!          Particle-particle interaction
!!
!!====================================================================

subroutine COLLISION_INTERACTION(TIME,NCYCLE,IG,NCLOSE,HOME,NEAR)

!!====================================================================
!!
!!
!!====================================================================

use PARTICLE_PARALLEL
use DNS_DIM               
use PARAM_PHYS
use GEOMETRIC_VARIABLE
use COLLISION_VARIABLE

implicit none

!!====================================================================
!! ARRAYS
!!====================================================================
!- Curent time
real(kind=8),               intent(in) :: TIME

!!- Cycle number
integer,                    intent(in) :: NCYCLE

!!- particle specy
integer,                    intent(in) :: IG

!!- Number of particle "close"
integer,                    intent(in) :: NCLOSE

!!- Index of particle "home"
integer, dimension(NCLOSE), intent(in) :: HOME

!!- Index of particle "near"
integer, dimension(NCLOSE), intent(in) :: NEAR

!!====================================================================
!- Impact parameter
real(kind=8) :: VKX, VKY, VKZ
real(kind=8) :: NRM_VK

!- Relative velocity
real(kind=8) :: WX, WY, WZ
real(kind=8) :: NRM_W

!- Scalar product between k and w
real(kind=8) :: WK
real(kind=8) :: THETA, THETAM, THETA2M

real(kind=8):: NCOL_LOC, NOVER_LOC, MEANOVER_LOC
real(kind=8):: NCOL,     NOVER,     MEANOVER
real(kind=8):: DQPART
real(kind=8):: DQPART_LOC

real(kind=8) :: TCOL, FCOL

real(kind=8):: M1, M2, EC, LI, BI, DELTAT

!!--------------------------------------------------------------------
integer :: NCYCLESTAT

!- Statistic pdf mesh
real(kind=8), dimension(NPDFMAX,5) :: PDF_TH
real(kind=8), dimension(        5) :: AREA_PDF
real(kind=8), dimension(        5) :: DPDF
integer :: IPDF
real(kind=8) :: X
real, dimension(15) :: MOMENT

!!- Index
integer :: N, NP1, NP2
!!====================================================================



EC = 1.0

!!- Initiation of colision number
NCOL_LOC = ZERO

!!- Number of overlapped particles
NOVER_LOC = ZERO

MEANOVER_LOC = ZERO

DQPART_LOC = ZERO

THETAM = ZERO
THETA2M = ZERO


!!- Cycle de debut de stat
NCYCLESTAT = NCYCLEMAX /2


!!-----------------------!!-----------------------
DPDF(1) = .5*PPI / real(NPDFMAX)
DPDF(2) =  3.    / real(NPDFMAX)
DPDF(3) =  3.    / real(NPDFMAX)
if(NCYCLE<NCYCLESTAT) then
 PDF(:,:) = ZERO
 MOMENT(:) = ZERO
 NMOM = ZERO
end if
!!-----------------------!!-----------------------




!- Direct method
do N = 1, NCLOSE

 NP1 = HOME(N)
 NP2 = NEAR(N)


!!====================================================================
!! 1. Compute the impact vector
!!====================================================================
!! The impact vector is the vector linkig the centres of the two
!! neighboring particles.
!!--------------------------------------------------------------------

!- Impact parameter
 VKX = PART(NP1,IG)%XP -  PART(NP2,IG)%XP
 VKY = PART(NP1,IG)%YP -  PART(NP2,IG)%YP
 VKZ = PART(NP1,IG)%ZP -  PART(NP2,IG)%ZP


!- Periodicity only in x-direction
 if(abs(VKX) > LXMAX/2.) VKX = VKX - SIGN(LXMAX,VKX)

!- Periodicity in y- and z-direction only if scalar computation
 if(NPROC==1) then
   
!- Periodicity only in y-direction
  if(abs(VKY) > LYMAX/2.) VKY = VKY - SIGN(LYMAX,VKY)

!- Periodicity only in z-direction
  if(abs(VKZ) > LZMAX/2.) VKZ = VKZ - SIGN(LZMAX,VKZ)

 end if



!- modulus of impact parameter
 NRM_VK = sqrt(VKX*VKX + VKY*VKY + VKZ*VKZ)

!write(*,*)'NRM_VK=',NRM_VK
!write(*,*)'    dp=',DPART(IG)

!!====================================================================
!! 2. Overlapped particles
!!====================================================================
 if(NRM_VK <= DPART(IG)) then

   NOVER_LOC = NOVER_LOC + 1.0

   MEANOVER_LOC = MEANOVER_LOC + NRM_VK

!- particle masses
   M1 = RHOP(IG)*PPI*DPART(IG)**3/6.D0
   M2 = RHOP(IG)*PPI*DPART(IG)**3/6.D0

!- Relative velocity
   WX = PART(NP2,IG)%UP - PART(NP1,IG)%UP
   WY = PART(NP2,IG)%VP - PART(NP1,IG)%VP
   WZ = PART(NP2,IG)%WP - PART(NP1,IG)%WP

   WK = VKX*WX + VKY*WY + VKZ*WZ

   NRM_W = sqrt(WX*WX + WY*WY + WZ*WZ)

!write(*,*)'Overlapping:'
!write(*,*)' |k|/dp =',NRM_VK/DPART(IG)
!write(*,*)'     kx =',VKX
!write(*,*)'     ky =',VKY
!write(*,*)'     kz =',VKZ
!write(*,*)'     xh =',PART(NP1,IG)%XP
!write(*,*)'     yh =',PART(NP1,IG)%YP
!write(*,*)'     zh =',PART(NP1,IG)%ZP
!write(*,*)'     xn =',PART(NP2,IG)%XP
!write(*,*)'     yn =',PART(NP2,IG)%YP
!write(*,*)'     zn =',PART(NP2,IG)%ZP
!write(*,*)'     wk =',WK
!pause

   if(WK>ZERO) then ! Particles are approaching !!


!write(400+MYID,*)'Collision:'
!write(400+MYID,*)' NP1:',NP1
!write(400+MYID,*)' NP2:',NP2
!write(400+MYID,*)' Before particle displacement'
!write(400+MYID,*)' |k|/dp =',NRM_VK/DPART(IG)
!write(400+MYID,*)'     kx =',VKX
!write(400+MYID,*)'     ky =',VKY
!write(400+MYID,*)'     kz =',VKZ
!write(400+MYID,*)'     zh =',PART(NP1,IG)%ZP
!write(400+MYID,*)'     zn =',PART(NP2,IG)%ZP
!write(400+MYID,*)'     wk =',WK

      NCOL_LOC = NCOL_LOC + 1.0

!!====================================================================
!! 3. Computation of the impact vector for contacting particle
!!====================================================================
!! Here the impact vector is re-computed for contacting particles
!! instead of overlapping particles. Geometric rules are used.
!!--------------------------------------------------------------------
      LI = abs(WK) / NRM_W
       
      BI = NRM_VK*NRM_VK - LI*LI
       
      if(BI<0.) then
        write(*,*)'Pb in particle displacement'
	write(*,*)'NP1 =', NP1,' NP2=',NP2
	write(*,*)'NRM_VK=',NRM_VK,' LI=',LI,' BI=',BI
      end if
       
      DELTAT = (sqrt(DPART(IG)**2-BI)-LI) / NRM_W

!- New impact vector
      VKX = VKX + WX*DELTAT
      VKY = VKY + WY*DELTAT
      VKZ = VKZ + WZ*DELTAT

!- New scalar product 
      WK = VKX*WX + VKY*WY + VKZ*WZ

!- New modulus of k
      NRM_VK = sqrt(VKX*VKX + VKY*VKY + VKZ*VKZ)
      
!- Normalization
      VKX = VKX / NRM_VK
      VKY = VKY / NRM_VK
      VKZ = VKZ / NRM_VK

      WK = WK / NRM_VK


      THETA = acos(WK/NRM_W)

      THETAM = THETAM + THETA
      THETA2M = THETA2M + THETA**2


if(NCYCLE >=NCYCLESTAT) then
!- Statistics on colliding particle

 IPDF = int(THETA/DPDF(1)) + 1  

 if(IPDF > NPDFMAX) IPDF = NPDFMAX
 if(IPDF < 1      ) IPDF = 1
 PDF(IPDF,1) = PDF(IPDF,1) + 1.0

 IPDF = int(WK/DPDF(2)) + 1
 if(IPDF > NPDFMAX) IPDF = NPDFMAX
 if(IPDF < 1      ) IPDF = 1
 PDF(IPDF,2) = PDF(IPDF,2) + 1.0

 IPDF = int(NRM_W/DPDF(3)) + 1
 if(IPDF > NPDFMAX) IPDF = NPDFMAX
 if(IPDF < 1      ) IPDF = 1
 PDF(IPDF,3) = PDF(IPDF,3) + 1.0

end if


!!====================================================================
!!4. Updating of particle velocities
!!====================================================================


!write(400+MYID,*)'Collision:'
!write(400+MYID,*)' Before particle displacement'
!write(400+MYID,*)' |k|/dp =',NRM_VK/DPART(IG)
!write(400+MYID,*)'     kx =',VKX
!write(400+MYID,*)'     ky =',VKY
!write(400+MYID,*)'     kz =',VKZ
!write(400+MYID,*)'     zh =',PART(NP1,IG)%ZP
!write(400+MYID,*)'     zn =',PART(NP2,IG)%ZP
!write(400+MYID,*)'     wk =',WK
!write(400+MYID,*)'Velocities before Collision:'
!write(400+MYID,*)'     u1 =',PART(NP1,IG)%UP
!write(400+MYID,*)'     v1 =',PART(NP1,IG)%VP
!write(400+MYID,*)'     w1 =',PART(NP1,IG)%WP
!write(400+MYID,*)'     u2 =',PART(NP2,IG)%UP
!write(400+MYID,*)'     v2 =',PART(NP2,IG)%VP
!write(400+MYID,*)'     w2 =',PART(NP2,IG)%WP

      DQPART_LOC = DQPART_LOC &
              + PART(NP1,IG)%UP**2+PART(NP1,IG)%VP**2+PART(NP1,IG)%WP**2 &
              + PART(NP2,IG)%UP**2+PART(NP2,IG)%VP**2+PART(NP2,IG)%WP**2


      PART(NP1,IG)%COLOR = PART(NP1,IG)%COLOR + 1.0
      PART(NP2,IG)%COLOR = PART(NP2,IG)%COLOR + 1.0

!- Hard sphere model
      PART(NP1,IG)%UP = PART(NP1,IG)%UP + M2/(M1+M2)*(1.0+EC)*WK*VKX
      PART(NP1,IG)%VP = PART(NP1,IG)%VP + M2/(M1+M2)*(1.0+EC)*WK*VKY
      PART(NP1,IG)%WP = PART(NP1,IG)%WP + M2/(M1+M2)*(1.0+EC)*WK*VKZ

      PART(NP2,IG)%UP = PART(NP2,IG)%UP - M1/(M1+M2)*(1.0+EC)*WK*VKX
      PART(NP2,IG)%VP = PART(NP2,IG)%VP - M1/(M1+M2)*(1.0+EC)*WK*VKY
      PART(NP2,IG)%WP = PART(NP2,IG)%WP - M1/(M1+M2)*(1.0+EC)*WK*VKZ


!write(400+MYID,*)'Velocities after Collision:'
!write(400+MYID,*)'     u1 =',PART(NP1,IG)%UP
!write(400+MYID,*)'     v1 =',PART(NP1,IG)%VP
!write(400+MYID,*)'     w1 =',PART(NP1,IG)%WP
!write(400+MYID,*)'     u2 =',PART(NP2,IG)%UP
!write(400+MYID,*)'     v2 =',PART(NP2,IG)%VP
!write(400+MYID,*)'     w2 =',PART(NP2,IG)%WP

      !- Energy conservation
      DQPART_LOC = DQPART_LOC &
             - (PART(NP1,IG)%UP**2+PART(NP1,IG)%VP**2+PART(NP1,IG)%WP**2) &
             - (PART(NP2,IG)%UP**2+PART(NP2,IG)%VP**2+PART(NP2,IG)%WP**2)

    end if !- end if( WK>0.)

  end if !- end NRM_VK < DPART(IG))

end do !- end loop do N = 1, NCLOSE




!!- Statistics for the whole domain
call RSUMCPU(NOVER_LOC,NOVER)
call RSUMCPU( NCOL_LOC, NCOL)

if(NOVER_LOC>ZERO) then
 MEANOVER_LOC = MEANOVER_LOC/NOVER_LOC
else
 MEANOVER_LOC = ZERO
end if
call RSUMCPU( MEANOVER_LOC, MEANOVER)


if(NCOL_LOC>ZERO) then
 DQPART_LOC = DQPART_LOC/NCOL_LOC
 THETAM = THETAM/NCOL_LOC 
 THETA2M = THETA2M /NCOL_LOC
else
 DQPART_LOC = ZERO
 THETAM = ZERO
 THETA2M = ZERO
end if
call RSUMCPU(DQPART_LOC,DQPART)

!!- Collision time scale
FCOL = NCOL*2./DTIME/NPART_FULL


if(MYID==0) write(580+IG,10000) &
                         TIME,  &
                 real(NCLOSE),  &
                        NOVER,  &
                         NCOL,  &
                         FCOL,  &
                       DQPART,  &
         1.-MEANOVER/DPART(IG), &
                          THETAM, THETA2M


if(NCOL>ZERO) then
write(*,*)' Overlaping particles=',NOVER
write(*,*)' Collision  particles=',NCOL
write(*,*)'  Mean overlap (% dp)=',100.*(1.-MEANOVER/DPART(IG))
write(*,*)' Modulation of particle agitation: <dqp>=',DQPART
end if
!pause





if(NCYCLE>=NCYCLESTAT) then
MOMENT(1) = MOMENT(1) + FCOL
do N = 1, NPART_LOC(IG)
 MOMENT(2) = MOMENT(2) + 0.5*(PART(N,IG)%UP**2+PART(N,IG)%VP**2+PART(N,IG)%WP**2)
end do
NMOM = NMOM + 1.0
end if




if(NCYCLE==NCYCLEMAX-1) then

MOMENT(:) = MOMENT(:) / NMOM/ NPART_LOC(IG)

!- Initiation
AREA_PDF(:) = ZERO

!- Surface of the pdf
do N = 1, NPDFMAX
 AREA_PDF(1) = AREA_PDF(1) + PDF(N,1)*DPDF(1) 
 AREA_PDF(2) = AREA_PDF(2) + PDF(N,2)*DPDF(2)  
 AREA_PDF(3) = AREA_PDF(3) + PDF(N,3)*DPDF(3)  
end do
!- Normalization
PDF(:,1) = PDF(:,1) / AREA_PDF(1)
PDF(:,2) = PDF(:,2) / AREA_PDF(2)
PDF(:,3) = PDF(:,3) / AREA_PDF(3)


!- Theoretical pdf
do N = 1, NPDFMAX
  !- Colliding angle
  X = (N - 0.5)*DPDF(1)
  PDF_TH(N,1) = sin(2.*X)

  !- Radial relative velocity
  X = (N - 0.5)*DPDF(2)
  PDF_TH(N,2) = 1./(2.*(4.*MOMENT(2)/3.)**2)*X**3*exp(-X**2/(8.*MOMENT(2)/3.))

  !- Radial relative velocity
  X = (N - 0.5)*DPDF(3)
  PDF_TH(N,3) = 1./(2.*(4.*MOMENT(2)/3.)**2)*X**3*exp(-X**2/(8.*MOMENT(2)/3.))

end do

do N = 1, NPDFMAX
 write(800,10000) (N-0.5)*DPDF(1)/PPI, PDF(N,1), PDF_TH(N,1), &
                  (N-0.5)*DPDF(2)    , PDF(N,2), PDF_TH(N,2), &
                  (N-0.5)*DPDF(3)    , PDF(N,3), PDF_TH(N,3)
end do



write(*,*)' tc =',1/MOMENT(1)
write(*,*)'qp2 =',MOMENT(2)
write(*,*)'NMOM=',NMOM

end if


!!--------------------------------------------------------------------
10000 format (30(e17.7))
10703 format (2x,' Part. in Class ',i3,' :  Max=',i7,' Min=',i7)

end subroutine COLLISION_INTERACTION

