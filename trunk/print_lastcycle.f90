!!====================================================================
!!
!!
!!====================================================================

subroutine PRINT_LASTCYCLE(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================
use DNS_DIM	       !- Dimension
use STATISTICS         !- Statistics
use PARAM_PHYS

implicit none

!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
integer, intent(in) :: NCYCLE

!!--------------------------------------------------------------------
!!- Local arrays
!!--------------------------------------------------------------------

!- Index
integer :: I, J, LGR, M ,N
!---------------------------------------------------------------------




!!====================================================================
!! Print in "lastcycle.info"
!!====================================================================


open(unit=200, file='lastcycle.info')

write(200,*)
write(200,*)'====================================================================='
write(200,*)'LAST CYCLE STATISTICS '
write(200,*)'====================================================================='
write(200,*)
write(200,10600)'Cycle Number =',real(NCYCLE)
write(200,*)
if(LEVEL0_STFLU) then
write(200,*)'---------------------------------------------------------------------'
write(200,*)' FLUID'
write(200,*)'---------------------------------------------------------------------'
write(200,*)
end if  !-> LEVEL0_STFLU
if(LEVEL1_STFLU) then
write(200,*)
write(200,*)'Mean fluid Velocity'
write(200,*)'-------------------'
write(200,10601)'<uf>f = ',MEAN_FLUID(1),' [m/s]'
write(200,10601)'<vf>f = ',MEAN_FLUID(2),' [m/s]'
write(200,10601)'<wf>f = ',MEAN_FLUID(3),' [m/s]'
write(200,*)
write(200,*)'Reynolds stress tensor'
write(200,*)'----------------------'
write(200,10601)'<uf.uf>f = ',MEAN_FLUID(4),' [m2/s2]'
write(200,10601)'<vf.vf>f = ',MEAN_FLUID(5),' [m2/s2]'
write(200,10601)'<wf.wf>f = ',MEAN_FLUID(6),' [m2/s2]'
write(200,*)
write(200,10601)'    <kf> = ',MEAN_FLUID(10),' [m2/s2]'
write(200,10601)'  <epsf> = ',MEAN_FLUID(11),' [m2/s3]'
write(200,*)
write(200,*    )'   u_rms=(2/3*kf)**0.5'
write(200,10601)'   u_rms = ',sqrt(2./3.*MEAN_FLUID(10)),' [m/s]'
write(200,*)
write(200,10601)'<uf.vf>f = ',MEAN_FLUID(7),' [m2/s2]'
write(200,10601)'<uf.wf>f = ',MEAN_FLUID(8),' [m2/s2]'
write(200,10601)'<vf.wf>f = ',MEAN_FLUID(9),' [m2/s2]'
write(200,*)
write(200,10601)'   eta_K = ',MEAN_FLUID(11),' [m]'
write(200,10601)'   tau_K = ',MEAN_FLUID(12),' [s]'
write(200,10601)'     v_K = ',MEAN_FLUID(17),' [m/s]'
write(200,*)
write(200,10601)'Taylor      lg = ',MEAN_FLUID(25),' [m]'
write(200,*)
write(200,*)'Reynolds number'
write(200,*)'---------------'
write(200,10602)' Re_Tay = ',sqrt(2.*MEAN_FLUID(10)/3.)*MEAN_FLUID(25)/VISC
write(200,*)
write(200,*)
end if  !-> LEVEL1_STFLU

if(LEVEL2_STFLU) then
write(200,*)
write(200,*)'Divergence'
write(200,*)'----------'
write(200,10601)'   < div(U)  > = ',MEAN_FLUID(15),' [1/s]'
write(200,10601)'   < div(U)^2> = ',MEAN_FLUID(16),' [1/s2]'
write(200,10601)'max(|div(U)| ) = ',MEAN_FLUID(14),' [1/s]'
write(200,*)
end if  !-> LEVEL2_STFLU

if(LEVEL3_STFLU) then
write(200,*)
write(200,*)'Velocity gradients'
write(200,*)'------------------'
write(200,10601)' <(dui/dxi)^2> = ',MEAN_FLUID(21),' [1/s2]'
write(200,10601)' <(dui/dxj)^2> = ',MEAN_FLUID(31),' [1/s2]'
write(200,*)
write(200,10602)'Skewness:   Sk = ',MEAN_FLUID(18)
write(200,10602)'Flatness:   Tk = ',MEAN_FLUID(19)
write(200,*)
write(200,10601)'Taylor      lf = ',MEAN_FLUID(24),' [m]'
write(200,10601)'Taylor      lg = ',MEAN_FLUID(25),' [m]'
write(200,10602)'   lf/lg/2^0.5 = ',MEAN_FLUID(24)/MEAN_FLUID(25)/sqrt(2.)
write(200,*)
end if  !-> LEVEL3_STFLU
write(200,*)
write(200,*)



if(LEVEL0_STSCL) then
write(200,*)'---------------------------------------------------------------------'
write(200,*)' SCALAR'
write(200,*)'---------------------------------------------------------------------'
write(200,*)
end if !-> LEVEL0_STSCL

if(LEVEL1_STSCL) then
write(200,*)
write(200,10601)'  <theta>f = ',MEAN_SCL(1),' [K]'
write(200,10601)'<theta^2>f = ',MEAN_SCL(2),' [K2]'
write(200,*)
write(200,10601)'<uf.theta>f = ',MEAN_SCL(3),' [K.m/s]'
write(200,10601)'<vf.theta>f = ',MEAN_SCL(4),' [K.m/s]'
write(200,10601)'<wf.theta>f = ',MEAN_SCL(5),' [K.m/s]'
end if !-> LEVEL1_STSCL

if(LEVEL2_STSCL) then
write(200,*)
write(200,10601)' Scalar dissipation = ',MEAN_SCL(18),' [K2/s]'
write(200,10601)'         production = ',MEAN_SCL(19),' [K2/s]'
write(200,*)
write(200,10601)'       Skewness: Sx = ',MEAN_SCL(20), ' [-]'
write(200,10601)'                 Sy = ',MEAN_SCL(21), ' [-]'
write(200,10601)'                 Sz = ',MEAN_SCL(22), ' [-]'
write(200,*)
write(200,10601)'       Flatness: Kx = ',MEAN_SCL(23), ' [-]'
write(200,10601)'                 Ky = ',MEAN_SCL(24), ' [-]'
write(200,10601)'                 Kz = ',MEAN_SCL(25), ' [-]'
end if !-> LEVEL2_STSCL



if(LEVEL0_STPAR) then
write(200,*)'---------------------------------------------------------------------'
write(200,*)' PARTICLE'
write(200,*)'---------------------------------------------------------------------'
write(200,*)
end if !->

if(LEVEL1_STPAR) then
do J = 1, NIG
write(200,*)'---------------------------------------------------------------------'
if(PARTDEF(J) == 0) then
write(200,10603)'Particle:',J,'  ==> Motionless'
elseif(PARTDEF(J) == 1) then
write(200,10603)'Particle:',J,'  ==> Fluid element'
elseif(PARTDEF(J) == 2) then
write(200,10603)'Particle:',J,'  ==> Solid particle'
end if
write(200,*)'---------------------------------------------------------------------'
write(200,*)
write(200,10601)'Diameter=',DPART(J),'  [m]'
write(200,10601)'Density=',RHOP(J),'  [kg/m3]'
write(200,*)
write(200,*)'Mean Velocity'
write(200,*)'-------------'
write(200,10601)'<up>p = ',MEAN_PART(1,J),' [m/s]'
write(200,10601)'<vp>p = ',MEAN_PART(2,J),' [m/s]'
write(200,10601)'<wp>p = ',MEAN_PART(3,J),' [m/s]'
write(200,*)
write(200,10601)'<uf>p = ',MEAN_PART(11,J),' [m/s]'
write(200,10601)'<vf>p = ',MEAN_PART(12,J),' [m/s]'
write(200,10601)'<wf>p = ',MEAN_PART(13,J),' [m/s]'
write(200,*)
write(200,*)'Fluctuating motion'
write(200,*)'------------------'
write(200,10601)'  qp = ',MEAN_PART(10,J),' [m2/s2]'
write(200,10601)'qf@p = ',MEAN_PART(20,J),' [m2/s2]'
write(200,10601)' qfp = ',MEAN_PART(27,J),' [m2/s2]'
write(200,*)
write(200,*)
write(200,*)'Particle-Particle kinetic stress tensor'
write(200,*)'--------------------------------------'
write(200,10601)'<upup>p = ',MEAN_PART(4,J),' [m2/s2]'
write(200,10601)'<vpvp>p = ',MEAN_PART(5,J),' [m2/s2]'
write(200,10601)'<wpwp>p = ',MEAN_PART(6,J),' [m2/s2]'
write(200,*)
write(200,10601)'<upvp>p = ',MEAN_PART(7,J),' [m2/s2]'
write(200,10601)'<upwp>p = ',MEAN_PART(8,J),' [m2/s2]'
write(200,10601)'<vpwp>p = ',MEAN_PART(9,J),' [m2/s2]'
write(200,*)
write(200,10601)'<upup>p/(2/3qp) = ',&
               MEAN_PART(4,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,10601)'<vpvp>p/(2/3qp) = ',&
               MEAN_PART(5,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,10601)'<wpwp>p/(2/3qp) = ',&
               MEAN_PART(6,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,*)
write(200,10601)'<upvp>p/(2/3qp) = ',&
               MEAN_PART(7,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,10601)'<upwp>p/(2/3qp) = ',&
               MEAN_PART(8,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,10601)'<vpwp>p/(2/3qp) = ',&
               MEAN_PART(9,J)/(2./3.*MEAN_PART(10,J)),' [-]'
write(200,*)
write(200,*)
write(200,*)'Fluid-fluid @p kinetic stress tensor'
write(200,*)'------------------------------------'
write(200,10601)'<ufuf>p = ',MEAN_PART(14,J),' [m2/s2]'
write(200,10601)'<vfvf>p = ',MEAN_PART(15,J),' [m2/s2]'
write(200,10601)'<wfwf>p = ',MEAN_PART(16,J),' [m2/s2]'
write(200,*)
write(200,10601)'<ufvf>p = ',MEAN_PART(17,J),' [m2/s2]'
write(200,10601)'<ufwf>p = ',MEAN_PART(18,J),' [m2/s2]'
write(200,10601)'<vfwf>p = ',MEAN_PART(19,J),' [m2/s2]'
write(200,*)
write(200,10601)'<ufuf>p/(2/3qf@p) = ',&
               MEAN_PART(14,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,10601)'<vfvf>p/(2/3qf@p) = ',&
               MEAN_PART(15,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,10601)'<wfwf>p/(2/3qf@p) = ',&
               MEAN_PART(16,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,*)
write(200,10601)'<ufvf>p/(2/3qf@p) = ',&
               MEAN_PART(17,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,10601)'<ufwf>p/(2/3qf@p) = ',&
               MEAN_PART(18,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,10601)'<vfwf>p/(2/3qf@p) = ',&
               MEAN_PART(19,J)/(2./3.*MEAN_PART(20,J)),' [-]'
write(200,*)
write(200,*)
write(200,*)'Fluid-particle covariance tensor'
write(200,*)'--------------------------------'
write(200,10601)'<ufup>p = ',MEAN_PART(21,J),' [m2/s2]'
write(200,10601)'<vfvp>p = ',MEAN_PART(22,J),' [m2/s2]'
write(200,10601)'<wfwp>p = ',MEAN_PART(23,J),' [m2/s2]'
write(200,*)
write(200,10601)'<ufvp>p = ',MEAN_PART(24,J),' [m2/s2]'
write(200,10601)'<ufwp>p = ',MEAN_PART(25,J),' [m2/s2]'
write(200,10601)'<vfwp>p = ',MEAN_PART(26,J),' [m2/s2]'
write(200,*)
write(200,10601)'<ufup>p/(1/3qfp)=',&
               MEAN_PART(21,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,10601)'<vfvp>p/(1/3qfp)=',&
               MEAN_PART(22,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,10601)'<wfwp>p/(1/3qfp)=',&
               MEAN_PART(23,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,*)
write(200,10601)'<ufvp>p/(1/3qfp)=',&
               MEAN_PART(24,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,10601)'<ufwp>p/(1/3qfp)=',&
               MEAN_PART(25,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,10601)'<vfwp>p/(1/3qfp)=',&
               MEAN_PART(26,J)/(1./3.*MEAN_PART(27,J)),' [-]'
write(200,*)
if(PARTDEF(J)>1) then
write(200,10601)'<tfp_F> = ',1./MEAN_PART(28,J),' [s]'
write(200,10601)'  <Rep> = ',MEAN_PART(31,J),' [-]'
write(200,10601)'  <Vfp> = ',MEAN_PART(29,J),' [m/s]'
write(200,10601)'   <Cd> = ',MEAN_PART(30,J),' [-]'
write(200,*)
end if !-> PARTDEF(J)>1

!!- Spatial distribution statistics
if(LEVEL3_STPAR) then
write(200,*)
write(200,*)' Particle spatial distribution'
write(200,*)'--------------------------------'
write(200,10601)'   <C> = ',MEAN_PART(61,J),' [-]'
write(200,10601)' <C^2> = ',MEAN_PART(63,J),' [-]'
write(200,10601)'<Cp^2> = ',MEAN_PART(62,J),' [-]'
write(200,*)
end if


if(FILTERING) then
write(200,*)
write(200,*)' Subgrid turbulence statistics'
write(200,*)'--------------------------------'
write(200,10601)' <duf>p = ',MEAN_PART(70,J),' [m/s]'
write(200,10601)' <dvf>p = ',MEAN_PART(71,J),' [m/s]'
write(200,10601)' <dwf>p = ',MEAN_PART(72,J),' [m/s]'
write(200,*)
write(200,*)
write(200,10601)' <duf.duf>p = ',MEAN_PART(73,J),' [m2/s2]'
write(200,10601)' <dvf.dvf>p = ',MEAN_PART(74,J),' [m2/s2]'
write(200,10601)' <dwf.dwf>p = ',MEAN_PART(75,J),' [m2/s2]'
write(200,*)
write(200,10601)' <duf.dvf>p = ',MEAN_PART(76,J),' [m2/s2]'
write(200,10601)' <duf.dwf>p = ',MEAN_PART(77,J),' [m2/s2]'
write(200,10601)' <dvf.dwf>p = ',MEAN_PART(78,J),' [m2/s2]'
write(200,*)
write(200,10601)'      dqfap = ',MEAN_PART(79,J),' [m2/s2]'
write(200,*)
write(200,10601)' <duf.dvf>p/(2/3qf@p) = ',MEAN_PART(76,J)/(2./3.*MEAN_PART(79,J)),' [-]'
write(200,10601)' <duf.dwf>p/(2/3qf@p) = ',MEAN_PART(77,J)/(2./3.*MEAN_PART(79,J)),' [-]'
write(200,10601)' <dvf.dwf>p/(2/3qf@p) = ',MEAN_PART(78,J)/(2./3.*MEAN_PART(79,J)),' [-]'
write(200,*)
write(200,*)
write(200,10601)' <duf.uf>p = ',MEAN_PART(80,J),' [m2/s2]'
write(200,10601)' <dvf.vf>p = ',MEAN_PART(81,J),' [m2/s2]'
write(200,10601)' <dwf.wf>p = ',MEAN_PART(82,J),' [m2/s2]'
write(200,*)
write(200,*)
write(200,10601)' <duf.up>p = ',MEAN_PART(83,J),' [m2/s2]'
write(200,10601)' <dvf.vp>p = ',MEAN_PART(84,J),' [m2/s2]'
write(200,10601)' <dwf.wp>p = ',MEAN_PART(85,J),' [m2/s2]'
end if

if(SOLVE_SCALAR) then
write(200,*)
write(200,*)
write(200,*)'Statistics on temperature'
write(200,*)'-------------------------'
write(200,10601)'  <Tp>p = ',MEAN_PART(35,J),' [K]'
write(200,10601)'<Tf@p>p = ',MEAN_PART(36,J),' [K]'
write(200,*)
write(200,10601)'  <Tp^2>p = ',MEAN_PART(37,J),' [K2]'
write(200,10601)'<Tf@p^2>p = ',MEAN_PART(38,J),' [K2]'
write(200,10601)'<Tp.Tf@p>p = ',MEAN_PART(53,J),' [K2]'
write(200,*)
write(200,10601)'<Tp.up>p = ',MEAN_PART(39,J),' [K.m/s]'
write(200,10601)'<Tp.vp>p = ',MEAN_PART(40,J),' [K.m/s]'
write(200,10601)'<Tp.wp>p = ',MEAN_PART(41,J),' [K.m/s]'
write(200,*)
write(200,10601)'<Tp.uf@p>p = ',MEAN_PART(42,J),' [K.m/s]'
write(200,10601)'<Tp.vf@p>p = ',MEAN_PART(43,J),' [K.m/s]'
write(200,10601)'<Tp.wf@p>p = ',MEAN_PART(44,J),' [K.m/s]'
write(200,*)
write(200,10601)'<Tf@p.up>p = ',MEAN_PART(45,J),' [K.m/s]'
write(200,10601)'<Tf@p.vp>p = ',MEAN_PART(46,J),' [K.m/s]'
write(200,10601)'<Tf@p.wp>p = ',MEAN_PART(47,J),' [K.m/s]'
write(200,*)
write(200,10601)'<Tf@p.uf@p>p = ',MEAN_PART(48,J),' [K.m/s]'
write(200,10601)'<Tf@p.vf@p>p = ',MEAN_PART(49,J),' [K.m/s]'
write(200,10601)'<Tf@p.wf@p>p = ',MEAN_PART(50,J),' [K.m/s]'
if(PARTDEF(J)>1) then
write(200,*)
write(200,10601)'<tp_theta> = ',1./MEAN_PART(51,J),' [s]'
write(200,10601)'     <Nup> = ',MEAN_PART(52,J),' [-]'
write(200,*)
end if !-> PARTDEF(J)>1

end if !-> SOLVE_SCALAR


end do !-> J = 1, NIG

end if !-> LEVEL1_STPAR

close(200)



!!====================================================================
!! Print statistic for postprocessing
!!====================================================================
open(unit=900, file='data.last', status='replace')
write(900,*)'#   Particle class 1 2 3 ....'
if(LEVEL1_STPAR) then

write(900,10001)(MEAN_PART(10,J),J=1,NIG) !- qp
write(900,10001)(MEAN_PART(20,J),J=1,NIG) !- qfap
write(900,10001)(MEAN_PART(27,J),J=1,NIG) !- qfp
write(900,10001)(MEAN_PART(28,J),J=1,NIG) !- <1/tp>

if(SOLVE_SCALAR)then
write(900,10001)(MEAN_PART(37,J),J=1,NIG) !- <Tp^2>p
write(900,10001)(MEAN_PART(38,J),J=1,NIG) !- <Tf@p^2>p
write(900,10001)(MEAN_PART(51,J),J=1,NIG) !- <1/tp_theta>
write(900,10001)(MEAN_PART(53,J),J=1,NIG) !- <Tp.Tf@p>p
write(900,10001)(MEAN_PART(40,J),J=1,NIG) !- <Tp.vp>p
write(900,10001)(MEAN_PART(43,J),J=1,NIG) !- <Tp.vf@p>p
write(900,10001)(MEAN_PART(46,J),J=1,NIG) !- <Tf@p.vp>p
write(900,10001)(MEAN_PART(49,J),J=1,NIG) !- <Tf@p.vf@p>p
end if !-> SOLVE_SCALAR

end if !-> end if LEVEL1_STPAR
close(900)






!!====================================================================
10000 format (15(e17.7))
10001 format (15(2x,e17.7))


10600 format (2x,A,f8.2)
10601 format (2x,A,E13.6,A)
10602 format (2x,A,E13.6)
10603 format (2x,A,I2,A)


end subroutine PRINT_LASTCYCLE
