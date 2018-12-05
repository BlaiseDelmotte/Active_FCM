!!====================================================================
!!
!!
!!====================================================================

subroutine PRINT_TIMESTAT(NCYCLE)

!!====================================================================
!!
!!
!!====================================================================
use DNS_DIM            !- Dimension
use STATISTICS         !- Statistics
use PARAM_PHYS
use GEOMETRIC_VARIABLE

implicit none

!!--------------------------------------------------------------------
!! ARRAYS STATEMENT
!!--------------------------------------------------------------------
integer, intent(in) :: NCYCLE

!!--------------------------------------------------------------------
!!- Local arrays
!!--------------------------------------------------------------------
real(kind=8) :: DYSTAT
real(kind=8), dimension(DIMSCOR) :: RUX, RVX, RWX

real(kind=8), dimension(NPDFCP) :: PDFCP
!- File name
character(len=40) :: FILENAME

!- Index
integer :: I, J, LGR, M ,N
!---------------------------------------------------------------------


!!====================================================================
!! 1. Fluid statistics normalization
!!====================================================================
if(LEVEL0_STFLU) then

!!--------------------------------------------------------------------
!! 1.1. Time-averaged statistics
!!--------------------------------------------------------------------
 MEAN_TIME_FLUID = MEAN_TIME_FLUID / real(NEVEN)

!!--------------------------------------------------------------------
!! 1.2. Spatial Correlation functions
!!--------------------------------------------------------------------
if(LEVEL4_STFLU) then
 call MPI_ALLREDUCE(MEAN_RUXLOC,RUX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(MEAN_RVXLOC,RVX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(MEAN_RWXLOC,RWX, &
                    DIMSCOR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
!
!
 if(MYID==0) then
  do I = 1, DIMSCOR 
   write(304,10000)XMESH(I),    RUX(I)                         &
                            ,0.5*(RVX(I)+RWX(I))               &
                            ,     RUX(I)        / real(NEVEN)  &
                            ,0.5*(RVX(I)+RWX(I))/ real(NEVEN)  &
                            ,     RUX(I)        /RUX(1)        &
                            ,(RVX(I)+RWX(I))/(RVX(1)+RWX(1))
  end do
 end if

end if !- endif(LEVEL4_STFLU)

end if !- end if(LEVEL0_STFLU)

!!====================================================================
!! 2. Scalar statistics normalization
!!====================================================================
if(LEVEL0_STSCL) then

!!--------------------------------------------------------------------
!! 2.1. Time-averaged statistics
!!--------------------------------------------------------------------
 MEAN_TIME_SCL = MEAN_TIME_SCL / real(NEVEN)

end if



!!====================================================================
!! 3. Particle statistics normalization
!!====================================================================
if(LEVEL0_STPAR) then

!!--------------------------------------------------------------------
!! 3.1. One point statistics
!!--------------------------------------------------------------------
 MEAN_TIME_PART = MEAN_TIME_PART / real(NEVEN)


!!--------------------------------------------------------------------
!! 3.2. Spatial distribution
!!--------------------------------------------------------------------
 if(LEVEL3_STPAR) then
 
  do J = 1, NIG
   call MPI_ALLREDUCE(PDFCP_LOC(:,J),PDFCP,NPDFCP,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
   write(FILENAME,10100)'part_l2_pdfcp_p',J,'.stat'
   open(unit=600, file=trim(FILENAME), status='replace')
   do I=1, NPDFCP
    write(600,10000)real(I-1),PDFCP(I)
   end do
   close(600)
  end do
 end if

end if


!!====================================================================
!! 3. Print in "stat.scilab" for matlab or scilab post-processing
!!====================================================================
!!call PRINT_SCILAB


!!====================================================================
!! 4. Print in "run.info"
!!====================================================================
if(MYID==0) then 

I = 1
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)'TIME AVERAGED STATISTICS '
write(UNIT_INFO(I),*)'====================================================================='
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10600)'Number of computed cycle      =',real(NCYCLE)
write(UNIT_INFO(I),10600)'Number of event for statistics=',real(NEVEN)
write(UNIT_INFO(I),*)
if(LEVEL0_STFLU) then
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)' FLUID'
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
end if
if(LEVEL1_STFLU) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Mean fluid Velocity'
write(UNIT_INFO(I),*)'-------------------'
write(UNIT_INFO(I),10601)'<uf>f = ',MEAN_TIME_FLUID(1),' [m/s]'
write(UNIT_INFO(I),10601)'<vf>f = ',MEAN_TIME_FLUID(2),' [m/s]'
write(UNIT_INFO(I),10601)'<wf>f = ',MEAN_TIME_FLUID(3),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Reynolds stress tensor'
write(UNIT_INFO(I),*)'----------------------'
write(UNIT_INFO(I),10601)'<uf.uf>f = ',MEAN_TIME_FLUID(4),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vf.vf>f = ',MEAN_TIME_FLUID(5),' [m2/s2]'
write(UNIT_INFO(I),10601)'<wf.wf>f = ',MEAN_TIME_FLUID(6),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'    <kf> = ',MEAN_TIME_FLUID(10),' [m2/s2]'
write(UNIT_INFO(I),10601)'  <epsf> = ',MEAN_TIME_FLUID(11),' [m2/s3]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*    )'   u_rms=(2/3*kf)**0.5'
write(UNIT_INFO(I),10601)'   u_rms = ',sqrt(2./3.*MEAN_TIME_FLUID(10)),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<uf.vf>f = ',MEAN_TIME_FLUID(7),' [m2/s2]'
write(UNIT_INFO(I),10601)'<uf.wf>f = ',MEAN_TIME_FLUID(8),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vf.wf>f = ',MEAN_TIME_FLUID(9),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'   eta_K = ',MEAN_TIME_FLUID(12),' [m]'
write(UNIT_INFO(I),10601)'   tau_K = ',MEAN_TIME_FLUID(13),' [s]'
write(UNIT_INFO(I),10601)'     v_K = ',MEAN_TIME_FLUID(17),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'Taylor      lg = ',MEAN_TIME_FLUID(25),' [m]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Reynolds number'
write(UNIT_INFO(I),*)'---------------'
write(UNIT_INFO(I),10602)' Re_Tay = ',sqrt(2.*MEAN_TIME_FLUID(10)/3.)*MEAN_TIME_FLUID(25)/VISC
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
end if
if(LEVEL2_STFLU) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Divergence'
write(UNIT_INFO(I),*)'----------'
write(UNIT_INFO(I),10601)'   < div(U)  > = ',MEAN_TIME_FLUID(15),' [1/s]'
write(UNIT_INFO(I),10601)'   < div(U)^2> = ',MEAN_TIME_FLUID(16),' [1/s2]'
write(UNIT_INFO(I),10601)'max(|div(U)| ) = ',MEAN_TIME_FLUID(14),' [1/s]'
write(UNIT_INFO(I),*)
end if

if(LEVEL3_STFLU) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Velocity gradients'
write(UNIT_INFO(I),*)'------------------'
write(UNIT_INFO(I),10601)' <(dui/dxi)^2> = ',MEAN_TIME_FLUID(21),' [1/s2]'
write(UNIT_INFO(I),10601)' <(dui/dxj)^2> = ',MEAN_TIME_FLUID(31),' [1/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10602)'Skewness:   Sk = ',MEAN_TIME_FLUID(18)
write(UNIT_INFO(I),10602)'Flatness:   Tk = ',MEAN_TIME_FLUID(19)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'Taylor      lf = ',MEAN_TIME_FLUID(24),' [m]'
write(UNIT_INFO(I),10601)'Taylor      lg = ',MEAN_TIME_FLUID(25),' [m]'
write(UNIT_INFO(I),10602)'   lf/lg/2^0.5 = ',MEAN_TIME_FLUID(24)/MEAN_TIME_FLUID(25)/sqrt(2.)
write(UNIT_INFO(I),*)
end if
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)


if(LEVEL0_STSCL) then
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)' SCALAR'
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
end if
if(LEVEL1_STSCL) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'  <Tf>f = ',MEAN_TIME_SCL(1),' [K]'
write(UNIT_INFO(I),10601)'<Tf^2>f = ',MEAN_TIME_SCL(2),' [K2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<uf.Tf>f = ',MEAN_TIME_SCL(3),' [K.m/s]'
write(UNIT_INFO(I),10601)'<vf.Tf>f = ',MEAN_TIME_SCL(4),' [K.m/s]'
write(UNIT_INFO(I),10601)'<wf.Tf>f = ',MEAN_TIME_SCL(5),' [K.m/s]'
end if
if(LEVEL2_STSCL) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' Scalar dissipation = ',MEAN_TIME_SCL(18),' [K2/s]'
write(UNIT_INFO(I),10601)'         production = ',MEAN_TIME_SCL(19),' [K2/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'       Skewness: Sx = ',MEAN_TIME_SCL(20), ' [-]'
write(UNIT_INFO(I),10601)'                 Sy = ',MEAN_TIME_SCL(21), ' [-]'
write(UNIT_INFO(I),10601)'                 Sz = ',MEAN_TIME_SCL(22), ' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'       Flatness: Kx = ',MEAN_TIME_SCL(23), ' [-]'
write(UNIT_INFO(I),10601)'                 Ky = ',MEAN_TIME_SCL(24), ' [-]'
write(UNIT_INFO(I),10601)'                 Kz = ',MEAN_TIME_SCL(25), ' [-]'
end if



if(LEVEL0_STPAR) then
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)' PARTICLE'
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
end if

if(LEVEL1_STPAR) then
do J = 1, NIG
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
if(PARTDEF(J) == 0) then
write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Motionless'
elseif(PARTDEF(J) == 1) then
write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Fluid element'
elseif(PARTDEF(J) == 2) then
write(UNIT_INFO(I),10603)'Particle:',J,'  ==> Solid particle'
end if
write(UNIT_INFO(I),*)'---------------------------------------------------------------------'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'Diameter=',DPART(J),'  [m]'
write(UNIT_INFO(I),10601)'Density=',RHOP(J),'  [kg/m3]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Mean Velocity'
write(UNIT_INFO(I),*)'-------------'
write(UNIT_INFO(I),10601)'<up>p = ',MEAN_TIME_PART(1,J),' [m/s]'
write(UNIT_INFO(I),10601)'<vp>p = ',MEAN_TIME_PART(2,J),' [m/s]'
write(UNIT_INFO(I),10601)'<wp>p = ',MEAN_TIME_PART(3,J),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<uf>p = ',MEAN_TIME_PART(11,J),' [m/s]'
write(UNIT_INFO(I),10601)'<vf>p = ',MEAN_TIME_PART(12,J),' [m/s]'
write(UNIT_INFO(I),10601)'<wf>p = ',MEAN_TIME_PART(13,J),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Fluctuating motion'
write(UNIT_INFO(I),*)'------------------'
write(UNIT_INFO(I),10601)'  qp = ',MEAN_TIME_PART(10,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'qf@p = ',MEAN_TIME_PART(20,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' qfp = ',MEAN_TIME_PART(27,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Particle-Particle kinetic stress tensor'
write(UNIT_INFO(I),*)'--------------------------------------'
write(UNIT_INFO(I),10601)'<upup>p = ',MEAN_TIME_PART(4,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vpvp>p = ',MEAN_TIME_PART(5,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<wpwp>p = ',MEAN_TIME_PART(6,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<upvp>p = ',MEAN_TIME_PART(7,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<upwp>p = ',MEAN_TIME_PART(8,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vpwp>p = ',MEAN_TIME_PART(9,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<upup>p/(2/3qp) = ',&
               MEAN_TIME_PART(4,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),10601)'<vpvp>p/(2/3qp) = ',&
               MEAN_TIME_PART(5,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),10601)'<wpwp>p/(2/3qp) = ',&
               MEAN_TIME_PART(6,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<upvp>p/(2/3qp) = ',&
               MEAN_TIME_PART(7,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),10601)'<upwp>p/(2/3qp) = ',&
               MEAN_TIME_PART(8,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),10601)'<vpwp>p/(2/3qp) = ',&
               MEAN_TIME_PART(9,J)/(2./3.*MEAN_TIME_PART(10,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Fluid-fluid @p kinetic stress tensor'
write(UNIT_INFO(I),*)'------------------------------------'
write(UNIT_INFO(I),10601)'<ufuf>p = ',MEAN_TIME_PART(14,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vfvf>p = ',MEAN_TIME_PART(15,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<wfwf>p = ',MEAN_TIME_PART(16,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufvf>p = ',MEAN_TIME_PART(17,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<ufwf>p = ',MEAN_TIME_PART(18,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vfwf>p = ',MEAN_TIME_PART(19,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufuf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(14,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),10601)'<vfvf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(15,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),10601)'<wfwf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(16,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufvf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(17,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),10601)'<ufwf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(18,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),10601)'<vfwf>p/(2/3qf@p) = ',&
               MEAN_TIME_PART(19,J)/(2./3.*MEAN_TIME_PART(20,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Fluid-particle covariance tensor'
write(UNIT_INFO(I),*)'--------------------------------'
write(UNIT_INFO(I),10601)'<ufup>p = ',MEAN_TIME_PART(21,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vfvp>p = ',MEAN_TIME_PART(22,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<wfwp>p = ',MEAN_TIME_PART(23,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufvp>p = ',MEAN_TIME_PART(24,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<ufwp>p = ',MEAN_TIME_PART(25,J),' [m2/s2]'
write(UNIT_INFO(I),10601)'<vfwp>p = ',MEAN_TIME_PART(26,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufup>p/(1/3qfp)=',&
               MEAN_TIME_PART(21,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),10601)'<vfvp>p/(1/3qfp)=',&
               MEAN_TIME_PART(22,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),10601)'<wfwp>p/(1/3qfp)=',&
               MEAN_TIME_PART(23,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<ufvp>p/(1/3qfp)=',&
               MEAN_TIME_PART(24,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),10601)'<ufwp>p/(1/3qfp)=',&
               MEAN_TIME_PART(25,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),10601)'<vfwp>p/(1/3qfp)=',&
               MEAN_TIME_PART(26,J)/(1./3.*MEAN_TIME_PART(27,J)),' [-]'
write(UNIT_INFO(I),*)

if(PARTDEF(J)>1) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<tfp_F>p = ',1./MEAN_TIME_PART(28,J),' [s]'
write(UNIT_INFO(I),10601)'  <Rep>p = ',MEAN_TIME_PART(31,J),' [-]'
write(UNIT_INFO(I),10601)'  <Vfp>p = ',MEAN_TIME_PART(29,J),' [m/s]'
write(UNIT_INFO(I),10601)'   <Cd>p = ',MEAN_TIME_PART(30,J),' [-]'
write(UNIT_INFO(I),*)
end if

!!- Spatial distribution statistics
if(LEVEL3_STPAR) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)' Particle spatial distribution'
write(UNIT_INFO(I),*)'--------------------------------'
write(UNIT_INFO(I),10601)'   <C> = ',MEAN_TIME_PART(61,J),' [-]'
write(UNIT_INFO(I),10601)' <C^2> = ',MEAN_TIME_PART(63,J),' [-]'
write(UNIT_INFO(I),10601)'<Cp^2> = ',MEAN_TIME_PART(62,J),' [-]'
write(UNIT_INFO(I),10601)' sigma =',(sqrt(MEAN_TIME_PART(62,J))-sqrt(MEAN_TIME_PART(61,J)))/MEAN_TIME_PART(61,J),' [-]'
write(UNIT_INFO(I),*)
end if



if(FILTERING) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)' Subgrid turbulence statistics'
write(UNIT_INFO(I),*)'--------------------------------'
write(UNIT_INFO(I),10601)' <duf>p = ',MEAN_TIME_PART(70,J),' [m/s]'
write(UNIT_INFO(I),10601)' <dvf>p = ',MEAN_TIME_PART(71,J),' [m/s]'
write(UNIT_INFO(I),10601)' <dwf>p = ',MEAN_TIME_PART(72,J),' [m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' <duf.duf>p = ',MEAN_TIME_PART(73,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dvf.dvf>p = ',MEAN_TIME_PART(74,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dwf.dwf>p = ',MEAN_TIME_PART(75,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' <duf.dvf>p = ',MEAN_TIME_PART(76,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <duf.dwf>p = ',MEAN_TIME_PART(77,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dvf.dwf>p = ',MEAN_TIME_PART(78,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'      dqfap = ',MEAN_TIME_PART(79,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' <duf.dvf>p/(2/3qf@p) = ',MEAN_TIME_PART(76,J)/(2./3.*MEAN_TIME_PART(79,J)),' [-]'
write(UNIT_INFO(I),10601)' <duf.dwf>p/(2/3qf@p) = ',MEAN_TIME_PART(77,J)/(2./3.*MEAN_TIME_PART(79,J)),' [-]'
write(UNIT_INFO(I),10601)' <dvf.dwf>p/(2/3qf@p) = ',MEAN_TIME_PART(78,J)/(2./3.*MEAN_TIME_PART(79,J)),' [-]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' <duf.uf>p = ',MEAN_TIME_PART(80,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dvf.vf>p = ',MEAN_TIME_PART(81,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dwf.wf>p = ',MEAN_TIME_PART(82,J),' [m2/s2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)' <duf.up>p = ',MEAN_TIME_PART(83,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dvf.vp>p = ',MEAN_TIME_PART(84,J),' [m2/s2]'
write(UNIT_INFO(I),10601)' <dwf.wp>p = ',MEAN_TIME_PART(85,J),' [m2/s2]'
end if




if(SOLVE_SCALAR) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),*)'Statistics on temperature'
write(UNIT_INFO(I),*)'-------------------------'
write(UNIT_INFO(I),10601)'     <Tp>p = ',MEAN_TIME_PART(35,J),' [K]'
write(UNIT_INFO(I),10601)'   <Tf@p>p = ',MEAN_TIME_PART(36,J),' [K]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'   <Tp^2>p = ',MEAN_TIME_PART(37,J),' [K2]'
write(UNIT_INFO(I),10601)' <Tf@p^2>p = ',MEAN_TIME_PART(38,J),' [K2]'
write(UNIT_INFO(I),10601)'<Tf@p.Tp>p = ',MEAN_TIME_PART(53,J),' [K2]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'  <Tp.up>p = ',MEAN_TIME_PART(39,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'  <Tp.vp>p = ',MEAN_TIME_PART(40,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'  <Tp.wp>p = ',MEAN_TIME_PART(41,J),' [K.m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<Tp.uf@p>p = ',MEAN_TIME_PART(42,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tp.vf@p>p = ',MEAN_TIME_PART(43,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tp.wf@p>p = ',MEAN_TIME_PART(44,J),' [K.m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<Tf@p.up>p = ',MEAN_TIME_PART(45,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tf@p.vp>p = ',MEAN_TIME_PART(46,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tf@p.wp>p = ',MEAN_TIME_PART(47,J),' [K.m/s]'
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<Tf@p.uf@p>p = ',MEAN_TIME_PART(48,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tf@p.vf@p>p = ',MEAN_TIME_PART(49,J),' [K.m/s]'
write(UNIT_INFO(I),10601)'<Tf@p.wf@p>p = ',MEAN_TIME_PART(50,J),' [K.m/s]'
if(PARTDEF(J)>1) then
write(UNIT_INFO(I),*)
write(UNIT_INFO(I),10601)'<tp_theta> = ',1./MEAN_TIME_PART(51,J),' [s]'
write(UNIT_INFO(I),10601)'     <Nup> = ',MEAN_TIME_PART(52,J),' [-]'
write(UNIT_INFO(I),*)
end if

end if !!- if(SOLVE_SCALAR) 

end do !!- end loop J = 1, NIG

end if





!!====================================================================
!! Print statistic for postprocessing
!!====================================================================
open(unit=900, file='data.mean', status='replace')
write(900,*)'#   Particle class 1 2 3 ....'
if(LEVEL1_STPAR) then
write(900,10001)(MEAN_TIME_PART( 1,J),J=1,NIG) !- <up>p
write(900,10001)(MEAN_TIME_PART( 2,J),J=1,NIG) !- <vp>p
write(900,10001)(MEAN_TIME_PART( 3,J),J=1,NIG) !- <wp>p
write(900,10001)(MEAN_TIME_PART(11,J),J=1,NIG) !- <uf>p
write(900,10001)(MEAN_TIME_PART(12,J),J=1,NIG) !- <vf>p
write(900,10001)(MEAN_TIME_PART(13,J),J=1,NIG) !- <wf>p
write(900,10001)(MEAN_TIME_PART(10,J),J=1,NIG) !-  qp
write(900,10001)(MEAN_TIME_PART(20,J),J=1,NIG) !-  qf@p
write(900,10001)(MEAN_TIME_PART(27,J),J=1,NIG) !-  qfp
write(900,10001)(MEAN_TIME_PART( 4,J),J=1,NIG) !- <upup>p
write(900,10001)(MEAN_TIME_PART( 5,J),J=1,NIG) !- <vpvp>p
write(900,10001)(MEAN_TIME_PART( 6,J),J=1,NIG) !- <wpwp>p
write(900,10001)(MEAN_TIME_PART( 7,J),J=1,NIG) !- <upvp>p
write(900,10001)(MEAN_TIME_PART( 8,J),J=1,NIG) !- <upwp>p
write(900,10001)(MEAN_TIME_PART( 9,J),J=1,NIG) !- <vpwp>p
write(900,10001)(MEAN_TIME_PART(14,J),J=1,NIG) !- <ufuf>p
write(900,10001)(MEAN_TIME_PART(15,J),J=1,NIG) !- <vfvf>p
write(900,10001)(MEAN_TIME_PART(16,J),J=1,NIG) !- <wfwf>p
write(900,10001)(MEAN_TIME_PART(17,J),J=1,NIG) !- <ufvf>p
write(900,10001)(MEAN_TIME_PART(18,J),J=1,NIG) !- <ufwf>p
write(900,10001)(MEAN_TIME_PART(19,J),J=1,NIG) !- <vfwf>p
write(900,10001)(MEAN_TIME_PART(21,J),J=1,NIG) !- <ufup>p
write(900,10001)(MEAN_TIME_PART(22,J),J=1,NIG) !- <vfvp>p
write(900,10001)(MEAN_TIME_PART(23,J),J=1,NIG) !- <wfwp>p
write(900,10001)(MEAN_TIME_PART(24,J),J=1,NIG) !- <ufvp>p
write(900,10001)(MEAN_TIME_PART(25,J),J=1,NIG) !- <ufwp>p
write(900,10001)(MEAN_TIME_PART(26,J),J=1,NIG) !- <vfwp>p
write(900,10001)(MEAN_TIME_PART(28,J),J=1,NIG) !- <1/tfp_F>p
write(900,10001)(MEAN_TIME_PART(31,J),J=1,NIG) !-  <Rep>p
write(900,10001)(MEAN_TIME_PART(29,J),J=1,NIG) !-  <Vfp>p
write(900,10001)(MEAN_TIME_PART(30,J),J=1,NIG) !-   <Cd>p

if(LEVEL3_STPAR) then
write(900,10001)(MEAN_TIME_PART(61,J),J=1,NIG) !- <C>
write(900,10001)(MEAN_TIME_PART(63,J),J=1,NIG) !-<C^2> 
write(900,10001)(MEAN_TIME_PART(62,J),J=1,NIG) !-<Cp^2>
write(900,10001)((sqrt(MEAN_TIME_PART(62,J))-sqrt(MEAN_TIME_PART(61,J)))/MEAN_TIME_PART(61,J),J=1,NIG) !- sigma
end if


if(SOLVE_SCALAR)then
write(900,10001)(MEAN_TIME_PART(37,J),J=1,NIG) !- <Tp^2>p
write(900,10001)(MEAN_TIME_PART(38,J),J=1,NIG) !- <Tf@p^2>p
write(900,10001)(MEAN_TIME_PART(51,J),J=1,NIG) !- <1/tp_theta>
write(900,10001)(MEAN_TIME_PART(53,J),J=1,NIG) !- <Tp.Tf@p>p
write(900,10001)(MEAN_TIME_PART(40,J),J=1,NIG) !- <Tp.vp>p
write(900,10001)(MEAN_TIME_PART(43,J),J=1,NIG) !- <Tp.vf@p>p
write(900,10001)(MEAN_TIME_PART(46,J),J=1,NIG) !- <Tf@p.vp>p
write(900,10001)(MEAN_TIME_PART(49,J),J=1,NIG) !- <Tf@p.vf@p>p
end if

if(FILTERING) then
write(900,10001)(MEAN_TIME_PART(70,J),J=1,NIG) !- <duf>p
write(900,10001)(MEAN_TIME_PART(71,J),J=1,NIG) !- <dvf>p
write(900,10001)(MEAN_TIME_PART(72,J),J=1,NIG) !- <dwf>p
write(900,10001)(MEAN_TIME_PART(73,J),J=1,NIG) !- <duf.duf>p
write(900,10001)(MEAN_TIME_PART(74,J),J=1,NIG) !- <dvf.dvf>p
write(900,10001)(MEAN_TIME_PART(75,J),J=1,NIG) !- <dwf.dwf>p
write(900,10001)(MEAN_TIME_PART(76,J),J=1,NIG) !- <duf.dvf>p
write(900,10001)(MEAN_TIME_PART(77,J),J=1,NIG) !- <duf.dwf>p
write(900,10001)(MEAN_TIME_PART(78,J),J=1,NIG) !- <dvf.dwf>p
write(900,10001)(MEAN_TIME_PART(79,J),J=1,NIG) !-      dqfap
write(900,10001)(MEAN_TIME_PART(80,J),J=1,NIG) !- <duf.uf>p
write(900,10001)(MEAN_TIME_PART(81,J),J=1,NIG) !- <dvf.vf>p
write(900,10001)(MEAN_TIME_PART(82,J),J=1,NIG) !- <dwf.wf>p
write(900,10001)(MEAN_TIME_PART(83,J),J=1,NIG) !- <duf.up>p
write(900,10001)(MEAN_TIME_PART(84,J),J=1,NIG) !- <dvf.vp>p
write(900,10001)(MEAN_TIME_PART(85,J),J=1,NIG) !- <dwf.wp>p
end if

end if !!- end if LEVEL1_STPAR
close(900)




end if !!- MyID==0









!!====================================================================
10000 format (40(e17.7))
10001 format (15(2x,e17.7))

10100 format (A,I2.2,A)


10600 format (2x,A,f8.2)
10601 format (2x,A,E13.6,A)
10602 format (2x,A,E13.6)
10603 format (2x,A,I2,A)



end subroutine PRINT_TIMESTAT
