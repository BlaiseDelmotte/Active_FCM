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

close(200)






!!====================================================================
10000 format (15(e17.7))
10001 format (15(2x,e17.7))


10600 format (2x,A,f8.2)
10601 format (2x,A,E13.6,A)
10602 format (2x,A,E13.6)
10603 format (2x,A,I2,A)


end subroutine PRINT_LASTCYCLE
