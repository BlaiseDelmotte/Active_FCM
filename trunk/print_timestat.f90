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
