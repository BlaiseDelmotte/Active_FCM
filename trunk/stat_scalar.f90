!!====================================================================
!!
!! This fortran file compute all statistics on scalar field
!!
!!====================================================================

subroutine STAT_SCALAR(NCYCLE,TIME)

!!====================================================================
!!
!! LEVEL1_STSCL: + Mean velocity
!!               + Variance
!!
!!--------------------------------------------------------------------
!! Time-averaged statistics are performed using a macro array
!! called MEAN_TIME_SCL. The averaging is done as
!!   --> MEAN_TIME_SCL = MEAN_TIME_SCL + MEAN_SCL
!!
!!--------------------------------------------------------------------
!! MEAN_SCL( 1): <theta>
!!           2 : <theta^2>
!!           3 : <uf.theta>
!!           4 : <vf.theta>
!!           5 : <wf.theta>
!!           6 : <(d theta/dx)>
!!           7 : <(d theta/dy)>
!!           8 : <(d theta/dz)>
!!           9 : <(d theta/dx)^2>
!!          10 : <(d theta/dy)^2>
!!          11 : <(d theta/dz)^2>
!!          12 : <(d theta/dx)^3>
!!          13 : <(d theta/dy)^3>
!!          14 : <(d theta/dz)^3>
!!          15 : <(d theta/dx)^4>
!!          16 : <(d theta/dy)^4>
!!          17 : <(d theta/dz)^4>
!!          18 : dissipation
!!          19 : production = -<v.theta>*gradient
!!          20 : skewness x
!!          21 : skewness y
!!          22 : skewness z
!!          23 : flatness x
!!          24 : flatness x
!!          25 : flatness x
!!====================================================================

use DNS_DIM	       !- Dimension
use PARAM_PHYS         !- Physical & numerical parameters
use FLUID_VARIABLE     !- Fluid velocity
use SCALAR_VARIABLE    !- Fluid velocity
use GEOMETRIC_VARIABLE !- 
use STATISTICS         !- Statistics
use CHECK_CPU	       !- CPU time checks
use WORK_ARRAYS

use P3DFFT

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
!- 
real(kind=8) :: RDUMMY

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

integer :: IFLAG1, I,J, K
!---------------------------------------------------------------------

!!- Check CPU time
!if(MYID == 0) then
! TIME_START=MPI_WTIME()
!end if



!!- Initiation
MEAN_SCL(:) = 0.


!!====================================================================
!! 1. 1st level of  statistic
!!====================================================================
if(LEVEL1_STSCL) then

!!--------------------------------------------------------------------
!! 1.1. Space-Averaging
!!--------------------------------------------------------------------
!!- <theta>
 call RSUMCPU(SUM(THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_SCL(1))

!!- <theta^2>
 call RSUMCPU(SUM(THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))**2),MEAN_SCL(2))

!!- <uf*theta>
 call RSUMCPU(SUM( UFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3))&
                 *THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_SCL(3))

!!- <vf*theta>
 call RSUMCPU(SUM( VFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)) &
                 *THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_SCL(4))

!!- <wf*theta>
 call RSUMCPU(SUM( WFLU(:,ISTART(2):IEND(2),ISTART(3):IEND(3)) &
                 *THETA(:,ISTART(2):IEND(2),ISTART(3):IEND(3))),MEAN_SCL(5))


!- Normalization by the full number of points
 MEAN_SCL = MEAN_SCL / NGLOB


!!--------------------------------------------------------------------
!! 1.2. Print in file
!!--------------------------------------------------------------------
 if(MYID==0) write(401,10000)TIME,       &
                   MEAN_SCL(1 ), & !- <theta>
                   MEAN_SCL(2 ), & !- <theta^2>
                   MEAN_SCL(3 ), & !- <uf.theta>
                   MEAN_SCL(4 ), & !- <vf.theta>
                   MEAN_SCL(5 )    !- <wf.theta>

end if


!!====================================================================
!! 2. 2nd level of  statistic
!!====================================================================
if(LEVEL2_STSCL) then

!!--------------------------------------------------------------------
!! 2.1 Scalar gradient in x-direction
!!--------------------------------------------------------------------
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  TMPFOU(I,J,K) = ICMPL*KX(I)*THETAFOU(I,J,K)
  end do
 end do
end do

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)       


!- < (d theta/dx)>
call RSUMCPU(SUM(TMPPHY   ),MEAN_SCL(6))

!- < (d theta/dx)^2>
call RSUMCPU(SUM(TMPPHY**2),MEAN_SCL(9))

!- < (d theta/dx)^3>
call RSUMCPU(SUM(TMPPHY**3),MEAN_SCL(12))

!- < (d theta/dx)^4>
call RSUMCPU(SUM(TMPPHY**4),MEAN_SCL(15))


!!--------------------------------------------------------------------
!! 2.2 Scalar gradient in y-direction
!!--------------------------------------------------------------------
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  TMPFOU(I,J,K) = ICMPL*KY(J)*THETAFOU(I,J,K)
  end do
 end do
end do

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)       


!- < (d theta/dy)>
call RSUMCPU(SUM(TMPPHY   ),MEAN_SCL(7))

!- < (d theta/dy)^2>
call RSUMCPU(SUM(TMPPHY**2),MEAN_SCL(10))

!- < (d theta/dy)^3>
call RSUMCPU(SUM(TMPPHY**3),MEAN_SCL(13))

!- < (d theta/dy)^4>
call RSUMCPU(SUM(TMPPHY**4),MEAN_SCL(16))



!!--------------------------------------------------------------------
!! 2.3 Scalar gradient in z-direction
!!--------------------------------------------------------------------
do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
  TMPFOU(I,J,K) = ICMPL*KZ(K)*THETAFOU(I,J,K)
  end do
 end do
end do

call P3DFFT_BTRAN_C2R(TMPFOU,TMPPHY,FFTFLAG)       


!- < (d theta/dz)>
call RSUMCPU(SUM(TMPPHY   ),MEAN_SCL( 8))

!- < (d theta/dz)^2>
call RSUMCPU(SUM(TMPPHY**2),MEAN_SCL(11))

!- < (d theta/dz)^3>
call RSUMCPU(SUM(TMPPHY**3),MEAN_SCL(14))

!- < (d theta/dz)^4>
call RSUMCPU(SUM(TMPPHY**4),MEAN_SCL(17))


!!--------------------------------------------------------------------
!! 2.4 Compute statistics based gradient of scalar
!!--------------------------------------------------------------------

!- Normalization by the full number of points
MEAN_SCL(6:17) = MEAN_SCL(6:17)  / NGLOB


!- Temperature dissipation
MEAN_SCL(18) = DIFF_SCL*(MEAN_SCL(9) + MEAN_SCL(10) + MEAN_SCL(11))

!- Production by the mean imposed gradient
MEAN_SCL(19) = -MEAN_SCL(4)*GRAD_SCL


!- Skewness in x-direction
MEAN_SCL(20) = MEAN_SCL(12)/MEAN_SCL( 9)**1.5

!- Skewness in y-direction
MEAN_SCL(21) = MEAN_SCL(13)/MEAN_SCL(10)**1.5

!- Skewness in z-direction
MEAN_SCL(22) = MEAN_SCL(14)/MEAN_SCL(11)**1.5



!- Flatness in x-direction
MEAN_SCL(23) = MEAN_SCL(15)/ MEAN_SCL( 9)**2

!- Flatness in y-direction
MEAN_SCL(24) = MEAN_SCL(16)/ MEAN_SCL(10)**2

!- Flatness in z-direction
MEAN_SCL(25) = MEAN_SCL(17)/ MEAN_SCL(11)**2




!!--------------------------------------------------------------------
!! 2.5 Print in file
!!--------------------------------------------------------------------
 if(MYID==0) write(402,10000)TIME,       &
                   MEAN_SCL( 2), & !- <theta^2>
                   MEAN_SCL(18), & !- dissipation
                   MEAN_SCL(19), & !- production
                   MEAN_SCL(20), & !- Sx
                   MEAN_SCL(21), & !- Sy
                   MEAN_SCL(22), & !- Sz
                   MEAN_SCL(23), & !- Kx
                   MEAN_SCL(24), & !- Ky
                   MEAN_SCL(25)    !- Kz

end if

!!======================================================================
!! 5. Time averaging if so
!!======================================================================
if(STAT_TIME) then
 MEAN_TIME_SCL = MEAN_TIME_SCL + MEAN_SCL
end if 

!!======================================================================
!! 6. Print statistics in "info" files
!!======================================================================

!- Print in file "stat.info"
if((mod(NCYCLE,FOUT0) == 0).and.(MYID==0)) then
IFLAG1 = 2
write(UNIT_INFO(IFLAG1),*)
write(UNIT_INFO(IFLAG1),*)'====================================='
write(UNIT_INFO(IFLAG1),*)'Statistic for cycle = ',NCYCLE
write(UNIT_INFO(IFLAG1),*)'====================================='
write(UNIT_INFO(IFLAG1),10601)'             <theta> = ',MEAN_SCL(1)
write(UNIT_INFO(IFLAG1),10601)'           <theta^2> = ',MEAN_SCL(2)
write(UNIT_INFO(IFLAG1),*)
!write(UNIT_INFO(IFLAG1),10601)'          <uf.theta> = ',MEAN_SCL(3)
!write(UNIT_INFO(IFLAG1),10601)'          <vf.theta> = ',MEAN_SCL(4)
!write(UNIT_INFO(IFLAG1),10601)'          <wf.theta> = ',MEAN_SCL(5)
!write(UNIT_INFO(IFLAG1),*)
end if





!!- CPU check
!if(MYID == 0) then
! TIME_END=MPI_WTIME()
! CPU_FLUID(5) = TIME_END - TIME_START
!end if


!!----------------------------------------------------------------------
10000 format (20(e17.7))
10601 format (2x,A,E13.6)
10602 format (2x,A,E13.6,1x,A,E13.6)
10603 format (2x,A,E13.6,2x,A,E13.6,2x,A,E13.6)

end subroutine STAT_SCALAR
