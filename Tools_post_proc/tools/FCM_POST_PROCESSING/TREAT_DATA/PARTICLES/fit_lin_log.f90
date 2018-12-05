subroutine FIT_LIN_LOG(NR, RANGE_R, CORREL_INTEGRAL, SLOPE)

implicit none

!-  Arguments
integer, intent(in) :: NR
real(kind=8), dimension(NR), intent(in) :: RANGE_R
real(kind=8), dimension(NR), intent(in) :: CORREL_INTEGRAL
real(kind=8), intent(out) :: SLOPE

!- Local parameters
real(kind=8) :: SUMX,SUMY,SUMSQX,SUMXY,DENO
integer :: L         

SUMX=0.0
SUMY=0.0
SUMSQX=0.0
SUMXY=0.0

! Linear Least Square fit of data 
do L=1,NR
 if(CORREL_INTEGRAL(L).gt.0)then
  SUMX = SUMX + log(RANGE_R(L))	
  SUMY = SUMY + log(CORREL_INTEGRAL(L))
  SUMSQX = SUMSQX + log(RANGE_R(L))*log(RANGE_R(L))
  SUMXY = SUMXY + log(RANGE_R(L))*log(CORREL_INTEGRAL(L))
 endif
enddo

DENO=NR*SUMSQX-SUMX*SUMX

SLOPE=(NR*SUMXY-SUMX*SUMY)/DENO
!~ B=(SUMSQX*SUMY-SUMX*SUMXY)/DENO


end subroutine FIT_LIN_LOG
