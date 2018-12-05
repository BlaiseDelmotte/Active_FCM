!!====================================================================
!!
!!
!!====================================================================
!!
subroutine PRINT_LAGFUNCTION
!!
!!====================================================================
!! The normalized autocorrelation function is computed as
!!
!!                 <u(t0).u(t0+t)>
!! R(t) = ---------------------------------
!!        sqrt(<u(t0)^2>).sqrt(<u(t0+t)^2>)
!!
!!====================================================================

use dns_dim
use param_phys
use particle_parallel
use statistics

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!- Local arrays
!---------------------------------------------------------------------
!- Lagrangian function
real(kind=8), dimension(DIMLGR) ::   RPX,   RPY,   RPZ
real(kind=8), dimension(DIMLGR) :: RFAPX, RFAPY, RFAPZ

real(kind=8), dimension(DIMLGR) :: RTP, RTFAP
real(kind=8), dimension(DIMLGR) :: RTFAPVFAP, RVFAPTFAP

!!- filtering statistics
real(kind=8), dimension(DIMLGR) :: RDUFAPDUFAP
real(kind=8), dimension(DIMLGR) :: RDVFAPDVFAP
real(kind=8), dimension(DIMLGR) :: RDWFAPDWFAP
real(kind=8), dimension(DIMLGR) :: RDUFAPUFAP
real(kind=8), dimension(DIMLGR) :: RDVFAPVFAP
real(kind=8), dimension(DIMLGR) :: RDWFAPWFAP

real(kind=8), dimension(DIMLGR) :: RDUFAPUP
real(kind=8), dimension(DIMLGR) :: RDVFAPVP
real(kind=8), dimension(DIMLGR) :: RDWFAPWP


!- True time step (accouting for frequency)
real(kind=8) :: DTIME0

!- File name
character(len=40) :: FILENAME

integer :: NUNIT

!- Index
integer :: I, J, LGR, K
!---------------------------------------------------------------------

if(MYID==0)write(*,*)' Lag. Func.: NBLGRMAX =',NBLGRMAX


DTIME0 = DTIME*FOUT2

!!- loop on lagrangian function
do K =1, NBLGRMAX


if(MYID==0)write(*,10700)' Lag. Func. #',K,' starts at N=',NT0(K)*FOUT2,' Time =',NT0(K)*DTIME0


 !!- Loop on particle kind
 do J = 1, NIG


!!======================================================================
!! 1. Add all Lagrangian function from each CPU
!!======================================================================
 call MPI_ALLREDUCE(RPX_LOC(:,J,K),RPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPY_LOC(:,J,K),RPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RPZ_LOC(:,J,K),RPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

 call MPI_ALLREDUCE(RFAPX_LOC(:,J,K),RFAPX,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RFAPY_LOC(:,J,K),RFAPY,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
 call MPI_ALLREDUCE(RFAPZ_LOC(:,J,K),RFAPZ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!- Normalization
 RPX = RPX / real(NPROC)
 RPY = RPY / real(NPROC)
 RPZ = RPZ / real(NPROC)

 RFAPX = RFAPX / real(NPROC)
 RFAPY = RFAPY / real(NPROC)
 RFAPZ = RFAPZ / real(NPROC)


 if(SOLVE_SCALAR) then
 
  call MPI_ALLREDUCE(      RTP_LOC(:,J,K),RTP      ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(    RTFAP_LOC(:,J,K),RTFAP    ,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RTFAPVFAP_LOC(:,J,K),RTFAPVFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RVFAPTFAP_LOC(:,J,K),RVFAPTFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

  RTP = RTP / real(NPROC)
  RTFAP = RTFAP / real(NPROC)

  RTFAPVFAP = RTFAPVFAP / real(NPROC)
  RVFAPTFAP = RVFAPTFAP / real(NPROC)
  
 end if !- SOLVE_SCALAR
 
 
 if(FILTERING) then
  call MPI_ALLREDUCE(RDUFAPDUFAP_LOC(:,J,K),RDUFAPDUFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDVFAPDVFAP_LOC(:,J,K),RDVFAPDVFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDWFAPDWFAP_LOC(:,J,K),RDWFAPDWFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

  call MPI_ALLREDUCE(RDUFAPUFAP_LOC(:,J,K),RDUFAPUFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDVFAPVFAP_LOC(:,J,K),RDVFAPVFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDWFAPWFAP_LOC(:,J,K),RDWFAPWFAP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

  call MPI_ALLREDUCE(RDUFAPUP_LOC(:,J,K),RDUFAPUP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDVFAPVP_LOC(:,J,K),RDVFAPVP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
  call MPI_ALLREDUCE(RDWFAPWP_LOC(:,J,K),RDWFAPWP,DIMLGR,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

  RDUFAPDUFAP = RDUFAPDUFAP / real(NPROC)
  RDVFAPDVFAP = RDVFAPDVFAP / real(NPROC)
  RDWFAPDWFAP = RDWFAPDWFAP / real(NPROC)
 
  RDUFAPUFAP = RDUFAPUFAP / real(NPROC)
  RDVFAPVFAP = RDVFAPVFAP / real(NPROC)
  RDWFAPWFAP = RDWFAPWFAP / real(NPROC)
 
  RDUFAPUP = RDUFAPUP / real(NPROC)
  RDVFAPVP = RDVFAPVP / real(NPROC)
  RDWFAPWP = RDWFAPWP / real(NPROC)
 
 end if
 
 

!!======================================================================
!! 2. Print in file
!!======================================================================
 if(MYID==0) then

!!----------------------------------------------------------------------
!!- 3.1. Define filename
!!----------------------------------------------------------------------
   write(FILENAME,10600)'part_l2_Rp_p',J,'_f',K,'.stat'
   open(unit=600, file=trim(FILENAME), status='replace')

   if(SOLVE_SCALAR) then 
    write(600,20106)
   elseif(FILTERING) then
    write(600,20107)
   else
    write(600,20104)
   end if 


!!----------------------------------------------------------------------
!!- 3.2. Print in file
!!----------------------------------------------------------------------
 if(SOLVE_SCALAR) then

 do I = 1,DIMLGR-1
  write(600,10000)  (I-1)*DTIME0, &  
                            RPX(I), &
		            RPY(I), &
	                    RPZ(I), &
        (RPX(I)+RPY(I)+RPZ(I))/3.0D0, &
                            RFAPX(I), &
		            RFAPY(I), &
	  		    RFAPZ(I), &
  (RFAPX(I)+RFAPY(I)+RFAPZ(I))/3.0D0, &
                              RTP(I), &
			    RTFAP(I), &
		        RTFAPVFAP(I), &
		        RVFAPTFAP(I) 
	                    
  end do 

 elseif(FILTERING) then
 
 do I = 1,DIMLGR-1
  write(600,10000)  (I-1)*DTIME0, &  
                            RPX(I), &
		            RPY(I), &
	                    RPZ(I), &
        (RPX(I)+RPY(I)+RPZ(I))/3.0D0, &
                            RFAPX(I), &
		            RFAPY(I), &
	  		    RFAPZ(I), &
  (RFAPX(I)+RFAPY(I)+RFAPZ(I))/3.0D0, &
                      RDUFAPDUFAP(I), &
                      RDVFAPDVFAP(I), &
	  	      RDWFAPDWFAP(I), &
  (RDUFAPDUFAP(I)+RDVFAPDVFAP(I)+RDWFAPDWFAP(I))/3.0D0, &
                      RDUFAPUFAP(I), &
                      RDVFAPVFAP(I), &
	  	      RDWFAPWFAP(I), &
  (RDUFAPUFAP(I)+RDVFAPVFAP(I)+RDWFAPWFAP(I))/3.0D0, &
                      RDUFAPUP(I), &
                      RDVFAPVP(I), &
	  	      RDWFAPWP(I), &
  (RDUFAPUP(I)+RDVFAPVP(I)+RDWFAPWP(I))/3.0D0                          
	                    
  end do 
 
 
 else

 do I = 1,DIMLGR-1
   write(600,10000) (I-1)*DTIME0,  &
                               RPX(I), &
			       RPY(I), &
			       RPZ(I), &
         (RPX(I)+RPY(I)+RPZ(I))/3.0D0, &
                             RFAPX(I), &
			     RFAPY(I), &
			     RFAPZ(I), &
        (RFAPX(I)+RFAPY(I)+RFAPZ(I))/3.0D0
 end do 


 end if !- SOLVE_SCALAR
 
 close(600)

 
end if !- If: (MYID==0)


end do  !!- Loop: J = 1, NIG

end do !!- Loop: K= 1, NBLGR




!!--------------------------------------------------------------------
10600 format (A,I2.2,A,I2.2,A)
10700 format (A,I2.2,A,I6.6,A,F12.3)
20104 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p')
20106 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p, RTp, RTf@p, RTfvf, RvfTf')
20107 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p, dRfapx, dRfapy, dRfapz, &
& Rdufufx,dRfap, Rdufufy, Rdufufz,Rdufuf, Rdufupx, Rdufupy, Rdufupz, Rdufup')


10000 format (30(e17.7))
10001 format (50(e17.7))

end subroutine PRINT_LAGFUNCTION 
