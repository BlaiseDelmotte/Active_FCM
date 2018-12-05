!!====================================================================
!!
!! 
!! This routine build the Right Hand Side of Navier-Stokes equations
!!
!!
!!====================================================================

subroutine BUILD_RHS

!!====================================================================
!!
!! This FORTRAN file does not use the vorticity formulation
!! for computing the non-linear terms.
!!
!!====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE, only: KX, KY, KZ
use RHS_VARIABLES
use PARAM_PHYS
use WORK_ARRAYS

use P3DFFT

implicit none


!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- squared wavenumber
real(kind=8) :: KAPPA2

!- Index
integer :: I, J, K, IJK
!------------------------------------------------------------------

!!=================================================================
!! 1.  Term: U*U - V*V
!!=================================================================

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = UFLU(I,J,K)*UFLU(I,J,K) - VFLU(I,J,K)*VFLU(I,J,K)
  end do
 end do
end do


!- Synchronize all the process
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR


do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   RHS_UFOU(I,J,K,TN) = -ICMPL*KX(I)*TMPFOU(I,K,J)
  end do
 end do
end do


!!====================================================================
!! 2. Term: U*V
!!====================================================================
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = UFLU(I,J,K)*VFLU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   !- x-momentum equation
   RHS_UFOU(I,J,K,TN) = RHS_UFOU(I,J,K,TN) - ICMPL*KY(J)*TMPFOU(I,J,K)

   !- y-momentum equation
   RHS_VFOU(I,J,K,TN) = RHS_VFOU(I,J,K,TN) - ICMPL*KX(I)*TMPFOU(I,J,K)
  end do
 end do
end do


!!====================================================================
!! 3. Term: U*W
!!====================================================================
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = UFLU(I,J,K)*WFLU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   !- x-momentum equation
   RHS_UFOU(I,J,K,TN) = RHS_UFOU(I,J,K,TN) - ICMPL*KZ(K)*TMPFOU(I,J,K)

   !- z-momentum equation
   RHS_WFOU(I,J,K,TN) = RHS_WFOU(I,J,K,TN) - ICMPL*KX(I)*TMPFOU(I,J,K)
  end do
 end do
end do

!!====================================================================
!! 4. Term: V*W
!!====================================================================
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = VFLU(I,J,K)*WFLU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   !- y-momentum equation
   RHS_VFOU(I,J,K,TN) = RHS_VFOU(I,J,K,TN) - ICMPL*KZ(K)*TMPFOU(I,J,K)

   !- z-momentum equation
   RHS_WFOU(I,J,K,TN) = RHS_WFOU(I,J,K,TN) - ICMPL*KY(J)*TMPFOU(I,J,K)
  end do
 end do
end do



!!====================================================================
!! 5. Term: W*W - V*V
!!====================================================================
do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
   TMPPHY(I,J,K) = WFLU(I,J,K)*WFLU(I,J,K)-VFLU(I,J,K)*VFLU(I,J,K)
  end do
 end do
end do

!- Synchronize all the process
 call MPI_BARRIER(MPI_COMM_WORLD,IERR)

!- Back in physical space
 call P3DFFT_FTRAN_R2C(TMPPHY,TMPFOU,FFTFLAG)       
 TMPFOU(:,:,:) = TMPFOU(:,:,:)*FACTOR

do K = FSTART(3), FEND(3)
 do J = FSTART(2), FEND(2)
  do I = FSTART(1), FEND(1)
   !- z-momentum equation
   RHS_WFOU(I,J,K,TN) = RHS_WFOU(I,J,K,TN) - ICMPL*KZ(K)*TMPFOU(I,J,K)
  end do
 end do
end do


end subroutine BUILD_RHS
