 !!====================================================================
!!
!! 
!!> @brief
!!> Routine filtering the fluid velocity field with the Gaussian monopole
!! enveloppe to get the particle translational velocities
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_VELOCITY_FILTER

!!====================================================================
!! Here thehe fluid velocity field is fieltered with the Gaussian monopole
!! enveloppe to get the particle translational velocities
!!====================================================================
!! Velocity filtering: 
!!------------------------------
!! WARNING: We use temporary variables before performing the sum of each
!! proc contribution
!!------------------------------
!!====================================================================

use WORK_ARRAYS
use FLUID_VARIABLE
use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use P3DFFT

implicit none

! Indices for loops
integer :: IP, K, J, I, IND

! Indices for position
integer :: IX, IY, IZ

! Temporary variables
! Rotation matrix
real(kind=8), dimension(3,3) :: Q
! Position from particle center
real(kind=8), dimension(3,1) :: XLOC
real(kind=8), dimension(3,1) :: XGLOB
real(kind=8) :: FCM_COEFF1, FCM_COEFF2, FCM_COEFF3
real(kind=8) :: FCM_DXDYDZ

! Zeros temp variables
FCM_UP_TEMP(:) = 0.0
FCM_VP_TEMP(:) = 0.0
FCM_WP_TEMP(:) = 0.0

!- dx*dy*dz for numerical Riemann sum
FCM_DXDYDZ = DX*DY*DZ

! Distribution in physical space

!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_NSPHERE

 do K = 1, FCM_NGD(IP)
 
  IZ = FCM_SPHERE_IZP(IP,K)
 
  ! Only compute the velocity average for the piece of gaussian within the local block
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then 
  
   FCM_COEFF1 = FCM_DXDYDZ * FCM_SPHERE_GAUSS3(IP,K)
	
   do J = 1, FCM_NGD(IP)
    
    IY = FCM_SPHERE_IYP(IP,J)
   		
    ! Only compute the velocity average for the piece of gaussian within the local block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
			
     FCM_COEFF2 = FCM_COEFF1 * FCM_SPHERE_GAUSS2(IP,J)
      
     do I = 1, FCM_NGD(IP)
     
      IX = FCM_SPHERE_IXP(IP,I)
      FCM_COEFF3 = FCM_COEFF2 * FCM_SPHERE_GAUSS1(IP,I)
     
      FCM_UP_TEMP(IP) = FCM_UP_TEMP(IP) + FCM_COEFF3 * UFLU(IX,IY,IZ)
      FCM_VP_TEMP(IP) = FCM_VP_TEMP(IP) + FCM_COEFF3 * VFLU(IX,IY,IZ)
      FCM_WP_TEMP(IP) = FCM_WP_TEMP(IP) + FCM_COEFF3 * WFLU(IX,IY,IZ) 									
     end do		
     		
    endif			
   end do		
  endif	
 end do  
end do


!!====================================================================
!! 2. ELLIPSOID GAUSSIAN ENVELOPPE
!!====================================================================

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND +1
  
 !- Transformation matrix calculated from Nikravesh 1985 Eq 7
 Q = FCM_ROT_MAT(IP,1:3,1:3)

 do K = 1, FCM_ELLIPSOID_NGD
 
  IZ = FCM_ELLIPSOID_IZP(IND,K)
  
  ! Only compute the force projection for the piece of gaussian within the block  
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then   
  
   FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH + K      
   XLOC(3,1) = dble(FCM_JZ)*DZ-FCM_ZP(IP)
   
   do J = 1, FCM_ELLIPSOID_NGD
   
    IY = FCM_ELLIPSOID_IYP(IND,J)
    
    ! Only compute the force projection for the piece of gaussian within the block   
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
   
     FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH + J  
     XLOC(2,1) = dble(FCM_JY)*DY-FCM_YP(IP)
       
     do I = 1, FCM_ELLIPSOID_NGD
     
      IX = FCM_ELLIPSOID_IXP(IND,I)
      FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH + I     
      XLOC(1,1) = dble(FCM_JX)*DX-FCM_XP(IP)
      
      XGLOB = matmul(Q,XLOC)
      
      call FCM_COMPUTE_GAUSSIAN_ENV_FUNCTION(IND, XGLOB, FCM_COEFF3)
      FCM_COEFF3 = FCM_DXDYDZ *FCM_COEFF3
    
      FCM_UP_TEMP(IP) = FCM_UP_TEMP(IP) + FCM_COEFF3 * UFLU(IX,IY,IZ)
      FCM_VP_TEMP(IP) = FCM_VP_TEMP(IP) + FCM_COEFF3 * VFLU(IX,IY,IZ)
      FCM_WP_TEMP(IP) = FCM_WP_TEMP(IP) + FCM_COEFF3 * WFLU(IX,IY,IZ) 									
     end do		
      
    endif			
   end do	
  end if	
 end do  
end do



! Simple addition of the contribution of each processor to the velocity average
call MPI_ALLREDUCE(FCM_UP_TEMP,FCM_UP(:,1),NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(FCM_VP_TEMP,FCM_VP(:,1),NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(FCM_WP_TEMP,FCM_WP(:,1),NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!!====================================================================
!! 3. ADD SHEAR INTENSITY IF ANY
!!====================================================================
FCM_UP(:,1) = FCM_UP(:,1) + FCM_SHEAR * FCM_ZP(:)


!!====================================================================
!! 4. ADD INTRINSIC SWIMMING SPEED IF ACTIVATED
!!====================================================================
!- 
if (FCM_NSWIM(1)>0) then
    if (FCM_NSWIM(3)>0) then
     !do IP = 1, NPART_FULL
        FCM_UP(1:FCM_NSWIM(3),1) = FCM_UP(1:FCM_NSWIM(3),1) + FCM_PSWIM(1:FCM_NSWIM(3),1)*FCM_VSW_TIME
        FCM_VP(1:FCM_NSWIM(3),1) = FCM_VP(1:FCM_NSWIM(3),1) + FCM_PSWIM(1:FCM_NSWIM(3),2)*FCM_VSW_TIME
        FCM_WP(1:FCM_NSWIM(3),1) = FCM_WP(1:FCM_NSWIM(3),1) + FCM_PSWIM(1:FCM_NSWIM(3),3)*FCM_VSW_TIME
     !end do
    else
        FCM_UP(1:FCM_NSWIM(1),1) = FCM_UP(1:FCM_NSWIM(1),1) + FCM_PSWIM(1:FCM_NSWIM(1),1)*FCM_VSW
        FCM_VP(1:FCM_NSWIM(1),1) = FCM_VP(1:FCM_NSWIM(1),1) + FCM_PSWIM(1:FCM_NSWIM(1),2)*FCM_VSW
        FCM_WP(1:FCM_NSWIM(1),1) = FCM_WP(1:FCM_NSWIM(1),1) + FCM_PSWIM(1:FCM_NSWIM(1),3)*FCM_VSW      
    end if
end if
   
end subroutine FCM_VELOCITY_FILTER
