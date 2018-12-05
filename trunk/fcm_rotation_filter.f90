
 !!====================================================================
!!
!! 
!!> @brief
!!> Routine filtering the fluid vorticity field with the Gaussian dipole
!! enveloppe to get the particlerotational velocities
!!
!! Date :  17/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ROTATION_FILTER

!!====================================================================
!! Here thehe fluid vorticity field is computed with VORTICITY
!! and  filtered with the Gaussian dipole enveloppe to get the particle
!! rotational velocities
!!====================================================================
!! Rotationfiltering: 
!!------------------------------
!!
!! TO DO : 
!!        1) If ellipsoid
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

!------------------------------------------------------------------
! ARRAYS STATEMENT
!------------------------------------------------------------------
!- Vorticity

!- Physical Space

!~ real(kind=8),   &
!~    dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAX
!~ real(kind=8),    &
!~    dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAY
!~ real(kind=8),    &
!~    dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: ZETAZ
   

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
real(kind=8) :: XX, YY, ZZ
real(kind=8) :: FCM_DXDYDZ

! Zeros temp variables
FCM_OMPX_TEMP(:) = 0d0
FCM_OMPY_TEMP(:) = 0d0
FCM_OMPZ_TEMP(:) = 0d0


!~ ! Compute vorticity in physical space from UFOU, VFOU, WFOU
!~ call VORTICITY(ZETAX,ZETAY,ZETAZ)

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
  
   FCM_COEFF1 = FCM_DXDYDZ * FCM_SPHERE_DIP_GAUSS3(IP,K)
   ZZ = FCM_SPHERE_GRAD_DIP_GAUSS3(IP,K)
	
   do J = 1, FCM_NGD(IP)
    
    IY = FCM_SPHERE_IYP(IP,J)
   		
    ! Only compute the velocity average for the piece of gaussian within the local block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
			
     FCM_COEFF2 = FCM_COEFF1 * FCM_SPHERE_DIP_GAUSS2(IP,J)
     YY = FCM_SPHERE_GRAD_DIP_GAUSS2(IP,J)
      
     do I = 1, FCM_NGD(IP)
     
      IX = FCM_SPHERE_IXP(IP,I)
      FCM_COEFF3 = FCM_COEFF2 * FCM_SPHERE_DIP_GAUSS1(IP,I)
      XX = FCM_SPHERE_GRAD_DIP_GAUSS1(IP,I)
      
     
      FCM_OMPX_TEMP(IP) = FCM_OMPX_TEMP(IP) &
                        - FCM_COEFF3 * (WFLU(IX,IY,IZ)*YY - VFLU(IX,IY,IZ)*ZZ)
      FCM_OMPY_TEMP(IP) = FCM_OMPY_TEMP(IP) &
                        - FCM_COEFF3 * (UFLU(IX,IY,IZ)*ZZ - WFLU(IX,IY,IZ)*XX)
      FCM_OMPZ_TEMP(IP) = FCM_OMPZ_TEMP(IP) &
                        - FCM_COEFF3 * (VFLU(IX,IY,IZ)*XX - UFLU(IX,IY,IZ)*YY)			
                       
      
     
!~       print*,'K,J,I = ',K,J,I
!~       print*,'FCM_COEFF3 = ',FCM_COEFF3
!~       print*,'XX, YY, ZZ = ',XX, YY, ZZ
!~       print*,'WFLU(IX,IY,IZ) = ',WFLU(IX,IY,IZ)
!~       print*,'VFLU(IX,IY,IZ) = ',VFLU(IX,IY,IZ)
!~       print*,'UFLU(IX,IY,IZ) = ',UFLU(IX,IY,IZ)
!~       print*,'FCM_OMPZ_TEMP(IP) = ',FCM_OMPZ_TEMP(IP)
!~       read(*,*)
     end do		
     		
    endif			
   end do		
  endif	
 end do  
 
 FCM_OMPX_TEMP(IP) = FCM_OMPX_TEMP(IP)  /2.0
 FCM_OMPY_TEMP(IP) = FCM_OMPY_TEMP(IP)  /2.0
 FCM_OMPZ_TEMP(IP) = FCM_OMPZ_TEMP(IP)  /2.0
end do



!!====================================================================
!! 2. ELLIPSOID GAUSSIAN ENVELOPPE
!!====================================================================

IND = 0

do IP = FCM_NSPHERE + 1, NPART_FULL

 IND = IND + 1

 !- Transformation matrix calculated from Nikravesh 1985 Eq 7
 Q = FCM_ROT_MAT(IP,1:3,1:3)

 do K = 1, FCM_ELLIPSOID_NGD
 
  IZ = FCM_ELLIPSOID_IZP(IND,K)
 
  ! Only compute the rotation average for the piece of gaussian within the local block
  if ( (IZ.GE.ISTART(3)) .AND. (IZ.LE.(IEND(3))) ) then 
   FCM_JZ = FCM_LHNODE(IP,3) - FCM_ELLIPSOID_NGDH + K      
   XLOC(3,1) = dble(FCM_JZ)*DZ-FCM_ZP(IP)
  
	
   do J = 1, FCM_ELLIPSOID_NGD
    
    IY = FCM_ELLIPSOID_IYP(IND,J)
   		
    ! Only compute the rotation average for the piece of gaussian within the local block
    if ( (IY.GE.ISTART(2)).AND.(IY.LE.(IEND(2))) ) then 
    
     FCM_JY = FCM_LHNODE(IP,2) - FCM_ELLIPSOID_NGDH + J  
     XLOC(2,1) = dble(FCM_JY)*DY-FCM_YP(IP)
      
     do I = 1, FCM_ELLIPSOID_NGD
     
      IX = FCM_ELLIPSOID_IXP(IND,I)
      FCM_JX = FCM_LHNODE(IP,1) - FCM_ELLIPSOID_NGDH + I     
      XLOC(1,1) = dble(FCM_JX)*DX-FCM_XP(IP)   
      
      XGLOB = matmul(Q,XLOC)
      
      call FCM_COMPUTE_GAUSSIAN_DIP_ENV_FUNCTION(IND, &
                                                 XGLOB, &
                                                 FCM_COEFF3 )
      FCM_COEFF3 = FCM_COEFF3 * FCM_DXDYDZ 
      
      call FCM_COMPUTE_GAUSSIAN_GRAD_DIP_ENV_FUNCTION(IND, &
                                                      XGLOB, &
                                                      Q, &                                                
                                                      XX, &
                                                      YY, &
                                                      ZZ  ) 
      

     
      FCM_OMPX_TEMP(IP) = FCM_OMPX_TEMP(IP) &
                        - FCM_COEFF3 * (WFLU(IX,IY,IZ)*YY - VFLU(IX,IY,IZ)*ZZ)
      FCM_OMPY_TEMP(IP) = FCM_OMPY_TEMP(IP) &
                        - FCM_COEFF3 * (UFLU(IX,IY,IZ)*ZZ - WFLU(IX,IY,IZ)*XX)
      FCM_OMPZ_TEMP(IP) = FCM_OMPZ_TEMP(IP) &
                        - FCM_COEFF3 * (VFLU(IX,IY,IZ)*XX - UFLU(IX,IY,IZ)*YY)		
                        
!~       print*,'K,J,I = ',K,J,I
!~       print*,'FCM_COEFF3 = ',FCM_COEFF3
!~       print*,'XX, YY, ZZ = ',XX, YY, ZZ
!~       print*,'WFLU(IX,IY,IZ) = ',WFLU(IX,IY,IZ)
!~       print*,'VFLU(IX,IY,IZ) = ',VFLU(IX,IY,IZ)
!~       print*,'UFLU(IX,IY,IZ) = ',UFLU(IX,IY,IZ)
!~       print*,'FCM_OMPZ_TEMP(IP) = ',FCM_OMPZ_TEMP(IP)
!~       read(*,*)
     end do		
     		
    endif			
   end do		
  endif	
 end do  
 
 FCM_OMPX_TEMP(IP) = FCM_OMPX_TEMP(IP)  /2.0
 FCM_OMPY_TEMP(IP) = FCM_OMPY_TEMP(IP)  /2.0
 FCM_OMPZ_TEMP(IP) = FCM_OMPZ_TEMP(IP)  /2.0
end do



! Simple addition of the contribution of each processor to the velocity average
call MPI_ALLREDUCE(FCM_OMPX_TEMP,FCM_OMPX,NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(FCM_OMPY_TEMP,FCM_OMPY,NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(FCM_OMPZ_TEMP,FCM_OMPZ,NPART_FULL,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)


!- Add rotation due to shear ux = shear*z
FCM_OMPY = FCM_OMPY + FCM_SHEAR /2.0


end subroutine FCM_ROTATION_FILTER
