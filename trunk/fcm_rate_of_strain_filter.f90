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

subroutine FCM_RATE_OF_STRAIN_FILTER

!!====================================================================
!! Here thehe fluid vorticity field is computed with VORTICITY
!! and  filtered with the Gaussian dipole enveloppe to get the particle
!! rotational velocities
!!====================================================================
!! Rotationfiltering: 
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
real(kind=8) :: FCM_E13_FLUID_SHEAR

! Zeros temporary variable
FCM_EIJ_TEMP(:,:) = 0.0

!- dx*dy*dz for numerical Riemann sum
FCM_DXDYDZ = DX*DY*DZ

! Distribution in physical space
if (FCM_CONSIDER_ROS_SHEAR==1) then
 !- Add shear intensity 
 FCM_E13_FLUID_SHEAR =  0.5 * FCM_SHEAR 
 
else
 FCM_E13_FLUID_SHEAR = 0.0
end if

!!====================================================================
!! 1. SPHERE GAUSSIAN ENVELOPPE
!!====================================================================
do IP = 1, FCM_ACTIVATE_STRESSLET

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
     
      FCM_EIJ_TEMP(IP,1) = FCM_EIJ_TEMP(IP,1) &
                         - FCM_COEFF3 * UFLU(IX,IY,IZ) * XX
      FCM_EIJ_TEMP(IP,2) = FCM_EIJ_TEMP(IP,2) &
                         - 0.5*FCM_COEFF3 * ( UFLU(IX,IY,IZ) * YY &
                                            + VFLU(IX,IY,IZ) * XX )
      FCM_EIJ_TEMP(IP,3) = FCM_EIJ_TEMP(IP,3) &
                         - 0.5*FCM_COEFF3 * ( UFLU(IX,IY,IZ) * ZZ &
                                            + WFLU(IX,IY,IZ) * XX &
                                            - 2.0*FCM_E13_FLUID_SHEAR )	
      FCM_EIJ_TEMP(IP,4) = FCM_EIJ_TEMP(IP,4) &
                         - FCM_COEFF3 * VFLU(IX,IY,IZ) * YY
      FCM_EIJ_TEMP(IP,5) = FCM_EIJ_TEMP(IP,5) &
                         - 0.5*FCM_COEFF3 * ( VFLU(IX,IY,IZ) * ZZ &
                                            + WFLU(IX,IY,IZ) * YY )
      							
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
                                                 
      call FCM_COMPUTE_GAUSSIAN_GRAD_DIP_ENV_FUNCTION(IND, &
                                                      XGLOB, &
                                                      Q, &                                                
                                                      XX, &
                                                      YY, &
                                                      ZZ  )                                                  
                                                 
      FCM_COEFF3 = FCM_COEFF3 * FCM_DXDYDZ           
      
      FCM_EIJ_TEMP(IP,1) = FCM_EIJ_TEMP(IP,1) &
                         - FCM_COEFF3 * UFLU(IX,IY,IZ) * XX
      FCM_EIJ_TEMP(IP,2) = FCM_EIJ_TEMP(IP,2) &
                         - 0.5*FCM_COEFF3 * ( UFLU(IX,IY,IZ) * YY &
                                            + VFLU(IX,IY,IZ) * XX )
      FCM_EIJ_TEMP(IP,3) = FCM_EIJ_TEMP(IP,3) &
                         - 0.5*FCM_COEFF3 * ( UFLU(IX,IY,IZ) * ZZ &
                                            + WFLU(IX,IY,IZ) * XX &
                                            - 2.0* FCM_E13_FLUID_SHEAR )	
      FCM_EIJ_TEMP(IP,4) = FCM_EIJ_TEMP(IP,4) &
                         - FCM_COEFF3 * VFLU(IX,IY,IZ) * YY
      FCM_EIJ_TEMP(IP,5) = FCM_EIJ_TEMP(IP,5) &
                         - 0.5*FCM_COEFF3 * ( VFLU(IX,IY,IZ) * ZZ &
                                            + WFLU(IX,IY,IZ) * YY )		
      					
     end do		
     		
    endif			
   end do		
  endif	
 end do  
 
end do
!

! Simple addition of the contribution of each processor to the velocity average
call MPI_ALLREDUCE(FCM_EIJ_TEMP,FCM_EIJ,FCM_ACTIVATE_STRESSLET*5,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)


   
end subroutine FCM_RATE_OF_STRAIN_FILTER
