!!=====================================================================
!!
!!   Arrays allocation
!!
!!=====================================================================
subroutine ALLOCATE_ARRAYS
!!=====================================================================
!!
!!
!!=====================================================================

use DNS_DIM
use FLUID_VARIABLE
use GEOMETRIC_VARIABLE
use RHS_VARIABLES
use PARAM_PHYS 
use FORCING
use STATISTICS
use WORK_ARRAYS
use SCALAR_VARIABLE
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none


!---------------------------------------------------------------------
integer :: K
!---------------------------------------------------------------------

!!=====================================================================
!! 1. Allocate arrays for the Fluid DNS
!!=====================================================================

!!---------------------------------------------------------------------
!! 1.1. Mesh and Wavenumber
!!---------------------------------------------------------------------

!- Allocate arrays for mesh
allocate(XMESH(ISTART(1):IEND(1)))
allocate(YMESH(ISTART(2):IEND(2)))
allocate(ZMESH(ISTART(3):IEND(3)))


if(SOLVE_FLUID >0) then
 !- Wavenumbers
 allocate(KX(FSTART(1):FEND(1)))
 allocate(KY(FSTART(2):FEND(2)))
 allocate(KZ(FSTART(3):FEND(3)))
end if

!!---------------------------------------------------------------------
!! 1.2. Arrays for the fluid
!!---------------------------------------------------------------------
if(SOLVE_FLUID >0) then

!- Fluid velocity in real space
allocate(UFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

allocate(VFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

allocate(WFLU(ISTART(1)         :IEND(1)             &
             ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
             ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL)   )

!- Fluid velocity in Fourier space
allocate( UFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)),  &
          VFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)),  &
          WFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3))   )


allocate(TMPPHY(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
allocate(TMPFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))



!!- Case of Navier-Stokes equations
if(SOLVE_FLUID==1) then

!- Integrating factor
 allocate(FILTER(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))

!- Right-Hand-Side --> Three times are stored
 allocate(RHS_UFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),3),  &
          RHS_VFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),3),  &
          RHS_WFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),3)   )


end if !- end if: SOLVE_FLUID==1

!!- Case of FCM with Stokes equations
if (SOLVE_FLUID==2)  then

!- Right-Hand-Side --> 1 time is stored
 allocate(RHS_UFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),1),  &
          RHS_VFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),1),  &
          RHS_WFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),1)   )

allocate(TMPPHY2(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
allocate(TMPPHY3(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))

allocate(FCM_FORCING_X(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
allocate(FCM_FORCING_Y(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
allocate(FCM_FORCING_Z(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))

 if (FCM_FLOW_TYPE==2) then
  allocate(FCM_RS_XX(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
  allocate(FCM_RS_YY(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
  allocate(FCM_RS_ZZ(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
  allocate(FCM_RS_XY(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
  allocate(FCM_RS_XZ(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))
  allocate(FCM_RS_YZ(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)))  
  
  allocate(FCM_RS_XX_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))
  allocate(FCM_RS_YY_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))
  allocate(FCM_RS_ZZ_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))
  allocate(FCM_RS_XY_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))
  allocate(FCM_RS_XZ_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))
  allocate(FCM_RS_YZ_FOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3)))  
  
 end if

end if !- end if: SOLVE_FLUID==2




end if !- end if: SOLVE_FLUID >0



!!---------------------------------------------------------------------
!! 1.4. Arrays for statistics
!!---------------------------------------------------------------------
if(LEVEL0_STFLU) then
 allocate(MEAN_FLUID(NSTAT))
 if(STAT_TIME) allocate(MEAN_TIME_FLUID(NSTAT))
end if

if(LEVEL0_STSCL) then
 allocate(MEAN_SCL(NSTAT))
 if(STAT_TIME) allocate(MEAN_TIME_SCL(NSTAT))
end if

!!--------------------------------------------------------------------
!! Spatial correlation
!!--------------------------------------------------------------------
if(STAT_TIME.and.LEVEL4_STFLU) then
 DIMSCOR = ISIZE(1)/2+1
 
 allocate(MEAN_RUXLOC(DIMSCOR))
 allocate(MEAN_RVXLOC(DIMSCOR))
 allocate(MEAN_RWXLOC(DIMSCOR))

 MEAN_RUXLOC(:) = ZERO
 MEAN_RVXLOC(:) = ZERO
 MEAN_RWXLOC(:) = ZERO
 
end if

!!---------------------------------------------------------------------
!! 2.2. Arrays for FCM particles
!!---------------------------------------------------------------------
if (SOLVE_FLUID ==2)  then


 !! 2.2.1. Common Gaussian parameters
 allocate(FCM_LHNODE(NPART_FULL,3)) 
  
 
 !! 2.2.1a. Sphere Gaussian parameters
 if (FCM_NSPHERE.gt.0) then
  
  allocate(FCM_SPHERE_ANORM(FCM_NSPHERE), &
           FCM_SPHERE_ANORM_DIP(FCM_NSPHERE)) 
  
  allocate(FCM_SPHERE_IXP(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_IYP(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_IZP(FCM_NSPHERE,FCM_NGD_MAX))  
           
  allocate(FCM_SPHERE_SIG2SQ(FCM_NSPHERE), &
           FCM_SPHERE_SIGSQ_DIP(FCM_NSPHERE), &
           FCM_SPHERE_SIG2SQ_DIP(FCM_NSPHERE)) 
           
  allocate(FCM_SPHERE_RADP(FCM_NSPHERE), &
           FCM_SPHERE_SIGMA(FCM_NSPHERE), &
           FCM_SPHERE_SIGMA_DIP(FCM_NSPHERE)) 
           
  allocate(FCM_SPHERE_GAUSS1(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_GAUSS2(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_GAUSS3(FCM_NSPHERE,FCM_NGD_MAX)) 
           
  allocate(FCM_SPHERE_DIP_GAUSS1(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_DIP_GAUSS2(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_DIP_GAUSS3(FCM_NSPHERE,FCM_NGD_MAX))  
           
  allocate(FCM_SPHERE_GRAD_DIP_GAUSS1(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_GRAD_DIP_GAUSS2(FCM_NSPHERE,FCM_NGD_MAX), &
           FCM_SPHERE_GRAD_DIP_GAUSS3(FCM_NSPHERE,FCM_NGD_MAX)) 
  
  !! 2.2.1b. Sphere Squirming Gaussian parameters
  
  if (FCM_NSWIM(2)>0) then 
  
   allocate(FCM_SPHERE_SIGSQ(FCM_NSWIM(2)))
   
   allocate(FCM_SPHERE_ANORM_POTDIP(FCM_NSWIM(2)))
    
            
   allocate(FCM_SPHERE_GRAD_GAUSS1_SQ(FCM_NSWIM(2),FCM_NGD_MAX), &
            FCM_SPHERE_GRAD_GAUSS2_SQ(FCM_NSWIM(2),FCM_NGD_MAX), &
            FCM_SPHERE_GRAD_GAUSS3_SQ(FCM_NSWIM(2),FCM_NGD_MAX))  
            
            
   allocate(FCM_COEFF_SPHERE_POTDIP_GAUSS1_SQ(FCM_NSWIM(2),FCM_NGD_MAX), &
            FCM_COEFF_SPHERE_POTDIP_GAUSS2_SQ(FCM_NSWIM(2),FCM_NGD_MAX), &
            FCM_COEFF_SPHERE_POTDIP_GAUSS3_SQ(FCM_NSWIM(2),FCM_NGD_MAX))
  end if 
  
 end if
 
 !! 2.2.1c. Ellipsoid Gaussian parameters
 if (FCM_NELLIPSOID.gt.0) then
  allocate(FCM_ELLIPSOID_ANORM(FCM_NELLIPSOID), &
           FCM_ELLIPSOID_ANORM_DIP(FCM_NELLIPSOID)) 
           
  allocate(FCM_ELLIPSOID_IXP(FCM_NELLIPSOID,FCM_ELLIPSOID_NGD), &
           FCM_ELLIPSOID_IYP(FCM_NELLIPSOID,FCM_ELLIPSOID_NGD), &
           FCM_ELLIPSOID_IZP(FCM_NELLIPSOID,FCM_ELLIPSOID_NGD))  
           
  allocate(FCM_ELLIPSOID_SIG2SQ(FCM_NELLIPSOID,3), &
           FCM_ELLIPSOID_SIGSQ_DIP(FCM_NELLIPSOID,3), &
           FCM_ELLIPSOID_SIG2SQ_DIP(FCM_NELLIPSOID,3)) 
           
  allocate(FCM_ELLIPSOID_RADP(FCM_NELLIPSOID,3), &
           FCM_ELLIPSOID_SIGMA(FCM_NELLIPSOID,3), &
           FCM_ELLIPSOID_SIGMA_DIP(FCM_NELLIPSOID,3)) 
           
  
  if (FCM_NSWIM(2)>0) then 
  
   allocate(FCM_ELLIPSOID_SIGSQ(FCM_NSWIM(2),3))
   
 
                                                 
                                                 
  end if                                                                                             
  
 end if
 
 !! 2.2.2. Particle positions
 allocate(FCM_XP(NPART_FULL), FCM_YP(NPART_FULL), FCM_ZP(NPART_FULL))
 allocate(FCM_XP_NOPER(NPART_FULL), FCM_YP_NOPER(NPART_FULL), FCM_ZP_NOPER(NPART_FULL))

 if (FCM_FLOW_TYPE==2) then
  allocate(FCM_XP0(NPART_FULL), FCM_YP0(NPART_FULL), FCM_ZP0(NPART_FULL))
  allocate(FCM_PSWIM0(NPART_FULL,3))
 end if

 !! 2.2.3. Particle orientations using quaternions or unit vectors, 3 times are stored
 allocate(FCM_QUAT(NPART_FULL,4))
 allocate(FCM_PSWIM(NPART_FULL,3))

 if (FCM_USE_QUAT == 1) then
  allocate(FCM_OM_VEC_Q(NPART_FULL,4,3))
 else  
  allocate(FCM_P2(NPART_FULL,3))
  allocate(FCM_P3(NPART_FULL,3))
  
  ! Test with AB4
  allocate(FCM_OM_VEC_PSWIM(NPART_FULL,3,4))
 end if

 
 if (FCM_NELLIPSOID>0) then
  allocate(FCM_ROT_MAT(NPART_FULL,3,3))
 end if
 
 !! 2.2.4. Particle velocities, 3 times are stored
 
 allocate(FCM_UP(NPART_FULL,4), FCM_VP(NPART_FULL,4), FCM_WP(NPART_FULL,4)) 
 allocate(FCM_OMPX(NPART_FULL), FCM_OMPY(NPART_FULL), FCM_OMPZ(NPART_FULL)) 
 
 if (FCM_FLOW_TYPE==2) then
  allocate(FCM_UP0(NPART_FULL), FCM_VP0(NPART_FULL), FCM_WP0(NPART_FULL)) 
  allocate(FCM_OMPX0(NPART_FULL), FCM_OMPY0(NPART_FULL), FCM_OMPZ0(NPART_FULL)) 
 end if
 
 allocate(FCM_UP_TEMP(NPART_FULL), FCM_VP_TEMP(NPART_FULL), FCM_WP_TEMP(NPART_FULL)) 
 allocate(FCM_OMPX_TEMP(NPART_FULL), FCM_OMPY_TEMP(NPART_FULL), FCM_OMPZ_TEMP(NPART_FULL)) 
 
  !! 2.2.5. Particle swimming velocities + stresslets if time-dependent swimming
  if (FCM_NSWIM(1)>0) then
   allocate(FCM_VSW(FCM_NSWIM(1)))
  end if
  
  if (FCM_NSWIM(3)>0) then
   allocate(FCM_VSW_TIME(FCM_NSWIM(3)))
   allocate(FCM_B2_TIME(FCM_NSWIM(3)))
   allocate(FCM_PHASE_SHIFT(FCM_NSWIM(3)))
  end if
 
 !! 2.2.6. Particle rate of strains + stresslets coeffs, only 5 coeff needed
 if (FCM_ACTIVATE_STRESSLET>0) then
  
  allocate(FCM_EIJ(FCM_ACTIVATE_STRESSLET,5))   
  allocate(FCM_EIJ_TEMP(FCM_ACTIVATE_STRESSLET,5))   
  allocate(FCM_SIJ(FCM_ACTIVATE_STRESSLET,5)) 
  allocate(FCM_DIR_MIN(FCM_ACTIVATE_STRESSLET,5)) 
  allocate(FCM_RES_ROS(FCM_ACTIVATE_STRESSLET,5)) 
  
 end if
 
 if (FCM_NSWIM(1)>0) then  
  allocate(FCM_SPIJ(FCM_NSWIM(1),5))    
   
  allocate(FCM_INT_R11(FCM_NSWIM(1),3,3), &
           FCM_INT_R12(FCM_NSWIM(1),3,3), &
           FCM_INT_R13(FCM_NSWIM(1),3,3), &
           FCM_INT_R22(FCM_NSWIM(1),3,3), &
           FCM_INT_R23(FCM_NSWIM(1),3,3), &
           FCM_INT_R33(FCM_NSWIM(1),3,3)  )
           
   if (FCM_NELLIPSOID.gt.0) then
    allocate(FCM_A1111(FCM_NSWIM(1)))
    allocate(FCM_A2222(FCM_NSWIM(1)))
    allocate(FCM_A1313(FCM_NSWIM(1)))
    allocate(FCM_A2323(FCM_NSWIM(1)))   
   end if
  
  if (FCM_NSWIM(2)>0) then
   allocate(FCM_INT_P(FCM_NSWIM(2),3,3))
   allocate(FCM_HI(FCM_NSWIM(2),3))
   
   if (FCM_NELLIPSOID.gt.0) then
    allocate(FCM_P11(FCM_NSWIM(2)))
    allocate(FCM_P22(FCM_NSWIM(2)))
   end if
  end if    
 end if


   
 !! 2.2.6. Force monopoles
 allocate(FCM_FORCE(NPART_FULL,3))
 allocate(FCM_FORCE_TEMP(NPART_FULL,3)) 
 if (FCM_BC==2) then
  allocate(FCM_FORCE_TEMP2(NPART_FULL,3)) 
 end if
 if (FCM_FLOW_TYPE==2) then
  allocate(FCM_RFD_RAND_FORCE(NPART_FULL,3))
 end if
   
 !! 2.2.6. Torques
 allocate(FCM_TORQUE(NPART_FULL,3))   
 allocate(FCM_AIJ(NPART_FULL,5)) 


end if

!!=====================================================================
!! 3. Allocate arrays for scalar
!!=====================================================================
if(SOLVE_SCALAR.and.(SOLVE_FLUID==1)) then

!- Scalar field in real space
allocate(THETA(ISTART(1)         :IEND(1)             &
              ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL    &
              ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL )  )

!- Scalar in Fourier space
allocate(THETAFOU(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3))   )

!- Right-Hand-Side of Scalar transport equation
allocate(RHS_SCL(FSTART(1):FEND(1),FSTART(2):FEND(2),FSTART(3):FEND(3),3))



end if


!!=====================================================================
if(MYID==0) write(*,*)'Arrays allocation --> OK'

end subroutine ALLOCATE_ARRAYS
