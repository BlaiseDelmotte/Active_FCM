 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating the particle positions and orientation from elli
!!> psoids algorithm for both spheres and ellipsoids
!!
!! Date :  10/07/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PARTICLE_POS_ORIENT

!!====================================================================
!! Particles are initiated 
!! according to the prescribed size of the particle
!!====================================================================
!! Particle initiation: 
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE


implicit none


real(kind=8), dimension(NPART_FULL,3) :: FCM_RADII

! Indices for particle initiation
integer :: NPART_START, NPART_END

! Indices for loops
integer :: IP, IND

!================ START ROUTINE ========================================


IND = 0

do IP = 1,NPART_FULL
 
!-  Save all radii in a temporary variable
 if (IP.le.FCM_NSPHERE) then
  FCM_RADII(IP,1) = FCM_SPHERE_RADP(IP)
  FCM_RADII(IP,2) = FCM_SPHERE_RADP(IP)
  FCM_RADII(IP,3) = FCM_SPHERE_RADP(IP)
 else
  IND = IND + 1 
  FCM_RADII(IP,1:3) = FCM_ELLIPSOID_RADP(IND,1:3)
 end if

end do


if (FCM_NSPHERE.gt.0) then
 NPART_START = 1
 NPART_END = FCM_NSPHERE

 if (FCM_INIT_PART_POS<4) then
  call INIT_ELLIPSOID( NPART_START,&                !Index of 1st particle to initiate
                       NPART_END,&                  !Index of last particle to initiate
                       FCM_RADII(1:NPART_FULL,1),&  !Radii in 1st principal direction
                       FCM_RADII(1:NPART_FULL,2),&  !Radii in 2nd principal direction
                       FCM_RADII(1:NPART_FULL,3),&  !Radii in 3rd principal direction
                       FCM_XP(1:NPART_FULL),&       !Particle x-positions
                       FCM_YP(1:NPART_FULL),&       !Particle y-positions
                       FCM_ZP(1:NPART_FULL),&       !Particle z-positions
                       FCM_QUAT(1:NPART_FULL,1),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,2),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,3),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,4) )   !Particle orientations

 else
                       
  call INIT_ELLIPSOID_2D( NPART_START,&                !Index of 1st particle to initiate
                          NPART_END,&                  !Index of last particle to initiate
                          FCM_RADII(1:NPART_FULL,1),&  !Radii in 1st principal direction
                          FCM_RADII(1:NPART_FULL,2),&  !Radii in 2nd principal direction
                          FCM_RADII(1:NPART_FULL,3),&  !Radii in 3rd principal direction
                          FCM_XP(1:NPART_FULL),&       !Particle x-positions
                          FCM_YP(1:NPART_FULL),&       !Particle y-positions
                          FCM_ZP(1:NPART_FULL),&       !Particle z-positions
                          FCM_QUAT(1:NPART_FULL,1),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,2),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,3),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,4) )   !Particle orientations                
 end if
end if  
                   
if (FCM_NELLIPSOID.gt.0) then   

               
 NPART_START = FCM_NSPHERE + 1
 NPART_END = NPART_FULL                 

 if ((FCM_INIT_PART_POS<4).AND.(FCM_INIT_PART_ORIENT<4)) then
 
  call INIT_ELLIPSOID( NPART_START,&                !Index of 1st particle to initiate
                       NPART_END,&                  !Index of last particle to initiate
                       FCM_RADII(1:NPART_FULL,1),&  !Radii in 1st principal direction
                       FCM_RADII(1:NPART_FULL,2),&  !Radii in 2nd principal direction
                       FCM_RADII(1:NPART_FULL,3),&  !Radii in 3rd principal direction
                       FCM_XP(1:NPART_FULL),&       !Particle x-positions
                       FCM_YP(1:NPART_FULL),&       !Particle y-positions
                       FCM_ZP(1:NPART_FULL),&       !Particle z-positions
                       FCM_QUAT(1:NPART_FULL,1),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,2),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,3),&   !Particle orientations
                       FCM_QUAT(1:NPART_FULL,4) )   !Particle orientations    
 else
                      
  call INIT_ELLIPSOID_2D( NPART_START,&                !Index of 1st particle to initiate
                          NPART_END,&                  !Index of last particle to initiate
                          FCM_RADII(1:NPART_FULL,1),&  !Radii in 1st principal direction
                          FCM_RADII(1:NPART_FULL,2),&  !Radii in 2nd principal direction
                          FCM_RADII(1:NPART_FULL,3),&  !Radii in 3rd principal direction
                          FCM_XP(1:NPART_FULL),&       !Particle x-positions
                          FCM_YP(1:NPART_FULL),&       !Particle y-positions
                          FCM_ZP(1:NPART_FULL),&       !Particle z-positions
                          FCM_QUAT(1:NPART_FULL,1),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,2),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,3),&   !Particle orientations
                          FCM_QUAT(1:NPART_FULL,4) )   !Particle orientations                
 end if                     
end if                 


! Initial non-periodic posistions = initial positions
FCM_XP_NOPER(1:NPART_FULL) = FCM_XP(1:NPART_FULL)
FCM_YP_NOPER(1:NPART_FULL) = FCM_YP(1:NPART_FULL)
FCM_ZP_NOPER(1:NPART_FULL) = FCM_ZP(1:NPART_FULL)


 ! Choose orientation depending on the eccentricity of the ellipsoids (if any)
 if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then
  do IP = 1, NPART_FULL
   FCM_PSWIM(IP,1) = 2.0*( FCM_QUAT(IP,2) * FCM_QUAT(IP,4) &
                         - FCM_QUAT(IP,1) * FCM_QUAT(IP,3) )
                           
   FCM_PSWIM(IP,2) = 2.0*( FCM_QUAT(IP,3) * FCM_QUAT(IP,4) &
                         + FCM_QUAT(IP,1) * FCM_QUAT(IP,2) )

   FCM_PSWIM(IP,3) = 2.0*( FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,4)**2 - 0.5 )
  end do 
 else
  do IP = 1, NPART_FULL
   FCM_PSWIM(IP,1) = 2.0*( FCM_QUAT(IP,1)**2 + &
                           FCM_QUAT(IP,2)**2 - 0.5d0 )
   FCM_PSWIM(IP,2) = 2.0*( FCM_QUAT(IP,2)*FCM_QUAT(IP,3) - &
                           FCM_QUAT(IP,1)*FCM_QUAT(IP,4) )
   FCM_PSWIM(IP,3) = 2.0*( FCM_QUAT(IP,2)*FCM_QUAT(IP,4) + &
                           FCM_QUAT(IP,1)*FCM_QUAT(IP,3) )
   if (FCM_USE_QUAT == 0) then                         
    FCM_P2(IP,1) = 2.0* ( FCM_QUAT(IP,2)*FCM_QUAT(IP,3) + FCM_QUAT(IP,1)*FCM_QUAT(IP,4) )
    FCM_P2(IP,2) = 2.0* ( FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,3)**2 -0.5  )
    FCM_P2(IP,3) = 2.0* ( FCM_QUAT(IP,3)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,2) )  
    
    FCM_P3(IP,1) = 2.0* ( FCM_QUAT(IP,2)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,3) )
    FCM_P3(IP,2) = 2.0* ( FCM_QUAT(IP,3)*FCM_QUAT(IP,4) + FCM_QUAT(IP,1)*FCM_QUAT(IP,2) ) 
    FCM_P3(IP,3) = 2.0* ( FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,4)**2 -0.5  ) 
   end if
  end do
 end if

if (FCM_NELLIPSOID>0) then

 do IP = 1, NPART_FULL
 !- Transformation matrix calculated from Nikravesh 1985 Eq 7
   FCM_ROT_MAT(IP,1,1) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,2)**2 -0.5 
   FCM_ROT_MAT(IP,1,2) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,3) - FCM_QUAT(IP,1)*FCM_QUAT(IP,4) 
   FCM_ROT_MAT(IP,1,3) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,4) + FCM_QUAT(IP,1)*FCM_QUAT(IP,3) 
   
   FCM_ROT_MAT(IP,2,1) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,3) + FCM_QUAT(IP,1)*FCM_QUAT(IP,4)  
   FCM_ROT_MAT(IP,2,2) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,3)**2 -0.5 
   FCM_ROT_MAT(IP,2,3) =  FCM_QUAT(IP,3)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,2)  
   
   FCM_ROT_MAT(IP,3,1) =  FCM_QUAT(IP,2)*FCM_QUAT(IP,4) - FCM_QUAT(IP,1)*FCM_QUAT(IP,3) 
   FCM_ROT_MAT(IP,3,2) =  FCM_QUAT(IP,3)*FCM_QUAT(IP,4) + FCM_QUAT(IP,1)*FCM_QUAT(IP,2) 
   FCM_ROT_MAT(IP,3,3) =  FCM_QUAT(IP,1)**2 + FCM_QUAT(IP,4)**2 -0.5 
  end do
  
  FCM_ROT_MAT = 2.0*FCM_ROT_MAT
  
end if
 

   
end subroutine FCM_INITIATION_PARTICLE_POS_ORIENT
