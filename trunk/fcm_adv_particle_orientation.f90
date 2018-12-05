
!!====================================================================
!!
!! 
!!> @brief
!!> Time-advancing particle orientation
!!
!! Date :  03/04/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADV_PARTICLE_ORIENTATION

!!=====================================================================
!! The numerical scheme used for time-advancing is a 2nd order
!! Adams-Bashforth coupled with integrating factor.
!! 
!!  df(t)                      d[f(t).exp(at)]
!! ------ + af(t) + b = 0  <=> --------------- + b.exp(at) = 0
!!   dt                               dt
!!             
!! 2nd order Adams-Bashforth is written as:
!!
!!  df(t)                f^n+1 - f^n     3            1
!! ------ = phi(t)  ==> ------------- = ---.phi^n  -.---.phi^n-1 
!!   dt                       Dt         2            2
!!
!!---------------------------------------------------------------------
!! Warning: The variable FLAG_TSCHEME is modified in the main
!!          
!!---------------------------------------------------------------------
!! Notes for further development:
!!------------------------------
!!
!!=====================================================================

use GEOMETRIC_VARIABLE
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE, only:FCM_SWIMMING, FCM_NSWIM
use DNS_DIM
use PARAM_PHYS, only: DTIME 

implicit none


!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Time step
real(kind=8) ::DTIME2, DTIME12

!- Index
integer :: IP, IC
!------------------------------------------------------------------



!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if


DTIME2 = DTIME / 2.
DTIME12 = DTIME / 12.

do IP =1, NPART_FULL 


!- RHS of quaternion integration Eq 59 Nikravesh 1985
! Add "-" to get the same rotation as vectors!!! PROBLEM TO BE SOLVED!!!
 FCM_OM_VEC_Q(IP,1,1) =  -0.5d0*( - FCM_OMPX(IP)*FCM_QUAT(IP,2) & 
                                  - FCM_OMPY(IP)*FCM_QUAT(IP,3) & 
                                  - FCM_OMPZ(IP)*FCM_QUAT(IP,4) )
                              
 FCM_OM_VEC_Q(IP,2,1) =  -0.5d0*(   FCM_OMPX(IP)*FCM_QUAT(IP,1) &
                                  - FCM_OMPZ(IP)*FCM_QUAT(IP,3) &
                                  + FCM_OMPY(IP)*FCM_QUAT(IP,4) )
                                
 FCM_OM_VEC_Q(IP,3,1) =  -0.5d0*(   FCM_OMPY(IP)*FCM_QUAT(IP,1) &
                                  + FCM_OMPZ(IP)*FCM_QUAT(IP,2) &
                                  - FCM_OMPX(IP)*FCM_QUAT(IP,4) )
                                
 FCM_OM_VEC_Q(IP,4,1) =  -0.5d0*(   FCM_OMPZ(IP)*FCM_QUAT(IP,1) &
                                  - FCM_OMPY(IP)*FCM_QUAT(IP,2) &
                                  + FCM_OMPX(IP)*FCM_QUAT(IP,3) )
                                

                             
                                                                
!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
 if(FLAG_TSCHEME == 1) then
  do IC = 1, 4  
   FCM_QUAT(IP,IC) = FCM_QUAT(IP,IC)  + DTIME*FCM_OM_VEC_Q(IP,IC,1)   
  enddo
  

!!====================================================================
!! 2. 2nd Order Adam-Bashford time integration
!!====================================================================
 elseif(FLAG_TSCHEME == 2) then
 
  do IC = 1, 4  
   FCM_QUAT(IP,IC) = FCM_QUAT(IP,IC) + &
                     DTIME2*( 3.*FCM_OM_VEC_Q(IP,IC,1) & 
                               - FCM_OM_VEC_Q(IP,IC,2) )   
  enddo


!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================

 elseif(FLAG_TSCHEME == 3 ) then
 
  do IC = 1, 4  
   FCM_QUAT(IP,IC) = FCM_QUAT(IP,IC) + &
                     DTIME12*( 23.*FCM_OM_VEC_Q(IP,IC,1) &
                              -16.*FCM_OM_VEC_Q(IP,IC,2) &   
                              + 5.*FCM_OM_VEC_Q(IP,IC,3) )
  enddo

 end if !!- If: FLAG_TSCHEME
 
 
 
 do IC = 1, 4 
  !- Update particle rotation history
  FCM_OM_VEC_Q(IP,IC,3) = FCM_OM_VEC_Q(IP,IC,2)
  FCM_OM_VEC_Q(IP,IC,2) = FCM_OM_VEC_Q(IP,IC,1)
 enddo
 
!~  !------------ONLY FOR 2D-----------!
 if (FCM_INIT_PART_ORIENT==4) then
  FCM_QUAT(IP,2)=0.0
  FCM_QUAT(IP,4)=0.0
 end if

 !- Renormalization

 FCM_QUAT(IP,:) = FCM_QUAT(IP,:)/ dsqrt( FCM_QUAT(IP,1)**2 &
                                        +FCM_QUAT(IP,2)**2 &
                                        +FCM_QUAT(IP,3)**2 &
                                        +FCM_QUAT(IP,4)**2 )
                                        

 !- Swimming Direction vector
 !- First row of transformation matrix between global and body-fixed coordinates 
  ! Change in signes 2014/02/05 +/-/+ ->  +/+/- 
   ! Change in signes 2014/02/10  back to +/-/+ as it is the right way 

 
end do !!- Loop: 1, NPART_FULL

if (FCM_NSWIM(1)>0) then
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
                           
  ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
  ! LAST CHANGE : 2014/02/10: back to (-)   
   FCM_PSWIM(IP,2) = 2.0*( FCM_QUAT(IP,2)*FCM_QUAT(IP,3) - &
                           FCM_QUAT(IP,1)*FCM_QUAT(IP,4) )
  ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
  ! LAST CHANGE : 2014/02/10: back to (+)  
   FCM_PSWIM(IP,3) = 2.0*( FCM_QUAT(IP,2)*FCM_QUAT(IP,4) + &
                           FCM_QUAT(IP,1)*FCM_QUAT(IP,3) )
  end do
 end if
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


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()

end if




end subroutine FCM_ADV_PARTICLE_ORIENTATION
