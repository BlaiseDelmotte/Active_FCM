
!!====================================================================
!!
!! 
!!> @brief
!!> Time-advancing particle position
!!
!! Date :  16/01/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADV_PARTICLE_POSITION

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
use FCM_FORCING_VARIABLE, only: FCM_SHEAR
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
!- Test with AB4
real(kind=8) ::DTIME24

!- Index
integer :: IP
!------------------------------------------------------------------



!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if






DTIME2 = DTIME / 2.
DTIME12 = DTIME / 12.
DTIME24 = DTIME / 24.

do IP =1, NPART_FULL

!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------

 if(FLAG_TSCHEME == 1) then

  FCM_XP(IP) = FCM_XP(IP)  + DTIME*FCM_UP(IP,1) 
  FCM_YP(IP) = FCM_YP(IP)  + DTIME*FCM_VP(IP,1) 
  FCM_ZP(IP) = FCM_ZP(IP)  + DTIME*FCM_WP(IP,1) 
  
  
  FCM_XP_NOPER(IP) = FCM_XP_NOPER(IP)  + DTIME*FCM_UP(IP,1) 
  FCM_YP_NOPER(IP) = FCM_YP_NOPER(IP)  + DTIME*FCM_VP(IP,1) 
  FCM_ZP_NOPER(IP) = FCM_ZP_NOPER(IP)  + DTIME*FCM_WP(IP,1) 
  
  
  

!!====================================================================
!! 2. 2nd Order Adam-Bashford time integration
!!====================================================================
 elseif(FLAG_TSCHEME == 2) then

  FCM_XP(IP) = FCM_XP(IP) &
               + DTIME2*(3.*FCM_UP(IP,1) -FCM_UP(IP,2))
  FCM_YP(IP) = FCM_YP(IP) &
               + DTIME2*(3.*FCM_VP(IP,1) -FCM_VP(IP,2))
  FCM_ZP(IP) = FCM_ZP(IP) &
               + DTIME2*(3.*FCM_WP(IP,1) -FCM_WP(IP,2))
               

  FCM_XP_NOPER(IP) = FCM_XP_NOPER(IP) &
               + DTIME2*(3.*FCM_UP(IP,1) -FCM_UP(IP,2))
  FCM_YP_NOPER(IP) = FCM_YP_NOPER(IP) &
               + DTIME2*(3.*FCM_VP(IP,1) -FCM_VP(IP,2))
  FCM_ZP_NOPER(IP) = FCM_ZP_NOPER(IP) &
               + DTIME2*(3.*FCM_WP(IP,1) -FCM_WP(IP,2))

!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================

 elseif(FLAG_TSCHEME == 3 ) then

  FCM_XP(IP) = FCM_XP(IP)               &
         +  DTIME12*( 23.* FCM_UP(IP,1)	&
                     -16.* FCM_UP(IP,2) &
                     + 5.* FCM_UP(IP,3) )
  FCM_YP(IP) = FCM_YP(IP)               &
         +  DTIME12*( 23.* FCM_VP(IP,1)	&
                     -16.* FCM_VP(IP,2) &
                     + 5.* FCM_VP(IP,3) )
  FCM_ZP(IP) = FCM_ZP(IP)               &
         +  DTIME12*( 23.* FCM_WP(IP,1)	&
                     -16.* FCM_WP(IP,2) &
                     + 5.* FCM_WP(IP,3) ) 


  FCM_XP_NOPER(IP) = FCM_XP_NOPER(IP)   &
         +  DTIME12*( 23.* FCM_UP(IP,1)	&
                     -16.* FCM_UP(IP,2) &
                     + 5.* FCM_UP(IP,3) )
  FCM_YP_NOPER(IP) = FCM_YP_NOPER(IP)   &
         +  DTIME12*( 23.* FCM_VP(IP,1)	&
                     -16.* FCM_VP(IP,2) &
                     + 5.* FCM_VP(IP,3) )
  FCM_ZP_NOPER(IP) = FCM_ZP_NOPER(IP)   &
         +  DTIME12*( 23.* FCM_WP(IP,1)	&
                     -16.* FCM_WP(IP,2) &
                     + 5.* FCM_WP(IP,3) )


!!====================================================================
!! 4. 4th Order Adam-Bashforth time integration
!!====================================================================

 elseif(FLAG_TSCHEME == 4 ) then

  FCM_XP(IP) = FCM_XP(IP)               &
         +  DTIME24*( 55.* FCM_UP(IP,1)	&
                     -59.* FCM_UP(IP,2) &
                     +37.* FCM_UP(IP,3) &
                     - 9.* FCM_UP(IP,4) )
                     
  FCM_YP(IP) = FCM_YP(IP)               &
         +  DTIME24*( 55.* FCM_VP(IP,1)	&
                     -59.* FCM_VP(IP,2) &
                     +37.* FCM_VP(IP,3) &
                     - 9.* FCM_VP(IP,4) )
                     
  FCM_ZP(IP) = FCM_ZP(IP)               &
         +  DTIME24*( 55.* FCM_WP(IP,1)	&
                     -59.* FCM_WP(IP,2) &
                     +37.* FCM_WP(IP,3) &
                     - 9.* FCM_WP(IP,4) )


  FCM_XP_NOPER(IP) = FCM_XP_NOPER(IP)   &
         +  DTIME24*( 55.* FCM_UP(IP,1)	&
                     -59.* FCM_UP(IP,2) &
                     +37.* FCM_UP(IP,3) &
                     - 9.* FCM_UP(IP,4) )
                     
  FCM_YP_NOPER(IP) = FCM_YP_NOPER(IP)   &
         +  DTIME24*( 55.* FCM_VP(IP,1)	&
                     -59.* FCM_VP(IP,2) &
                     +37.* FCM_VP(IP,3) &
                     - 9.* FCM_VP(IP,4) )
                     
  FCM_ZP_NOPER(IP) = FCM_ZP_NOPER(IP)   &
         +  DTIME24*( 55.* FCM_WP(IP,1)	&
                     -59.* FCM_WP(IP,2) &
                     +37.* FCM_WP(IP,3) &
                     - 9.* FCM_WP(IP,4) )

 end if !!- If: FLAG_TSCHEME



 !- Periodicity 
 !- WOULD BE BETTER TO USE KEAVENY CODE FORMULA
 if (FCM_XP(IP).gt.LXMAX) then 
  FCM_XP(IP) = FCM_XP(IP) - LXMAX  
 end if
 
 if (FCM_XP(IP).lt.0.0) then 
  FCM_XP(IP) = FCM_XP(IP) + LXMAX  
 end if
  
 if (FCM_YP(IP).gt.LYMAX) then 
  FCM_YP(IP) = FCM_YP(IP) - LYMAX    
 end if
 
 if (FCM_YP(IP).lt.0.0) then 
  FCM_YP(IP) = FCM_YP(IP) + LYMAX    
 end if 
 
 if (FCM_ZP(IP).gt.LZMAX) then 
  FCM_ZP(IP) = FCM_ZP(IP) - LZMAX
  FCM_UP(IP,1) = FCM_UP(IP,1) - FCM_SHEAR*LZMAX  
  FCM_UP(IP,2) = FCM_UP(IP,2) - FCM_SHEAR*LZMAX  
  FCM_UP(IP,3) = FCM_UP(IP,3) - FCM_SHEAR*LZMAX  
  FCM_UP(IP,4) = FCM_UP(IP,4) - FCM_SHEAR*LZMAX  
 end if 
 
 if (FCM_ZP(IP).lt.0.0) then 
  FCM_ZP(IP) = FCM_ZP(IP) + LZMAX
  FCM_UP(IP,1) =  FCM_UP(IP,1) + FCM_SHEAR*LZMAX  
  FCM_UP(IP,2) =  FCM_UP(IP,2) + FCM_SHEAR*LZMAX
  FCM_UP(IP,3) =  FCM_UP(IP,3) + FCM_SHEAR*LZMAX
  FCM_UP(IP,4) =  FCM_UP(IP,4) + FCM_SHEAR*LZMAX
 end if 
 
 
!~  !------------ONLY FOR 2D-----------!
if (FCM_INIT_PART_POS==4) then
 FCM_YP(IP) = LYMAX /2.0
end if
 

end do !!- Loop: 1, NPART_FULL



!- Update particle velocity history
FCM_UP(:,4) = FCM_UP(:,3)
FCM_UP(:,3) = FCM_UP(:,2)
FCM_UP(:,2) = FCM_UP(:,1)

FCM_VP(:,4) = FCM_VP(:,3)
FCM_VP(:,3) = FCM_VP(:,2)
FCM_VP(:,2) = FCM_VP(:,1)

FCM_WP(:,4) = FCM_WP(:,3)
FCM_WP(:,3) = FCM_WP(:,2)
FCM_WP(:,2) = FCM_WP(:,1)







!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()
end if




end subroutine FCM_ADV_PARTICLE_POSITION
