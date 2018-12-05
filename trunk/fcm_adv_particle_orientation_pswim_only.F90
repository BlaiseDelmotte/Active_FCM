
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

subroutine FCM_ADV_PARTICLE_ORIENTATION_PSWIM_ONLY(DT_ARG)

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
use FCM_FORCING_VARIABLE
use DNS_DIM
!use PARAM_PHYS, only: DTIME 

implicit none

#define SIMPLE_SPRNG	
#define USE_MPI 
#include "sprng_f.h"

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
real(kind=8), intent(in) :: DT_ARG

!- Arrays for run_tumble
real(kind=8), dimension(NPART_FULL) :: NO_TUMBLE
real(kind=8), dimension(NPART_FULL,3) :: PSWIM_TUMBLE
!- Random number for run and tumble
real(kind=8) :: RANDNUM, RAND1, RAND2, RAND3

!- Time control variable
real(kind=8) :: TIME_START, TIME_END

!- Time step
real(kind=8) ::DTIME2, DTIME12
!- Test with AB4
real(kind=8) ::DTIME24

!- Index
integer :: IP, IC
!------------------------------------------------------------------



!!- CPU check
if(MYID == 0) then
 TIME_START=MPI_WTIME()
end if

NO_TUMBLE = 1.0
PSWIM_TUMBLE = 0.0

if (FCM_RUNTUMBLE==1) then
 if (MYID ==0) then
  do IP = 1, NPART_FULL

   RANDNUM = sprng()

   if (RANDNUM.lt.(DT_ARG/FCM_TAU_RUN)) then
    NO_TUMBLE(IP) = 0.0
    RAND1 = 1.0-2.0*sprng()
    RAND2 = TWOPI*sprng()
    RAND3 = dsqrt(1.0-RAND1**2)
    PSWIM_TUMBLE(IP,1) = RAND3*dcos(RAND2)
    PSWIM_TUMBLE(IP,2) = RAND3*dsin(RAND2)
    PSWIM_TUMBLE(IP,3) = RAND1
   end if
  end do
 end if
 call MPI_BCAST(NO_TUMBLE,NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 call MPI_BCAST(PSWIM_TUMBLE,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
end if


DTIME2 = DT_ARG/ 2.
DTIME12 = DT_ARG / 12.
DTIME24 = DT_ARG / 24.

do IP = 1, NPART_FULL 


!- RHS of pswim
 FCM_OM_VEC_PSWIM(IP,1,1) =   FCM_OMPY(IP)*FCM_PSWIM(IP,3) & 
                            - FCM_OMPZ(IP)*FCM_PSWIM(IP,2)                               
 FCM_OM_VEC_PSWIM(IP,2,1) =   FCM_OMPZ(IP)*FCM_PSWIM(IP,1) &
                            - FCM_OMPX(IP)*FCM_PSWIM(IP,3)                                
 FCM_OM_VEC_PSWIM(IP,3,1) =   FCM_OMPX(IP)*FCM_PSWIM(IP,2) &
                            - FCM_OMPY(IP)*FCM_PSWIM(IP,1) 
                                    
! FCM_OM_VEC_P2(IP,1,1) =  FCM_OMPY(IP)*FCM_P2(IP,3) & 
!                        - FCM_OMPZ(IP)*FCM_P2(IP,2)                               
! FCM_OM_VEC_P2(IP,2,1) =  FCM_OMPZ(IP)*FCM_P2(IP,1) &
!                        - FCM_OMPX(IP)*FCM_P2(IP,3)                                
! FCM_OM_VEC_P2(IP,3,1) =  FCM_OMPX(IP)*FCM_P2(IP,2) &
!                        - FCM_OMPY(IP)*FCM_P2(IP,1) 
!                                    
! FCM_OM_VEC_P3(IP,1,1) =  FCM_OMPY(IP)*FCM_P3(IP,3) & 
!                        - FCM_OMPZ(IP)*FCM_P3(IP,2)                               
! FCM_OM_VEC_P3(IP,2,1) =  FCM_OMPZ(IP)*FCM_P3(IP,1) &
!                        - FCM_OMPX(IP)*FCM_P3(IP,3)                                 
! FCM_OM_VEC_P3(IP,3,1) =  FCM_OMPX(IP)*FCM_P3(IP,2) &
!                        - FCM_OMPY(IP)*FCM_P3(IP,1) 

                                                                
!!====================================================================
!! 1. First time step 
!!====================================================================
!! The Right-Hand Side is filled in order that the Adam-Bashfort 
!! scheme becomes an Euler explicite scheme. The user should be awared
!! that the first time step must be small enough for stability.
!!--------------------------------------------------------------------
 if(FLAG_TSCHEME == 1) then
  do IC = 1, 3  
   FCM_PSWIM(IP,IC) = (FCM_PSWIM(IP,IC)  + DT_ARG*FCM_OM_VEC_PSWIM(IP,IC,1))*NO_TUMBLE(IP) &
                     + PSWIM_TUMBLE(IP,IC)
!   FCM_P2(IP,IC) = FCM_P2(IP,IC)  + DTIME*FCM_OM_VEC_P2(IP,IC,1) 
!   FCM_P3(IP,IC) = FCM_P3(IP,IC)  + DTIME*FCM_OM_VEC_P3(IP,IC,1) 

  enddo
  

!!====================================================================
!! 2. 2nd Order Adam-Bashford time integration
!!====================================================================
 elseif(FLAG_TSCHEME == 2) then
 
  do IC = 1, 3  
   FCM_PSWIM(IP,IC) = FCM_PSWIM(IP,IC) + &
                     DTIME2*( 3.*FCM_OM_VEC_PSWIM(IP,IC,1) & 
                               - FCM_OM_VEC_PSWIM(IP,IC,2) )  
                               
!   FCM_P2(IP,IC) = FCM_P2(IP,IC) + &
!                     DTIME2*( 3.*FCM_OM_VEC_P2(IP,IC,1) & 
!                               - FCM_OM_VEC_P2(IP,IC,2) ) 
!                               
!   FCM_P3(IP,IC) = FCM_P3(IP,IC) + &
!                     DTIME2*( 3.*FCM_OM_VEC_P3(IP,IC,1) & 
!                               - FCM_OM_VEC_P3(IP,IC,2) ) 
        
                            
  enddo


!!====================================================================
!! 3. 3rd Order Adam-Bashforth time integration
!!====================================================================

 elseif(FLAG_TSCHEME == 3 ) then
 
  do IC = 1, 3  
   FCM_PSWIM(IP,IC) = FCM_PSWIM(IP,IC) + &
                     DTIME12*( 23.*FCM_OM_VEC_PSWIM(IP,IC,1) &
                              -16.*FCM_OM_VEC_PSWIM(IP,IC,2) &   
                              + 5.*FCM_OM_VEC_PSWIM(IP,IC,3) )
                              
!   FCM_P2(IP,IC) = FCM_P2(IP,IC) + &
!                     DTIME12*( 23.*FCM_OM_VEC_P2(IP,IC,1) &
!                              -16.*FCM_OM_VEC_P2(IP,IC,2) &   
!                              + 5.*FCM_OM_VEC_P2(IP,IC,3) )
!                              
!   FCM_P3(IP,IC) = FCM_P3(IP,IC) + &
!                     DTIME12*( 23.*FCM_OM_VEC_P3(IP,IC,1) &
!                              -16.*FCM_OM_VEC_P3(IP,IC,2) &   
!                              + 5.*FCM_OM_VEC_P3(IP,IC,3) )                    

                             
  enddo
  
!!====================================================================
!! 4. 4th Order Adam-Bashforth time integration
!!====================================================================  
 elseif(FLAG_TSCHEME == 4 ) then
 
  do IC = 1, 3  
   FCM_PSWIM(IP,IC) = FCM_PSWIM(IP,IC) + &
                     DTIME24*( 55.*FCM_OM_VEC_PSWIM(IP,IC,1) &
                              -59.*FCM_OM_VEC_PSWIM(IP,IC,2) &   
                              +37.*FCM_OM_VEC_PSWIM(IP,IC,3) &
                              - 9.*FCM_OM_VEC_PSWIM(IP,IC,4))
                              
!   FCM_P2(IP,IC) = FCM_P2(IP,IC) + &
!                     DTIME24*( 55.*FCM_OM_VEC_P2(IP,IC,1) &
!                              -59.*FCM_OM_VEC_P2(IP,IC,2) &   
!                              +37.*FCM_OM_VEC_P2(IP,IC,3) &
!                              - 9.*FCM_OM_VEC_P2(IP,IC,4) )
!                              
!   FCM_P3(IP,IC) = FCM_P3(IP,IC) + &
!                     DTIME24*( 55.*FCM_OM_VEC_P3(IP,IC,1) &
!                              -59.*FCM_OM_VEC_P3(IP,IC,2) & 
!                              +37.*FCM_OM_VEC_P3(IP,IC,3) &  
!                              - 9.*FCM_OM_VEC_P3(IP,IC,4) )                    

                             
  enddo
 end if !!- If: FLAG_TSCHEME
 
 
 
 do IC = 1, 3 
  !- Update particle rotation history
  FCM_OM_VEC_PSWIM(IP,IC,4) = FCM_OM_VEC_PSWIM(IP,IC,3)
  FCM_OM_VEC_PSWIM(IP,IC,3) = FCM_OM_VEC_PSWIM(IP,IC,2)
  FCM_OM_VEC_PSWIM(IP,IC,2) = FCM_OM_VEC_PSWIM(IP,IC,1)
  
!  FCM_OM_VEC_P2(IP,IC,4) = FCM_OM_VEC_P2(IP,IC,3)
!  FCM_OM_VEC_P2(IP,IC,3) = FCM_OM_VEC_P2(IP,IC,2)
!  FCM_OM_VEC_P2(IP,IC,2) = FCM_OM_VEC_P2(IP,IC,1)
!  
!  FCM_OM_VEC_P3(IP,IC,4) = FCM_OM_VEC_P3(IP,IC,3)
!  FCM_OM_VEC_P3(IP,IC,3) = FCM_OM_VEC_P3(IP,IC,2)
!  FCM_OM_VEC_P3(IP,IC,2) = FCM_OM_VEC_P3(IP,IC,1)
 enddo
 
 

!~  !------------ONLY FOR 2D-----------!
 if (FCM_INIT_PART_POS==4) then
   FCM_PSWIM(IP,2)=0.0
 end if

 
 !- Renormalization
 FCM_PSWIM(IP,:) = FCM_PSWIM(IP,:)/ dsqrt( FCM_PSWIM(IP,1)**2 &
                                          +FCM_PSWIM(IP,2)**2 &
                                          +FCM_PSWIM(IP,3)**2 )
                                          
! FCM_P2(IP,:) = FCM_P2(IP,:)/ dsqrt( FCM_P2(IP,1)**2 &
!                                    +FCM_P2(IP,2)**2 &
!                                    +FCM_P2(IP,3)**2 )
!                                          
! FCM_P3(IP,:) = FCM_P3(IP,:)/ dsqrt( FCM_P3(IP,1)**2 &
!                                    +FCM_P3(IP,2)**2 &
!                                    +FCM_P3(IP,3)**2 )

 
end do !!- Loop: 1, NPART_FULL

!if (FCM_NELLIPSOID>0) then
!
! do IP = 1, NPART_FULL
! !- Transformation matrix 
!   FCM_ROT_MAT(IP,1,:) = FCM_PSWIM(IP,:)
!   
!   FCM_ROT_MAT(IP,2,:) = FCM_P2(IP,:)
!   
!   FCM_ROT_MAT(IP,3,:) = FCM_P3(IP,:)
!  end do
!  
!end if


!!- CPU check
if(MYID == 0) then
 TIME_END = MPI_WTIME()

end if




end subroutine FCM_ADV_PARTICLE_ORIENTATION_PSWIM_ONLY
