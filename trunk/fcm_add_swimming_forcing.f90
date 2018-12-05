  !!====================================================================
!!
!! 
!!> @author 
!!> Blaise Delmotte
!!
!! DESCRIPTION: 
!!> @brief
!!> Compute the swimming stresslet from the particle intrinsic velocity
!!> and FCM_PSWIM + potential dipole for squirming
!!
!!
!!
!! Date :  21/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_ADD_SWIMMING_FORCING

!!====================================================================
!! Here the swimming stresslet forcing is set according to the formula 
!! of Pedley Kessler JFM 1990
!!====================================================================
!! Forcing: 
!!------------------------------
!! TO DO : 
!!        1) If ellipse ...
!! !!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use PARAM_PHYS,  only: VISC

implicit none

! Squirming variables deduced from BETA and VSW
real(kind=8) :: B1, B2 

! Indices for loops
integer :: IP, IND



if ((FCM_NSWIM(1)>0).and.(FCM_NSWIM(2).eq.0)) then

 do IP = 1, FCM_NSWIM(1)
 

  FCM_S0 =  6.0 * PPI * FCM_SPHERE_RADP(IP)**2 * FCM_BETA * FCM_VSW(IP)   ! S0 : magnitude of the dipole

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = FCM_PSWIM(IP,1)**2 - 1.0/3.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = FCM_PSWIM(IP,2)**2 - 1.0/3.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)
  
!~   FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0    ! Add unit tensor on the diagonal terms
!~   FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
!~   FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
!~   FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0   ! Add unit tensor on the diagonal terms
!~   FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)
  
  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
   
 end do

 IND = 0

 do IP = FCM_NSPHERE + 1, NPART_FULL
 

  IND = IND + 1

 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!
  FCM_S0 = PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**2 * FCM_BETA * FCM_VSW(IP)   ! S0 : magnitude of the dipole
 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = FCM_PSWIM(IP,1)**2 - 1.0/3.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = FCM_PSWIM(IP,2)**2 - 1.0/3.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)
  
!~   FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0    ! Add unit tensor on the diagonal terms
!~   FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
!~   FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
!~   FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0   ! Add unit tensor on the diagonal terms
!~   FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)

  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
   
 end do
 
else if ((FCM_NSWIM(2)>0).and.(FCM_NSWIM(3).eq.0)) then


  
 do IP = 1, FCM_NSWIM(2)
  
  B1 = FCM_VSW(IP) * 3.0/2.0 
  B2 = FCM_BETA *B1 

  ! Potential dipole = -4/3*pi*a^3*B1*psw = -2*pi*a^3*VSW
  FCM_H0 =  -4.0/3.0 * PPI * FCM_SPHERE_RADP(IP)**3 * B1
  FCM_S0 =   4.0/3.0 * PPI * FCM_SPHERE_RADP(IP)**2 * B2   ! S0 : magnitude of the dipole

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0   ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)
  
  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
  
  FCM_HI(IP,:) = FCM_PSWIM(IP,:) * FCM_H0
   
 end do

 IND = 0
 
 B1 = FCM_VSW_REF * 3.0/2.0 
 B2 = FCM_BETA * B1

 do IP = FCM_NSPHERE + 1, NPART_FULL
 
  IND = IND + 1
  
  !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!
  FCM_H0 = -4.0/3.0 * PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**3 * B1  ! S0 : magnitude of the dipole
 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!

 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!
  FCM_S0 = 4.0/3.0 * PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**2 * B2   ! S0 : magnitude of the dipole
 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0   ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)

  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
  FCM_HI(IP,:) = FCM_PSWIM(IP,:) * FCM_H0
   
 end do



else if (FCM_NSWIM(3)>0) then


  
 do IP = 1, FCM_NSWIM(3)
  
   B1 = FCM_VSW_TIME(IP) * 3.0/2.0 
   B2 = FCM_B2_TIME(IP)

  ! Potential dipole = -4/3*pi*a^3*B1*psw = -2*pi*a^3*VSW
  FCM_H0 =  -4.0/3.0 * PPI * FCM_SPHERE_RADP(IP)**3 * B1
  FCM_S0 =   4.0/3.0 * PPI * FCM_SPHERE_RADP(IP)**2 * B2   ! S0 : magnitude of the dipole

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0   ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)
  
  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
  
  FCM_HI(IP,:) = FCM_PSWIM(IP,:) * FCM_H0
   
 end do

 IND = 0
 


 do IP = FCM_NSPHERE + 1, NPART_FULL
 
  IND = IND + 1
  
  B1 = FCM_VSW_TIME(IP) * 3.0/2.0 
  B2 = FCM_B2_TIME(IP)
  
  !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!
  FCM_H0 = -4.0/3.0 * PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**3 * B1  ! S0 : magnitude of the dipole
 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!

 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!
  FCM_S0 = 4.0/3.0 * PPI * maxval(FCM_ELLIPSOID_RADP(IND,:))**2 * B2   ! S0 : magnitude of the dipole
 !!======= IT HAS TO BE MODIFIED AS THE DRAG ON A SPHEROID IS NOT 6*pi*mu*amax======!!

  ! Compute the swimming stresslet
  FCM_SPIJ(IP,1) = 3.0*FCM_PSWIM(IP,1)**2 - 1.0   ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,2) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,2)
  FCM_SPIJ(IP,3) = 3.0*FCM_PSWIM(IP,1)*FCM_PSWIM(IP,3)
  FCM_SPIJ(IP,4) = 3.0*FCM_PSWIM(IP,2)**2 - 1.0    ! Add unit tensor on the diagonal terms
  FCM_SPIJ(IP,5) = 3.0*FCM_PSWIM(IP,2)*FCM_PSWIM(IP,3)

  FCM_SPIJ(IP,:) = FCM_SPIJ(IP,:) * FCM_S0
  FCM_HI(IP,:) = FCM_PSWIM(IP,:) * FCM_H0
   
 end do


end if

!~ if (MYID==0) then
!~ 
! print*,'B1, FCM_VSW, B2, FCM_H0, FCM_S0 = ',B1, FCM_VSW, B2, FCM_H0, FCM_S0
! read(*,*)
!~ 
!~ end if


end subroutine FCM_ADD_SWIMMING_FORCING
