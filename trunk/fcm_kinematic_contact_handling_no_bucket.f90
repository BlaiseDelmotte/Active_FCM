 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the admissible for particles subject to enter in contact
!!> with no bucket strategy first
!! Date :  24/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_KINEMATIC_CONTACT_HANDLING_NO_BUCKET_OPTIM

!!====================================================================
!! POS_REL is an antisymmetric matrix
!!====================================================================
!! Forcing: 
!!------------------------------
!! TO DO : 
!!        1) check if necessary to save POS_REL as a vector instead
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use PARAM_PHYS



implicit none


! Temp variables
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2,3) :: EIJ_VEC
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2,2) :: LIST_CONTACT
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2)   :: DIJ_VEC
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2)   :: LAMBDA
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2)   :: DIJ_NEXT

real(kind=8),dimension(NPART_FULL,3)  :: V0_UZA
real(kind=8),dimension(NPART_FULL,3)  :: BLAMBDA



real(kind=8)     :: XI, YI, ZI
real(kind=8)     :: XIJ, YIJ, ZIJ
real(kind=8)     :: RIJ
real(kind=8)     :: TWOA



! Indices for Uzawa algorithm
integer :: ITE_UZA, ITE_UZA_MAX, NCONTACT
real(kind=8)     :: EPS_UZA
real(kind=8)     :: RHO_UZA

! Indices for loops
integer :: I, J, IND

IND = 0

EPS_UZA = 10.0**(-10.0)
RHO_UZA =2000.0  !DT = 0.001
!~ RHO_UZA = 20.0  !DT = 0.1
!~ RHO_UZA = 100.0 !DT = 0.05

ITE_UZA_MAX = 20000

LAMBDA = 0.0

NCONTACT = NPART_FULL*(NPART_FULL-1)/2


do I = 1,NPART_FULL
  V0_UZA(I,1) = FCM_UP(I,1)
  V0_UZA(I,2) = FCM_VP(I,1)
  V0_UZA(I,3) = FCM_WP(I,1)
end do
 
do I = 1,NPART_FULL
 
 XI = FCM_XP(I)
 YI = FCM_YP(I)
 ZI = FCM_ZP(I)
  
  do J = I+1, NPART_FULL
  
   IND = IND +1
   
   LIST_CONTACT(IND,1) = I 
   LIST_CONTACT(IND,2) = J
  
   XIJ = XI - FCM_XP(J)
   YIJ = YI - FCM_YP(J)
   ZIJ = ZI - FCM_ZP(J)
   
   XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
   YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
   ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
   
   RIJ = dsqrt(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ);   
      
   TWOA = FCM_SPHERE_RADP(I) + FCM_SPHERE_RADP(J)
   
   
   DIJ_VEC(IND) = RIJ - TWOA
   
   
   
   EIJ_VEC(IND,1) = -XIJ/RIJ
   EIJ_VEC(IND,2) = -YIJ/RIJ
   EIJ_VEC(IND,3) = -ZIJ/RIJ  

   DIJ_NEXT(IND) = DIJ_VEC(IND) &
                 - DTIME*( EIJ_VEC(IND,1)*(V0_UZA(I,1)-V0_UZA(J,1)) &
                         + EIJ_VEC(IND,2)*(V0_UZA(I,2)-V0_UZA(J,2)) &
                         + EIJ_VEC(IND,3)*(V0_UZA(I,3)-V0_UZA(J,3)) )
  
 end do
 
end do


!~ print*,'DIJ_VEC = ', DIJ_VEC
!~ read(*,*)
!~ print*,'maxval(-DIJ_NEXT) = ',maxval(-DIJ_NEXT)
!~ read(*,*)



ITE_UZA = 0



do while ((maxval(-DIJ_NEXT)>EPS_UZA).and.(ITE_UZA<ITE_UZA_MAX))

 ITE_UZA = ITE_UZA +1
 
!~  print*, 'ITE_UZA = ', ITE_UZA
!~  read(*,*)

 
 BLAMBDA = 0.0
 
 if (ITE_UZA==ITE_UZA_MAX) then
   print*, 'MAX UZA ITERATIONS REACHED : ITE_UZA = ', ITE_UZA_MAX
 end if
 
 
 
 do IND = 1, NCONTACT  
  
  I = LIST_CONTACT(IND,1)
  J = LIST_CONTACT(IND,2)   
    
  BLAMBDA(I,1) = BLAMBDA(I,1) + LAMBDA(IND)*EIJ_VEC(IND,1)
  BLAMBDA(I,2) = BLAMBDA(I,2) + LAMBDA(IND)*EIJ_VEC(IND,2)
  BLAMBDA(I,3) = BLAMBDA(I,3) + LAMBDA(IND)*EIJ_VEC(IND,3)
  
  BLAMBDA(J,1) = BLAMBDA(J,1) - LAMBDA(IND)*EIJ_VEC(IND,1)
  BLAMBDA(J,2) = BLAMBDA(J,2) - LAMBDA(IND)*EIJ_VEC(IND,2)
  BLAMBDA(J,3) = BLAMBDA(J,3) - LAMBDA(IND)*EIJ_VEC(IND,3)

      
 end do
 

 V0_UZA(:,1) = FCM_UP(:,1) - DTIME*BLAMBDA(:,1)
 V0_UZA(:,2) = FCM_VP(:,1) - DTIME*BLAMBDA(:,2)
 V0_UZA(:,3) = FCM_WP(:,1) - DTIME*BLAMBDA(:,3)
 

!~  read(*,*)
 
 do IND = 1, NCONTACT
 
  I = LIST_CONTACT(IND,1)
  J = LIST_CONTACT(IND,2)
     
  DIJ_NEXT(IND) = DIJ_VEC(IND) &
                - DTIME*( EIJ_VEC(IND,1)*(V0_UZA(I,1)-V0_UZA(J,1)) &
                        + EIJ_VEC(IND,2)*(V0_UZA(I,2)-V0_UZA(J,2)) &
                        + EIJ_VEC(IND,3)*(V0_UZA(I,3)-V0_UZA(J,3)) )
  
  LAMBDA(IND) =  LAMBDA(IND) - RHO_UZA*DIJ_NEXT(IND)
  
  if (LAMBDA(IND)<=0) then
   LAMBDA(IND)=0
  end if
   
 end do
 
!~  print*,'DIJ_NEXT = ', DIJ_NEXT
!~  read(*,*)
!~  
!~  print*,'LAMBDA = ', LAMBDA
!~  read(*,*)
  
end do

if (ITE_UZA>0) then
 FCM_UP(:,1) = V0_UZA(:,1)
 FCM_VP(:,1) = V0_UZA(:,2)
 FCM_WP(:,1) = V0_UZA(:,3)
end if
 
print*, 'ITE_UZA = ', ITE_UZA


   
end subroutine FCM_KINEMATIC_CONTACT_HANDLING_NO_BUCKET_OPTIM
