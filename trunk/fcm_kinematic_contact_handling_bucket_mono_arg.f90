 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the admissible for particles subject to enter in contact
!!> with bucket strategy 
!! Date :  27/01/2014
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_KINEMATIC_CONTACT_HANDLING_BUCKET_MONO_ARG

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
use FCM_BUCKET_VARIABLE
use PARAM_PHYS



implicit none

! Known size
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2,3) :: EIJ_VEC
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2)   :: DIJ_VEC
real(kind=8),dimension(NPART_FULL*(NPART_FULL-1)/2)   :: DIJ_NEXT
integer,dimension(NPART_FULL*(NPART_FULL-1)/2,2) :: LIST_CONTACT

integer, dimension(NPART_FULL)  :: NEW
integer, dimension(NPART_FULL)  :: INVOLVED
integer, dimension(NPART_FULL)  :: INVOLVED_REVERSE


!Contact Parameters
real(kind=8)     :: XI, YI, ZI
real(kind=8)     :: XIJ, YIJ, ZIJ
real(kind=8)     :: RIJ
real(kind=8)     :: RAD_I
real(kind=8)     :: RAD_J
real(kind=8)     :: TWOA
real(kind=8)     :: RREF
integer          :: NCONTACT, NINVOLVED



! Indices for loops
integer :: I, J, IND, ICELL, NABOR
integer :: IND_INV1, IND_INV2

!~ ! Index which indicates contacts belonging to him
!~ integer :: IND_MINE

! Indices which indicate current bucket
integer :: JCELL, JCELLO


IND = 0
IND_INV1 = 0
IND_INV2 = 0


NEW = 0
INVOLVED = -1
INVOLVED_REVERSE = -1


 
!~ do ICELL = FCM_LOC_BUCKET_START, FCM_LOC_BUCKET_STOP 
do ICELL = 1, FCM_BUCKET_NB_TOT
 
 I = FCM_BUCKET_HEAD(ICELL)
 
 do while(I>-1)
 
  XI = FCM_XP(I)
  YI = FCM_YP(I)
  ZI = FCM_ZP(I)
  
  if (I>FCM_NSPHERE) then
   RAD_I = FCM_ELLIPSOID_RADP(I,1)
  else
   RAD_I = FCM_SPHERE_RADP(I)
  end if
   
  J = FCM_BUCKET_PART_LIST(I)
  
    
  do while(J>-1)
   
   XIJ = XI - FCM_XP(J)
   YIJ = YI - FCM_YP(J)
   ZIJ = ZI - FCM_ZP(J)
   
   if (J>FCM_NSPHERE) then
    RAD_J = FCM_ELLIPSOID_RADP(J,1)
   else
    RAD_J = FCM_SPHERE_RADP(J)
   end if
   
   XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
   YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
   ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
   
   RIJ = dsqrt(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ);   
     
      
   TWOA = RAD_I + RAD_J
   
   RREF = FCM_FBRANGE  * TWOA 

   ! If contact possible then construct Uzawa system before solving
   if (RIJ-TWOA<RREF) then
   
    IND = IND +1
    
    LIST_CONTACT(IND,1) = I 
    LIST_CONTACT(IND,2) = J
    
    DIJ_VEC(IND) = RIJ - TWOA 
    
    EIJ_VEC(IND,1) = -XIJ/RIJ
    EIJ_VEC(IND,2) = -YIJ/RIJ
    EIJ_VEC(IND,3) = -ZIJ/RIJ  

    DIJ_NEXT(IND) = DIJ_VEC(IND) &
                  - DTIME*( EIJ_VEC(IND,1)*(FCM_UP(I,1)-FCM_UP(J,1)) &
                          + EIJ_VEC(IND,2)*(FCM_VP(I,1)-FCM_VP(J,1)) &
                          + EIJ_VEC(IND,3)*(FCM_WP(I,1)-FCM_WP(J,1)) )
    
     if (NEW(I)==0) then
    
      IND_INV1 = IND_INV1 + 1
      INVOLVED(IND_INV1) = I
      INVOLVED_REVERSE(I) = IND_INV1
      
      NEW(I) = NEW(I) + 1
      
     end if
    
     if (NEW(J)==0) then
    
      IND_INV1 = IND_INV1 + 1
      INVOLVED(IND_INV1) = J
      INVOLVED_REVERSE(J) = IND_INV1
      
      NEW(J) = NEW(J) + 1
      
     end if
                          
   end if
   
   J = FCM_BUCKET_PART_LIST(J)
   
  end do
  
  JCELLO = 13*(ICELL-1)
  
  do NABOR = 1, 13
  
   JCELL = FCM_BUCKET_MAPLIST(JCELLO + NABOR)
   
   if ((JCELL>=FCM_LOC_BUCKET_START).and.(JCELL<=FCM_LOC_BUCKET_STOP)) then
   
    J = FCM_BUCKET_HEAD(JCELL)
    
    do while(J>-1)     
      
     XIJ = XI - FCM_XP(J)
     YIJ = YI - FCM_YP(J)
     ZIJ = ZI - FCM_ZP(J)
     
     if (J>FCM_NSPHERE) then
      RAD_J = FCM_ELLIPSOID_RADP(J,1)
     else
      RAD_J = FCM_SPHERE_RADP(J)
     end if
     
     XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
     YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
     ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
     
     RIJ = dsqrt(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ);   
        
     TWOA = RAD_I + RAD_J
    
     RREF = FCM_FBRANGE  * TWOA 

     ! If contact possible then construct Uzawa system before solving
     if (RIJ-TWOA<RREF) then
     
      IND = IND +1
      
      LIST_CONTACT(IND,1) = I 
      LIST_CONTACT(IND,2) = J
      
      DIJ_VEC(IND) = RIJ - TWOA 
    
      EIJ_VEC(IND,1) = -XIJ/RIJ
      EIJ_VEC(IND,2) = -YIJ/RIJ
      EIJ_VEC(IND,3) = -ZIJ/RIJ  

      DIJ_NEXT(IND) = DIJ_VEC(IND) &
                    - DTIME*( EIJ_VEC(IND,1)*(FCM_UP(I,1)-FCM_UP(J,1)) &
                            + EIJ_VEC(IND,2)*(FCM_VP(I,1)-FCM_VP(J,1)) &
                            + EIJ_VEC(IND,3)*(FCM_WP(I,1)-FCM_WP(J,1)) )
              
      if (NEW(I)==0) then
     
       IND_INV1 = IND_INV1 + 1
       INVOLVED(IND_INV1) = I
       INVOLVED_REVERSE(I) = IND_INV1
       
       NEW(I) = NEW(I) + 1
       
      end if
     
      if (NEW(J)==0) then
     
       IND_INV1 = IND_INV1 + 1
       INVOLVED(IND_INV1) = J
       INVOLVED_REVERSE(J) = IND_INV1
       
       NEW(J) = NEW(J) + 1
       
      end if
                            
     end if
     
     J = FCM_BUCKET_PART_LIST(J)
    
    end do 
   
   else
   
    J = FCM_BUCKET_HEAD(JCELL)
    
    do while(J>-1)
    
            
     XIJ = XI - FCM_XP(J)
     YIJ = YI - FCM_YP(J)
     ZIJ = ZI - FCM_ZP(J)
     
     if (J>FCM_NSPHERE) then
      RAD_J = FCM_ELLIPSOID_RADP(J,1)
     else
      RAD_J = FCM_SPHERE_RADP(J)
     end if
     
     XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
     YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
     ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
     
     RIJ = dsqrt(XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ);   
        
     TWOA = RAD_I + RAD_J
    
     RREF = FCM_FBRANGE  * TWOA 

     ! If contact possible then construct Uzawa system before solving
     if (RIJ-TWOA<RREF) then
     
      IND = IND +1
      
      LIST_CONTACT(IND,1) = I 
      LIST_CONTACT(IND,2) = J
      
      DIJ_VEC(IND) = RIJ - TWOA 
    
      EIJ_VEC(IND,1) = -XIJ/RIJ
      EIJ_VEC(IND,2) = -YIJ/RIJ
      EIJ_VEC(IND,3) = -ZIJ/RIJ  

      DIJ_NEXT(IND) = DIJ_VEC(IND) &
                    - DTIME*( EIJ_VEC(IND,1)*(FCM_UP(I,1)-FCM_UP(J,1)) &
                            + EIJ_VEC(IND,2)*(FCM_VP(I,1)-FCM_VP(J,1)) &
                            + EIJ_VEC(IND,3)*(FCM_WP(I,1)-FCM_WP(J,1)) )
        
      if (NEW(I)==0) then
     
       IND_INV1 = IND_INV1 + 1
       INVOLVED(IND_INV1) = I
       INVOLVED_REVERSE(I) = IND_INV1
       
       NEW(I) = NEW(I) + 1
       
      end if
     
      if (NEW(J)==0) then
     
       IND_INV1 = IND_INV1 + 1
       INVOLVED(IND_INV1) = J
       INVOLVED_REVERSE(J) = IND_INV1
       
       NEW(J) = NEW(J) + 1
       
      end if
                            
     end if
     
     J = FCM_BUCKET_PART_LIST(J)
    
    end do 
    
   end if
   
  end do

  I = FCM_BUCKET_PART_LIST(I)
  
 end do
 
end do
  



NCONTACT = IND
NINVOLVED = IND_INV1

!~ print*,'LIST_CONTACT(1:NCONTACT,1:2) = ',LIST_CONTACT(1:NCONTACT,1:2)

call FCM_UZAWA_MONO(NCONTACT,                     &
                    NINVOLVED,                    &
                    LIST_CONTACT(1:NCONTACT,1:2), &
                    INVOLVED,                     &
                    INVOLVED_REVERSE,             &
                    EIJ_VEC(1:NCONTACT,1:3),      &
                    DIJ_VEC(1:NCONTACT),          &
                    DIJ_NEXT(1:NCONTACT)          ) 

end subroutine FCM_KINEMATIC_CONTACT_HANDLING_BUCKET_MONO_ARG
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
!--------------------------------------------------------------------
subroutine FCM_UZAWA_MONO(NCONTACT,                   &
                          NINVOLVED,                  &
                          LIST_CONTACT,               &
                          INVOLVED,                   &
                          INVOLVED_REVERSE,           &
                          EIJ_VEC,                    &
                          DIJ_VEC,                    &
                          DIJ_NEXT                    )


use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use PARAM_PHYS



implicit none
!--------------------------------------------------------------------
! ARRAY STATEMENTS
!--------------------------------------------------------------------
!- Arguments
integer, intent(in)                                :: NCONTACT
integer, intent(in)                                :: NINVOLVED
integer, dimension(NCONTACT,2), intent(in)         :: LIST_CONTACT
integer, dimension(NPART_FULL), intent(in)         :: INVOLVED
integer, dimension(NPART_FULL), intent(in)         :: INVOLVED_REVERSE
real(kind=8), dimension(NCONTACT,3), intent(in)    :: EIJ_VEC
real(kind=8), dimension(NCONTACT), intent(in)      :: DIJ_VEC
real(kind=8), dimension(NCONTACT), intent(inout)   :: DIJ_NEXT


!- Local parameters
real(kind=8), dimension(NCONTACT)     :: LAMBDA
real(kind=8), dimension(NINVOLVED,3)  :: V0_UZA
real(kind=8), dimension(NINVOLVED,3)  :: BLAMBDA


! Parameters for Uzawa algorithm
integer          :: ITE_UZA, ITE_UZA_MAX
real(kind=8)     :: EPS_UZA
real(kind=8)     :: RHO_UZA

! Indices for loops
integer :: I, J, IND, ICELL, NABOR
integer :: IND_INV1, IND_INV2

!~ 
!~  print*,MYID,'NEW = ', NEW
!~  print*,MYID,'NCONTACT = ', NCONTACT
!~  print*,MYID,'NINVOLVED = ', NINVOLVED
!~  print*,MYID,'NINVOLVED_MINE = ', NINVOLVED_MINE
!~  print*,MYID,'INVOLVED = ', INVOLVED
!~  print*,MYID,'INVOLVED_MINE = ', INVOLVED_MINE
!~  print*,MYID,'LIST_CONTACT(1:NCONTACT,1) = ',LIST_CONTACT(1:NCONTACT,1)
!~  print*,MYID,'LIST_CONTACT(1:NCONTACT,2) = ',LIST_CONTACT(1:NCONTACT,2)
!~ 
!~ print*,'NCONTACT, NINVOLVED = ',NCONTACT, NINVOLVED
!~ 
!~ print*,'LIST_CONTACT(1:NCONTACT,1:2) = ',LIST_CONTACT(1:NCONTACT,1:2)


EPS_UZA = 10.0**(-10.0)

if (DTIME<=0.005) then
 RHO_UZA = 15000.0  !DT = 0.0.005
else if (DTIME<=0.01) then
 RHO_UZA = 2000.0 !DT = 0.01
else if (DTIME<=0.05) then
 RHO_UZA = 100.0 !DT = 0.05
else if (DTIME<=0.1) then
 RHO_UZA = 20.0  !DT = 0.1
end if



!~ ITE_UZA_MAX = 2
ITE_UZA_MAX = 20000

LAMBDA = 0.0
BLAMBDA = 0.0
V0_UZA = 0.0


ITE_UZA = 0



do while ((maxval(-DIJ_NEXT)>EPS_UZA).and.(ITE_UZA<ITE_UZA_MAX))

 ITE_UZA = ITE_UZA +1
 
!~  print*,MYID, 'ITE_UZA = ', ITE_UZA
!~  read(*,*)
 
 BLAMBDA = 0.0
 
 if (ITE_UZA==ITE_UZA_MAX) then
   print*,MYID, 'MAX UZA ITERATIONS REACHED : ITE_UZA = ', ITE_UZA_MAX
 end if
 
 
 
 do IND = 1, NCONTACT
  
!~   print*,'IND = ', IND
  I = LIST_CONTACT(IND,1)
  J = LIST_CONTACT(IND,2)
!~   print*,'I, J = ', I, J 
  
  IND_INV1 = INVOLVED_REVERSE(I)
  IND_INV2 = INVOLVED_REVERSE(J)
  
!~   print*,'Hello435'

  BLAMBDA(IND_INV1,1) = BLAMBDA(IND_INV1,1) + LAMBDA(IND)*EIJ_VEC(IND,1)
  BLAMBDA(IND_INV1,2) = BLAMBDA(IND_INV1,2) + LAMBDA(IND)*EIJ_VEC(IND,2)
  BLAMBDA(IND_INV1,3) = BLAMBDA(IND_INV1,3) + LAMBDA(IND)*EIJ_VEC(IND,3)
  
!~   print*,'Hello441'
  
  BLAMBDA(IND_INV2,1) = BLAMBDA(IND_INV2,1) - LAMBDA(IND)*EIJ_VEC(IND,1)
  BLAMBDA(IND_INV2,2) = BLAMBDA(IND_INV2,2) - LAMBDA(IND)*EIJ_VEC(IND,2)
  BLAMBDA(IND_INV2,3) = BLAMBDA(IND_INV2,3) - LAMBDA(IND)*EIJ_VEC(IND,3)
  
!~   print*,'Hello447'
  
!~    if (ITE_UZA==ITE_UZA_MAX) then
!~     if ((I==1).or.(J==1)) then
!~      print*,MYID,'IND CONTACT = ', IND
!~      print*,MYID,'I,J = ', I, J
!~      print*,MYID,'BLAMBDA(IND_INV1,:) = ', BLAMBDA(IND_INV1,:)
!~      print*,MYID,'BLAMBDA(IND_INV2,:) = ', BLAMBDA(IND_INV2,:)
!~     !read(*,*)
!~      end if
!~    end if
   
      
 end do
 
!~   if (ITE_UZA==ITE_UZA_MAX) then
   !if (MYID==1) then
!~      do IND_INV1 = 1,NINVOLVED
!~         
!~         I = INVOLVED(IND_INV1)
!~         print*,MYID,'I = ', I
!~         print*,MYID,'BLAMBDA(IND_INV1,:) = ', (BLAMBDA(IND_INV1,J), J=1,3)
!~         
!~       end do      
     !end if
!~    end if
   



 do IND_INV1 = 1, NINVOLVED

  I = INVOLVED(IND_INV1)
  V0_UZA(IND_INV1,1) = FCM_UP(I,1) - DTIME*BLAMBDA(IND_INV1,1)
  V0_UZA(IND_INV1,2) = FCM_VP(I,1) - DTIME*BLAMBDA(IND_INV1,2)
  V0_UZA(IND_INV1,3) = FCM_WP(I,1) - DTIME*BLAMBDA(IND_INV1,3)  
  
end do

 

 
 do IND = 1, NCONTACT
 
  I = LIST_CONTACT(IND,1)
  J = LIST_CONTACT(IND,2)
  
  
  IND_INV1 = INVOLVED_REVERSE(I)
  IND_INV2 = INVOLVED_REVERSE(J)
  
!~   if ((I==11).or.(J==11)) then


!~  if ((ITE_UZA>1)) then
!~    print*,MYID,'I,J = ', I,J
!~    print*,MYID,'IND = ', IND
!~    print*,MYID,'BLAMBDA(IND_INV1,:) = ', BLAMBDA(IND_INV1,:)
   
!~   end if
  

     
  DIJ_NEXT(IND) = DIJ_VEC(IND) &
                - DTIME*( EIJ_VEC(IND,1)*(V0_UZA(IND_INV1,1)-V0_UZA(IND_INV2,1)) &
                        + EIJ_VEC(IND,2)*(V0_UZA(IND_INV1,2)-V0_UZA(IND_INV2,2)) &
                        + EIJ_VEC(IND,3)*(V0_UZA(IND_INV1,3)-V0_UZA(IND_INV2,3)) )
  
  LAMBDA(IND) =  LAMBDA(IND) - RHO_UZA*DIJ_NEXT(IND)
  
  if (LAMBDA(IND)<=0) then
   LAMBDA(IND)=0   
  end if
!~   if (ITE_UZA==ITE_UZA_MAX) then
!~    !print*,MYID,'LAMBDA(IND) = ', LAMBDA(IND)
!~    print*,MYID,'DIJ_NEXT(IND) = ', DIJ_NEXT(IND)
!~   end if
   
 end do
!~  if (ITE_UZA==ITE_UZA_MAX) then
!~   !print*,MYID,'ITE_UZA = ',ITE_UZA,'LAMBDA = ', LAMBDA
!~   print*,MYID,'ITE_UZA = ',ITE_UZA,'DIJ_NEXT = ', DIJ_NEXT
!~  end if
  
end do


!~ NEW = 0
!~ do IND = 1, NCONTACT
!~  
!~   I = LIST_CONTACT(IND,1)
!~   J = LIST_CONTACT(IND,2)
!~   
!~   if (LAMBDA(IND)>0) then
!~    if (NEW(I)==0) then
!~     NEW(I) = NEW(I) + 1
!~    end if   
!~    if (NEW(J)==0) then
!~     NEW(J) = NEW(J) + 1
!~    end if
!~   end if
!~   
!~ end do






if (ITE_UZA>0) then
 
 do IND_INV1 = 1,NINVOLVED
   
  I = INVOLVED(IND_INV1)
 
  FCM_UP(I,1) = V0_UZA(IND_INV1,1)   
  FCM_VP(I,1) = V0_UZA(IND_INV1,2) 
  FCM_WP(I,1) = V0_UZA(IND_INV1,3)
    
 end do 
 
 
end if



!~ print*,MYID,'ITE_UZA = ', ITE_UZA



!~ if (ITE_UZA>0) then
!~  do IND_MINE_INV1 = 1,NINVOLVED_MINE
!~ 
!~   I = INVOLVED_MINE(IND_MINE_INV1)
!~   IND_INV1 = INVOLVED_REVERSE(I)
!~   FCM_UP_TEMP(I) = V0_UZA(IND_INV1,1)
!~   FCM_VP_TEMP(I) = V0_UZA(IND_INV1,2)
!~   FCM_WP_TEMP(I) = V0_UZA(IND_INV1,3)
!~   
!~  end do
!~  
!~ else
!~ 
!~  do IND_MINE_INV1 = 1,NINVOLVED_MINE
!~ 
!~   I = INVOLVED_MINE(IND_MINE_INV1)
!~   FCM_UP_TEMP(I) = FCM_UP(I,1)
!~   FCM_VP_TEMP(I) = FCM_VP(I,1)
!~   FCM_WP_TEMP(I) = FCM_WP(I,1)
!~   
!~  end do
!~  
!~ end if





!~ print*,MYID, ITE_UZA
!~ print*,MYID,'BLAMBDA(INVOLVED_REVERSE(6),:) = ', BLAMBDA(INVOLVED_REVERSE(6),:)
!~ print*,MYID,'BLAMBDA(INVOLVED_REVERSE(11),:) = ', BLAMBDA(INVOLVED_REVERSE(11),:)


!~ FCM_VP_TEMP = FCM_VP_TEMP - real(NEW_TOT-1)/real(NPROC)*FCM_VP(I,1)
!~ FCM_WP_TEMP = FCM_WP_TEMP - real(NEW_TOT-1)/real(NPROC)*FCM_WP(I,1)






   
end subroutine FCM_UZAWA_MONO
