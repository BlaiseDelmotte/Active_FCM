 !!====================================================================
!!
!! 
!!> @brief
!!> Routine computing the repulsive shifted GB potential to avoid overlapping with 
!!> bucket sorting strategy (cf Keaveny listed list)
!! Date :  2014/02/07
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_GB_SHIFTED_UNIAXIAL

!!====================================================================
!! 
!!====================================================================
!! Forcing: 
!!------------------------------
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use FCM_BUCKET_VARIABLE



implicit none


!---- Temp variables
!- Geometric quantities
real(kind=8), dimension(3,3) :: MAT_ROT_I
real(kind=8), dimension(3,3) :: MAT_ROT_J
real(kind=8), dimension(3,3) :: MAT_RADSQ_I
real(kind=8), dimension(3,3) :: MAT_RADSQ_J
real(kind=8), dimension(3,3) :: MAT_ELL_I
real(kind=8), dimension(3,3) :: MAT_ELL_J
real(kind=8), dimension(3,3) :: MAT_H
real(kind=8), dimension(3,3) :: MAT_HINV

real(kind=8), dimension(3)   :: KAPPA

real(kind=8)     :: XI, YI, ZI
real(kind=8)     :: XIJ, YIJ, ZIJ
real(kind=8)     :: RIJSQ, RIJ
real(kind=8)     :: RREFSQ
real(kind=8)     :: TWOA
real(kind=8)     :: RHO
real(kind=8)     :: SIGMA


!- Dynamic quantities
real(kind=8), dimension(NPART_FULL,3) :: FCM_TORQUE_TEMP

real(kind=8), dimension(3)   :: FKAPPA

real(kind=8)     :: FPHI, FR
real(kind=8)     :: FXI, FYI, FZI
real(kind=8)     :: FXIJ, FYIJ, FZIJ
real(kind=8)     :: TXI, TYI, TZI
real(kind=8)     :: TXIJ, TYIJ, TZIJ
real(kind=8)     :: TXJI, TYJI, TZJI

!- Temp variables
real(kind=8), dimension(3)   :: TEMP_VEC
real(kind=8)     :: TEMP


! Indices which indicate current bucket
integer :: JCELL, JCELLO

! Indices for loops
integer :: I, J, ICELL, NABOR

!- Zeros necessary variables
FXI = 0.0
FYI = 0.0
FZI = 0.0

!- Zeros temp variable
FCM_FORCE_TEMP(:,:) = 0.0
FCM_TORQUE_TEMP(:,:) = 0.0


MAT_RADSQ_I = 0.0
MAT_RADSQ_J = 0.0

do ICELL = FCM_LOC_BUCKET_START, FCM_LOC_BUCKET_STOP 
 
 I = FCM_BUCKET_HEAD(ICELL)
 
 do while(I>-1)
 
  XI = FCM_XP(I)
  YI = FCM_YP(I)
  ZI = FCM_ZP(I)
  
  MAT_RADSQ_I(1,1) = FCM_ELLIPSOID_RADP(I,1)**2
  MAT_RADSQ_I(2,2) = FCM_ELLIPSOID_RADP(I,2)**2
  MAT_RADSQ_I(3,3) = FCM_ELLIPSOID_RADP(I,3)**2
  
!~   MAT_ROT_I =  transpose(FCM_ROT_MAT(I,1:3,1:3))
  MAT_ROT_I = FCM_ROT_MAT(I,1:3,1:3)
  
!~   print*,'I = ', I
!~   print*,'FCM_QUAT(I,:)= ', FCM_QUAT(I,:)
!~   read(*,*)
!~   print*,'MAT_ROT_I(1,:)= ', MAT_ROT_I(1,:)
!~   print*,'MAT_ROT_I(2,:)= ', MAT_ROT_I(2,:)
!~   print*,'MAT_ROT_I(3,:)= ', MAT_ROT_I(3,:)
!~   read(*,*)



  MAT_ELL_I = matmul( MAT_ROT_I, MAT_RADSQ_I)
  MAT_ELL_I = matmul( MAT_ROT_I, transpose(MAT_ELL_I))

  
!~   print*,'MAT_ELL_I(1,:)= ', MAT_ELL_I(1,:)
!~   print*,'MAT_ELL_I(2,:)= ', MAT_ELL_I(2,:)
!~   print*,'MAT_ELL_I(3,:)= ', MAT_ELL_I(3,:)
!~   read(*,*)
  
  FXI = FCM_FORCE_TEMP(I,1)
  FYI = FCM_FORCE_TEMP(I,2)
  FZI = FCM_FORCE_TEMP(I,3)
  
  TXI = FCM_TORQUE_TEMP(I,1)
  TYI = FCM_TORQUE_TEMP(I,2)
  TZI = FCM_TORQUE_TEMP(I,3)
  
  J = FCM_BUCKET_PART_LIST(I)
  
  do while(J>-1)
  
!~    XIJ = -( XI - FCM_XP(J) )
!~    YIJ = -( YI - FCM_YP(J) )
!~    ZIJ = -( ZI - FCM_ZP(J) )
   
   XIJ = XI - FCM_XP(J)
   YIJ = YI - FCM_YP(J)
   ZIJ = ZI - FCM_ZP(J)
        
   XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
   YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
   ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
   
   RIJSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ
   RIJ = dsqrt(RIJSQ)
   
   TWOA = maxval(FCM_ELLIPSOID_RADP(I,:)) + maxval(FCM_ELLIPSOID_RADP(J,:))
   
   RREFSQ = (FCM_FBRANGE * TWOA)**2
   
   if (RIJSQ.lt.RREFSQ) then    
   
    MAT_RADSQ_J(1,1) = FCM_ELLIPSOID_RADP(J,1)**2
    MAT_RADSQ_J(2,2) = FCM_ELLIPSOID_RADP(J,2)**2
    MAT_RADSQ_J(3,3) = FCM_ELLIPSOID_RADP(J,3)**2
   
!~     MAT_ROT_J = transpose(FCM_ROT_MAT(J,1:3,1:3))
    MAT_ROT_J = FCM_ROT_MAT(J,1:3,1:3)    
    
!~     print*,'J = ', J
!~     print*,'FCM_QUAT(J,:)= ', FCM_QUAT(J,:)
!~     read(*,*)
!~     print*,'MAT_ROT_J(1,:)= ', MAT_ROT_J(1,:)
!~     print*,'MAT_ROT_J(2,:)= ', MAT_ROT_J(2,:)
!~     print*,'MAT_ROT_J(3,:)= ', MAT_ROT_J(3,:)
!~     read(*,*)



    MAT_ELL_J = matmul( MAT_ROT_J, MAT_RADSQ_I)
    MAT_ELL_J = matmul( MAT_ROT_J, transpose(MAT_ELL_J))

    
!~     print*,'MAT_ELL_J(1,:)= ', MAT_ELL_J(1,:)
!~     print*,'MAT_ELL_J(2,:)= ', MAT_ELL_J(2,:)
!~     print*,'MAT_ELL_J(3,:)= ', MAT_ELL_J(3,:)
!~     read(*,*)
    
    MAT_H = MAT_ELL_I + MAT_ELL_J
        
    call FCM_INV_3_3(MAT_H,MAT_HINV)
    
    KAPPA(1) = MAT_HINV(1,1)*XIJ &
             + MAT_HINV(1,2)*YIJ & 
             + MAT_HINV(1,3)*ZIJ
             
    KAPPA(2) = MAT_HINV(2,1)*XIJ &
             + MAT_HINV(2,2)*YIJ & 
             + MAT_HINV(2,3)*ZIJ
             
    KAPPA(3) = MAT_HINV(3,1)*XIJ &
             + MAT_HINV(3,2)*YIJ & 
             + MAT_HINV(3,3)*ZIJ
             
    TEMP = KAPPA(1)*XIJ + KAPPA(2)*YIJ + KAPPA(3)*ZIJ
    
    SIGMA = RIJ*dsqrt(2.0/TEMP)
    
    RHO = ( RIJ - SIGMA + SIGMA_MIN)/SIGMA_MIN
    
    
    
!~     print*,'RIJ = ', RIJ 
!~      print*,'RIJ - SIGMA = ', RIJ - SIGMA
!~      print*,'RHO  = ', RHO
!~      read(*,*)
    
    if ((RHO**(POW_RHO)<THRES_RHO).AND.(RHO>1.0)) then
    
    
     FR = FCM_FBLEVEL*4.0*(POW_REP-1)*RHO**(-POW_REP)/SIGMA_MIN
     
     FPHI = FR * SIGMA**3/2.0 
     
     
     
     FXIJ = FR * XIJ/RIJ &
          - FPHI/RIJSQ*(KAPPA(1) - TEMP*XIJ/RIJSQ)
     FYIJ = FR * YIJ/RIJ &
          - FPHI/RIJSQ*(KAPPA(2) - TEMP*YIJ/RIJSQ)
     FZIJ = FR * ZIJ/RIJ &
          - FPHI/RIJSQ*(KAPPA(3) - TEMP*ZIJ/RIJSQ)
          
     FXI = FXI + FXIJ
     FYI = FYI + FYIJ
     FZI = FZI + FZIJ
     
     FCM_FORCE_TEMP(J,1) = FCM_FORCE_TEMP(J,1) - FXIJ
     FCM_FORCE_TEMP(J,2) = FCM_FORCE_TEMP(J,2) - FYIJ
     FCM_FORCE_TEMP(J,3) = FCM_FORCE_TEMP(J,3) - FZIJ
          
     FKAPPA = -FPHI/RIJSQ*KAPPA
     
     TEMP_VEC(1) =  MAT_ELL_I(1,1)*KAPPA(1) &
                  + MAT_ELL_I(1,2)*KAPPA(2) & 
                  + MAT_ELL_I(1,3)*KAPPA(3)
                  
     TEMP_VEC(2) =  MAT_ELL_I(2,1)*KAPPA(1) &
                  + MAT_ELL_I(2,2)*KAPPA(2) & 
                  + MAT_ELL_I(2,3)*KAPPA(3)   
                  
     TEMP_VEC(3) =  MAT_ELL_I(3,1)*KAPPA(1) &
                  + MAT_ELL_I(3,2)*KAPPA(2) & 
                  + MAT_ELL_I(3,3)*KAPPA(3) 
                  
     TXIJ = -(TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
     TYIJ = -(TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
     TZIJ = -(TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))
!~      
!~      TXIJ = (TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
!~      TYIJ = (TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
!~      TZIJ = (TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))
     
     TXI = TXI + TXIJ
     TYI = TYI + TYIJ
     TZI = TZI + TZIJ
     
     TEMP_VEC(1) =  MAT_ELL_J(1,1)*KAPPA(1) &
                  + MAT_ELL_J(1,2)*KAPPA(2) & 
                  + MAT_ELL_J(1,3)*KAPPA(3)
                  
     TEMP_VEC(2) =  MAT_ELL_J(2,1)*KAPPA(1) &
                  + MAT_ELL_J(2,2)*KAPPA(2) & 
                  + MAT_ELL_J(2,3)*KAPPA(3)   
                  
     TEMP_VEC(3) =  MAT_ELL_J(3,1)*KAPPA(1) &
                  + MAT_ELL_J(3,2)*KAPPA(2) & 
                  + MAT_ELL_J(3,3)*KAPPA(3) 
                  
     TXJI = -(TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
     TYJI = -(TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
     TZJI = -(TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))   
!~      
!~      TXJI = (TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
!~      TYJI = (TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
!~      TZJI = (TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1)) 
           
     
     FCM_TORQUE_TEMP(J,1) = FCM_TORQUE_TEMP(J,1) + TXJI
     FCM_TORQUE_TEMP(J,2) = FCM_TORQUE_TEMP(J,2) + TYJI
     FCM_TORQUE_TEMP(J,3) = FCM_TORQUE_TEMP(J,3) + TZJI
    
    end if
       
   end if
   

   
   J = FCM_BUCKET_PART_LIST(J)
   
  end do
  
  JCELLO = 13*(ICELL-1)
  
  do NABOR = 1, 13
  
   JCELL = FCM_BUCKET_MAPLIST(JCELLO + NABOR)
   J = FCM_BUCKET_HEAD(JCELL)
   
   do while(J>-1)
  
!~     XIJ = -( XI - FCM_XP(J) )
!~     YIJ = -( YI - FCM_YP(J) )
!~     ZIJ = -( ZI - FCM_ZP(J) )
!~     
    XIJ =  XI - FCM_XP(J)
    YIJ =  YI - FCM_YP(J)
    ZIJ =  ZI - FCM_ZP(J)
         
    XIJ = XIJ - LXMAX* real(int(XIJ/(0.5*LXMAX)))
    YIJ = YIJ - LYMAX* real(int(YIJ/(0.5*LYMAX)))
    ZIJ = ZIJ - LZMAX* real(int(ZIJ/(0.5*LZMAX)))
    
    RIJSQ = XIJ*XIJ + YIJ*YIJ + ZIJ*ZIJ
    RIJ = dsqrt(RIJSQ)
    
    TWOA = maxval(FCM_ELLIPSOID_RADP(I,:)) + maxval(FCM_ELLIPSOID_RADP(J,:))
    
    RREFSQ = (FCM_FBRANGE * TWOA)**2
    
    if (RIJSQ.lt.RREFSQ) then    
    
     MAT_RADSQ_J(1,1) = FCM_ELLIPSOID_RADP(J,1)**2
     MAT_RADSQ_J(2,2) = FCM_ELLIPSOID_RADP(J,2)**2
     MAT_RADSQ_J(3,3) = FCM_ELLIPSOID_RADP(J,3)**2
    
!~      MAT_ROT_J = transpose(FCM_ROT_MAT(J,1:3,1:3))
     MAT_ROT_J = FCM_ROT_MAT(J,1:3,1:3)

!~      print*,'J = ', J
!~      print*,'FCM_QUAT(J,:)= ', FCM_QUAT(J,:)
!~      read(*,*)
!~      print*,'MAT_ROT_J(1,:)= ', MAT_ROT_J(1,:)
!~      print*,'MAT_ROT_J(2,:)= ', MAT_ROT_J(2,:)
!~      print*,'MAT_ROT_J(3,:)= ', MAT_ROT_J(3,:)
!~      read(*,*)



     MAT_ELL_J = matmul( MAT_ROT_J, MAT_RADSQ_I)
     MAT_ELL_J = matmul( MAT_ROT_J, transpose(MAT_ELL_J))

     
!~      print*,'MAT_ELL_J(1,:)= ', MAT_ELL_J(1,:)
!~      print*,'MAT_ELL_J(2,:)= ', MAT_ELL_J(2,:)
!~      print*,'MAT_ELL_J(3,:)= ', MAT_ELL_J(3,:)
!~      read(*,*)
     
     MAT_H = MAT_ELL_I + MAT_ELL_J
         
     call FCM_INV_3_3(MAT_H,MAT_HINV)
     
     KAPPA(1) = MAT_HINV(1,1)*XIJ &
              + MAT_HINV(1,2)*YIJ & 
              + MAT_HINV(1,3)*ZIJ
              
     KAPPA(2) = MAT_HINV(2,1)*XIJ &
              + MAT_HINV(2,2)*YIJ & 
              + MAT_HINV(2,3)*ZIJ
              
     KAPPA(3) = MAT_HINV(3,1)*XIJ &
              + MAT_HINV(3,2)*YIJ & 
              + MAT_HINV(3,3)*ZIJ
              
     TEMP = KAPPA(1)*XIJ + KAPPA(2)*YIJ + KAPPA(3)*ZIJ
     
     SIGMA = RIJ*dsqrt(2.0/TEMP)
     
     
     RHO = ( RIJ - SIGMA + SIGMA_MIN)/SIGMA_MIN
     
!~      print*,'RIJ = ', RIJ 
!~      print*,'RIJ - SIGMA = ', RIJ - SIGMA
!~      print*,'RHO  = ', RHO
!~      read(*,*)
     
     if ((RHO**(POW_RHO)<THRES_RHO).AND.(RHO>1.0)) then
     
      FR = FCM_FBLEVEL*4.0*(POW_REP-1)*RHO**(-POW_REP)/SIGMA_MIN
      
      FPHI = FR * SIGMA**3/2.0 
      
      
      
      FXIJ = FR * XIJ/RIJ &
           - FPHI/RIJSQ*(KAPPA(1) - TEMP*XIJ/RIJSQ)
      FYIJ = FR * YIJ/RIJ &
           - FPHI/RIJSQ*(KAPPA(2) - TEMP*YIJ/RIJSQ)
      FZIJ = FR * ZIJ/RIJ &
           - FPHI/RIJSQ*(KAPPA(3) - TEMP*ZIJ/RIJSQ)
           
      FXI = FXI + FXIJ
      FYI = FYI + FYIJ
      FZI = FZI + FZIJ
      
      FCM_FORCE_TEMP(J,1) = FCM_FORCE_TEMP(J,1) - FXIJ
      FCM_FORCE_TEMP(J,2) = FCM_FORCE_TEMP(J,2) - FYIJ
      FCM_FORCE_TEMP(J,3) = FCM_FORCE_TEMP(J,3) - FZIJ
           
      FKAPPA = -FPHI/RIJSQ*KAPPA
      
      TEMP_VEC(1) =  MAT_ELL_I(1,1)*KAPPA(1) &
                   + MAT_ELL_I(1,2)*KAPPA(2) & 
                   + MAT_ELL_I(1,3)*KAPPA(3)
                   
      TEMP_VEC(2) =  MAT_ELL_I(2,1)*KAPPA(1) &
                   + MAT_ELL_I(2,2)*KAPPA(2) & 
                   + MAT_ELL_I(2,3)*KAPPA(3)   
                   
      TEMP_VEC(3) =  MAT_ELL_I(3,1)*KAPPA(1) &
                   + MAT_ELL_I(3,2)*KAPPA(2) & 
                   + MAT_ELL_I(3,3)*KAPPA(3) 
                   
      
      TXIJ = -(TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
      TYIJ = -(TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
      TZIJ = -(TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))
!~       
!~       TXIJ = (TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
!~       TYIJ = (TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
!~       TZIJ = (TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))
      
      TXI = TXI + TXIJ
      TYI = TYI + TYIJ
      TZI = TZI + TZIJ
      
      TEMP_VEC(1) =  MAT_ELL_J(1,1)*KAPPA(1) &
                   + MAT_ELL_J(1,2)*KAPPA(2) & 
                   + MAT_ELL_J(1,3)*KAPPA(3)
                   
      TEMP_VEC(2) =  MAT_ELL_J(2,1)*KAPPA(1) &
                   + MAT_ELL_J(2,2)*KAPPA(2) & 
                   + MAT_ELL_J(2,3)*KAPPA(3)   
                   
      TEMP_VEC(3) =  MAT_ELL_J(3,1)*KAPPA(1) &
                   + MAT_ELL_J(3,2)*KAPPA(2) & 
                   + MAT_ELL_J(3,3)*KAPPA(3) 
                   
      TXJI = -(TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
      TYJI = -(TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
      TZJI = -(TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))   
!~       
!~       TXJI = (TEMP_VEC(2)*FKAPPA(3) - TEMP_VEC(3)*FKAPPA(2))
!~       TYJI = (TEMP_VEC(3)*FKAPPA(1) - TEMP_VEC(1)*FKAPPA(3))
!~       TZJI = (TEMP_VEC(1)*FKAPPA(2) - TEMP_VEC(2)*FKAPPA(1))      
 
      
      FCM_TORQUE_TEMP(J,1) = FCM_TORQUE_TEMP(J,1) + TXJI
      FCM_TORQUE_TEMP(J,2) = FCM_TORQUE_TEMP(J,2) + TYJI
      FCM_TORQUE_TEMP(J,3) = FCM_TORQUE_TEMP(J,3) + TZJI
     
     end if
        
    end if
    
!~     print*,'Neighbouring Bucket'
!~     print*, 'I = ', I
!~     print*, 'J = ', J
!~     print*, 'FCM_POS_REL_X(I,J) = ', FCM_POS_REL_X(I,J)
!~     print*, 'FCM_POS_REL_Y(I,J) = ', FCM_POS_REL_Y(I,J)
!~     print*, 'FCM_POS_REL_Z(I,J)  = ',FCM_POS_REL_Z(I,J) 
!~     read(*,*)
       
   
    J = FCM_BUCKET_PART_LIST(J)
   
   end do 
  
  end do
  
  FCM_FORCE_TEMP(I,1) = FXI
  FCM_FORCE_TEMP(I,2) = FYI
  FCM_FORCE_TEMP(I,3) = FZI 
  
  FCM_TORQUE_TEMP(I,1) = TXI
  FCM_TORQUE_TEMP(I,2) = TYI
  FCM_TORQUE_TEMP(I,3) = TZI 
  
  I = FCM_BUCKET_PART_LIST(I)
  
 end do
 
end do



!~ print*,'FCM_PSWIM = ',FCM_PSWIM 
!~ print*,'XIJ, YIJ, ZIJ = ', XIJ, YIJ, ZIJ
!~ read(*,*) 



!~ print*,'MAT_HINV = ', MAT_HINV
!~ read(*,*) 
!~ 
!~ print*,'KAPPA = ', KAPPA
!~ print*,'SIGMA = ', SIGMA
!~ 
!~ read(*,*)


! Simple addition of the contribution of each processor to the force monopole and torque
call MPI_ALLREDUCE(FCM_FORCE_TEMP,FCM_FORCE,NPART_FULL*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)
call MPI_ALLREDUCE(FCM_TORQUE_TEMP,FCM_TORQUE,NPART_FULL*3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERR)

!~ print*,'FCM_FORCE = ', FCM_FORCE 
!~ print*,'FCM_TORQUE = ', FCM_TORQUE  
!~ 
!~ print*,'FCM_TORQUE(1,1) + FCM_TORQUE(2,1) - (YIJ*FZI-ZIJ*FYI) = ',FCM_TORQUE(1,1) + FCM_TORQUE(2,1) - (YIJ*FZI-ZIJ*FYI)
!~ print*,'FCM_TORQUE(1,2) + FCM_TORQUE(2,2) - (ZIJ*FXI-XIJ*FZI) = ',FCM_TORQUE(1,2) + FCM_TORQUE(2,2) - (ZIJ*FXI-XIJ*FZI) 
!~ print*,'FCM_TORQUE(1,3) + FCM_TORQUE(2,3) - (XIJ*FYI-YIJ*FXI) = ',FCM_TORQUE(1,3) + FCM_TORQUE(2,3) - (XIJ*FYI-YIJ*FXI) 
!~ read(*,*)


end subroutine FCM_GB_SHIFTED_UNIAXIAL
