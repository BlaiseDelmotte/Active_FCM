!!====================================================================
!!
!!
!!====================================================================

subroutine POLAR_ORDER(NSAVES, &
                       PART_START, &
                       PART_END, &
                       NBOX_DIR, &
                       LX, &
                       LY, &
                       LZ, &
                       TIME_VEC, &
                       PSWIM, &
                       P2, &
                       P3, &       
                       POSI, &                
                       PSWIM_SPH, &
                       MEAN_PSWIM , &
                       PSWIM_SPH_PMEAN  )
                       

!!====================================================================
!!
!!
!!====================================================================

use MPI

implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer, intent(in) :: NSAVES
! Number of swimmers
integer, intent(in) :: PART_START, PART_END
! Number of boxes in one direction to compute local PO
integer, intent(in) :: NBOX_DIR
! Physical Time
real(kind=8), intent(in) :: LX, LY, LZ
! Physical Time
real(kind=8), dimension(NSAVES), intent(in) :: TIME_VEC
! Swimmers orientation
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: PSWIM
! Swimmers orientation
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: P2
! Swimmers orientation
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: P3
! Swimmers position
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,3), intent(in) :: POSI
! Swimmers orientation along time in spherical coordinates (theta,phi)
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,2), intent(out) :: PSWIM_SPH
! Mean Swimmers orientation along time
real(kind=8), dimension(NSAVES,3), intent(out) :: MEAN_PSWIM
! Fluctuations of Swimmers orientation around Pmean along time in spherical coordinates (theta,phi)
real(kind=8), dimension(NSAVES,PART_END-PART_START+1,2), intent(out) :: PSWIM_SPH_PMEAN


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
real(kind=8), dimension(NSAVES,2) :: MEAN_PSWIM_SPH
real(kind=8), dimension(NSAVES,NBOX_DIR*NBOX_DIR*NBOX_DIR,3) :: LOCAL_MEAN_PSWIM
real(kind=8), dimension(NSAVES,NBOX_DIR*NBOX_DIR*NBOX_DIR) :: COUNT_BOX
real(kind=8), dimension(NSAVES,3) :: MEAN_PSWIM_NON_NORM
real(kind=8), dimension(NSAVES,3) :: MEAN_P2
real(kind=8), dimension(NSAVES,3) :: MEAN_P3
real(kind=8), dimension(NSAVES,3) :: MEAN_P2_ORTH
real(kind=8), dimension(NSAVES,3) :: MEAN_P3_ORTH

!- Polar order
real(kind=8), dimension(NSAVES) :: POLAR_ORDER_TIME

!- Polar order
real(kind=8), dimension(NSAVES,NBOX_DIR*NBOX_DIR*NBOX_DIR) :: LOCAL_POLAR_ORDER_TIME


!- Module of Mean Pswim for one time
real(kind=8) :: MOD_MEAN_PSWIM, &
                LOCAL_MOD_MEAN_PSWIM, &
                MOD_MEAN_P2, &
                MOD_MEAN_P3, &
                MOD_VECTOR, &
                MOD_VECTOR2 

!- 
real(kind=8) :: SCAL_MEAN_P2_MEAN_PSWIM

! - Box size for local PO
real(kind=8) :: DBOX_X, DBOX_Y, DBOX_Z

!- Pi
real(kind=8) :: PPI

!- File name 
character(len=40) :: FILENAME
character(len=20) :: MYFMT

!- Total number of boxes
integer :: NBOX_TOT

!- Box Index
integer :: IND_BOX

!- Index
integer :: I, J, K, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_POLAR



PPI = 4.0*datan(1.0)

!---------------------------------------------------------------------
!=====================================================================
! 1. COMPUTE POLAR ORDER ALONG TIME
!=====================================================================
print*,' '
print*,'-------------------START POLAR ORDER--------------------------- '
print*,' '

NBOX_TOT = NBOX_DIR*NBOX_DIR*NBOX_DIR

DBOX_X = LX/real(NBOX_DIR)
DBOX_Y = LY/real(NBOX_DIR)
DBOX_Z = LZ/real(NBOX_DIR)

print*,'DBOX_X, DBOX_Y, DBOX_Z = ', DBOX_X, DBOX_Y, DBOX_Z

PSWIM_SPH = 0.0
!~ PSWIM_SPH_PMEAN = 0.0
LOCAL_POLAR_ORDER_TIME = 0.0

do IND = 1, NSAVES

 LOCAL_MEAN_PSWIM(IND,1:NBOX_TOT,1:3) = 0.0
 COUNT_BOX(IND,1:NBOX_TOT) = 0.0
 MEAN_PSWIM(IND,1:3) = 0.0
 MEAN_PSWIM_NON_NORM(IND,1:3) = 0.0
 MEAN_P2(IND,1:3) = 0.0
 MEAN_P3(IND,1:3) = 0.0
 
 MEAN_PSWIM_SPH(IND,1:2) = 0.0
 
 
 
 do I = 1, PART_END-PART_START+1 
 
  MEAN_PSWIM(IND,1:3) = MEAN_PSWIM(IND,1:3) + PSWIM(IND,I,1:3)
  MEAN_P2(IND,1:3) = MEAN_P2(IND,1:3) + P2(IND,I,1:3)
  MEAN_P3(IND,1:3) = MEAN_P3(IND,1:3) + P3(IND,I,1:3)
  
  IND_BOX = (floor( POSI(IND,I,1)/DBOX_X)) &
          + (floor( POSI(IND,I,2)/DBOX_Y))  * NBOX_DIR &
          + (floor( POSI(IND,I,3)/DBOX_Z))  * NBOX_DIR*NBOX_DIR  +1
          
  LOCAL_MEAN_PSWIM(IND,IND_BOX,1:3) = LOCAL_MEAN_PSWIM(IND,IND_BOX,1:3) + PSWIM(IND,I,1:3)
  COUNT_BOX(IND,IND_BOX) = COUNT_BOX(IND,IND_BOX) + 1.0
  
!~   if (IND==1) then
!~    if (I==1) then
!~     print*,'I = ', I
!~     print*, 'POSI(IND,I,1:3) = ', POSI(IND,I,1:3)
!~     print*,'IND_BOX = ', IND_BOX
!~     print*,'PSWIM(IND,I,:) = ', PSWIM(IND,I,:)
!~    print*,'P2(IND,I,:) = ', P2(IND,I,:)
!~    print*,'P3(IND,I,:) = ', P3(IND,I,:)
!~    read(*,*)
!~    end if
!~   end if
  ! In spherical coordinates (theta,phi)
  PSWIM_SPH(IND,I,1) = dacos( PSWIM(IND,I,3) ) 
  
  PSWIM_SPH(IND,I,2) = datan2( PSWIM(IND,I,2), PSWIM(IND,I,1) )
  if (PSWIM_SPH(IND,I,2)<0) then
   PSWIM_SPH(IND,I,2) = PSWIM_SPH(IND,I,2) + 2.0*PPI
  end if
  
  MEAN_PSWIM_SPH(IND,1:2) = MEAN_PSWIM_SPH(IND,1:2) + PSWIM_SPH(IND,I,1:2)
  
 end do
 

 MOD_MEAN_PSWIM  = dsqrt(MEAN_PSWIM(IND,1)**2 &
                      +  MEAN_PSWIM(IND,2)**2 &      
                      +  MEAN_PSWIM(IND,3)**2 ) 
 
 do J = 1,NBOX_TOT  
  if (COUNT_BOX(IND,J)>0.0) then
   LOCAL_MOD_MEAN_PSWIM  = dsqrt(LOCAL_MEAN_PSWIM(IND,J,1)**2 &
                           +  LOCAL_MEAN_PSWIM(IND,J,2)**2 &      
                           +  LOCAL_MEAN_PSWIM(IND,J,3)**2 ) 
                           
   LOCAL_MEAN_PSWIM(IND,J,1:3) = &
               LOCAL_MEAN_PSWIM(IND,J,1:3)/COUNT_BOX(IND,J)   
               
   LOCAL_POLAR_ORDER_TIME(IND,J) = LOCAL_MOD_MEAN_PSWIM/COUNT_BOX(IND,J)
  else
   LOCAL_POLAR_ORDER_TIME(IND,J) = -1.0
  end if
 end do  
 
            
    
 MOD_MEAN_P2  = dsqrt(MEAN_P2(IND,1)**2 &
                   +  MEAN_P2(IND,2)**2 &      
                   +  MEAN_P2(IND,3)**2 )    
                       
 MOD_MEAN_P3  = dsqrt(MEAN_P3(IND,1)**2 &
                   +  MEAN_P3(IND,2)**2 &      
                   +  MEAN_P3(IND,3)**2 ) 
                                             
 POLAR_ORDER_TIME(IND) = MOD_MEAN_PSWIM / real(PART_END-PART_START+1)
 
 
 
!~  if (IND==1) then
!~   print*,'MOD_MEAN_PSWIM  = ', MOD_MEAN_PSWIM 
!~   print*,'POLAR_ORDER_TIME(IND) = ', POLAR_ORDER_TIME(IND)
!~   read(*,*)
!~  end if
               
 MEAN_PSWIM_NON_NORM(IND,1:3) =  MEAN_PSWIM(IND,1:3)/ real(PART_END-PART_START+1)
 MEAN_PSWIM(IND,1:3) = MEAN_PSWIM(IND,1:3) / MOD_MEAN_PSWIM
 MEAN_P2(IND,1:3) = MEAN_P2(IND,1:3) / MOD_MEAN_P2
 
 SCAL_MEAN_P2_MEAN_PSWIM = MEAN_P2(IND,1)*MEAN_PSWIM(IND,1) &
                         + MEAN_P2(IND,2)*MEAN_PSWIM(IND,2) &
                         + MEAN_P2(IND,3)*MEAN_PSWIM(IND,3) 
                         
 MEAN_P2_ORTH(IND,1:3) = MEAN_P2(IND,1:3) &
                       - SCAL_MEAN_P2_MEAN_PSWIM * MEAN_PSWIM(IND,1:3)

 MOD_MEAN_P2 = dsqrt(MEAN_P2_ORTH(IND,1)**2 &
                  +  MEAN_P2_ORTH(IND,2)**2 &      
                  +  MEAN_P2_ORTH(IND,3)**2 )
                                         
 MEAN_P2_ORTH(IND,1:3) = MEAN_P2_ORTH(IND,1:3) / MOD_MEAN_P2
                       
 MEAN_P3_ORTH(IND,1) = MEAN_PSWIM(IND,2)*MEAN_P2_ORTH(IND,3) &
                     - MEAN_PSWIM(IND,3)*MEAN_P2_ORTH(IND,2)
                     
 MEAN_P3_ORTH(IND,2) = MEAN_PSWIM(IND,3)*MEAN_P2_ORTH(IND,1) &
                     - MEAN_PSWIM(IND,1)*MEAN_P2_ORTH(IND,3)
                     
 MEAN_P3_ORTH(IND,3) = MEAN_PSWIM(IND,1)*MEAN_P2_ORTH(IND,2) &
                     - MEAN_PSWIM(IND,2)*MEAN_P2_ORTH(IND,1)
                  
 MEAN_P3(IND,1:3) = MEAN_P3(IND,1:3) / MOD_MEAN_P3
 
!~  read(*,*)
 MEAN_PSWIM_SPH(IND,1:2) = MEAN_PSWIM_SPH(IND,1:2)/ real(PART_END-PART_START+1)

 if (POLAR_ORDER_TIME(IND).ge.1.0) then
  print*,'IND = ', IND 
  print*,'POLAR_ORDER_TIME(IND) = ', POLAR_ORDER_TIME(IND)
  print*,'real(PART_END-PART_START+1) = ', real(PART_END-PART_START+1)
  print*,'PART_END,PART_START+1, = ', PART_END, PART_START+1
!~  print*,'MEAN_PSWIM_SPH(IND,:) = ', MEAN_PSWIM_SPH(IND,:)
  print*,'MEAN_PSWIM(IND,:) = ', MEAN_PSWIM(IND,:)
!~  print*,'MEAN_P2_ORTH(IND,:) = ', MEAN_P2_ORTH(IND,:)
!~  print*,'MEAN_P3_ORTH(IND,:) = ', MEAN_P3_ORTH(IND,:)
!~                                     
  read(*,*)
  end if
 
 
 do I = 1, PART_END-PART_START+1      
  ! In spherical coordinates (theta,phi)
!~   MOD_VECTOR = dsqrt((PSWIM(IND,I,1)-MEAN_PSWIM(IND,1))**2 &
!~                   +  (PSWIM(IND,I,2)-MEAN_PSWIM(IND,2))**2 &      
!~                   +  (PSWIM(IND,I,3)-MEAN_PSWIM(IND,3))**2 )    
!~                   
!~   PSWIM_SPH_PMEAN(IND,I,1) = dacos( -(PSWIM(IND,I,1)-MEAN_PSWIM(IND,1))*MEAN_PSWIM(IND,1)/MOD_VECTOR &
!~                                     -(PSWIM(IND,I,2)-MEAN_PSWIM(IND,2))*MEAN_PSWIM(IND,2)/MOD_VECTOR &
!~                                     -(PSWIM(IND,I,3)-MEAN_PSWIM(IND,3))*MEAN_PSWIM(IND,3)/MOD_VECTOR ) 
                          
  PSWIM_SPH_PMEAN(IND,I,1) = dacos( PSWIM(IND,I,1)*MEAN_PSWIM(IND,1) &
                                  + PSWIM(IND,I,2)*MEAN_PSWIM(IND,2) &
                                  + PSWIM(IND,I,3)*MEAN_PSWIM(IND,3) ) 

                                                                   
!~   PSWIM_SPH_PMEAN(IND,I,2) = datan2( -(PSWIM(IND,I,1)-MEAN_PSWIM(IND,1))*MEAN_P3_ORTH(IND,1)/MOD_VECTOR &
!~                                    - (PSWIM(IND,I,2)-MEAN_PSWIM(IND,2))*MEAN_P3_ORTH(IND,2)/MOD_VECTOR &
!~                                    - (PSWIM(IND,I,3)-MEAN_PSWIM(IND,3))*MEAN_P3_ORTH(IND,3)/MOD_VECTOR, &
!~                                      -(PSWIM(IND,I,1)-MEAN_PSWIM(IND,1))*MEAN_P2_ORTH(IND,1)/MOD_VECTOR &
!~                                    - (PSWIM(IND,I,2)-MEAN_PSWIM(IND,2))*MEAN_P2_ORTH(IND,2)/MOD_VECTOR &
!~                                    - (PSWIM(IND,I,3)-MEAN_PSWIM(IND,3))*MEAN_P2_ORTH(IND,3)/MOD_VECTOR)
                                   
  PSWIM_SPH_PMEAN(IND,I,2) = datan2( PSWIM(IND,I,1)*MEAN_P3_ORTH(IND,1) &
                                   + PSWIM(IND,I,2)*MEAN_P3_ORTH(IND,2) &
                                   + PSWIM(IND,I,3)*MEAN_P3_ORTH(IND,3), &
                                     PSWIM(IND,I,1)*MEAN_P2_ORTH(IND,1) &
                                   + PSWIM(IND,I,2)*MEAN_P2_ORTH(IND,2) &
                                   + PSWIM(IND,I,3)*MEAN_P2_ORTH(IND,3))                                 
  if (PSWIM_SPH_PMEAN(IND,I,2)<0) then
   PSWIM_SPH_PMEAN(IND,I,2) = PSWIM_SPH_PMEAN(IND,I,2) + 2.0*PPI
  end if
  
  if (PSWIM_SPH_PMEAN(IND,I,1)+10.0 == PSWIM_SPH_PMEAN(IND,I,1)) then
   print*,'IND = ', IND
   print*,'I = ', I
   print*,'PSWIM_SPH_PMEAN(IND,I,1) = ', PSWIM_SPH_PMEAN(IND,I,1)
   read(*,*)
  end if
  
  if (PSWIM_SPH_PMEAN(IND,I,2)+10.0 == PSWIM_SPH_PMEAN(IND,I,2)) then
   print*,'IND = ', IND
   print*,'I = ', I
   print*,'PSWIM_SPH_PMEAN(IND,I,2) = ', PSWIM_SPH_PMEAN(IND,I,2)
   read(*,*)
  end if


 end do

 
end do

print*,'COMPUTE POLAR_ORDER_TIME--->  OK '


!!-Print filename
write(FILENAME,10200)'POLAR_ORDER.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(301,'(2(e17.7))') TIME_VEC(IND), POLAR_ORDER_TIME(IND)
  if (POLAR_ORDER_TIME(IND).ge.1.0) then
   print*,'IND = ', IND 
   print*,'POLAR_ORDER_TIME(IND) = ', POLAR_ORDER_TIME(IND)
   print*,'real(PART_END-PART_START+1) = ', real(PART_END-PART_START+1)
   print*,'PART_END,PART_START+1, = ', PART_END, PART_START+1
 !~  print*,'MEAN_PSWIM_SPH(IND,:) = ', MEAN_PSWIM_SPH(IND,:)
   print*,'MEAN_PSWIM(IND,:) = ', MEAN_PSWIM(IND,:)
 !~  print*,'MEAN_P2_ORTH(IND,:) = ', MEAN_P2_ORTH(IND,:)
 !~  print*,'MEAN_P3_ORTH(IND,:) = ', MEAN_P3_ORTH(IND,:)
 !~                                     
   read(*,*)
  end if
end do

!- close file
close(301)

!!-Print filename
write(FILENAME,10204)'LOCAL_POLAR_ORDER_',NBOX_TOT,'.dat'

!!- ASCII
open(unit=301,file=trim(FILENAME))
!!- Ecriture ASCII

write(MYFMT, '("(",I0,"(e17.7))")') NBOX_TOT+1

do IND = 1, NSAVES
 write(301,fmt=MYFMT) TIME_VEC(IND), ( LOCAL_POLAR_ORDER_TIME(IND,J), J=1,NBOX_TOT )
end do

!- close file
close(301)


!!-Print filename
write(FILENAME,10200)'MEAN_PSWIM.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

do IND = 1, NSAVES
 write(302,'(4(e17.7))') TIME_VEC(IND), MEAN_PSWIM_NON_NORM(IND,1), MEAN_PSWIM_NON_NORM(IND,2), MEAN_PSWIM_NON_NORM(IND,3)
end do

!- close file
close(302)


!!-Print filename
write(FILENAME,10204)'LOCAL_MEAN_PSWIM_',NBOX_TOT,'.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII

write(MYFMT, '("(",I0,"(e17.7))")') NBOX_TOT*3+1

do IND = 1, NSAVES
 write(302,fmt=MYFMT) TIME_VEC(IND), ( LOCAL_MEAN_PSWIM(IND,J,1:3), J=1,NBOX_TOT )
end do

!- close file
close(302)


!!-Print filename
write(FILENAME,10204)'LOCAL_BOX_COUNT_',NBOX_TOT,'.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

write(MYFMT, '("(",I0,"(e17.7))")') NBOX_TOT+1

do IND = 1, NSAVES
 write(303,fmt=MYFMT) TIME_VEC(IND), ( COUNT_BOX(IND,J), J=1,NBOX_TOT )
end do

!- close file
close(303)

print*,'SAVE POLAR_ORDER_TIME--->  OK '

print*,' '
print*,'-------------------END POLAR ORDER------------------------------ '
print*,' '

!!====================================================================
1999 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "n", "p"')
2001 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',e12.5)

1998 format ('VARIABLES = "xp noper", "yp noper", "zp noper"')
2002 format ('VARIABLES = "xp", "yp", "zp"')
2003 format ('VARIABLES = "up", "vp", "wp"')
2004 format ('VARIABLES = "ompx", "ompy", "ompz"')
2005 format ('VARIABLES = "quat1", "quat2", "quat3", "quat4"')
2006 format ('VARIABLES = "pswimx", "pswimy", "pswimz"')
2007 format ('VARIABLES = "Sxx", "Sxy", "Sxz", "Syy", "Syz"')
2008 format ('VARIABLES = "a1", "a2", "a3"')

2010 format ('NPART_FULL = ', i4)
2011 format ('NPART_FULL = ')

10200 format (A)
10201 format (A,I1,A)
10202 format (A,I2,A)
10203 format (A,I3,A)
10204 format (A,I4.4,A)
10205 format (I8.8)
10101 format (A,A,A)

end subroutine POLAR_ORDER
