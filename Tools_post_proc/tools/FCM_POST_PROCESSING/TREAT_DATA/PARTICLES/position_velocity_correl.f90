!!====================================================================
!!
!!
!!====================================================================

subroutine POSITION_VELOCITY_CORREL(NSAVES, NPART, L, RAD, ROT_MAT, POSI, VEL)

!!====================================================================
!!
!!
!!====================================================================


implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer, intent(in) :: NSAVES
! Number of swimmers
integer, intent(in) :: NPART
! Domain size
real(kind=8), intent(in):: L
! MAximal particle radius
real(kind=8), intent(in):: RAD
! Particle positions, orientation and velocities
real(kind=8), dimension(NSAVES,NPART,3,3), intent(in) :: ROT_MAT
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: POSI
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: VEL


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- r-component
real(kind=8), allocatable, dimension(:) :: RANGE_R
!- theta-component
real(kind=8), allocatable, dimension(:) :: RANGE_THETA
!~ real(kind=8), allocatable, dimension(:,:) :: RANGE_THETA_IRREG
!- phi-component
real(kind=8), allocatable, dimension(:) :: RANGE_PHI

!- Cor_Vel discretization
real(kind=8), allocatable, dimension(:) :: RANGE_VELCOR


!- g(r,theta,phi)
real(kind=8), allocatable, dimension(:,:,:) :: RDF
!- g(r,theta)
real(kind=8), allocatable, dimension(:,:) :: RDF_R_THETA
!- g(r)
real(kind=8), allocatable, dimension(:) :: RDF_R

!~ !- velcor_pdf(range_cor_vel,r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:,:) :: VELCOR_PDF
!~ !- velcor_pdf(range_cor_vel,r,theta)
!~ real(kind=8), allocatable, dimension(:,:,:) :: VELCOR_R_THETA_PDF
!~ !- velcor_pdf(range_cor_vel,r)
!~ real(kind=8), allocatable, dimension(:,:) :: VELCOR_R_PDF

!- Iu(vel,r,theta,phi)
real(kind=8), allocatable, dimension(:,:,:) :: VELCOR
!- Iu(r,theta)
real(kind=8), allocatable, dimension(:,:) :: VELCOR_R_THETA
!- Iu(r)
real(kind=8), allocatable, dimension(:) :: VELCOR_R

!-Vol discretized sector
real(kind=8), allocatable, dimension(:,:,:) :: VOL_SECTOR
!-Vol of ring for discretized sector for g(r,theta)
real(kind=8), allocatable, dimension(:,:) :: VOL_RING
!-Vol of shell for discretized sector for g(r)
real(kind=8), allocatable, dimension(:) :: VOL_SHELL

!- Step size to discretize the interval of positions and vel corr
real(kind=8) ::  DR, DTHETA, DPHI, DVELCOR
!~ real(kind=8), allocatable, dimension(:) :: DTHETA_IRREG


real(kind=8) ::  PPI

!- Half domain size
real(kind=8) ::  L_2

!-2*rad
real(kind=8) ::  TWORAD


!-Maximal distance accepted for RDF computation
real(kind=8) ::  MAX_DIST

!- Pos, Distance,...
real(kind=8), dimension(3) :: POS_I, POS_J 
real(kind=8), dimension(3) :: POS_IJ, POS_IJ_I, POS_IJ_J
real(kind=8) :: RIJ
real(kind=8), dimension(3,3) :: ROT_I, ROT_J
real(kind=8) :: PHI_I, PHI_J
real(kind=8) :: THETA_I, THETA_J

real(kind=8) :: SCAL_VEL_IJ
real(kind=8) :: VEL_I_SQ
real(kind=8) :: VELCOR_IJ


!- Number of steps to discretize the interval of pos and velocities
integer :: NR, NTHETA, NPHI, NVELCOR

!- From when to start correlation computation
integer :: SAVE_START

!- Indices for corr functions
integer :: IND_VELCOR, IND_R, IND_THETA_I, IND_THETA_J, IND_PHI_I, IND_PHI_J

!- File name 
character(len=40) :: FILENAME



!- Index
integer :: I, J, K,M, IND

!- Integers to signalize problems when writing
integer :: ERR_FILE_RDF


PPI = 4.0*datan(1.0)

print*,' '
print*,'-------------------END POSITIONS, VELOCITIES CORREL--------------------------- '
print*,' '

!=====================================================================
! 2. DEFINE DISCRETIZATION OF PDF
!=====================================================================
print*, 'DISCRETIZE RANGE OF POSITIONS AND VELOCITY CORRELATION'
DR = 0.01*RAD
DPHI = 0.1
DTHETA = 0.1


L_2 = L/2.0 
TWORAD = 2.0*RAD
MAX_DIST = 2.5*TWORAD
print*,'MAX_DIST/RAD = ', MAX_DIST/RAD

NR = ceiling((MAX_DIST - TWORAD)/DR)
NTHETA = ceiling((PPI-DTHETA)/DTHETA)
NPHI = ceiling((2.0*PPI)/DPHI)

DVELCOR = 0.1
NVELCOR = ceiling(2.0/(DVELCOR))+1

print*,'NVELCOR, NR, NTHETA, NPHI = '
print*, NVELCOR, NR, NTHETA, NPHI
     

allocate(RANGE_R(NR))
allocate(RANGE_THETA(NTHETA))
!~ allocate(RANGE_THETA_IRREG(NR,NPHI))
allocate(RANGE_PHI(NPHI))
allocate(RANGE_VELCOR(NVELCOR))

allocate(RDF(NR,NTHETA,NPHI))
allocate(VOL_SECTOR(NR,NTHETA,NPHI))
allocate(RDF_R_THETA(NR,NTHETA))
allocate(VOL_RING(NR,NTHETA))
allocate(RDF_R(NR))
allocate(VOL_SHELL(NR))

!~ allocate(VELCOR_PDF(NVELCOR,NR,NTHETA,NPHI))
!~ allocate(VELCOR_R_THETA_PDF(NVELCOR,NR,NTHETA))
!~ allocate(VELCOR_R_PDF(NVELCOR,NR))

allocate(VELCOR(NR,NTHETA,NPHI))
allocate(VELCOR_R_THETA(NR,NTHETA))
allocate(VELCOR_R(NR))

RANGE_R(1) = TWORAD + DR
do J = 2, NR
 RANGE_R(J) = RANGE_R(J-1) + DR
end do

RANGE_THETA(1) = DTHETA 
do J = 2, NTHETA
 RANGE_THETA(J) = RANGE_THETA(J-1) + DTHETA
end do

RANGE_PHI(1) = -PPI + DPHI 
do J = 2, NPHI
 RANGE_PHI(J) = RANGE_PHI(J-1) + DPHI
end do

RANGE_VELCOR(1) = -1.0 
do J = 2, NVELCOR
 RANGE_VELCOR(J) = RANGE_VELCOR(J-1) + DVELCOR
end do



print*, 'DISCRETIZE RANGE OF POSITIONS---> OK'

!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE RDF AND VELCOR '

RDF = 0.0
RDF_R_THETA = 0.0
RDF_R = 0.0

VELCOR = 0.0
VELCOR_R_THETA = 0.0
VELCOR_R = 0.0

!~ VELCOR_PDF = 0.0
!~ VELCOR_R_THETA_PDF = 0.0
!~ VELCOR_R_PDF = 0.0

VOL_SECTOR = 1000000000000000000.0
VOL_RING = 1000000000000000000.0
VOL_SHELL = 1000000000000000000.0

SAVE_START = 1
!~ SAVE_START = floor(0.0*real(NSAVES)/1.0)

M = 0

do IND = SAVE_START, NSAVES 

 M = M+1
 if (mod(M,int((NSAVES-SAVE_START + 1)/10))==0) then
  print*,'TIMESAVE N° ', M, ' / ', NSAVES-SAVE_START + 1
 end if

 do I = 1, NPART

  POS_I = POSI(IND,I,1:3)
  ROT_I = ROT_MAT(IND,I,1:3,1:3)
  
  VEL_I_SQ = VEL(IND,I,1)**2 + VEL(IND,I,2)**2 + VEL(IND,I,3)**2
  
  do J = I+1,NPART
  
   POS_J = POSI(IND,J,1:3)
   ROT_J = ROT_MAT(IND,J,1:3,1:3)
   
   POS_IJ = POS_I - POS_J
   POS_IJ = POS_IJ - L* real(int(POS_IJ/(L_2)))   
   
   RIJ = dsqrt(POS_IJ(1)**2 + POS_IJ(2)**2 + POS_IJ(3)**2)
   
   if ((RIJ>=TWORAD).and.(RIJ<MAX_DIST)) then
   
    IND_R = ceiling( ( RIJ-RANGE_R(1) )/DR ) + 1
    
!~     if (IND_R>NR) then
!~      print*,'!!!!   ERROR -->  IND_R> NR  !!!!!!!!! '
!~      print*,'IND_R = ' ,IND_R 
!~      read(*,*)
!~     end if
    
    POS_IJ_I = matmul(ROT_I,-POS_IJ)
    PHI_I = real(atan2(POS_IJ_I(2),POS_IJ_I(1)))
    THETA_I = real(acos(POS_IJ_I(3)/RIJ))
    
    IND_PHI_I = ceiling( ( PHI_I-RANGE_PHI(1) )/DPHI ) + 1
    IND_THETA_I = ceiling( ( THETA_I-RANGE_THETA(1) )/DTHETA ) + 1
        
    POS_IJ_J = matmul(ROT_J,POS_IJ)
    PHI_J = real(atan2(POS_IJ_J(2),POS_IJ_J(1)))
    THETA_J = real(acos(POS_IJ_J(3)/RIJ))
    
    IND_PHI_J = ceiling( ( PHI_J-RANGE_PHI(1) )/DPHI ) + 1
    IND_THETA_J = ceiling( ( THETA_J-RANGE_THETA(1) )/DTHETA ) + 1
    
    SCAL_VEL_IJ = VEL(IND,I,1)*VEL(IND,J,1) &
                + VEL(IND,I,2)*VEL(IND,J,2) &
                + VEL(IND,I,3)*VEL(IND,J,3) 
                
    VELCOR_IJ = SCAL_VEL_IJ/VEL_I_SQ
    
    IND_VELCOR = int(VELCOR_IJ/DVELCOR)
    
    
!~     if (IND_PHI_I>NPHI) then
!~      print*,'!!!!   ERROR -->  IND_PHI_I> NPHI  !!!!!!!!! '
!~      print*,'IND_PHI_I = ' ,IND_PHI_I 
!~      print*,'PHI_I = ' ,PHI_I 
!~      print*,'maxval(RANGE_PHI) = ' ,maxval(RANGE_PHI)
!~      read(*,*)
!~     end if
    
!~     if (IND_PHI_J>NPHI) then
!~      print*,'!!!!   ERROR -->  IND_PHI_J> NPHI  !!!!!!!!! '
!~      print*,'IND_PHI_J = ' ,IND_PHI_J 
!~      print*,'PHI_J = ' ,PHI_J
!~      print*,'maxval(RANGE_PHI) = ' ,maxval(RANGE_PHI)
!~      read(*,*)
!~      read(*,*)
!~     end if
    
    if (IND_THETA_I>NTHETA) then    
!~      print*,'!!!!   ERROR -->  IND_THETA_I> NTHETA  !!!!!!!!! '
!~      print*,'IND_THETA_I = ' ,IND_THETA_I 
!~      print*,'THETA_I = ' ,THETA_I 
!~      print*,'maxval(RANGE_THETA) = ' ,maxval(RANGE_THETA)
     IND_THETA_I  = IND_THETA_I  - 1
!~      !read(*,*)
    end if
    
    if (IND_THETA_J>NTHETA) then
!~      print*,'!!!!   ERROR -->  IND_THETA_J> NTHETA  !!!!!!!!! '
!~      print*,'IND_THETA_J = ' ,IND_THETA_J 
!~      print*,'THETA_J = ' ,THETA_J
!~      print*,'maxval(RANGE_THETA) = ' ,maxval(RANGE_THETA)
     IND_THETA_J  = IND_THETA_J  - 1
!~      !read(*,*)
    end if
    
    RDF(IND_R, IND_THETA_I, IND_PHI_I) = RDF(IND_R, IND_THETA_I, IND_PHI_I) + 1
    RDF(IND_R, IND_THETA_J, IND_PHI_J) = RDF(IND_R, IND_THETA_J, IND_PHI_J) + 1
    
    RDF_R_THETA(IND_R, IND_THETA_I) = RDF_R_THETA(IND_R, IND_THETA_I) + 1
    RDF_R_THETA(IND_R, IND_THETA_J) = RDF_R_THETA(IND_R, IND_THETA_J) + 1
    
    RDF_R(IND_R) = RDF_R(IND_R) + 2    
     
!~     VELCOR_PDF(IND_VELCOR,IND_R,IND_THETA_I, IND_PHI_I) = &
!~         VELCOR_PDF(IND_VELCOR,IND_R,IND_THETA_I, IND_PHI_I) + 1        
!~     VELCOR_PDF(IND_VELCOR,IND_R,IND_THETA_J, IND_PHI_J) = &
!~         VELCOR_PDF(IND_VELCOR,IND_R,IND_THETA_J, IND_PHI_J) + 1
        
    VELCOR(IND_R,IND_THETA_I, IND_PHI_I) = &
        VELCOR(IND_R,IND_THETA_I, IND_PHI_I) + VELCOR_IJ
    VELCOR(IND_R,IND_THETA_J, IND_PHI_J) = &
        VELCOR(IND_R,IND_THETA_J, IND_PHI_J) + VELCOR_IJ
     
!~     VELCOR_R_THETA_PDF(IND_VELCOR,IND_R,IND_THETA_I) = &
!~         VELCOR_R_THETA_PDF(IND_VELCOR,IND_R,IND_THETA_I) + 1        
!~     VELCOR_R_THETA_PDF(IND_VELCOR,IND_R,IND_THETA_J) = &
!~         VELCOR_R_THETA_PDF(IND_VELCOR,IND_R,IND_THETA_J) + 1
        
        
    VELCOR_R_THETA(IND_R,IND_THETA_I) = &
        VELCOR_R_THETA(IND_R,IND_THETA_I) + VELCOR_IJ        
    VELCOR_R_THETA(IND_R,IND_THETA_J) = &
        VELCOR_R_THETA(IND_R,IND_THETA_J) + VELCOR_IJ
        
!~     VELCOR_R_PDF(IND_VELCOR,IND_R) = &
!~         VELCOR_R_PDF(IND_VELCOR,IND_R) + 2
        
    VELCOR_R(IND_R) = &
        VELCOR_R(IND_R) + 2.0*VELCOR_IJ
    
    VOL_SECTOR(IND_R, IND_THETA_I, IND_PHI_I) = &
        RANGE_R(IND_R)**(2)*sin(RANGE_THETA(IND_THETA_I))*DR*DTHETA*DPHI        
    VOL_SECTOR(IND_R, IND_THETA_J, IND_PHI_J) = &
        RANGE_R(IND_R)**(2)*sin(RANGE_THETA(IND_THETA_J))*DR*DTHETA*DPHI
        
    VOL_RING(IND_R, IND_THETA_I) = &
        RANGE_R(IND_R)**(2)*sin(RANGE_THETA(IND_THETA_I))*DR*DTHETA*2.0*PPI
    VOL_RING(IND_R, IND_THETA_J) = &
        RANGE_R(IND_R)**(2)*sin(RANGE_THETA(IND_THETA_J))*DR*DTHETA*2.0*PPI
        
    VOL_SHELL(IND_R) = &
        4.0/3.0*PPI*( (RANGE_R(IND_R)+DR)**3 - RANGE_R(IND_R)**3 )
        
!~         print*,'VOL_SHELL(IND_R) = ', VOL_SHELL(IND_R)
!~         print*,'RANGE_R(IND_R) = ', RANGE_R(IND_R)
!~         print*,'DR = ', DR
!~         read(*,*)
        
!~     print*,'J = ', J
!~     print*,'POS_IJ = ', POS_IJ
!~     print*,'RIJ = ', RIJ
!~     print*,'ROT_I = ', ROT_I
    
   end if
   
  end do  
 end do
end do

RDF = RDF/(NPART*(NSAVES-SAVE_START + 1))
RDF_R_THETA  = RDF_R_THETA /(NPART*(NSAVES-SAVE_START + 1))
RDF_R = RDF_R/(NPART*(NSAVES-SAVE_START + 1))

print*,'maxval(RDF ) = ',maxval(RDF )
print*,'maxval(RDF_R_THETA ) = ',maxval(RDF_R_THETA )
print*,'maxval(RDF_R ) = ',maxval(RDF_R )

RDF = RDF/(NPART/L**(3)*VOL_SECTOR)
RDF_R_THETA = RDF_R_THETA /(NPART/L**(3)*VOL_RING)
RDF_R = RDF_R /(NPART/L**(3)*VOL_SHELL)

print*,'maxval(RDF/(c*Vol_sector) ) = ',maxval(RDF )
print*,'maxval(RDF_R_THETA/(c*Vol_ring) ) = ',maxval(RDF_R_THETA )
print*,'maxval(RDF_R/(c*Vol_shell) ) = ',maxval(RDF_R )

VELCOR = VELCOR/(NPART*(NSAVES-SAVE_START + 1))
VELCOR_R_THETA  = VELCOR_R_THETA /(NPART*(NSAVES-SAVE_START + 1))
VELCOR_R = VELCOR_R/(NPART*(NSAVES-SAVE_START + 1))
print*,'maxval(VELCOR ) = ',maxval(VELCOR )
print*,'maxval(VELCOR_R_THETA ) = ',maxval(VELCOR_R_THETA )
print*,'maxval(VELCOR_R ) = ',maxval(VELCOR_R )
  
!~ 
!=====================================================================
! 3. WRITE RDF'S
!=====================================================================
!- Save only if <100MB
if (NR*NTHETA*NPHI*4.0*16.0/10.0**6<100.0) then
 print*,'SAVE RDF  '
 print*,NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
 print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '


 !---------- g(r,theta,phi)-----------
 !!-Print filename
 write(FILENAME,10200)'RDF.dat'

 !!- ASCII
 open(unit=301,file=trim(FILENAME))
 !!- Ecriture ASCII
 write(301,'(4(i17.7))') NR, NTHETA, NPHI, NR*NTHETA*NPHI
 do I = 1, NR
  do J = 1, NTHETA
   do K = 1, NPHI
    write(301,'(4(e17.7))') RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), RDF(I,J,K)
   end do
  end do
 end do

 !- close file
 close(301)
 print*,'SAVE RDF -->  OK '
else
 print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
end if

print*,'SAVE RDF_R_THETA  '
print*,NR*NTHETA*3.0,' VALUES TO SAVE  '
print*,NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- g(r,theta)-----------
!!-Print filename
write(FILENAME,10200)'RDF_R_THETA.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))') NR, NTHETA, NR*NTHETA
do I = 1, NR
 do J = 1, NTHETA
   write(302,'(3(e17.7))') RANGE_R(I), RANGE_THETA(J), RDF_R_THETA(I,J)
 end do
end do

!- close file
close(302)


print*,'SAVE RDF_R_THETA -->  OK '

print*,'SAVE RDF_R  '
print*,NR*2.0,' VALUES TO SAVE  '
print*,NR*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- g(r,theta,phi)-----------
!!-Print filename
write(FILENAME,10200)'RDF_R.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1, NR
   write(303,'(2(e17.7))') RANGE_R(I), RDF_R(I)
end do

!- close file
close(303)


print*,'SAVE RDF_R -->  OK '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Save only if <200MB
if (NR*NTHETA*NPHI*4.0*16.0/10.0**6<200.0) then
 print*,'SAVE VELCOR  '

 print*,NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
 print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
 !---------- velcor(r,theta,phi)-----------
 !!-Print filename
 write(FILENAME,10200)'VELCOR.dat'

 !!- ASCII
 open(unit=301,file=trim(FILENAME))
 !!- Ecriture ASCII
 write(301,'(4(i17.7))')  NR, NTHETA, NPHI, NR*NTHETA*NPHI


 do I = 1, NR
  do J = 1, NTHETA
   do K = 1, NPHI
    write(301,'(4(e17.7))') &
     RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), VELCOR(I,J,K)
   end do
  end do
 end do


 !- close file
 close(301)


 print*,'SAVE VELCOR -->  OK '
else
 print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
end if

print*,'SAVE VELCOR_R_THETA  '
print*,NR*NTHETA*3.0,' VALUES TO SAVE  '
print*,NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- velcor(r,theta)-----------
!!-Print filename
write(FILENAME,10200)'VELCOR_R_THETA.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))')  NR, NTHETA, NR*NTHETA


do I = 1, NR
 do J = 1, NTHETA
   write(302,'(3(e17.7))') &
    RANGE_R(I), RANGE_THETA(J), VELCOR_R_THETA(I,J)
 end do
end do


!- close file
close(302)


print*,'SAVE VELCOR_R_THETA -->  OK '

print*,'SAVE VELCOR_R  '
print*,NR*2.0,' VALUES TO SAVE  '
print*,NR*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- velcor(r)-----------
!!-Print filename
write(FILENAME,10200)'VELCOR_R.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NR
   write(303,'(2(e17.7))') &
    RANGE_R(I), VELCOR_R(I)
end do


!- close file
close(303)


print*,'SAVE VELCOR_R -->  OK '


!~ !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ !- Save only if <100MB
!~ if (NVELCOR*NR*NTHETA*NPHI*4.0*16.0/10.0**6<200.0) then
!~  print*,'SAVE VELCOR_PDF  '
!~ 
!~  print*,NVELCOR*NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NVELCOR*NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~  !---------- g(r,theta,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'VELCOR_PDF.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(5(i17.7))') NVELCOR, NR, NTHETA, NPHI, NVELCOR*NR*NTHETA*NPHI
!~ 
!~  do M = 1, NVELCOR
!~   do I = 1, NR
!~    do J = 1, NTHETA
!~     do K = 1, NPHI
!~      write(301,'(5(e17.7))') &
!~      RANGE_VELCOR(M), RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), VELCOR_PDF(M,I,J,K)
!~     end do
!~    end do
!~   end do
!~  end do
!~ 
!~  !- close file
!~  close(301)
!~ 
!~ 
!~  print*,'SAVE VELCOR_PDF -->  OK '
!~ else
!~  print*,NVELCOR*NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if
!~ 
!~ print*,'SAVE VELCOR_R_THETA_PDF  '
!~ print*,NVELCOR*NR*NTHETA*3.0,' VALUES TO SAVE  '
!~ print*,NVELCOR*NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!~ !---------- g(r,theta)-----------
!~ !!-Print filename
!~ write(FILENAME,10200)'VELCOR_R_THETA_PDF.dat'
!~ 
!~ !!- ASCII
!~ open(unit=302,file=trim(FILENAME))
!~ !!- Ecriture ASCII
!~ write(302,'(4(i17.7))') NVELCOR, NR, NTHETA, NVELCOR*NR*NTHETA
!~ 
!~ do M = 1, NVELCOR
!~  do I = 1, NR
!~   do J = 1, NTHETA
!~     write(302,'(4(e17.7))') &
!~     RANGE_VELCOR(M), RANGE_R(I), RANGE_THETA(J), VELCOR_R_THETA_PDF(M,I,J)
!~   end do
!~  end do
!~ end do
!~ 
!~ !- close file
!~ close(302)
!~ 
!~ 
!~ print*,'SAVE VELCOR_R_THETA_PDF -->  OK '
!~ 
!~ print*,'SAVE VELCOR_R_PDF  '
!~ print*,NVELCOR*NR*2.0,' VALUES TO SAVE  '
!~ print*,NVELCOR*NR*2.0*16.0/10.0**6,' MB TO STORE  '
!~ !---------- g(r,theta,phi)-----------
!~ !!-Print filename
!~ write(FILENAME,10200)'VELCOR_R_PDF.dat'
!~ 
!~ !!- ASCII
!~ open(unit=303,file=trim(FILENAME))
!~ !!- Ecriture ASCII
!~ write(303,'(3(i17.7))') NVELCOR, NR, NVELCOR*NR
!~ 
!~ do M = 1, NVELCOR
!~  do I = 1, NR
!~     write(303,'(3(e17.7))') &
!~     RANGE_VELCOR(M), RANGE_R(I), VELCOR_R_PDF(M,I)
!~  end do
!~ end do
!~ 
!~ !- close file
!~ close(303)
!~ 
!~ 
!~ print*,'SAVE VELCOR_R_PDF -->  OK '







print*,' '
print*,'-------------------END POSITIONS, VELOCITIES CORREL--------------------------- '
print*,' '

 
deallocate(RANGE_R)
deallocate(RANGE_THETA)
!~ allocate(RANGE_THETA_IRREG(NR,NPHI))
deallocate(RANGE_PHI)
deallocate(RANGE_VELCOR)

deallocate(RDF)
deallocate(VOL_SECTOR)
deallocate(RDF_R_THETA)
deallocate(VOL_RING)
deallocate(RDF_R)
deallocate(VOL_SHELL)

deallocate(VELCOR)
deallocate(VELCOR_R_THETA)
deallocate(VELCOR_R)

!~ deallocate(VELCOR_PDF(NVELCOR,NR,NTHETA,NPHI))
!~ deallocate(VELCOR_R_THETA_PDF(NVELCOR,NR,NTHETA))
!~ deallocate(VELCOR_R_PDF(NVELCOR,NR))



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
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10205 format (I8.8)
10101 format (A,A,A)

end subroutine POSITION_VELOCITY_CORREL
