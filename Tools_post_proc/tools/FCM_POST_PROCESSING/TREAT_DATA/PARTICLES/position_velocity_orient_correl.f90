!!====================================================================
!!
!!
!!====================================================================

subroutine POSITION_VELOCITY_ORIENT_CORREL(NSAVES, &
                                           SAVE_START, &
                                           NPART, &
                                           NPSTAT, &
                                           L, &
                                           RAD, &
                                           ROT_MAT, &
                                           POSI, &
                                           VEL        )

!!====================================================================
!! TO IMPROVE : COMPUTATION OF RDF_R_THETA WITH RINGS VOLUME!!!!!
!! TO IMPROVE : COMPUTATION OF RDF_R_THETA WITH RINGS VOLUME!!!!!
!! TO IMPROVE : COMPUTATION OF RDF_R_THETA WITH RINGS VOLUME!!!!!
!! TO IMPROVE : COMPUTATION OF RDF_R_THETA WITH RINGS VOLUME!!!!!
!!====================================================================


implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer, intent(in) :: NSAVES
! Starting snapshot to treat
integer, intent(in) :: SAVE_START
! Number of swimmers
integer, intent(in) :: NPART
! Number of swimmers to treat
integer, intent(in) :: NPSTAT
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

!! Per time step values
!~ !- g(r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: RDF
!- g(r,theta)
real(kind=8), allocatable, dimension(:,:) :: RDF_R_THETA
!- g(r)
real(kind=8), allocatable, dimension(:) :: RDF_R
!~ !- Iu(vel,r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: VELCOR
!- Iu(r,theta)
real(kind=8), allocatable, dimension(:,:) :: VELCOR_R_THETA
!- Iu(r)
real(kind=8), allocatable, dimension(:) :: VELCOR_R
!~ !- Ie(vel,r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: POLARCOR
!- Ie(r,theta)
real(kind=8), allocatable, dimension(:,:) :: POLARCOR_R_THETA
!- Ie(r)
real(kind=8), allocatable, dimension(:) :: POLARCOR_R

!!Mean values
!~ !- g(r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: RDF_MEAN
!- g(r,theta)
real(kind=8), allocatable, dimension(:,:) :: RDF_R_THETA_MEAN
!- g(r)
real(kind=8), allocatable, dimension(:) :: RDF_R_MEAN
!~ !- Iu(vel,r,theta,phi) mean
!~ real(kind=8), allocatable, dimension(:,:,:) :: VELCOR_MEAN
!- Iu(r,theta)
real(kind=8), allocatable, dimension(:,:) :: VELCOR_R_THETA_MEAN
!- Iu(r)
real(kind=8), allocatable, dimension(:) :: VELCOR_R_MEAN
!~ !- Ie(vel,r,theta,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: POLARCOR_MEAN
!- Ie(r,theta)
real(kind=8), allocatable, dimension(:,:) :: POLARCOR_R_THETA_MEAN
!- Ie(r)
real(kind=8), allocatable, dimension(:) :: POLARCOR_R_MEAN

!- Correlation integral along R for Npairs at a given time
real(kind=8), allocatable, dimension(:,:) :: CORREL_INTEGRAL_PART
!- Correlation integral averaged over pairs along time
real(kind=8), allocatable, dimension(:,:) :: CORREL_INTEGRAL_MEAN_TIME
!- Test to compare mean computations
real(kind=8), allocatable, dimension(:) :: TEST_MEAN

!~ !-Vol discretized sector
!~ real(kind=8), allocatable, dimension(:,:,:) :: VOL_SECTOR
!-Vol of ring for discretized sector for g(r,theta)
real(kind=8), allocatable, dimension(:,:) :: VOL_RING
!-Vol of shell for discretized sector for g(r)
real(kind=8), allocatable, dimension(:) :: VOL_SHELL

!- Vector to select randomly NPSTAT particles
real(kind=8), dimension(5*NPART) :: VEC_RAND
integer, dimension(5*NPART) :: IND_RAND
integer, dimension(5*NPART) :: IND_RAND_UNIQUE

!- Step size to discretize the interval of positions and vel corr
real(kind=8) ::  DR, DTHETA, DPHI
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

real(kind=8) :: POLARCOR_IJ


!- Correlation dimension
real(kind=8), dimension(NPART) :: CORREL_DIM_PART
real(kind=8), dimension(NSAVES-SAVE_START + 1) :: CORREL_DIM_MEAN_TIME

!- Number of steps to discretize the interval of pos and velocities
integer :: NR, NTHETA, NPHI


integer :: NUNIQUE

!- Indices for corr functions
integer :: IND_VELCOR, IND_R, IND_THETA_I, IND_THETA_J, IND_PHI_I, IND_PHI_J

!- File name 
character(len=40) :: FILENAME



!- Index
integer :: I, J, K,M, IND, IND_I, IND_J

!- Integers to signalize problems when writing
integer :: ERR_FILE_RDF


PPI = 4.0*datan(1.0d0)

print*,' '
print*,'-------------------END POSITIONS, VELOCITIES CORREL--------------------------- '
print*,' '

!=====================================================================
! 2. DEFINE DISCRETIZATION OF PDF
!=====================================================================
print*, 'DISCRETIZE RANGE OF POSITIONS AND VELOCITY CORRELATION'

L_2 = L/2.0 
TWORAD = 2.0*RAD
MAX_DIST = L_2
print*,'MAX_DIST/RAD = ', MAX_DIST/RAD

DR = 0.1*RAD
NR = ceiling((MAX_DIST - TWORAD)/DR)

NTHETA = 30
NPHI = 60

!~ DTHETA = (PPI)/(real(NTHETA)-1.0)
DTHETA = (PPI)/(real(NTHETA))
DPHI = 2.0*PPI/(real(NPHI) -1.0)





print*,'NR, NTHETA, NPHI = '
print*, NR, NTHETA, NPHI
     

allocate(RANGE_R(NR))
allocate(RANGE_THETA(NTHETA))
!~ allocate(RANGE_THETA_IRREG(NR,NPHI))
allocate(RANGE_PHI(NPHI))

!~ allocate(VOL_SECTOR(NR,NTHETA,NPHI))
allocate(VOL_RING(NR,NTHETA))
allocate(VOL_SHELL(NR))

!~ allocate(RDF(NR,NTHETA,NPHI))
allocate(RDF_R_THETA(NR,NTHETA))
allocate(RDF_R(NR))


!~ allocate(VELCOR(NR,NTHETA,NPHI))
allocate(VELCOR_R_THETA(NR,NTHETA))
allocate(VELCOR_R(NR))

!~ allocate(POLARCOR(NR,NTHETA,NPHI))
allocate(POLARCOR_R_THETA(NR,NTHETA))
allocate(POLARCOR_R(NR))

!~ allocate(RDF_MEAN(NR,NTHETA,NPHI))
allocate(RDF_R_THETA_MEAN(NR,NTHETA))
allocate(RDF_R_MEAN(NR))

!~ allocate(VELCOR_MEAN(NR,NTHETA,NPHI))
allocate(VELCOR_R_THETA_MEAN(NR,NTHETA))
allocate(VELCOR_R_MEAN(NR))

!~ allocate(POLARCOR_MEAN(NR,NTHETA,NPHI))
allocate(POLARCOR_R_THETA_MEAN(NR,NTHETA))
allocate(POLARCOR_R_MEAN(NR))

allocate(CORREL_INTEGRAL_PART(NPART,NR))
allocate(CORREL_INTEGRAL_MEAN_TIME(NSAVES-SAVE_START + 1,NR))
allocate(TEST_MEAN(NR))

RANGE_R(1) = TWORAD
do J = 2, NR
 RANGE_R(J) = RANGE_R(J-1) + DR
end do

RANGE_THETA(1) = 0.0
do J = 2, NTHETA
 RANGE_THETA(J) = RANGE_THETA(J-1) + DTHETA
end do


!~ RANGE_PHI(1) = -PPI + DPHI 
!~ do J = 2, NPHI
!~  RANGE_PHI(J) = RANGE_PHI(J-1) + DPHI
!~ end do



do I = 1, NR
 VOL_SHELL(I) = &
     4.0/3.0*PPI*( (RANGE_R(I)+DR)**3 - RANGE_R(I)**3 )
     
 do J = 1, NTHETA      
!~    VOL_SECTOR(I, J, 1:NPHI) = &
!~        RANGE_R(I)**(2)*sin(RANGE_THETA(J))*DR*DTHETA*DPHI

   VOL_RING(I, J) = &
       1.0/3.0*((RANGE_R(I)+DR)**3- RANGE_R(I)**3)*( dcos(RANGE_THETA(J)) &
                                                   - dcos(RANGE_THETA(J)+DTHETA) )*2.0*PPI        

 end do
end do
 



print*, 'DISCRETIZE RANGE OF POSITIONS---> OK'

!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE RDF AND VELCOR AND POLARCOR'

!~ RDF_MEAN = 0.0
RDF_R_THETA_MEAN = 0.0
RDF_R_MEAN = 0.0

!~ VELCOR_MEAN = 0.0
VELCOR_R_THETA_MEAN = 0.0
VELCOR_R_MEAN = 0.0

!~ POLARCOR_MEAN = 0.0
POLARCOR_R_THETA_MEAN = 0.0
POLARCOR_R_MEAN = 0.0

CORREL_INTEGRAL_PART = 0.0
CORREL_INTEGRAL_MEAN_TIME = 0.0
TEST_MEAN = 0.0
print*,'NB OF PART USED TO COMPUTE CORRELATIONS = ', NPSTAT


if (NPSTAT==NPART) then
 do I = 1, NPART
  IND_RAND_UNIQUE(I) = I
 end do
end if


print*,'START FROM TIME SAVE N° = ', SAVE_START

M = 0

do IND = SAVE_START, NSAVES 

 M = M + 1
 
 if (NPSTAT<NPART) then   
  call init_random_seed()
  call random_number(VEC_RAND)

  
  IND_RAND = int(VEC_RAND*10.0**9)
  
  IND_RAND = int(real(IND_RAND)*real(NPART)/10.0**9) +1
  
  call REMOVE_DUPS(5*NPART,IND_RAND,NUNIQUE,IND_RAND_UNIQUE)
  if (NUNIQUE<NPSTAT) then
   print*,'NUNIQUE = ', NUNIQUE
   print*,'Too much repeated elements!!!'
  end if 
 end if
 
 
 if (mod(M,int((NSAVES-SAVE_START + 1)/10))==0) then
  print*,'TIMESAVE N° ', M, ' / ', NSAVES-SAVE_START + 1
 end if
 
!~  RDF = 0.0
 RDF_R_THETA = 0.0
 RDF_R = 0.0
 
!~  VELCOR = 0.0
 VELCOR_R_THETA = 0.0
 VELCOR_R = 0.0

!~  POLARCOR = 0.0
 POLARCOR_R_THETA = 0.0
 POLARCOR_R = 0.0

 do IND_I = 1, NPSTAT
  
  
  I = IND_RAND_UNIQUE(IND_I)
  
  if ((I<1).or.(I>NPART)) then
   print*,'IND_I = ', IND_I
   print*,'NUNIQUE = ', NUNIQUE
   print*,'I = ', I
   read(*,*)
  end if
 

  POS_I = POSI(IND,I,1:3)
  ROT_I = ROT_MAT(IND,I,1:3,1:3)
  
  VEL_I_SQ = VEL(IND,I,1)**2 + VEL(IND,I,2)**2 + VEL(IND,I,3)**2
    
  do IND_J = IND_I+1, NPSTAT
  
   J = IND_RAND_UNIQUE(IND_J)   
  
   POS_J = POSI(IND,J,1:3)
   ROT_J = ROT_MAT(IND,J,1:3,1:3)
   
   POS_IJ = POS_I - POS_J
    
   
   POS_IJ = POS_IJ - L* real(int(POS_IJ/(L_2)))   
   
   RIJ = dsqrt(POS_IJ(1)**2 + POS_IJ(2)**2 + POS_IJ(3)**2)
   
   if ((RIJ>=TWORAD).and.(RIJ<MAX_DIST)) then
   
!~     IND_R = ceiling( ( RIJ-RANGE_R(1) )/DR ) + 1
    IND_R = ceiling( ( RIJ-RANGE_R(1) )/DR ) 
        
    POS_IJ_I = matmul(ROT_I,-POS_IJ)
    
!~    print*,'POS_I = ', POS_I
!~    print*,'POS_J = ', POS_J
!~    print*,'POS_IJ = '
!~    print*, POS_IJ
!~    read(*,*)
!~    
!~    print*,'ROT_I(1,:) = '
!~    print*, ROT_I(1,:)
!~    print*,'ROT_I(2,:) = '
!~    print*, ROT_I(2,:)
!~    print*,'ROT_I(3,:) = '
!~    print*, ROT_I(3,:)
!~    print*,'POS_IJ_I = '
!~    print*, POS_IJ_I
!~    read(*,*)
   
   
  
   
!~     PHI_I = real(atan2(POS_IJ_I(2),POS_IJ_I(1)))
    THETA_I = real(acos(POS_IJ_I(3)/RIJ))
    
!~     IND_PHI_I = ceiling( ( PHI_I-RANGE_PHI(1) )/DPHI ) + 1
    IND_THETA_I = ceiling( ( THETA_I )/DTHETA ) 
    
!~     print*,'RIJ = ', RIJ
!~     print*,'IND_R = ', IND_R
!~     print*,'RANGE_R(IND_R) = ', RANGE_R(IND_R)
!~     print*,'RANGE_R(IND_R+1) = ', RANGE_R(IND_R+1)
!~     
!~     print*,'THETA_I = ', THETA_I
!~     print*,'IND_THETA_I = ', IND_THETA_I
!~     print*,'RANGE_THETA(IND_THETA_I) = ', RANGE_THETA(IND_THETA_I)
!~     print*,'RANGE_THETA(IND_THETA_I+1) = ', RANGE_THETA(IND_THETA_I+1)
!~     read(*,*)
        
    POS_IJ_J = matmul(ROT_J,POS_IJ)
    
    
!~     PHI_J = real(atan2(POS_IJ_J(2),POS_IJ_J(1)))

!~     print*,'I = ', I
!~     print*,'J = ', J
!~     print*,'RIJ = ', RIJ
!~     print*,' matmul(ROT_I,transpose(ROT_I)) = ', matmul(ROT_I,transpose(ROT_I))
!~     print*,'POS_IJ_I = ', POS_IJ_I
!~     print*,'NORM(POS_IJ_I) = ', dsqrt(POS_IJ_I(1)**2 + POS_IJ_I(2)**2 + POS_IJ_I(3)**2 )
!~     print*,' matmul(ROT_J,transpose(ROT_J)) = ' , matmul(ROT_J,transpose(ROT_J))
!~     print*,'POS_IJ_J = ', POS_IJ_J
!~     print*,'NORM(POS_IJ_J) = ', dsqrt(POS_IJ_J(1)**2 + POS_IJ_J(2)**2 + POS_IJ_J(3)**2 )
!~     read(*,*)

    THETA_J = real(acos(POS_IJ_J(3)/RIJ))
    
!~     IND_PHI_J = ceiling( ( PHI_J-RANGE_PHI(1) )/DPHI ) + 1
    IND_THETA_J = ceiling( ( THETA_J )/DTHETA )
    
    SCAL_VEL_IJ = VEL(IND,I,1)*VEL(IND,J,1) &
                + VEL(IND,I,2)*VEL(IND,J,2) &
                + VEL(IND,I,3)*VEL(IND,J,3) 
                
    VELCOR_IJ = SCAL_VEL_IJ/VEL_I_SQ
    
!~     IND_VELCOR = int(VELCOR_IJ/DVELCOR)

    POLARCOR_IJ = ROT_I(3,1)*ROT_J(3,1) &
                + ROT_I(3,2)*ROT_J(3,2) &
                + ROT_I(3,3)*ROT_J(3,3)
                
!~     POLARCOR_IJ = ROT_I(1,1)*ROT_J(1,1) &
!~                 + ROT_I(1,2)*ROT_J(1,2) &
!~                 + ROT_I(1,3)*ROT_J(1,3)
      
    
!~     if ((IND_PHI_I>NPHI).or.(IND_PHI_I<1)) then
!~      print*,'!!!!   ERROR -->  IND_PHI_I> NPHI  !!!!!!!!! '
!~      print*,'IND_PHI_I = ' ,IND_PHI_I 
!~      print*,'PHI_I = ' ,PHI_I 
!~      print*,'maxval(RANGE_PHI) = ' ,maxval(RANGE_PHI)
!~      read(*,*)
!~     end if
!~     
!~     if ((IND_PHI_J>NPHI).or.(IND_PHI_J<1)) then
!~      print*,'!!!!   ERROR -->  IND_PHI_J> NPHI  !!!!!!!!! '
!~      print*,'IND_PHI_J = ' ,IND_PHI_J 
!~      print*,'PHI_J = ' ,PHI_J
!~      print*,'maxval(RANGE_PHI) = ' ,maxval(RANGE_PHI)
!~      read(*,*)
!~      read(*,*)
!~     end if

    
    if ((IND_THETA_I>NTHETA).or.(IND_THETA_I<1)) then    
     print*,'!!!!   ERROR -->  IND_THETA_I> NTHETA  !!!!!!!!! '
     print*,'IND_THETA_I = ' ,IND_THETA_I 
     print*,'THETA_I = ' ,THETA_I 
     print*,'maxval(RANGE_THETA) = ' ,maxval(RANGE_THETA)
!~      IND_THETA_I  = IND_THETA_I  - 1
     read(*,*)
    end if
    
    if ((IND_THETA_J>NTHETA).or.(IND_THETA_J<1)) then
     print*,'!!!!   ERROR -->  IND_THETA_J> NTHETA  !!!!!!!!! '
     print*,'IND_THETA_J = ' ,IND_THETA_J 
     print*,'THETA_J = ' ,THETA_J
     print*,'maxval(RANGE_THETA) = ' ,maxval(RANGE_THETA)
!~      IND_THETA_J  = IND_THETA_J  - 1
     read(*,*)
    end if

    
    RDF_R_THETA(IND_R, IND_THETA_I) = RDF_R_THETA(IND_R, IND_THETA_I) + 1.0
    RDF_R_THETA(IND_R, IND_THETA_J) = RDF_R_THETA(IND_R, IND_THETA_J) + 1.0
    
    RDF_R(IND_R) = RDF_R(IND_R) + 2.0 
!~     RDF_R(IND_R) = RDF_R(IND_R) + 1    
     
    CORREL_INTEGRAL_PART(IND_I,IND_R) = CORREL_INTEGRAL_PART(IND_I,IND_R) + 1.0

    VELCOR_R_THETA(IND_R,IND_THETA_I) = &
       VELCOR_R_THETA(IND_R,IND_THETA_I) + VELCOR_IJ       
    VELCOR_R_THETA(IND_R,IND_THETA_J) = &
        VELCOR_R_THETA(IND_R,IND_THETA_J) + VELCOR_IJ 
       
    VELCOR_R(IND_R) = &
        VELCOR_R(IND_R) + 2.0*VELCOR_IJ 
!~     VELCOR_R(IND_R) = &
!~         VELCOR_R(IND_R) + VELCOR_IJ
        

    POLARCOR_R_THETA(IND_R,IND_THETA_I) = &
       POLARCOR_R_THETA(IND_R,IND_THETA_I) + POLARCOR_IJ         
    POLARCOR_R_THETA(IND_R,IND_THETA_J) = &
        POLARCOR_R_THETA(IND_R,IND_THETA_J) + POLARCOR_IJ 
       
    POLARCOR_R(IND_R) = &
        POLARCOR_R(IND_R) + 2.0*POLARCOR_IJ  
!~     POLARCOR_R(IND_R) = &
!~         POLARCOR_R(IND_R) + POLARCOR_IJ
                
   end if
  
   
  end do   !do IND_J = IND_I+1, NPSTAT

 end do !do IND_I = 1, NPSTAT
 

 
 do I = 1,NR
   if (RDF_R(I)>0) then
   VELCOR_R(I) = VELCOR_R(I)/RDF_R(I)!*(real(NPART)/real(NPSTAT))
   POLARCOR_R(I) = POLARCOR_R(I)/RDF_R(I)!*(real(NPART)/real(NPSTAT)) 
  end if  
  do J = 1, NTHETA
   if (RDF_R_THETA(I,J)>0) then
    VELCOR_R_THETA(I,J) = VELCOR_R_THETA(I,J)/RDF_R_THETA(I,J)
    POLARCOR_R_THETA(I,J) = POLARCOR_R_THETA(I,J)/RDF_R_THETA(I,J)
   end if
  end do
  
  if (I>1) then
   CORREL_INTEGRAL_PART(IND_I,I) = CORREL_INTEGRAL_PART(IND_I,I) &
                                 + CORREL_INTEGRAL_PART(IND_I,I-1)                                            
   
   CORREL_INTEGRAL_MEAN_TIME(M,I) = RDF_R(I) 
   CORREL_INTEGRAL_MEAN_TIME(M,I) = CORREL_INTEGRAL_MEAN_TIME(M,I) &
                                    + CORREL_INTEGRAL_MEAN_TIME(M,I-1)  
  else
   CORREL_INTEGRAL_MEAN_TIME(M,I) = RDF_R(I) 
  end if
  
!~   do J = 1, NPSTAT
!~    TEST_MEAN(I) =  TEST_MEAN(I) +  CORREL_INTEGRAL_PART(J,I) 
!~   end do
!~   
!~   TEST_MEAN(I) = TEST_MEAN(I) / NPSTAT
 end do
 
!~  call FIT_LIN_LOG(NR, RANGE_R, CORREL_INTEGRAL_PART(15,:), CORREL_DIM_PART(15))  


 
 CORREL_INTEGRAL_MEAN_TIME(M,1:NR) = CORREL_INTEGRAL_MEAN_TIME(M,1:NR)/real(NPSTAT*(NPSTAT-1))
 
 call FIT_LIN_LOG(int(NR/2), RANGE_R(1:int(NR/2)), CORREL_INTEGRAL_MEAN_TIME(M,1:int(NR/2)), CORREL_DIM_MEAN_TIME(M))
 
!~   print*, 'IND = ', IND
!~  print*, 'CORREL_INTEGRAL_PART(1,1:15) = ', CORREL_INTEGRAL_PART(1,1:15)
!~  read(*,*)
!~  print*, 'TEST_MEAN(1:15) = ', TEST_MEAN(1:15)
!~  read(*,*)
!~  print*, 'CORREL_DIM_PART(15) = ', CORREL_DIM_PART(15)
!~  read(*,*)
!~  print*, 'CORREL_INTEGRAL_MEAN_TIME(IND,:) = ', CORREL_INTEGRAL_MEAN_TIME(IND,:) 
!~   OPEN(UNIT=12, FILE="output.txt", ACTION="write")
!~   DO i=1,NR
!~     WRITE(12,*) CORREL_INTEGRAL_MEAN_TIME(IND,I) 
!~   END DO
!~   close(12)
!~  
!~  print*, 'CORREL_DIM_MEAN_TIME(IND) = ', CORREL_DIM_MEAN_TIME(IND)
!~ read(*,*)
!~  

 

 RDF_R_THETA = RDF_R_THETA/(NPSTAT**(2)/L**(3)*VOL_RING)
 RDF_R = RDF_R/(real(NPSTAT)**(2)/L**(3)*VOL_SHELL)

 
  
 RDF_R_THETA_MEAN = RDF_R_THETA_MEAN + RDF_R_THETA!*(real(NPART)/real(NPSTAT))**2
 RDF_R_MEAN = RDF_R_MEAN + RDF_R!*(real(NPART)/real(NPSTAT))**2
 

 VELCOR_R_THETA_MEAN = VELCOR_R_THETA_MEAN + VELCOR_R_THETA!*(real(NPART)/real(NPSTAT))**2
 VELCOR_R_MEAN = VELCOR_R_MEAN + VELCOR_R!*(real(NPART)/real(NPSTAT))**2

 POLARCOR_R_THETA_MEAN = POLARCOR_R_THETA_MEAN + POLARCOR_R_THETA!*(real(NPART)/real(NPSTAT))**2
 POLARCOR_R_MEAN = POLARCOR_R_MEAN + POLARCOR_R!*(real(NPART)/real(NPSTAT))**2
 
!~ print*, 'RDF_R_THETA(1:10,1) second= ', RDF_R_THETA(1:10,1)
!~ print*, 'RDF_R_THETA_MEAN(1:10,1) = ', RDF_R_THETA_MEAN(1:10,1)
!~ read(*,*)
end do  ! IND = ....

!~ RDF_MEAN = RDF_MEAN/(NSAVES-SAVE_START + 1)
RDF_R_THETA_MEAN  = RDF_R_THETA_MEAN /real(NSAVES-SAVE_START + 1)
RDF_R_MEAN = RDF_R_MEAN/real(NSAVES-SAVE_START + 1)
!~ print*,'maxval(RDF_MEAN ) = ',maxval(RDF_MEAN )
print*,'maxval(abs(RDF_R_THETA_MEAN) ) = ',maxval(abs(RDF_R_THETA_MEAN) )
print*,'maxval(abs(RDF_R_MEAN) ) = ',maxval(abs(RDF_R_MEAN) )

VELCOR_R_THETA_MEAN  = VELCOR_R_THETA_MEAN /real(NSAVES-SAVE_START + 1)
VELCOR_R_MEAN = VELCOR_R_MEAN/real(NSAVES-SAVE_START + 1)
print*,'maxval(abs(VELCOR_R_THETA_MEAN) ) = ',maxval(abs(VELCOR_R_THETA_MEAN) )
print*,'maxval(abs(VELCOR_R_MEAN) ) = ',maxval(abs(VELCOR_R_MEAN) )
  
POLARCOR_R_THETA_MEAN  = POLARCOR_R_THETA_MEAN /real(NSAVES-SAVE_START + 1)
POLARCOR_R_MEAN = POLARCOR_R_MEAN/real(NSAVES-SAVE_START + 1)
print*,'maxval(abs(POLARCOR_R_THETA_MEAN) ) = ',maxval(abs(POLARCOR_R_THETA_MEAN) )
print*,'maxval(abs(POLARCOR_R_MEAN) ) = ',maxval(abs(POLARCOR_R_MEAN) )

!~ 
!=====================================================================
! 3. WRITE RDF'S
!=====================================================================
!~ !- Save only if <100MB
!~ if (NR*NTHETA*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE RDF  '
!~  print*,NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~ 
!~ 
!~  !---------- g(r,theta,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'RDF_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))') NR, NTHETA, NPHI, NR*NTHETA*NPHI
!~  do I = 1, NR
!~   do J = 1, NTHETA
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), RDF_MEAN (I,J,K)
!~    end do
!~   end do
!~  end do
!~ 
!~  !- close file
!~  close(301)
!~  print*,'SAVE RDF -->  OK '
!~ else
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE RDF_R_THETA  '
print*,NR*NTHETA*3.0,' VALUES TO SAVE  '
print*,NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- g(r,theta)-----------
!!-Print filename
write(FILENAME,10200)'RDF_R_THETA_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))') NR, NTHETA, NR*NTHETA
do I = 1, NR
 do J = 1, NTHETA
   write(302,'(3(e17.7))') RANGE_R(I), RANGE_THETA(J), RDF_R_THETA_MEAN(I,J)
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
write(FILENAME,10200)'RDF_R_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1, NR
   write(303,'(2(e17.7))') RANGE_R(I), RDF_R_MEAN(I)
end do

!- close file
close(303)


print*,'SAVE RDF_R -->  OK '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ !- Save only if <100MB
!~ if (NR*NTHETA*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE VELCOR  '
!~ 
!~  print*,NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~  !---------- velcor(r,theta,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'VELCOR_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))')  NR, NTHETA, NPHI, NR*NTHETA*NPHI
!~ 
!~ 
!~  do I = 1, NR
!~   do J = 1, NTHETA
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') &
!~      RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), VELCOR_MEAN(I,J,K)
!~    end do
!~   end do
!~  end do
!~ 
!~ 
!~  !- close file
!~  close(301)
!~ 
!~ 
!~  print*,'SAVE VELCOR -->  OK '
!~ else
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE VELCOR_R_THETA  '
print*,NR*NTHETA*3.0,' VALUES TO SAVE  '
print*,NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- velcor(r,theta)-----------
!!-Print filename
write(FILENAME,10200)'VELCOR_R_THETA_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))')  NR, NTHETA, NR*NTHETA


do I = 1, NR
 do J = 1, NTHETA
   write(302,'(3(e17.7))') &
    RANGE_R(I), RANGE_THETA(J), VELCOR_R_THETA_MEAN(I,J)
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
write(FILENAME,10200)'VELCOR_R_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NR
   write(303,'(2(e17.7))') &
    RANGE_R(I), VELCOR_R_MEAN(I)
end do


!- close file
close(303)


print*,'SAVE VELCOR_R -->  OK '



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Save only if <100MB
!~ if (NR*NTHETA*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE POLARCOR  '
!~ 
!~  print*,NR*NTHETA*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~  !---------- POLARcor(r,theta,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'POLARCOR_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))')  NR, NTHETA, NPHI, NR*NTHETA*NPHI
!~ 
!~ 
!~  do I = 1, NR
!~   do J = 1, NTHETA
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') &
!~      RANGE_R(I), RANGE_THETA(J), RANGE_PHI(K), POLARCOR_MEAN(I,J,K)
!~    end do
!~   end do
!~  end do
!~ 
!~ 
!~  !- close file
!~  close(301)
!~ 
!~ 
!~  print*,'SAVE POLARCOR -->  OK '
!~ else
!~  print*,NR*NTHETA*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE POLARCOR_R_THETA  '
print*,NR*NTHETA*3.0,' VALUES TO SAVE  '
print*,NR*NTHETA*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- POLARcor(r,theta)-----------
!!-Print filename
write(FILENAME,10200)'POLARCOR_R_THETA_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))')  NR, NTHETA, NR*NTHETA


do I = 1, NR
 do J = 1, NTHETA
   write(302,'(3(e17.7))') &
    RANGE_R(I), RANGE_THETA(J), POLARCOR_R_THETA_MEAN(I,J)
 end do
end do


!- close file
close(302)


print*,'SAVE POLARCOR_R_THETA -->  OK '

print*,'SAVE POLARCOR_R  '
print*,NR*2.0,' VALUES TO SAVE  '
print*,NR*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- POLARcor(r)-----------
!!-Print filename
write(FILENAME,10200)'POLARCOR_R_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NR
   write(303,'(2(e17.7))') &
    RANGE_R(I), POLARCOR_R_MEAN(I)
end do


!- close file
close(303)


print*,'SAVE POLARCOR_R -->  OK '


print*,'SAVE CORREL_DIM_MEAN_TIME  '
print*,NR*2.0,' VALUES TO SAVE  '
print*,NR*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- CORREL_DIM_MEAN_TIME(t)-----------
!!-Print filename
write(FILENAME,10200)'CORREL_DIM_MEAN_TIME.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NSAVES - SAVE_START + 1
   write(303,'(2(e17.7))') &
    I + SAVE_START, CORREL_DIM_MEAN_TIME(I)
end do


!- close file
close(303)


print*,'SAVE CORREL_DIM_MEAN_TIME -->  OK '
print*,' '
print*,'-------------------END POSITIONS, VELOCITIES, ORIENT CORREL--------------------------- '
print*,' '

 
deallocate(RANGE_R)
deallocate(RANGE_THETA)
!~ allocate(RANGE_THETA_IRREG(NR,NPHI))
deallocate(RANGE_PHI)

!~ deallocate(RDF)
!~ deallocate(VOL_SECTOR)

deallocate(VOL_RING)
deallocate(VOL_SHELL)

deallocate(RDF_R)
deallocate(RDF_R_THETA)

!~ deallocate(VELCOR)
deallocate(VELCOR_R_THETA)
deallocate(VELCOR_R)

!~ deallocate(POLARCOR)
deallocate(POLARCOR_R_THETA)
deallocate(POLARCOR_R)

deallocate(RDF_R_MEAN)
deallocate(RDF_R_THETA_MEAN)

!~ deallocate(VELCOR)
deallocate(VELCOR_R_THETA_MEAN)
deallocate(VELCOR_R_MEAN)

!~ deallocate(POLARCOR)
deallocate(POLARCOR_R_THETA_MEAN)
deallocate(POLARCOR_R_MEAN)




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

end subroutine POSITION_VELOCITY_ORIENT_CORREL

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine init_random_seed()
 implicit none
 integer, allocatable :: seed(:)
 integer :: i, n, un, istat, dt(8), pid, t(2), s
 integer(8) :: count, tms

 call random_seed(size = n)
 allocate(seed(n))
 ! First try if the OS provides a random number generator
 open(newunit=un, file="/dev/urandom", access="stream", &
      form="unformatted", action="read", status="old", iostat=istat)
 if (istat == 0) then
    read(un) seed
    close(un)
 else
    ! Fallback to XOR:ing the current time and pid. The PID is
    ! useful in case one launches multiple instances of the same
    ! program in parallel.
    call system_clock(count)
    if (count /= 0) then
       t = transfer(count, t)
    else
       call date_and_time(values=dt)
       tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
            + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
            + dt(3) * 24 * 60 * 60 * 60 * 1000 &
            + dt(5) * 60 * 60 * 1000 &
            + dt(6) * 60 * 1000 + dt(7) * 1000 &
            + dt(8)
       t = transfer(tms, t)
    end if
    s = ieor(t(1), t(2))
    pid = 1 + 1099279 ! Add a prime
    s = ieor(s, pid)
    if (n >= 3) then
       seed(1) = t(1) + 36269
       seed(2) = t(2) + 72551
       seed(3) = pid
       if (n > 3) then
          seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
       end if
    else
       seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
    end if
 end if
 call random_seed(put=seed)
end subroutine init_random_seed
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine REMOVE_DUPS(length,input,k,output)
  implicit none
  integer, intent(in) :: length
  integer, dimension(length), intent(in) :: input        ! The input
  integer, intent(out) :: k                   ! The number of unique elements
  integer, dimension(length), intent(out) :: output  ! The output

                  ! The number of unique elements
  integer :: i, j
 
  k = 1
  output(1) = input(1)
  outer: do i=2,length
           do j=1,k
              if (output(j) == input(i)) then
                 ! Found a match so start looking again
                 cycle outer
              end if
           end do
           ! No match found so add it to the output
           k = k + 1
           output(k) = input(i)
         end do outer
 
!~   write(*,advance='no',fmt='(a,i0,a)') 'Unique list has ',k,' elements: '
!~   write(*,*) output(1:k)
end subroutine REMOVE_DUPS

