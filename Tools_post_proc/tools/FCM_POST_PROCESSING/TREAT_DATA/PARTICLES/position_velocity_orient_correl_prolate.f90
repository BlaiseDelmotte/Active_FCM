!!====================================================================
!!
!!
!!====================================================================

subroutine POSITION_VELOCITY_ORIENT_CORREL_PROLATE(NSAVES, &
                                                   SAVE_START, &
                                                   NPART, &
                                                   NPSTAT, &
                                                   L, &
                                                   RADMIN, &
                                                   RADMAX, &
                                                   ROT_MAT, &
                                                   POSI, &
                                                   VEL         )

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
! Starting snapshot to treat
integer, intent(in) :: SAVE_START
! Number of swimmers
integer, intent(in) :: NPART
! Number of swimmers to treat
integer, intent(in) :: NPSTAT
! Domain size
real(kind=8), intent(in):: L
! Minimal particle radius
real(kind=8), intent(in):: RADMIN
! MAximal particle radius
real(kind=8), intent(in):: RADMAX
! Particle positions, orientation and velocities
real(kind=8), dimension(NSAVES,NPART,3,3), intent(in) :: ROT_MAT
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: POSI
real(kind=8), dimension(NSAVES,NPART,3), intent(in) :: VEL


!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- r-component
real(kind=8), allocatable, dimension(:) :: RANGE_SIGMA
!- TAU-component
real(kind=8), allocatable, dimension(:) :: RANGE_TAU
!~ real(kind=8), allocatable, dimension(:,:) :: RANGE_TAU_IRREG
!~ !- phi-component
!~ real(kind=8), allocatable, dimension(:) :: RANGE_PHI

!! Per time step values
!~ !- g(r,TAU,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: RDF
!- g(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: RDF_SIGMA_TAU
!- g(r)
real(kind=8), allocatable, dimension(:) :: RDF_SIGMA
!~ !- Iu(vel,r,TAU,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: VELCOR
!- Iu(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: VELCOR_SIGMA_TAU
!- Iu(r)
real(kind=8), allocatable, dimension(:) :: VELCOR_SIGMA
!~ !- Ie(vel,r,TAU,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: POLARCOR
!- Ie(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: POLARCOR_SIGMA_TAU
!- Ie(r)
real(kind=8), allocatable, dimension(:) :: POLARCOR_SIGMA

!!Mean values
!~ !- g(r,TAU,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: RDF_MEAN
!- g(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: RDF_SIGMA_TAU_MEAN
!- g(r)
real(kind=8), allocatable, dimension(:) :: RDF_SIGMA_MEAN
!~ !- Iu(vel,r,TAU,phi) mean
!~ real(kind=8), allocatable, dimension(:,:,:) :: VELCOR_MEAN
!- Iu(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: VELCOR_SIGMA_TAU_MEAN
!- Iu(r)
real(kind=8), allocatable, dimension(:) :: VELCOR_SIGMA_MEAN
!~ !- Ie(vel,r,TAU,phi)
!~ real(kind=8), allocatable, dimension(:,:,:) :: POLARCOR_MEAN
!- Ie(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: POLARCOR_SIGMA_TAU_MEAN
!- Ie(r)
real(kind=8), allocatable, dimension(:) :: POLARCOR_SIGMA_MEAN

!~ !-Vol discretized sector
!~ real(kind=8), allocatable, dimension(:,:,:) :: VOL_SECTOR
!-Vol of ring for discretized sector for g(r,TAU)
real(kind=8), allocatable, dimension(:,:) :: VOL_RING
!-Vol of shell for discretized sector for g(r)
real(kind=8), allocatable, dimension(:) :: VOL_SHELL

!- Vector to select randomly NPSTAT particles
real(kind=8), dimension(5*NPART) :: VEC_RAND
integer, dimension(5*NPART) :: IND_RAND
integer, dimension(5*NPART) :: IND_RAND_UNIQUE

!- Step size to discretize the interval of positions and vel corr
real(kind=8) ::  DSIGMA, DTAU, DPHI
!~ real(kind=8), allocatable, dimension(:) :: DTAU_IRREG


real(kind=8) ::  PPI

!- Half domain size
real(kind=8) ::  L_2

!-2*RADMIN
real(kind=8) ::  TWORADMIN

!-Aspect Ration
real(kind=8) ::  AR


!-Maximal distance accepted for RDF computation
real(kind=8) ::  MAX_DIST

real(kind=8) :: FOCAL_DIST

real(kind=8) :: SIGMA_MIN, SIGMA_MAX
real(kind=8) :: TAU_MIN, TAU_MAX

real(kind=8) :: RADMAX_LOC, RADMAX_LOC_NEXT
real(kind=8) :: RADMIN_LOC, RADMIN_LOC_NEXT


!- Pos, Distance,...
real(kind=8), dimension(3) :: POS_I, POS_J 
real(kind=8), dimension(3) :: POS_IJ, POS_IJ_I, POS_IJ_J
real(kind=8) :: RIJ
real(kind=8) :: TEMP_IJ_I, TEMP2_IJ_I
real(kind=8) :: TEMP_IJ_J, TEMP2_IJ_J
real(kind=8), dimension(3,3) :: ROT_I, ROT_J
real(kind=8) :: SIGMA_I, SIGMA_J
real(kind=8) :: TAU_I, TAU_J

real(kind=8) :: SCAL_VEL_IJ
real(kind=8) :: VEL_I_SQ
real(kind=8) :: VELCOR_IJ

real(kind=8) :: POLARCOR_IJ

!- Number of steps to discretize the interval of pos and velocities
integer :: NSIGMA, NTAU!, NPHI


integer :: NUNIQUE

!- Indices for corr functions
integer :: IND_VELCOR, IND_SIGMA, IND_SIGMA_I, IND_SIGMA_J, IND_TAU_I, IND_TAU_J, IND_PHI_I, IND_PHI_J

!- File name 
character(len=100) :: FILENAME

integer :: alloc_status

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

NSIGMA = 250
NTAU = 100
!~ NPHI = 100

L_2 = L/2.0 
TWORADMIN = 2.0*RADMIN
MAX_DIST = L_2*dsqrt(3.0d0)
print*,'MAX_DIST/RADMAX = ', MAX_DIST/RADMAX



AR = RADMAX/RADMIN

FOCAL_DIST = dsqrt(RADMAX**2 - RADMIN**2)

SIGMA_MIN = dcosh(datanh(1.0d0/AR))
SIGMA_MAX = 1.01*MAX_DIST/FOCAL_DIST

TAU_MIN = -0.99999999
TAU_MAX = 0.99999999

DSIGMA = (SIGMA_MAX - SIGMA_MIN)/real(NSIGMA)
DTAU = ( TAU_MAX-TAU_MIN )/real(NTAU-1)
!~ DPHI = (2.0*PPI)/real(NPHI)

!~ NSIGMA = ceiling((SIGMA_MAX - SIGMA_MIN)/DSIGMA)
!~ NTAU = ceiling((TAU_MAX-TAU_MIN )/DTAU)
!~ NPHI = ceiling((2.0*PPI)/DPHI)

print*,'NSIGMA, NTAU = '
print*, NSIGMA, NTAU!, NPHI
     

allocate(RANGE_SIGMA(NSIGMA))
allocate(RANGE_TAU(NTAU))
!~ allocate(RANGE_TAU_IRREG(NSIGMA,NPHI))
!~ allocate(RANGE_PHI(NPHI))

!~ allocate(VOL_SECTOR(NSIGMA,NTAU,NPHI))
allocate(VOL_RING(NSIGMA,NTAU))
allocate(VOL_SHELL(NSIGMA))

!~ allocate(RDF(NSIGMA,NTAU,NPHI))
allocate(RDF_SIGMA_TAU(NSIGMA,NTAU))
allocate(RDF_SIGMA(NSIGMA))


!~ allocate(VELCOR(NSIGMA,NTAU,NPHI))
allocate(VELCOR_SIGMA_TAU(NSIGMA,NTAU))
allocate(VELCOR_SIGMA(NSIGMA))

!~ allocate(POLARCOR(NSIGMA,NTAU,NPHI))
allocate(POLARCOR_SIGMA_TAU(NSIGMA,NTAU))
allocate(POLARCOR_SIGMA(NSIGMA))

!~ allocate(RDF_MEAN(NSIGMA,NTAU,NPHI))
allocate(RDF_SIGMA_TAU_MEAN(NSIGMA,NTAU))
allocate(RDF_SIGMA_MEAN(NSIGMA))

!~ allocate(VELCOR_MEAN(NSIGMA,NTAU,NPHI))
allocate(VELCOR_SIGMA_TAU_MEAN(NSIGMA,NTAU))
allocate(VELCOR_SIGMA_MEAN(NSIGMA))

!~ allocate(POLARCOR_MEAN(NSIGMA,NTAU,NPHI))
allocate(POLARCOR_SIGMA_TAU_MEAN(NSIGMA,NTAU))
allocate(POLARCOR_SIGMA_MEAN(NSIGMA))



RANGE_SIGMA(1) = SIGMA_MIN + DSIGMA
do J = 2, NSIGMA
 RANGE_SIGMA(J) = RANGE_SIGMA(J-1) + DSIGMA
end do

RANGE_TAU(1) = TAU_MIN 
do J = 2, NTAU
 RANGE_TAU(J) = RANGE_TAU(J-1) + DTAU
end do

!~ RANGE_PHI(1) = -PPI + DPHI 
!~ do J = 2, NPHI
!~  RANGE_PHI(J) = RANGE_PHI(J-1) + DPHI
!~ end do


VOL_SHELL = 0.0
VOL_RING = 0.0
do I = 1, NSIGMA
 
 RADMIN_LOC = FOCAL_DIST*dsqrt((RANGE_SIGMA(I)-DSIGMA)**2-1.0) 
 RADMAX_LOC = FOCAL_DIST*(RANGE_SIGMA(I)-DSIGMA)
  
 RADMIN_LOC_NEXT = FOCAL_DIST*dsqrt(RANGE_SIGMA(I)**2-1.0)
 RADMAX_LOC_NEXT = FOCAL_DIST*RANGE_SIGMA(I)

 
 VOL_SHELL(I) = 4.0/3.0*PPI*( RADMAX_LOC_NEXT*RADMIN_LOC_NEXT**2 &
                            - RADMAX_LOC*RADMIN_LOC**2           )
                                 
 do J = 1, NTAU            
   VOL_RING(I, J) =  &
       FOCAL_DIST**3*(RANGE_SIGMA(I)**2-RANGE_TAU(J)**2)*DTAU*DSIGMA *2.0*PPI
 end do
end do


print*, 'DISCRETIZE RANGE OF POSITIONS---> OK'

!=====================================================================
! 2. COMPUTE PDF
!=====================================================================
print*,'COMPUTE RDF AND VELCOR AND POLARCOR'

!~ RDF_MEAN = 0.0
RDF_SIGMA_TAU_MEAN = 0.0
RDF_SIGMA_MEAN = 0.0

!~ VELCOR_MEAN = 0.0
VELCOR_SIGMA_TAU_MEAN = 0.0
VELCOR_SIGMA_MEAN = 0.0

!~ POLARCOR_MEAN = 0.0
POLARCOR_SIGMA_TAU_MEAN = 0.0
POLARCOR_SIGMA_MEAN = 0.0



if (NPSTAT==NPART) then
 do I = 1, NPART
  IND_RAND_UNIQUE(I) = I
 end do
end if



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
 RDF_SIGMA_TAU = 0.0
 RDF_SIGMA = 0.0
 
!~  VELCOR = 0.0
 VELCOR_SIGMA_TAU = 0.0
 VELCOR_SIGMA = 0.0

!~  POLARCOR = 0.0
 POLARCOR_SIGMA_TAU = 0.0
 POLARCOR_SIGMA = 0.0

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
   POS_IJ = (POS_IJ - L* real(int(POS_IJ/(L_2))))  
   
   RIJ = dsqrt(POS_IJ(1)**2 + POS_IJ(2)**2 + POS_IJ(3)**2)
   
   if ((RIJ>=TWORADMIN).and.(RIJ<MAX_DIST)) then
   
    POS_IJ_I = matmul(ROT_I,-POS_IJ)  
    
    TEMP_IJ_I = dsqrt( POS_IJ_I(1)**2 &
                     + POS_IJ_I(2)**2 &
                     + (POS_IJ_I(3)+FOCAL_DIST)**2 )
                   
    TEMP2_IJ_I = dsqrt( POS_IJ_I(1)**2 &
                      + POS_IJ_I(2)**2 &
                      + (POS_IJ_I(3)-FOCAL_DIST)**2 )
                    
    SIGMA_I = 1.0/(2.0*FOCAL_DIST)* (TEMP_IJ_I + TEMP2_IJ_I)
    TAU_I = 1.0/(2.0*FOCAL_DIST)* (TEMP_IJ_I - TEMP2_IJ_I)  
      
   
!~     IND_SIGMA_I = ceiling( ( SIGMA_I-RANGE_SIGMA(1) )/DSIGMA ) + 1
    IND_TAU_I = ceiling( ( TAU_I-RANGE_TAU(1) )/DTAU ) !+ 1
    
        
    POS_IJ_J = matmul(ROT_J,POS_IJ)

    TEMP_IJ_J = dsqrt( POS_IJ_J(1)**2 &
                     + POS_IJ_J(2)**2 &
                     + (POS_IJ_J(3)+FOCAL_DIST)**2 )
                   
    TEMP2_IJ_J = dsqrt( POS_IJ_J(1)**2 &
                      + POS_IJ_J(2)**2 &
                      + (POS_IJ_J(3)-FOCAL_DIST)**2 )
                    
    SIGMA_J = 1.0/(2.0*FOCAL_DIST)* (TEMP_IJ_J + TEMP2_IJ_J)
    

    
    TAU_J = 1.0/(2.0*FOCAL_DIST)* (TEMP_IJ_I - TEMP2_IJ_J)  
    
    
!~       print*,'I = ', I
!~     print*,'J = ', J
!~     print*,'RIJ = ', RIJ
!~     print*,'POS_IJ_I = ', POS_IJ_I
!~     print*,'NORM(POS_IJ_I) = ', dsqrt(POS_IJ_I(1)**2 + POS_IJ_I(2)**2 + POS_IJ_I(3)**2 )
!~     print*,'POS_IJ_J = ', POS_IJ_J
!~     print*,'NORM(POS_IJ_!j) = ', dsqrt(POS_IJ_J(1)**2 + POS_IJ_J(2)**2 + POS_IJ_J(3)**2 )
!~     print*,'SIGMA_I = ', SIGMA_I
!~     print*,'SIGMA_J = ', SIGMA_J
!~     print*,'TAU_I = ', TAU_I
!~     print*,'TAU_J = ', TAU_J
!~     read(*,*)
   
    
!~     IND_SIGMA_J = ceiling( ( SIGMA_J-RANGE_SIGMA(1) )/DSIGMA ) + 1
    IND_TAU_J = ceiling( ( TAU_J-RANGE_TAU(1) )/DTAU ) !+ 1
    
    SCAL_VEL_IJ = VEL(IND,I,1)*VEL(IND,J,1) &
                + VEL(IND,I,2)*VEL(IND,J,2) &
                + VEL(IND,I,3)*VEL(IND,J,3) 
                
    VELCOR_IJ = SCAL_VEL_IJ/VEL_I_SQ
    
!~     IND_VELCOR = int(VELCOR_IJ/DVELCOR)

    POLARCOR_IJ = ROT_I(3,1)*ROT_J(3,1) &
                + ROT_I(3,2)*ROT_J(3,2) &
                + ROT_I(3,3)*ROT_J(3,3)
      
    
!~     if (IND_PHI_I>NPHI) then
!~      print*,'!!!!   ERROR -->  IND_PHI_I> NPHI  !!!!!!!!! '
!~      print*,'IND_PHI_I = ' ,IND_PHI_I 
!~      print*,'PHI_I = ' ,PHI_I 
!~      print*,'maxval(RANGE_PHI) = ' ,maxval(RANGE_PHI)
!~      read(*,*)
!~     end if
    
    ! As sigma_i and sigma_j should be equal , we compute the mean
    IND_SIGMA = ceiling( (min(SIGMA_J , SIGMA_I)-RANGE_SIGMA(1))/DSIGMA ) + 1
    
    if (IND_SIGMA>NSIGMA) then
     print*,'!!!!   ERROR -->  IND_SIGMA> NSIGMA  !!!!!!!!! '
     print*,'IND_SIGMA = ' ,IND_SIGMA 
     print*,'SIGMA_I = ' ,SIGMA_I
     print*,'SIGMA_J = ' ,SIGMA_J
     print*,'maxval(RANGE_SIGMA) = ' ,maxval(RANGE_SIGMA)
     
     print*,'RIJ = ', RIJ
     print*,'MAX_DIST = ', MAX_DIST 

     read(*,*)
    end if
    
    if (IND_TAU_I>NTAU) then
     if ( ((TAU_I-maxval(RANGE_TAU))/maxval(RANGE_TAU))>0.1 ) then
      print*,'!!!!   ERROR -->  IND_TAU_I> NTAU  !!!!!!!!! '
      print*,'IND_TAU_I = ' ,IND_TAU_I 
      print*,'TAU_I = ' ,TAU_I
      print*,'maxval(RANGE_TAU) = ' ,maxval(RANGE_TAU)
      IND_TAU_I  = NTAU
      read(*,*)
     else
      print*,'((TAU_I-maxval(RANGE_TAU))/maxval(RANGE_TAU)) = ', &
              ((TAU_I-maxval(RANGE_TAU))/maxval(RANGE_TAU))
      print*,'IND_TAU_I = ' ,IND_TAU_I 
      print*,'TAU_I = ' ,TAU_I
      print*,'maxval(RANGE_TAU) = ' ,maxval(RANGE_TAU)
      IND_TAU_I = NTAU
     end if
    end if
    
    if (IND_TAU_J>NTAU) then
     if ( ((TAU_J-maxval(RANGE_TAU))/maxval(RANGE_TAU))>0.1 ) then
      print*,'!!!!   ERROR -->  IND_TAU_J> NTAU  !!!!!!!!! '
      print*,'IND_TAU_J = ' ,IND_TAU_J 
      print*,'TAU_J = ' ,TAU_J
      print*,'maxval(RANGE_TAU) = ' ,maxval(RANGE_TAU)
      IND_TAU_J  = NTAU
      read(*,*)
     else
      print*,'((TAU_J-maxval(RANGE_TAU))/maxval(RANGE_TAU)) = ', &
              ((TAU_J-maxval(RANGE_TAU))/maxval(RANGE_TAU))
      print*,'IND_TAU_J = ' ,IND_TAU_J 
      print*,'TAU_J = ' ,TAU_J
      print*,'maxval(RANGE_TAU) = ' ,maxval(RANGE_TAU)
      IND_TAU_J = NTAU
     end if
    end if

    
    RDF_SIGMA_TAU(IND_SIGMA, IND_TAU_I) = RDF_SIGMA_TAU(IND_SIGMA, IND_TAU_I) + 1.0
    RDF_SIGMA_TAU(IND_SIGMA, IND_TAU_J) = RDF_SIGMA_TAU(IND_SIGMA, IND_TAU_J) + 1.0
    
    RDF_SIGMA(IND_SIGMA) = RDF_SIGMA(IND_SIGMA) + 2.0
!~     RDF_SIGMA(IND_SIGMA) = RDF_SIGMA(IND_SIGMA) + 1    
     
        

    VELCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_I) = &
        VELCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_I) + VELCOR_IJ       
    VELCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_J) = &
        VELCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_J) + VELCOR_IJ 
       
    VELCOR_SIGMA(IND_SIGMA) = &
        VELCOR_SIGMA(IND_SIGMA) + 2.0*VELCOR_IJ 

!~     VELCOR_SIGMA(IND_SIGMA) = &
!~         VELCOR_SIGMA(IND_SIGMA) + VELCOR_IJ
         

    POLARCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_I) = &
       POLARCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_I) + POLARCOR_IJ         
    POLARCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_J) = &
        POLARCOR_SIGMA_TAU(IND_SIGMA,IND_TAU_J) + POLARCOR_IJ 
       
    POLARCOR_SIGMA(IND_SIGMA) = &
        POLARCOR_SIGMA(IND_SIGMA) + 2.0*POLARCOR_IJ  
 
!~     POLARCOR_SIGMA(IND_SIGMA) = &
!~         POLARCOR_SIGMA(IND_SIGMA) + POLARCOR_IJ
                
   end if

   
  end do  
 end do
 

 
 do I = 1,NSIGMA
   if (RDF_SIGMA(I)>0) then
   VELCOR_SIGMA(I) = VELCOR_SIGMA(I)/RDF_SIGMA(I)!*(real(NPART)/real(NPSTAT))
   POLARCOR_SIGMA(I) = POLARCOR_SIGMA(I)/RDF_SIGMA(I)!*(real(NPART)/real(NPSTAT)) 
  end if  
  do J = 1, NTAU
   if (RDF_SIGMA_TAU(I,J)>0) then
    VELCOR_SIGMA_TAU(I,J) = VELCOR_SIGMA_TAU(I,J)/RDF_SIGMA_TAU(I,J)
    POLARCOR_SIGMA_TAU(I,J) = POLARCOR_SIGMA_TAU(I,J)/RDF_SIGMA_TAU(I,J)
   end if
  end do
 end do
 
 
 RDF_SIGMA_TAU = RDF_SIGMA_TAU/(NPSTAT**(2)/L**(3)*VOL_RING)
 RDF_SIGMA = RDF_SIGMA/(NPSTAT**(2)/L**(3)*VOL_SHELL)

 
  
 RDF_SIGMA_TAU_MEAN = RDF_SIGMA_TAU_MEAN + RDF_SIGMA_TAU!*(real(NPART)/real(NPSTAT))**2
 RDF_SIGMA_MEAN = RDF_SIGMA_MEAN + RDF_SIGMA!*(real(NPART)/real(NPSTAT))**2
 

 VELCOR_SIGMA_TAU_MEAN = VELCOR_SIGMA_TAU_MEAN + VELCOR_SIGMA_TAU!*(real(NPART)/real(NPSTAT))**2
 VELCOR_SIGMA_MEAN = VELCOR_SIGMA_MEAN + VELCOR_SIGMA!*(real(NPART)/real(NPSTAT))**2

 POLARCOR_SIGMA_TAU_MEAN = POLARCOR_SIGMA_TAU_MEAN + POLARCOR_SIGMA_TAU!*(real(NPART)/real(NPSTAT))**2
 POLARCOR_SIGMA_MEAN = POLARCOR_SIGMA_MEAN + POLARCOR_SIGMA!*(real(NPART)/real(NPSTAT))**2
 

 
end do

!~ RDF_MEAN = RDF_MEAN/(NSAVES-SAVE_START + 1)
RDF_SIGMA_TAU_MEAN  = RDF_SIGMA_TAU_MEAN /(NSAVES-SAVE_START + 1)
RDF_SIGMA_MEAN = RDF_SIGMA_MEAN/(NSAVES-SAVE_START + 1)
!~ print*,'maxval(RDF_MEAN ) = ',maxval(RDF_MEAN )
print*,'maxval(abs(RDF_SIGMA_TAU_MEAN) ) = ',maxval(abs(RDF_SIGMA_TAU_MEAN) )
print*,'maxval(abs(RDF_SIGMA_MEAN) ) = ',maxval(abs(RDF_SIGMA_MEAN) )

VELCOR_SIGMA_TAU_MEAN  = VELCOR_SIGMA_TAU_MEAN /(NSAVES-SAVE_START + 1)
VELCOR_SIGMA_MEAN = VELCOR_SIGMA_MEAN/(NSAVES-SAVE_START + 1)
print*,'maxval(abs(VELCOR_SIGMA_TAU_MEAN) ) = ',maxval(abs(VELCOR_SIGMA_TAU_MEAN) )
print*,'maxval(abs(VELCOR_SIGMA_MEAN) ) = ',maxval(abs(VELCOR_SIGMA_MEAN) )
  
POLARCOR_SIGMA_TAU_MEAN  = POLARCOR_SIGMA_TAU_MEAN /(NSAVES-SAVE_START + 1)
POLARCOR_SIGMA_MEAN = POLARCOR_SIGMA_MEAN/(NSAVES-SAVE_START + 1)
print*,'maxval(abs(POLARCOR_SIGMA_TAU_MEAN) ) = ',maxval(abs(POLARCOR_SIGMA_TAU_MEAN) )
print*,'maxval(abs(POLARCOR_SIGMA_MEAN) ) = ',maxval(abs(POLARCOR_SIGMA_MEAN) )

!~ 
!=====================================================================
! 3. WRITE RDF'S
!=====================================================================
!~ !- Save only if <100MB
!~ if (NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE RDF  '
!~  print*,NSIGMA*NTAU*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~ 
!~ 
!~  !---------- g(r,TAU,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'RDF_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))') NSIGMA, NTAU, NPHI, NSIGMA*NTAU*NPHI
!~  do I = 1, NSIGMA
!~   do J = 1, NTAU
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') RANGE_SIGMA(I), RANGE_TAU(J), RANGE_PHI(K), RDF_MEAN (I,J,K)
!~    end do
!~   end do
!~  end do
!~ 
!~  !- close file
!~  close(301)
!~  print*,'SAVE RDF -->  OK '
!~ else
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE RDF_SIGMA_TAU  '
print*,NSIGMA*NTAU*3.0,' VALUES TO SAVE  '
print*,NSIGMA*NTAU*3.0*16.0/10.0**6,' MB TO STORE  '

!---------- g(r,TAU)-----------
!!-Print filename
write(FILENAME,10200)'RDF_SIGMA_TAU_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))') NSIGMA, NTAU, NSIGMA*NTAU

do I = 1, NSIGMA
 do J = 1, NTAU
   write(302,'(3(e17.7))') RANGE_SIGMA(I), RANGE_TAU(J), RDF_SIGMA_TAU_MEAN(I,J)
 end do
end do

!- close file
close(302)


print*,'SAVE RDF_SIGMA_TAU -->  OK '

print*,'SAVE RDF_SIGMA  '
print*,NSIGMA*2.0,' VALUES TO SAVE  '
print*,NSIGMA*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- g(r,TAU,phi)-----------
!!-Print filename
write(FILENAME,10200)'RDF_SIGMA_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))
!!- Ecriture ASCII

do I = 1, NSIGMA
   write(303,'(2(e17.7))') RANGE_SIGMA(I), RDF_SIGMA_MEAN(I)
end do

!- close file
close(303)


print*,'SAVE RDF_SIGMA -->  OK '

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!~ !- Save only if <100MB
!~ if (NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE VELCOR  '
!~ 
!~  print*,NSIGMA*NTAU*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~  !---------- velcor(r,TAU,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'VELCOR_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))')  NSIGMA, NTAU, NPHI, NSIGMA*NTAU*NPHI
!~ 
!~ 
!~  do I = 1, NSIGMA
!~   do J = 1, NTAU
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') &
!~      RANGE_SIGMA(I), RANGE_TAU(J), RANGE_PHI(K), VELCOR_MEAN(I,J,K)
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
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE VELCOR_SIGMA_TAU  '
print*,NSIGMA*NTAU*3.0,' VALUES TO SAVE  '
print*,NSIGMA*NTAU*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- velcor(r,TAU)-----------
!!-Print filename
write(FILENAME,10200)'VELCOR_SIGMA_TAU_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))')  NSIGMA, NTAU, NSIGMA*NTAU


do I = 1, NSIGMA
 do J = 1, NTAU
   write(302,'(3(e17.7))') &
    RANGE_SIGMA(I), RANGE_TAU(J), VELCOR_SIGMA_TAU_MEAN(I,J)
 end do
end do


!- close file
close(302)


print*,'SAVE VELCOR_SIGMA_TAU -->  OK '

print*,'SAVE VELCOR_SIGMA  '
print*,NSIGMA*2.0,' VALUES TO SAVE  '
print*,NSIGMA*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- velcor(r)-----------
!!-Print filename
write(FILENAME,10200)'VELCOR_SIGMA_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NSIGMA
   write(303,'(2(e17.7))') &
    RANGE_SIGMA(I), VELCOR_SIGMA_MEAN(I)
end do


!- close file
close(303)


print*,'SAVE VELCOR_SIGMA -->  OK '



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!- Save only if <100MB
!~ if (NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6<100.0) then
!~  print*,'SAVE POLARCOR  '
!~ 
!~  print*,NSIGMA*NTAU*NPHI*4.0,' VALUES TO SAVE  '
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  '
!~  !---------- POLARcor(r,TAU,phi)-----------
!~  !!-Print filename
!~  write(FILENAME,10200)'POLARCOR_MEAN.dat'
!~ 
!~  !!- ASCII
!~  open(unit=301,file=trim(FILENAME))
!~  !!- Ecriture ASCII
!~  write(301,'(4(i17.7))')  NSIGMA, NTAU, NPHI, NSIGMA*NTAU*NPHI
!~ 
!~ 
!~  do I = 1, NSIGMA
!~   do J = 1, NTAU
!~    do K = 1, NPHI
!~     write(301,'(4(e17.7))') &
!~      RANGE_SIGMA(I), RANGE_TAU(J), RANGE_PHI(K), POLARCOR_MEAN(I,J,K)
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
!~  print*,NSIGMA*NTAU*NPHI*4.0*16.0/10.0**6,' MB TO STORE  , TOO BIG!!!'
!~ end if

print*,'SAVE POLARCOR_SIGMA_TAU  '
print*,NSIGMA*NTAU*3.0,' VALUES TO SAVE  '
print*,NSIGMA*NTAU*3.0*16.0/10.0**6,' MB TO STORE  '
!---------- POLARcor(r,TAU)-----------
!!-Print filename
write(FILENAME,10200)'POLARCOR_SIGMA_TAU_MEAN.dat'

!!- ASCII
open(unit=302,file=trim(FILENAME))
!!- Ecriture ASCII
write(302,'(3(i17.7))')  NSIGMA, NTAU, NSIGMA*NTAU


do I = 1, NSIGMA
 do J = 1, NTAU
   write(302,'(3(e17.7))') &
    RANGE_SIGMA(I), RANGE_TAU(J), POLARCOR_SIGMA_TAU_MEAN(I,J)
 end do
end do


!- close file
close(302)


print*,'SAVE POLARCOR_SIGMA_TAU -->  OK '

print*,'SAVE POLARCOR_SIGMA  '
print*,NSIGMA*2.0,' VALUES TO SAVE  '
print*,NSIGMA*2.0*16.0/10.0**6,' MB TO STORE  '
!---------- POLARcor(r)-----------
!!-Print filename
write(FILENAME,10200)'POLARCOR_SIGMA_MEAN.dat'

!!- ASCII
open(unit=303,file=trim(FILENAME))

do I = 1, NSIGMA
   write(303,'(2(e17.7))') &
    RANGE_SIGMA(I), POLARCOR_SIGMA_MEAN(I)
end do


!- close file
close(303)


print*,'SAVE POLARCOR_SIGMA -->  OK '

print*,' '
print*,'-------------------END POSITIONS, VELOCITIES, ORIENT CORREL--------------------------- '
print*,' '

 print*,'Hello846'
deallocate(RANGE_SIGMA, stat=alloc_status)
    if ( alloc_status /= 0 ) then
     print *, 'deallocate failed.  status = ', alloc_status
     stop
    endif
 print*,'Hello848'
deallocate(RANGE_TAU)
 print*,'Hello850'
!~ !allocate(RANGE_TAU_IRREG(NSIGMA,NPHI))
!~ !print*,'RANGE_PHI = ', RANGE_PHI
!~ !deallocate(RANGE_PHI)

!~ !deallocate(RDF)
!~ !deallocate(VOL_SECTOR)

 print*,'Hello855'
 
deallocate(VOL_RING)
deallocate(VOL_SHELL)

deallocate(RDF_SIGMA)
deallocate(RDF_SIGMA_TAU)

 print*,'Hello863'
!~ !deallocate(VELCOR)
deallocate(VELCOR_SIGMA_TAU)
deallocate(VELCOR_SIGMA)

!~ !deallocate(POLARCOR)
deallocate(POLARCOR_SIGMA_TAU)
deallocate(POLARCOR_SIGMA)

 print*,'Hello872'
 
deallocate(RDF_SIGMA_MEAN)
deallocate(RDF_SIGMA_TAU_MEAN)

!~ !deallocate(VELCOR)
deallocate(VELCOR_SIGMA_TAU_MEAN)
deallocate(VELCOR_SIGMA_MEAN)

!~ !deallocate(POLARCOR)
deallocate(POLARCOR_SIGMA_TAU_MEAN)
deallocate(POLARCOR_SIGMA_MEAN)


print*,'Hello881'

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

end subroutine POSITION_VELOCITY_ORIENT_CORREL_PROLATE
!~ 
!~ !-----------------------------------------------------------------------
!~ !-----------------------------------------------------------------------
!~ subroutine init_random_seed()
!~  implicit none
!~  integer, allocatable :: seed(:)
!~  integer :: i, n, un, istat, dt(8), pid, t(2), s
!~  integer(8) :: count, tms
!~ 
!~  call random_seed(size = n)
!~  allocate(seed(n))
!~  ! First try if the OS provides a random number generator
!~  open(newunit=un, file="/dev/urandom", access="stream", &
!~       form="unformatted", action="read", status="old", iostat=istat)
!~  if (istat == 0) then
!~     read(un) seed
!~     close(un)
!~  else
!~     ! Fallback to XOR:ing the current time and pid. The PID is
!~     ! useful in case one launches multiple instances of the same
!~     ! program in parallel.
!~     call system_clock(count)
!~     if (count /= 0) then
!~        t = transfer(count, t)
!~     else
!~        call date_and_time(values=dt)
!~        tms = (dt(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 &
!~             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
!~             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
!~             + dt(5) * 60 * 60 * 1000 &
!~             + dt(6) * 60 * 1000 + dt(7) * 1000 &
!~             + dt(8)
!~        t = transfer(tms, t)
!~     end if
!~     s = ieor(t(1), t(2))
!~     pid = 1 + 1099279 ! Add a prime
!~     s = ieor(s, pid)
!~     if (n >= 3) then
!~        seed(1) = t(1) + 36269
!~        seed(2) = t(2) + 72551
!~        seed(3) = pid
!~        if (n > 3) then
!~           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
!~        end if
!~     else
!~        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
!~     end if
!~  end if
!~  call random_seed(put=seed)
!~ end subroutine init_random_seed
!~ !-----------------------------------------------------------------------
!~ !-----------------------------------------------------------------------
!~ subroutine REMOVE_DUPS(length,input,k,output)
!~   implicit none
!~   integer, intent(in) :: length
!~   integer, dimension(length), intent(in) :: input        ! The input
!~   integer, intent(out) :: k                   ! The number of unique elements
!~   integer, dimension(length), intent(out) :: output  ! The output
!~ 
!~                   ! The number of unique elements
!~   integer :: i, j
!~  
!~   k = 1
!~   output(1) = input(1)
!~   outer: do i=2,length
!~            do j=1,k
!~               if (output(j) == input(i)) then
!~                  ! Found a match so start looking again
!~                  cycle outer
!~               end if
!~            end do
!~            ! No match found so add it to the output
!~            k = k + 1
!~            output(k) = input(i)
!~          end do outer
!~  
 !write(*,advance='no',fmt='(a,i0,a)') 'Unique list has ',k,' elements: '
 ! write(*,*) output(1:k)
!~ end subroutine REMOVE_DUPS
!~ 
!~ 
!~ !-----------------------------------------------------------------------
