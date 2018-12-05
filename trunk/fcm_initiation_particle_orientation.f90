
 !!====================================================================
!!
!! 
!!> @brief
!!> Routine initiating particle orientations with different ways
!!
!! Date :  03/04/2013
!!
!!
!!> @author 
!!> Blaise Delmotte
!!====================================================================

subroutine FCM_INITIATION_PARTICLE_ORIENTATION

!!====================================================================
!! Here the particle positions are initiated according
!! to the variable FCM_INIT_PART_ORIENT specified by user in the
!! parameter file: 'fcm_param.in'.
!!====================================================================
!! Particle orientation initiation: 
!!------------------------------
!! FCM_INIT_PART_ORIENT=
!!
!!  0: The particle orientation are specified in the Fortran file.
!!     The default orientation correspond to TO WRITE
!!
!!  1: Random particle orientation distribution in the box.
!!
!!  2: Particle orientation are read from the binary file 'FCM_ORIENT_PART.ini'
!!
!!  3: Particle oriented according to a routine  ?????? TO IMPLEMENT.
!!
!!====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE 
use FCM_PART_VARIABLE
use FCM_FORCING_VARIABLE
use MPI_STRUCTURES


implicit none

! Process information variables
integer :: IERROR

!- Integers to signalize a non-existing file
integer :: ERR_FILE_ORIENT

!- Local variables
real(kind=8), dimension(4) :: PPPP
integer :: IB, IP

!- File name 
character(len=40) :: FILENAME

!- To delete
real(kind=8) ::  theta_0


! Prescribed orientations in x-direction (NPART_FULL<=4)
if (FCM_INIT_PART_ORIENT ==0) then

 if (NPART_FULL ==1) then
  
  !- Rotation angle of local body frame
  PPPP(1) =0.0*PPI/3.0  !0.0*PPI/1.0
  
  !- Rotation axis of local body frame
  PPPP(2) = 0.0
  PPPP(3) = 1.0
  PPPP(4) = 0.0
  
  !- Quaternion from Euler parameters
  FCM_QUAT(NPART_FULL,1) = cos(PPPP(1)/2.0)
  FCM_QUAT(NPART_FULL,2:4) = PPPP(2:4)*sin(PPPP(1)/2.0)
  
  !- Quaternion normalization
  FCM_QUATNORM = dsqrt( FCM_QUAT(NPART_FULL,1)**2 + & 
                        FCM_QUAT(NPART_FULL,2)**2 + & 
                        FCM_QUAT(NPART_FULL,3)**2 + & 
                        FCM_QUAT(NPART_FULL,4)**2 )
                        
  FCM_QUAT(NPART_FULL,:) = FCM_QUAT(NPART_FULL,:)/FCM_QUATNORM
  
  !- Swimming Direction vector
  !- First row of transformation matrix between global and body-fixed coordinates 
  
!~   if (FCM_SWIMMING>=1) then
   if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,2) * FCM_QUAT(IB,4) &
                           - FCM_QUAT(IB,1) * FCM_QUAT(IB,3) )
                           
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,3) * FCM_QUAT(IB,4) &
                           + FCM_QUAT(IB,1) * FCM_QUAT(IB,2) )

     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 - 0.5 )
    end do 
   else
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,1)**2 + &
                             FCM_QUAT(IB,2)**2 - 0.5d0 )
                             
    ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
    ! LAST CHANGE : 2014/02/10: back to (-)   
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
    ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
    ! LAST CHANGE : 2014/02/10: back to (+)  
     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )

     if (FCM_USE_QUAT == 0) then                              
      FCM_P2(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) + FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
      FCM_P2(IB,2) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,3)**2 -0.5  )
      FCM_P2(IB,3) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,2) )  
      
      FCM_P3(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
      FCM_P3(IB,2) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,2) ) 
      FCM_P3(IB,3) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 -0.5  ) 
     end if
    end do
   end if
!~   end if
  
  
 else if (NPART_FULL ==2) then
 
  !- Rotation angle of local body frame
  PPPP(1) = 0.0  
  !- Rotation axis of local body frame
  PPPP(2) = 0.0
  PPPP(3) = 1.0
  PPPP(4) = 0.0    
  !- Quaternion from Euler parameters
  FCM_QUAT(1,1) = cos(PPPP(1)/2d0)
  FCM_QUAT(1,2:4) = PPPP(2:4)*sin(PPPP(1)/2d0) 
  
  theta_0 = 0.0!10.0*PPI/20.0
  !- Rotation angle of local body frame
  PPPP(1) = PPI + theta_0   
  !- Rotation axis of local body frame
  PPPP(2) = 0.0
  PPPP(3) = 1.0
  PPPP(4) = 0.0    
  !- Quaternion from Euler parameters
  FCM_QUAT(2,1) = cos(PPPP(1)/2d0)
  FCM_QUAT(2,2:4) = PPPP(2:4)*sin(PPPP(1)/2d0) 
  
  do IB = 1, NPART_FULL
 
   !- Quaternion normalization
   FCM_QUATNORM = dsqrt( FCM_QUAT(IB,1)**2 + & 
                         FCM_QUAT(IB,2)**2 + & 
                         FCM_QUAT(IB,3)**2 + & 
                         FCM_QUAT(IB,4)**2 )
                        
   FCM_QUAT(IB,:) = FCM_QUAT(IB,:)/FCM_QUATNORM
   
   !- Swimming Direction vector
   !- First row of transformation matrix between global and body-fixed coordinates 
   
!~    if (FCM_SWIMMING>=1) then
    if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then
    
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,2) * FCM_QUAT(IB,4) &
                           - FCM_QUAT(IB,1) * FCM_QUAT(IB,3) )
                           
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,3) * FCM_QUAT(IB,4) &
                           + FCM_QUAT(IB,1) * FCM_QUAT(IB,2) )

     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 - 0.5 )

    else
    
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,1)**2 + &
                             FCM_QUAT(IB,2)**2 - 0.5d0 )
                             
    ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
    ! LAST CHANGE : 2014/02/10: back to (-)   
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
    ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
    ! LAST CHANGE : 2014/02/10: back to (+)  
     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )                          
     if (FCM_USE_QUAT == 0) then                                                             
      FCM_P2(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) + FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
      FCM_P2(IB,2) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,3)**2 -0.5  )
      FCM_P2(IB,3) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,2) )  
      
      FCM_P3(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
      FCM_P3(IB,2) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,2) ) 
      FCM_P3(IB,3) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 -0.5  ) 
     end if
                             
    end if
!~    end if
  
  enddo
  
  print*,'FCM_PSWIM = ', FCM_PSWIM
!~   read(*,*)
 else if (NPART_FULL >=4) then
 

  
  !- Rotation axis of local body frame
  PPPP(2) = 0.0
  PPPP(3) = 0.0
  PPPP(4) = 0.0
  
  do IB = 1, NPART_FULL
  
  
!~    if (NPART_FULL == 4) then
!~    !- Rotation angle of local body frame
!~     if (IB==1) then    
!~      PPPP(1) = PPI/2.0
!~     else if (IB==2) then
!~      PPPP(1) = PPI
!~     else if (IB==3) then
!~      PPPP(1) = -PPI/2.0
!~     else if (IB==4) then
!~      PPPP(1) = 0.0
!~     end if
!~    else
    PPPP(1) = 0.0
!~    end if
  
   !- Quaternion from Euler parameters
   FCM_QUAT(IB,1) = cos(PPPP(1)/2d0)
   FCM_QUAT(IB,2:4) = PPPP(2:4)*sin(PPPP(1)/2d0)
  
   !- Quaternion normalization
   FCM_QUATNORM = dsqrt( FCM_QUAT(IB,1)**2 + & 
                         FCM_QUAT(IB,2)**2 + & 
                         FCM_QUAT(IB,3)**2 + & 
                         FCM_QUAT(IB,4)**2 )
                        
   FCM_QUAT(IB,:) = FCM_QUAT(IB,:)/FCM_QUATNORM
   
   !- Swimming Direction vector
   !- First row of transformation matrix between global and body-fixed coordinates 
   
!~    if (FCM_SWIMMING>=1) then
    if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then

     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,2) * FCM_QUAT(IB,4) &
                           - FCM_QUAT(IB,1) * FCM_QUAT(IB,3) )
                             
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,3) * FCM_QUAT(IB,4) &
                           + FCM_QUAT(IB,1) * FCM_QUAT(IB,2) )

     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 - 0.5 )

    else

     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,1)**2 + &
                             FCM_QUAT(IB,2)**2 - 0.5d0 )
                             
    ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
    ! LAST CHANGE : 2014/02/10: back to (-)   
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
    ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
    ! LAST CHANGE : 2014/02/10: back to (+)  
     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
     if (FCM_USE_QUAT == 0) then                                
      FCM_P2(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) + FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
      FCM_P2(IB,2) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,3)**2 -0.5  )
      FCM_P2(IB,3) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,2) )  
      
      FCM_P3(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
      FCM_P3(IB,2) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,2) ) 
      FCM_P3(IB,3) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 -0.5  ) 
     end if
    end if
!~    end if
  
  enddo
  
 end if

! Orientations from a file (Any value of NPART_FULL)
elseif (FCM_INIT_PART_ORIENT ==2) then
 
 if (FCM_USE_QUAT == 1) then
  FILENAME='FCM_PART_ORIENT.ini'
  call FCM_READ_VAR_MPIIO(NPART_FULL,4,1,FCM_QUAT,FILENAME,ERR_FILE_ORIENT)
  
  if (ERR_FILE_ORIENT.eq.1) then
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!                     ERROR                    !!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!'
   write(*,*)'!! file     : FCM_PART_ORIENT.ini' 
   write(*,*)'!!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop
  end if
  


  call MPI_BCAST(FCM_QUAT,4*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)
 

   if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,2) * FCM_QUAT(IB,4) &
                           - FCM_QUAT(IB,1) * FCM_QUAT(IB,3) )
                             
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,3) * FCM_QUAT(IB,4) &
                           + FCM_QUAT(IB,1) * FCM_QUAT(IB,2) )

     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 - 0.5 )
    end do 
   else
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,1)**2 + &
                             FCM_QUAT(IB,2)**2 - 0.5d0 )
                             
    ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
    ! LAST CHANGE : 2014/02/10: back to (-)   
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
    ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
    ! LAST CHANGE : 2014/02/10: back to (+)  
     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
    end do
   end if

 else if (FCM_USE_QUAT == 0) then

  FILENAME='FCM_PART_SWIM.ini'
  call FCM_READ_VAR_MPIIO(NPART_FULL,3,1,FCM_PSWIM,FILENAME,ERR_FILE_ORIENT)
  
  if (ERR_FILE_ORIENT.eq.1) then
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!                     ERROR                    !!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!'
   write(*,*)'!! file     : FCM_PART_SWIM.ini' 
   write(*,*)'!!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop
  end if

  call MPI_BCAST(FCM_PSWIM,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)

  DO IP=1,NPART_FULL
   FCM_P2(IP,1) = -FCM_PSWIM(IP,2) 
   FCM_P2(IP,2) = FCM_PSWIM(IP,1) 
   FCM_P2(IP,3) = 0.0

   FCM_P3(IP,1) = FCM_PSWIM(IP,2)*FCM_P2(IP,3) - FCM_PSWIM(IP,3)*FCM_P2(IP,2)     
   FCM_P3(IP,1) = -FCM_PSWIM(IP,1)*FCM_P2(IP,3) + FCM_PSWIM(IP,3)*FCM_P2(IP,1)     
   FCM_P3(IP,1) = FCM_PSWIM(IP,1)*FCM_P2(IP,2) - FCM_PSWIM(IP,2)*FCM_P2(IP,1)     
  END DO

  !FILENAME='FCM_PART_P2.ini'
  !call FCM_READ_VAR_MPIIO(NPART_FULL,3,1,FCM_P2,FILENAME,ERR_FILE_ORIENT)
  !if (ERR_FILE_ORIENT.eq.1) then
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! write(*,*)'!!                     ERROR                    !!'
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! write(*,*)'!!'
  ! write(*,*)'!! file     : FCM_PART_P2.ini' 
  ! write(*,*)'!!'
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! stop
  !end if
  !
  !call MPI_BCAST(FCM_P2,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)  
  !
  !FILENAME='FCM_PART_P3.ini'
  !call FCM_READ_VAR_MPIIO(NPART_FULL,3,1,FCM_P3,FILENAME,ERR_FILE_ORIENT)
  !if (ERR_FILE_ORIENT.eq.1) then
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! write(*,*)'!!                     ERROR                    !!'
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! write(*,*)'!!'
  ! write(*,*)'!! file     : FCM_PART_P3.ini' 
  ! write(*,*)'!!'
  ! write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  ! stop
  !end if
  !
  !call MPI_BCAST(FCM_P3,3*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)

 end if

 
! FCM_INIT_PART_ORIENT  
!     =3 : quaternion from ASCII file
elseif (FCM_INIT_PART_ORIENT ==3) then
  
 if (MYID==0) then
  FILENAME='FCM_PART_ORIENT.ini.ascii' 
  
  open(unit=200,file=FILENAME,status='old', iostat=ERR_FILE_ORIENT)
 
  if (ERR_FILE_ORIENT.eq.1) then
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!                     ERROR                    !!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   write(*,*)'!!'
   write(*,*)'!! file     : FCM_PART_ORIENT.ini.ascii'
   write(*,*)'!!'
   write(*,*)'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   stop
   
   
  else
  
   do IB = 1, NPART_FULL
    read(200,*) FCM_QUAT(IB,1:4)
   end do
   
  end if
  
  close(200)

 end if
  
 call MPI_BCAST(FCM_QUAT,4*NPART_FULL,MPI_DOUBLE,0,MPI_COMM_WORLD,IERR)

 
!~   if (FCM_SWIMMING>=1) then
   if ((FCM_NELLIPSOID>0).and.(FCM_ELLIPSOID_SIZE(1)/FCM_ELLIPSOID_SIZE(2)<1)) then
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,2) * FCM_QUAT(IB,4) &
                           - FCM_QUAT(IB,1) * FCM_QUAT(IB,3) )
                             
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,3) * FCM_QUAT(IB,4) &
                           + FCM_QUAT(IB,1) * FCM_QUAT(IB,2) )

     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 - 0.5 )
    end do 
   else
    do IB = 1, NPART_FULL
     FCM_PSWIM(IB,1) = 2.0*( FCM_QUAT(IB,1)**2 + &
                             FCM_QUAT(IB,2)**2 - 0.5d0 )
                             
    ! LAST CHANGE : 2013/12/4: Changed signs (-)->(+)       
    ! LAST CHANGE : 2014/02/10: back to (-)   
     FCM_PSWIM(IB,2) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
    ! LAST CHANGE : 2013/12/4: Changed signs (+)->(-) 
    ! LAST CHANGE : 2014/02/10: back to (+)  
     FCM_PSWIM(IB,3) = 2.0*( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + &
                             FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
                                  
     if (FCM_USE_QUAT == 0) then                            
      FCM_P2(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,3) + FCM_QUAT(IB,1)*FCM_QUAT(IB,4) )
      FCM_P2(IB,2) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,3)**2 -0.5  )
      FCM_P2(IB,3) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,2) )  
      
      FCM_P3(IB,1) = 2.0* ( FCM_QUAT(IB,2)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,3) )
      FCM_P3(IB,2) = 2.0* ( FCM_QUAT(IB,3)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,2) ) 
      FCM_P3(IB,3) = 2.0* ( FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 -0.5  ) 
     end if
    end do
   end if
!~   end if
 
 
else

 print *, "CPU -- ", MYID, ":: ERROR: Wrong initialization mode!"
 call MPI_FINALIZE(IERROR)	

endif ! if FCM_INIT_PART_ORIENT
   
   
if (FCM_NELLIPSOID>0) then

 
  do IB = 1, NPART_FULL
  !- Transformation matrix calculated from Nikravesh 1985 Eq 7
   FCM_ROT_MAT(IB,1,1) =  FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,2)**2 -0.5 
   FCM_ROT_MAT(IB,1,2) =  FCM_QUAT(IB,2)*FCM_QUAT(IB,3) - FCM_QUAT(IB,1)*FCM_QUAT(IB,4) 
   FCM_ROT_MAT(IB,1,3) =  FCM_QUAT(IB,2)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,3) 
   
   FCM_ROT_MAT(IB,2,1) =  FCM_QUAT(IB,2)*FCM_QUAT(IB,3) + FCM_QUAT(IB,1)*FCM_QUAT(IB,4)  
   FCM_ROT_MAT(IB,2,2) =  FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,3)**2 -0.5 
   FCM_ROT_MAT(IB,2,3) =  FCM_QUAT(IB,3)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,2) 
   
   FCM_ROT_MAT(IB,3,1) =  FCM_QUAT(IB,2)*FCM_QUAT(IB,4) - FCM_QUAT(IB,1)*FCM_QUAT(IB,3) 
   FCM_ROT_MAT(IB,3,2) =  FCM_QUAT(IB,3)*FCM_QUAT(IB,4) + FCM_QUAT(IB,1)*FCM_QUAT(IB,2) 
   FCM_ROT_MAT(IB,3,3) =  FCM_QUAT(IB,1)**2 + FCM_QUAT(IB,4)**2 -0.5 
  end do
 
  FCM_ROT_MAT = 2.0*FCM_ROT_MAT
  
end if   
   
end subroutine FCM_INITIATION_PARTICLE_ORIENTATION
