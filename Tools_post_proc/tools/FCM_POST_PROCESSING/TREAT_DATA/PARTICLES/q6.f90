!!====================================================================
!!
!!
!!====================================================================

subroutine Q6_NEIGHBOR(NSAVES, &
                       SAVE_START, &
                       PART_START, &
                       PART_END, &
                       LX, &
                       LY, &
                       LZ, &
                       RADMAX, &
                       POS, &
                       DTIME, &
                       FOUT3)

!!====================================================================
!!
!!
!!====================================================================


implicit none


!--------------------------------------------------------------------
! ARGUMENTS
!--------------------------------------------------------------------
! Number of snapshots
integer  :: NSAVES, FOUT3 !intent(in)
! Where to start
integer, intent(in) :: SAVE_START
! Number of swimmers
integer, intent(in) :: PART_START, PART_END
! Height of the domain
real(kind=8), intent(in) :: LX, LY, LZ
! Max radius
real(kind=8), intent(in) :: RADMAX
! Particle velocities
real(kind=8), dimension(NSAVES-SAVE_START,PART_END-PART_START+1,3), intent(in) :: POS
! Time step
real(kind=8), intent(in) :: DTIME
!--------------------------------------------------------------------
! ARRAYS STATEMENT
!--------------------------------------------------------------------
!- Edge time distrib and its standars deviation for error bars
real ,dimension(:,:), allocatable :: LIST_PART_HAUT, LIST_PART_BAS
real ,dimension(:,:,:), allocatable :: LIST_VOISIN_HAUT, LIST_VOISIN_BAS
real ,dimension(:,:), allocatable :: Q6_HAUT, Q6_BAS
real ,dimension(:,:), allocatable :: Q_ALL

integer,dimension(1) :: MAX_POS_1

!TYPE int_real
!  INTEGER   :: ints
!  REAL      :: floats
!ENDTYPE int_real
!
!TYPE(int_real), allocatable, dimension(:,:,:,:) :: LIST_PART_HAUT, LIST_PART_BAS

!- Step size to discretize the interval
real(kind=8) :: start, finish,norme,  R_range,theta,y1,y2,z1,z2,y_relatif,z_relatif
complex(kind=4) :: q
!- File name 
character(len=100) :: FILENAME

!- String for saves
character(len=10) :: FILE_NUM

!- Index
integer :: IND,i,loc,NB_BAS,NB_HAUT,dt, ID, IND_HAUT_PART_1, IND_HAUT_PART_2, IND_HAUT, COMPT_PART_HAUT, K, L
integer :: IND_BAS_PART_1, IND_BAS_PART_2, IND_BAS, compt_dt, voisin, NPART_FULL, COMPT_PART_BAS

character(len=40) :: FCM_FILENAME_Q



print*,' '
print*,'-------------------START Q6 ------------------------ '
print*,' '
!---------------------------------------------------------------------

 
!=====================================================================
! 1. RUN Q6
!=====================================================================
call CPU_TIME(start)

NPART_FULL = PART_END - PART_START + 1
R_range = 1.1*RADMAX
allocate(Q_ALL(NSAVES,NPART_FULL))
Q_ALL(:,1:NPART_FULL) = 0

print *, ' NSAVES = ',NSAVES,' SAVE_START = ',SAVE_START

print *, ' NPART_FULL = ', NPART_FULL
print *, ' NSAVES-SAVE_START = ', NSAVES-SAVE_START

do dt  = 1,NSAVES-SAVE_START 


    NB_HAUT = 0;
    NB_BAS = 0;
        
    do ID = 1, NPART_FULL
        if (POS(dt, ID, 1) > ((LX/2)-R_range)) then
                NB_HAUT = NB_HAUT + 1
        endif
        if (POS(dt, ID, 1) < R_range) then        !parois du bas
                NB_BAS = NB_BAS + 1
        endif

    enddo

    !print *, 'dt = ', dt , 'NB_HAUT = ' , NB_HAUT , 'NB_BAS = ' , NB_BAS
    
    allocate(LIST_PART_HAUT(NB_HAUT,4))
    allocate(LIST_PART_BAS(NB_BAS,4))
    allocate(LIST_VOISIN_HAUT(NB_HAUT,6,3))
    allocate(LIST_VOISIN_BAS(NB_BAS,6,3))
    allocate(Q6_HAUT(NB_HAUT, 2))
    allocate(Q6_BAS(NB_BAS, 2))
        
    LIST_VOISIN_HAUT(:,:,1) = 0
    LIST_VOISIN_BAS(:,:,1) = 0        
    LIST_VOISIN_HAUT(:,:,2) = LX
    LIST_VOISIN_BAS(:,:,2) = LX
    LIST_VOISIN_HAUT(:,:,3) = -1
    LIST_VOISIN_BAS(:,:,3) = -1

    IND_HAUT = 1;
    IND_BAS = 1;

    do ID = 1, NPART_FULL
        if (POS(dt, ID, 1) > ((LX/2)-R_range)) then
                
                LIST_PART_HAUT(IND_HAUT,1) = ID
                LIST_PART_HAUT(IND_HAUT,2:4) = POS(dt, ID, 1:3)
                IND_HAUT = IND_HAUT + 1
        endif
        if (POS(dt, ID, 1) < R_range) then        !parois du bas
                
                LIST_PART_BAS(IND_BAS,1) = ID
                LIST_PART_BAS(IND_BAS,2:4) = POS(dt, ID, 1:3)
                IND_BAS = IND_BAS + 1
        endif
        
    enddo
    
    !print *, ' '
    !print *, ' LIST HAUT ET BAS OK'
    !print *, ' '

    COMPT_PART_HAUT = 0
    if (size(LIST_PART_HAUT,1) > 1) then
            
        do IND_HAUT_PART_1 = 1 , size(LIST_PART_HAUT,1) -1

                y1 = LIST_PART_HAUT(IND_HAUT_PART_1,3)
                z1 = LIST_PART_HAUT(IND_HAUT_PART_1,4)
        
                do IND_HAUT_PART_2 =  IND_HAUT_PART_1 + 1 , size(LIST_PART_HAUT,1)
            
                        y2 = LIST_PART_HAUT(IND_HAUT_PART_2,3)
                        z2 = LIST_PART_HAUT(IND_HAUT_PART_2,4)
                 
                        y_relatif = y2 - y1
                        z_relatif = z2 - z1
            
                        if (y_relatif >= 0) then
                                y_relatif = y_relatif - LY*real(floor(y_relatif/(0.5*LY)))
                        endif
                        if (z_relatif >= 0) then
                                z_relatif = z_relatif - LZ*real(floor(z_relatif/(0.5*LZ)))
                        endif
                        if (y_relatif <= 0) then
                                y_relatif = y_relatif - LY*real(ceiling(y_relatif/(0.5*LY)))
                        endif
                        if (z_relatif <= 0) then
                                z_relatif = z_relatif - LZ*real(ceiling(z_relatif/(0.5*LZ)))
                        endif

            
                        norme = sqrt((y_relatif**2) + (z_relatif**2))                !calcule de la norme            
                        
                        if (any(LIST_VOISIN_HAUT(IND_HAUT_PART_1,:,1) /= LIST_PART_HAUT(IND_HAUT_PART_2,1))) then          
                        
                                MAX_POS_1 = maxloc(LIST_VOISIN_HAUT(IND_HAUT_PART_1,:,2))
                                                                                                
                                if (norme < LIST_VOISIN_HAUT(IND_HAUT_PART_1,MAX_POS_1(1),2))  then
                                                       
                                        theta = atan2(z_relatif,y_relatif) + 4*datan(1.D0)
                    
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_1,MAX_POS_1(1),1) = LIST_PART_HAUT(IND_HAUT_PART_2,1)
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_1,MAX_POS_1(1),2) = norme
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_1,MAX_POS_1(1),3) = theta !LIST_VOISIN : Dim 1(Nb iteration temporeel) : pas de temps , Dim 2(Nb particule parois) : ID particule dans LIST_PART_BAS , Dim 3(6) : Liste des voisins, Dim 4(3) : 1-> ID particule dans POS 2-> norme 3-> angle
                                endif
                        endif
            
                        if (any(LIST_VOISIN_HAUT(IND_HAUT_PART_2,:,1) /= LIST_PART_HAUT(IND_HAUT_PART_1,1))) then             
              
                                MAX_POS_1 = maxloc(LIST_VOISIN_HAUT(IND_HAUT_PART_2,:,2))
                
                                if (norme < LIST_VOISIN_HAUT(IND_HAUT_PART_2,MAX_POS_1(1),2))  then
                                                       
                                        theta = atan2(z_relatif,y_relatif) + 4*datan(1.D0)    
                    
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_2,MAX_POS_1(1),1) = LIST_PART_HAUT(IND_HAUT_PART_1,1)
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_2,MAX_POS_1(1),2) = norme
                                        LIST_VOISIN_HAUT(IND_HAUT_PART_2,MAX_POS_1(1),3) = theta                          
                                endif
                        endif
                enddo
                COMPT_PART_HAUT = COMPT_PART_HAUT +1;
        enddo
    endif
    


    !print *, ' '
    !print *, ' LIST VOISIN HAUT OK'
    !print *, ' '

    COMPT_PART_BAS = 0
    if (size(LIST_PART_BAS,1) > 1) then
        do IND_BAS_PART_1 = 1 , size(LIST_PART_BAS,1) -1

                y1 = LIST_PART_BAS(IND_BAS_PART_1,3)
                z1 = LIST_PART_BAS(IND_BAS_PART_1,4)

                do IND_BAS_PART_2 =  IND_BAS_PART_1 + 1 , size(LIST_PART_BAS,1)

                        y2 = LIST_PART_BAS(IND_BAS_PART_2,3)
                        z2 = LIST_PART_BAS(IND_BAS_PART_2,4)

                        y_relatif = y2 - y1
                        z_relatif = z2 - z1

                        if (y_relatif >= 0) then
                                y_relatif = y_relatif - LY*real(floor(y_relatif/(0.5*LY)))
                        endif
                        if (z_relatif >= 0) then
                                z_relatif = z_relatif - LZ*real(floor(z_relatif/(0.5*LZ)))
                        endif
                        if (y_relatif <= 0) then
                                y_relatif = y_relatif - LY*real(ceiling(y_relatif/(0.5*LY)))
                        endif
                        if (z_relatif <= 0) then
                                z_relatif = z_relatif - LZ*real(ceiling(z_relatif/(0.5*LZ)))
                        endif

                        norme = sqrt((y_relatif**2) + (z_relatif**2))                !calcule de la norme            

                        if (any(LIST_VOISIN_BAS(IND_BAS_PART_1,:,1) /= LIST_PART_BAS(IND_BAS_PART_2,1))) then

                                MAX_POS_1 = maxloc(LIST_VOISIN_BAS(IND_BAS_PART_1,:,2))

                                if (norme < LIST_VOISIN_BAS(IND_BAS_PART_1,MAX_POS_1(1),2)) then

                                        theta = atan2(z_relatif,y_relatif) + 4*datan(1.D0)     

                                        LIST_VOISIN_BAS(IND_BAS_PART_1,MAX_POS_1(1),1) = LIST_PART_BAS(IND_BAS_PART_2,1)
                                        LIST_VOISIN_BAS(IND_BAS_PART_1,MAX_POS_1(1),2) = norme
                                        LIST_VOISIN_BAS(IND_BAS_PART_1,MAX_POS_1(1),3) = theta !LIST_VOISIN : Dim 1(Nb iteration temporeel) : pas de temps , Dim 2(Nb particule parois) : ID particule dans LIST_PART_BAS , Dim 3(6) : Liste des voisins, Dim 4(3) : 1-> ID particule dans POS 2-> norme 3-> angle


                                endif
                        endif
                        if (any(LIST_VOISIN_BAS(IND_BAS_PART_2,:,1) /= LIST_PART_BAS(IND_BAS_PART_1,1))) then

                                MAX_POS_1 = maxloc(LIST_VOISIN_BAS(IND_BAS_PART_2,:,2))
                                
                                if (norme < LIST_VOISIN_BAS(IND_BAS_PART_2,MAX_POS_1(1),2)) then

                                        theta = atan2(z_relatif,y_relatif) + 4*datan(1.D0)

                                        LIST_VOISIN_BAS(IND_BAS_PART_2,MAX_POS_1(1),1) = LIST_PART_BAS(IND_BAS_PART_1,1)
                                        LIST_VOISIN_BAS(IND_BAS_PART_2,MAX_POS_1(1),2) = norme
                                        LIST_VOISIN_BAS(IND_BAS_PART_2,MAX_POS_1(1),3) = theta
                                endif
                        endif
                enddo
                COMPT_PART_BAS = COMPT_PART_BAS +1;
        enddo

    endif

    
        do IND_HAUT = 1 , size(LIST_VOISIN_HAUT,1)
                q = 0
                        do voisin = 1 , 6
                                if (LIST_VOISIN_HAUT( IND_HAUT, voisin, 3) /= -1) then
                                        q = cexp(6*CMPLX(0.,1)*LIST_VOISIN_HAUT( IND_HAUT, voisin, 3)) + q
                                else if (LIST_VOISIN_HAUT( IND_HAUT, voisin, 3) == -1) then
                                        q = 0 + q
                                endif
                        enddo

                        Q6_HAUT(IND_HAUT,2)= cabs(q)
                        Q6_HAUT(IND_HAUT,2)=((1./6)*Q6_HAUT(IND_HAUT,2))**2

                        if (Q6_HAUT(IND_HAUT,2) > 1) then
                                print *, Q6_HAUT(IND_HAUT,2)
                                print *, '--- ERROR HAUT ---'
                                stop
                        end if

                        Q6_HAUT(IND_HAUT,1)= LIST_PART_HAUT(IND_HAUT,1)
        enddo
    

    
        do IND_BAS = 1 , size(LIST_VOISIN_BAS,1)      
                q = 0
                        do voisin = 1 , 6            
                                if (LIST_VOISIN_BAS(IND_BAS, voisin, 3) /= -1) then
                                        q = cexp(6*CMPLX(0.,1)*LIST_VOISIN_BAS(IND_BAS, voisin, 3)) + q
                                else if (LIST_VOISIN_BAS( IND_BAS, voisin, 3) == -1) then
                                        q = 0 + q
                                endif
                        enddo
                        
                        Q6_BAS(IND_BAS,2)= cabs(q)
                        Q6_BAS(IND_BAS,2)=((1./6)*Q6_BAS(IND_BAS,2))**2

                        if (Q6_BAS(IND_BAS,2) > 1) then
                                print *, '--- ERROR BAS ---'
                                stop
                        end if

                        Q6_BAS(IND_BAS,1)= LIST_PART_BAS(IND_BAS,1)
        enddo
    

    do IND = 1, NPART_FULL 
                do i = 1, size(Q6_HAUT,1)
                        if (Q6_HAUT(i,1) == IND) then
                                Q_ALL(dt,IND) = Q6_HAUT(i,2)
                        endif
                enddo
                do i = 1, size(Q6_BAS,1)
                        if (Q6_BAS(i,1) == IND) then
                                Q_ALL(dt,IND) = Q6_BAS(i,2)
                        endif
                enddo            
    enddo


    deallocate(LIST_PART_HAUT)
    deallocate(LIST_PART_BAS)
    deallocate(LIST_VOISIN_HAUT)
    deallocate(LIST_VOISIN_BAS)
    deallocate(Q6_HAUT)
    deallocate(Q6_BAS)
    
    compt_dt = compt_dt + 1

 enddo
 
print*,' '
print*,'  Q6 --> OK  '
print*,' '

!=====================================================================
! 2. WRITE Q6 .dat end .vtu
!=====================================================================
do K = 1, NSAVES-SAVE_START 
write(FILE_NUM,10205) (K-1)*FOUT3 + 1
    if(K <  NSAVES-SAVE_START+1) then
            write(FCM_FILENAME_Q,10101) 'FCM_PART_Q_t',trim(FILE_NUM),'.dat'
    else
            FCM_FILENAME_Q = 'FCM_PART_Q.end'
    end if

    open(unit=302,file=trim(FCM_FILENAME_Q))

    do L = 1, NPART_FULL
        write(302,'(3(e17.7))') Q_ALL(K,L)
       !print *, 'Q_ALL = ',Q_ALL(K,L) 
    end do
    close(302)
    write(*,*) 'Creation du Fichier ', FCM_FILENAME_Q
enddo



print*,' '
print*,'  WRITE Q6 --> OK '
print*,' '
call cpu_time(finish)
print '("Time = ",f6.3," seconds.")',finish-start
print*,' '
print*,'-------------------END Q6 ------------------------ '
print*,' '

deallocate(Q_ALL)
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
23 format (1x,'<Piece NumberOfPoints=" ',i12,' " NumberOfCells="0">')
10100 format (A,I8.8,A)

10200 format (A)
10201 format (A,I1)
10202 format (A,I2)
10203 format (A,I3)
10204 format (A,I4)
10205 format (I8.8)
10101 format (A,A,A)
10102 format (A,A,A,A,A)
10103 format (A,A,A,A,A,A,A,A,A)
end subroutine Q6_NEIGHBOR
