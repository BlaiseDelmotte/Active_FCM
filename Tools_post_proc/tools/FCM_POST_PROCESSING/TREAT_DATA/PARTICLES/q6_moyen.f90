subroutine MEAN_Q6( NSAVES, &
                    SAVE_START, &
                    PART_START, &
                    PART_END, &
                    LX, &
                    RADMAX, &
                    POS, &
                    FOUT3 )

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
! Where to start
integer, intent(in) :: SAVE_START
! Number of swimmers
integer, intent(in) :: PART_START, PART_END
! Height of the domain
real(kind=8), intent(in) :: LX
! Max radius
real(kind=8), intent(in) :: RADMAX
! Particle velocities
real(kind=8), dimension(NSAVES-SAVE_START,PART_END-PART_START+1,3), intent(in) :: POS
! Screenshot
integer, intent(in) :: FOUT3


integer :: K, L, IND
integer :: NPART_FULL
character(len=40) :: FILENAME, FCM_FILENAME_Q, FILE_NUM
real(kind=8) :: Q, Q_T, R_RANGE
real(kind=8), dimension(NSAVES-SAVE_START) :: Q_TOT

R_RANGE = 1.1*RADMAX
NPART_FULL = PART_END-PART_START+1

print *, ''
print *, '------- MEAN Q6 START --------'
print *, ''

do K = 1, NSAVES-SAVE_START
        write(FILE_NUM,10205) (K-1)*FOUT3 + 1

        if(K <  NSAVES-SAVE_START + 1) then
                write(FCM_FILENAME_Q,10101) 'FCM_PART_Q_t',trim(FILE_NUM),'.dat'
        else
                FCM_FILENAME_Q = 'FCM_PART_Q.end'
        end if
        open(unit=302,file=trim(FCM_FILENAME_Q),status = 'old')

        Q_T = 0
        IND = 0
        do L = 1, NPART_FULL 
                read(302,*) Q
                if ( (POS(K,L,1) < R_RANGE) .OR. (POS(K,L,1) > ((LX/2) - R_RANGE )) ) then
                        Q_T = Q_T + Q
                        IND = IND + 1
                end if
        end do

        close(302)
        if (IND /= 0) then
                Q_TOT(K) = Q_T/IND
        else 
                Q_TOT(K) = 0
        end if

enddo

open(unit=314,file='FCM_Q_TEMPS.dat')

do K = 1,NSAVES-SAVE_START
        write(314,*) K,Q_TOT(K)
enddo

print *, ' FCM_Q_TEMPS.dat --> OK '

print *, ''
print *, '------- MEAN Q6 END --------'
print *, ''


close(314)
10205 format (I8.8)
10101 format (A,A,A)

end subroutine MEAN_Q6
