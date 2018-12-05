!!=====================================================================
!!
!!   Print_particle
!!
!!=====================================================================

subroutine PRINT_TECPLOT(TIME)

!!=====================================================================
!!
!!
!!=====================================================================

use DNS_DIM
use GEOMETRIC_VARIABLE !- 
use FLUID_VARIABLE     !- Fluid velocity
use PARTICLE_PARALLEL


implicit none



!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time
real(kind=8), intent(in) :: TIME

!- File name
character(len=30) :: FILENAME


!!--------------------------------------------------------------------
!! Tecplot's variable
!!--------------------------------------------------------------------
integer :: TECINI110, TECZNE110, TECDAT110, TECEND110 ! Tecplot functions
Integer VISDOUBLE
Integer STRANDID,PARENTZN
CHARACTER*1 NULCHAR
INTEGER   IDEBUG,III,DISDOUBLE
INTEGER    IMAX, JMAX, KMAX
integer,POINTER  ::   NullPtr => null()

real(kind=8), dimension(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3)) :: LOC 
!---------------------------------------------------------------------


!- Index
integer :: I, J, K, IJK
!---------------------------------------------------------------------




IDEBUG    = 1 !- 0: no informations, 1: info, 2: info+data
VISDOUBLE = 0
DISDOUBLE = 1
NULCHAR   = CHAR(0)


!!====================================================================
!! 1. Print fluid velocity field
!!====================================================================
!~ !!-Print filename
!~ write(FILENAME,10101)'uf_t',NFILEOUT,trim(FILE_EXT),'.dat'
!~ !!- ASCII
!~ open(unit=300,file=trim(FILENAME))
!~ !!- Ecriture ASCII
!~ write(300,1999)
!~ write(300,2001)ISIZE(1),ISIZE(2),ISIZE(3),TIME
!~ do K = ISTART(3), IEND(3)
!~  do J = ISTART(2), IEND(2)
!~   do I = ISTART(1), IEND(1)
!~    write(300,'(6(e17.7))')XMESH(I), YMESH(J), ZMESH(K), &
!~                                    UFLU(I,J,K), &
!~                                    VFLU(I,J,K), &
!~                                    WFLU(I,J,K)
!~                                    
!~   end do
!~  end do
!~ end do
!~ !- close file
!~ close(300)


write(FILENAME,10101)'uf_t',NFILEOUT,trim(FILE_EXT),'.plt'

!!- Open file for tecplot format
III = TECINI110('Fluid'//NULCHAR,'X Y Z U V W'//NULCHAR,trim(FILENAME)//NULCHAR, &
              '.'//NULCHAR,IDEBUG,VISDOUBLE)


STRANDID = 1
PARENTZN = 0
III =  TECZNE110('Fluid'//NULCHAR,0,ISIZE(1),ISIZE(2),ISIZE(3),&
                   0,0,0,TIME,STRANDID,PARENTZN,1,& 
                   0,0,NULLPTR,NULLPTR,NULLPTR,0)


do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
  LOC(I,J,K) = XMESH(I) 
  end do
 end do
end do

IMAX=ISIZE(1)*ISIZE(2)*ISIZE(3)



III = TECDAT110(IMAX,(LOC(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
  LOC(I,J,K) = YMESH(J) 
  end do
 end do
end do

IMAX=ISIZE(1)*ISIZE(2)*ISIZE(3)
III = TECDAT110(IMAX,(LOC(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)

do K = ISTART(3), IEND(3)
 do J = ISTART(2), IEND(2)
  do I = ISTART(1), IEND(1)
  LOC(I,J,K) = ZMESH(K) 
  end do
 end do
end do

IMAX=ISIZE(1)*ISIZE(2)*ISIZE(3)
III = TECDAT110(IMAX,(LOC(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)




IMAX=ISIZE(1)*ISIZE(2)*ISIZE(3)
III = TECDAT110(IMAX,(UFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)
III = TECDAT110(IMAX,(VFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)
III = TECDAT110(IMAX,(WFLU(ISTART(1):IEND(1),ISTART(2):IEND(2),ISTART(3):IEND(3))),DISDOUBLE)


III = TECEND110()


if(MYID==0) write(*,*) 'Print tecplot: fluid velocity field --> OK'



!!- Print particle trajectory
if(SOLVE_PART) then

do J = 1, NIG

write(FILENAME,10102)'p',J,'_t',NFILEOUT,'.p',MYID,'.plt'

!!- Open file for tecplot format
III = TECINI110('Part'//NULCHAR,'X Y Z U V W'//NULCHAR,trim(FILENAME)//NULCHAR, &
              '.'//NULCHAR,IDEBUG,VISDOUBLE)

STRANDID = 1
PARENTZN = 0
III =  TECZNE110('Part'//NULCHAR,0,NPART_LOC(J),1,1,&
                   0,0,0,TIME,STRANDID,PARENTZN,1,   & 
                   0,0,NULLPTR,NULLPTR,NULLPTR,0)

III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%XP),DISDOUBLE)
III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%YP),DISDOUBLE)
III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%ZP),DISDOUBLE)
III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%UP),DISDOUBLE)
III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%VP),DISDOUBLE)
III = TECDAT110(NPART_LOC(J),sngl(PART(1:NPART_LOC(J),J)%WP),DISDOUBLE)


III = TECEND110()

end do


if(MYID==0) write(*,*) 'Print tecplot: particle velocity field --> OK'

end if



!pause

!===========================================================================
1998 format ('VARIABLES = "x", "y", "z"')
1999 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "n", "p"')
2001 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',e12.5)

10100 format (A,A,A)
10101 format (A,I2.2,A,A,A)
10102 format (A,I2.2,A,I2.2,A,I2.2,A)
10205 format(A,I4.4)

end subroutine PRINT_TECPLOT
