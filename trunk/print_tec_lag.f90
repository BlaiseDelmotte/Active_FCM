!!=====================================================================
!!
!!   print Lagrangian field with TECPLOT
!!
!!=====================================================================

subroutine PRINT_TEC_LAG(TIME,NP,XP,YP,ZP,UP,VP,WP,SCL,NAME,BINFLAG)

!!=====================================================================
!!
!!=====================================================================

implicit none

!!=====================================================================
!!
!!
!!=====================================================================

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
real(kind=8),                intent(in) :: TIME
integer,                     intent(in) :: NP
real(kind=8), dimension(NP), intent(in) :: XP
real(kind=8), dimension(NP), intent(in) :: YP
real(kind=8), dimension(NP), intent(in) :: ZP
real(kind=8), dimension(NP), intent(in) :: UP
real(kind=8), dimension(NP), intent(in) :: VP
real(kind=8), dimension(NP), intent(in) :: WP
real(kind=8), dimension(NP), intent(in) :: SCL
character(len=40),           intent(in) :: NAME
logical,                     intent(in) :: BINFLAG


!!--------------------------------------------------------------------
!! Tecplot's variable
!!--------------------------------------------------------------------
integer :: TECINI110, TECZNE110, TECDAT110, TECEND110 ! Tecplot functions
Integer VIsDouble
Integer StrandID,ParentZn
CHARACTER*1 NULCHAR
INTEGER   IDebug,III,DIsDouble
INTEGER    IMAX, JMAX, KMAX
integer,POINTER  ::   NullPtr => null()


!---------------------------------------------------------------------
character(len=40) :: FILENAME


!- Index
integer :: I, J, K
!---------------------------------------------------------------------

IDEBUG    = 0 !- 0: no informations, 1: info, 2: info+data
VISDOUBLE = 0
DISDOUBLE = 1
NULCHAR   = CHAR(0)


if(BINFLAG) then

!!--------------------------------------------------------------------
!! Tecplot's binary case
!!--------------------------------------------------------------------
 write(FILENAME,10000)trim(NAME),'.plt'

!!- Open file for tecplot format
III = TECINI110('Part'//NULCHAR,'X Y Z U V W COLOR'//NULCHAR,trim(FILENAME)//NULCHAR, &
              '.'//NULCHAR,IDebug,VIsDouble)

StrandID = 1
ParentZn = 0
III =  TECZNE110('Part'//NULCHAR,0,NP,1,1,        &
                   0,0,0,TIME,StrandID,ParentZn,1,   & 
                   0,0,NullPtr,NullPtr,NullPtr,0)

III = TECDAT110(NP,sngl(XP),DISDOUBLE)
III = TECDAT110(NP,sngl(YP),DISDOUBLE)
III = TECDAT110(NP,sngl(ZP),DISDOUBLE)
III = TECDAT110(NP,sngl(UP),DISDOUBLE)
III = TECDAT110(NP,sngl(VP),DISDOUBLE)
III = TECDAT110(NP,sngl(WP),DISDOUBLE)
III = TECDAT110(NP,sngl(SCL),DISDOUBLE)


III = TECEND110()





else

!!--------------------------------------------------------------------
!! ASCII case
!!--------------------------------------------------------------------
 write(FILENAME,10000)trim(NAME),'.dat'

!- Open descriptor
 open(unit=300,file=trim(FILENAME))

 write(300,20000)
 write(300,20001)NP,1,1,TIME
 do I = 1,NP
  write(300,'(8(e17.7))')XP(I), YP(I), ZP(I), &
                         UP(I), VP(I), WP(I)
 end do
 close(300)

end if






10000 format (A,A)
20000 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
20001 format ('ZONE F=POINT I=',i6,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

end subroutine PRINT_TEC_LAG
