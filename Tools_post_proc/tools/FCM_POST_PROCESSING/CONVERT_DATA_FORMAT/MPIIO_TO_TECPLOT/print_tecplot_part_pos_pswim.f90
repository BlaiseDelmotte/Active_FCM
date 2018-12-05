!!=====================================================================
!!
!!   print Lagrangian field with TECPLOT
!!
!!=====================================================================

subroutine PRINT_TECPLOT_PART_POS_PSWIM(TIME,NP,XP,YP,ZP,UP,VP,WP,BINFLAG)

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
integer,                intent(in) :: TIME
integer,                     intent(in) :: NP
real(kind=8), dimension(NP), intent(in) :: XP
real(kind=8), dimension(NP), intent(in) :: YP
real(kind=8), dimension(NP), intent(in) :: ZP
real(kind=8), dimension(NP), intent(in) :: UP
real(kind=8), dimension(NP), intent(in) :: VP
real(kind=8), dimension(NP), intent(in) :: WP
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
!- String for time
character(len=10) :: FILE_EXT


!- Index
integer :: I, J, K
!---------------------------------------------------------------------

IDEBUG    = 1 !- 0: no informations, 1: info, 2: info+data
VISDOUBLE = 0
DISDOUBLE = 1
NULCHAR   = CHAR(0)


if(BINFLAG) then
!!--------------------------------------------------------------------
!! Tecplot's binary case
!!--------------------------------------------------------------------
write(FILE_EXT,10205) TIME

write(FILENAME,10100)'PART_POS_PSWIM_t_',trim(FILE_EXT),'.plt'

!!- Open file for tecplot format
III = TECINI110('Part'//NULCHAR,'XP YP ZP PX PY PZ'//NULCHAR,trim(FILENAME)//NULCHAR, &
              '.'//NULCHAR,IDebug,VIsDouble)

StrandID = 1
ParentZn = 0
III =  TECZNE110('Part'//NULCHAR,0,NP,1,1,        &
                   0,0,0,TIME,StrandID,ParentZn,1,   & 
                   0,0,NullPtr,NullPtr,NullPtr,0)


print*,'Hello81'
III = TECDAT110(NP,sngl(XP),DISDOUBLE)
print*,'Hello83'
III = TECDAT110(NP,sngl(YP),DISDOUBLE)
III = TECDAT110(NP,sngl(ZP),DISDOUBLE)
III = TECDAT110(NP,sngl(UP),DISDOUBLE)
III = TECDAT110(NP,sngl(VP),DISDOUBLE)
III = TECDAT110(NP,sngl(WP),DISDOUBLE)


III = TECEND110()
 write(*,*) 'Print tecplot: Particles pos and orientations --> OK'

else

!!--------------------------------------------------------------------
!! ASCII case
!!--------------------------------------------------------------------
!~ write(FILE_EXT,10205) TIME

write(FILENAME,10100)'PART_POS_PSWIM_t_',TIME,'.dat'

print*,'TIME = ', TIME
print*,'FILENAME = ', FILENAME

!- Open descriptor
 open(unit=300,file=trim(FILENAME))

 write(300,20000)
 write(300,20001)NP,1,1,TIME
 do I = 1,NP
  write(300,'(8(e17.7))')XP(I), YP(I), ZP(I), &
                         UP(I), VP(I), WP(I)
 end do
 close(300)
  write(*,*) 'Print ASCII: Particles pos and orientations --> OK'

end if






10000 format (A,A)
20000 format ('VARIABLES = "X", "Y", "Z", "px", "py", "pz"')
20001 format ('ZONE F=POINT I=',i6,' J=',i4,' K=',i4,' SOLUTIONTIME=',E13.6)

10100 format (A,I8.8,A)
10101 format (A,I2.2,A,A,A)
10102 format (A,I2.2,A,I2.2,A,I2.2,A)
10205 format (I8.8)

end subroutine PRINT_TECPLOT_PART_POS_PSWIM
