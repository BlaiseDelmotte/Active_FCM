!!=====================================================================
!!
!!   Print_field + mesh in tecplot format
!!
!!=====================================================================

subroutine PRINT_TECPLOT(TIME,NX,NY,NZ,X,Y,Z,U,V,W)

!!=====================================================================
!!
!!
!!=====================================================================


implicit none



!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------
!- Time
integer,                            intent(in) :: TIME
integer,                            intent(in) :: NX, NY, NZ
real(kind=8),  dimension(NX),       intent(in) :: X
real(kind=8),  dimension(NY),       intent(in) :: Y
real(kind=8),  dimension(NZ),       intent(in) :: Z
real(kind=8),  dimension(NX,NY,NZ), intent(in) :: U, V, W

!- File name
character(len=30) :: FILENAME

!- String for time
character(len=10) :: FILE_EXT

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

real(kind=8), dimension(1:NX, 1:NY, 1:NZ) :: LOC 
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
!~ write(300,2001)NX,NY,NZ,TIME
!~ do K = 1, NZ
!~  do J = 1, NY
!~   do I = 1, NX
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

write(FILE_EXT,10205) TIME

write(FILENAME,10100)'uf_t_',trim(FILE_EXT),'.plt'

!!- Open file for tecplot format
III = TECINI110('Fluid'//NULCHAR,'X Y Z U V W'//NULCHAR,trim(FILENAME)//NULCHAR, &
              '.'//NULCHAR,IDEBUG,VISDOUBLE)


STRANDID = 1
PARENTZN = 0
III =  TECZNE110('Fluid'//NULCHAR,0,NX,NY,NZ,&
                   0,0,0,TIME,STRANDID,PARENTZN,1,& 
                   0,0,NULLPTR,NULLPTR,NULLPTR,0)


!!- Mesh
do K = 1, NZ
 do J = 1, NY
  do I = 1, NX
  LOC(I,J,K) = X(I) 
  end do
 end do
end do


IMAX=NX*NY*NZ

III = TECDAT110(IMAX,(LOC(1:NX, 1:NY, 1:NZ)),DISDOUBLE)

do K = 1, NZ
 do J = 1, NY
  do I = 1, NX
  LOC(I,J,K) = Y(J) 
  end do
 end do
end do


III = TECDAT110(IMAX,(LOC(1:NX, 1:NY, 1:NZ)),DISDOUBLE)

do K = 1, NZ
 do J = 1, NY
  do I = 1, NX
  LOC(I,J,K) = Z(K) 
  end do
 end do
end do

III = TECDAT110(IMAX,(LOC(1:NX, 1:NY, 1:NZ)),DISDOUBLE)



!!- Fields
III = TECDAT110(IMAX,U(1:NX, 1:NY, 1:NZ),DISDOUBLE)
III = TECDAT110(IMAX,V(1:NX, 1:NY, 1:NZ),DISDOUBLE)
III = TECDAT110(IMAX,W(1:NX, 1:NY, 1:NZ),DISDOUBLE)


III = TECEND110()


 write(*,*) 'Print tecplot: fluid velocity field --> OK'


!===========================================================================
1998 format ('VARIABLES = "x", "y", "z"')
1999 format ('VARIABLES = "x", "y", "z", "u", "v", "w"')
2000 format ('VARIABLES = "x", "y", "z", "u", "v", "w", "n", "p"')
2001 format ('ZONE F=POINT I=',i4,' J=',i4,' K=',i4,' SOLUTIONTIME=',e12.5)

10100 format (A,A,A)
10101 format (A,I2.2,A,A,A)
10102 format (A,I2.2,A,I2.2,A,I2.2,A)
10205 format (I8.8)

end subroutine PRINT_TECPLOT
