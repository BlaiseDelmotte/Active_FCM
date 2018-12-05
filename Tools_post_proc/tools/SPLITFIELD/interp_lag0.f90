!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     
!                   Lagrangian Polynomial INTERPOLATION               
!                                                                     
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                     
! The velocity field "UAIN" located on the mesh "XAIN,YAIN,ZAIN" is   
! interpolated on the mesh "XINT,YINT,ZINT" and the result is "UINT". 
!                                                                     
!=====================================================================
subroutine INTERP_LAG0(XAIN,YAIN,ZAIN, &
                           UAIN,           &
                           XINT,YINT,ZINT, &
                           UINT            )
!=====================================================================
!
!
!                  8--------7
!     k+1         /|       /|
!                5--------6 |
!                | |      | |
!                | 4------|-3
!                |/       |/
!      k         1--------2
!
!---------------------------------------------------------------------
!                             P. FEDE     --  I.M.F.T. --  31/03/2011
!---------------------------------------------------------------------

use dns_dim

implicit none

!---------------------------------------------------------------------
! ARRAYS STATEMENT
!---------------------------------------------------------------------

!=====================================================================
! Input Arrays
!=============
!- Mesh of data for interpolation
real(kind=8), dimension(ISTART(1):IEND(1)), intent(in) :: XAIN
real(kind=8), dimension(ISTART(2):IEND(2)), intent(in) :: YAIN
real(kind=8), dimension(ISTART(3):IEND(3)), intent(in) :: ZAIN

  
!- Data field for interpolation
real(kind=8),                                         &
 dimension(ISTART(1)         :IEND(1)               &
              ,ISTART(2)-NGHTCELL:IEND(2)+NGHTCELL      &
              ,ISTART(3)-NGHTCELL:IEND(3)+NGHTCELL ) ,  &
 intent(in) :: UAIN
  
  
  
!- Positions for interpolation 
real(kind=8), dimension(ISTARTI(1):IENDI(1)), intent(in) :: XINT
real(kind=8), dimension(ISTARTI(2):IENDI(2)), intent(in) :: YINT
real(kind=8), dimension(ISTARTI(3):IENDI(3)), intent(in) :: ZINT


!- Interpolated velocity field
real(kind=8),                                         &
 dimension(ISTARTI(1)         :IENDI(1)             &
              ,ISTARTI(2)-NGHTCELL:IENDI(2)+NGHTCELL    &
              ,ISTARTI(3)-NGHTCELL:IENDI(3)+NGHTCELL ), intent(out) :: UINT
  
                               
!!- Normalized distance to the reference node
real(kind=8) :: ALFA, BETA, GAMA


!!- Lagrange's Polynom
real(kind=8) :: FX0, FX1
real(kind=8) :: FY0, FY1
real(kind=8) :: FZ0, FZ1

!!- Mesh step
real(kind=8) :: DXAINT
real(kind=8) :: DYAINT
real(kind=8) :: DZAINT

!!- Index
integer :: I, J ,K
integer :: I0,  J0,  K0  
integer :: IP1, JP1, KP1


real(kind=8) :: EPSILON


intrinsic int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

EPSILON = 1E-8

DXAINT = XAIN(ISTART(1)+1) - XAIN(ISTART(1))
DYAINT = YAIN(ISTART(2)+1) - YAIN(ISTART(2))
DZAINT = ZAIN(ISTART(3)+1) - ZAIN(ISTART(3))

EPSILON = DXAINT/1000.



!- Debug time
!
!do I = ISTART(1),IEND(1)
!write(*,*)' I=',real(I),'x=',XAIN(I),' y=',YAIN(I),' z=',ZAIN(I)
!end do
!write(*,*)
!do I = ISTARTI(1),IENDI(1)
!write(*,*)' I=',real(I),'xi=',XINT(I),' yi=',YINT(I),' zi=',ZINT(I)
!end do
!
!- End debug time



do I = ISTARTI(1),IENDI(1)
 do J = ISTARTI(2),IENDI(2)
  do K = ISTARTI(3),IENDI(3)


!!====================================================================
!! 1. Interpolation nodes location
!!====================================================================
 I0 = int((XINT(I)+EPSILON)/ DXAINT) + 1
 J0 = int((YINT(J)+EPSILON)/ DYAINT) + 1
 K0 = int((ZINT(K)+EPSILON)/ DZAINT) + 1


! x-direction --> periodic boundary condition
!----------------------------------------------
 IP1 = I0 + 1
 if(I0 == IEND(1)) IP1 = 1


! y-direction
!----------------
 JP1 = J0 + 1

! z-direction
!----------------
 KP1 = K0 + 1


!!====================================================================
!! 2. Distance to the first node location
!!====================================================================
 ALFA = (XINT(I)-XAIN(I0)) / DXAINT
 BETA = (YINT(J)-YAIN(J0)) / DYAINT
 GAMA = (ZINT(K)-ZAIN(K0)) / DZAINT


!!====================================================================
!! 3. Lagrange polynomial functions
!!====================================================================


!!====================================================================
!! 4.Interpolation
!!====================================================================

   UINT(I,J,K) = UAIN(I0,J0,K0)

  end do
 end do
end do 




if(MYID==0) write(*,*)'Interpolation Lagrange 1st order --> OK'
if(MYID==0) write(*,*)'                         EPSILON --> ',EPSILON


10000 format (8(A,f6.2,4x))


end subroutine INTERP_LAG0
