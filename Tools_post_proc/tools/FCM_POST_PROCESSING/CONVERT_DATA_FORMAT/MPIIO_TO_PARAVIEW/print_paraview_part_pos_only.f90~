!***** JADIM_TOOLS/traj2vtu
! NAME
!   write_vtu
! SYNOPSIS
! subroutine write_vtu
! AUTHOR
!   A. Chouippe 
! CREATION DATE
!  13/04/2010
! DESCRIPTION
! Cree des fichiers au format vtu pour la visualisation des particules
! Subroutine d'écriture 
!!! MODIFICATION HISTORY
! PARENTS
! CHILDREN
!!! ARGUMENTS
! SOURCE
subroutine PRINT_PARAVIEW_PART_POS_ONLY(TIME,&
                               NP, &
                               POSI, &                               
                               RAD)


implicit none


integer,                     intent(in) :: TIME
integer,                     intent(in) :: NP
real(kind=8), dimension(NP,3), intent(in) :: POSI
real(kind=8),intent(in) :: RAD

!---------------------------------------------------------------------
character(len=40) :: FILENAME

real(kind=8) :: MINRAD
integer :: IP
!
!.. Executable Statements ..



write(FILENAME,10100)'PART_KINEMATICS_t_',TIME,'.vtu'
!
open(1984,file=trim(FILENAME))
write(1984,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(1984,'(A)') '<UnstructuredGrid>' 
write(1984,23) NP
! Sortie des positions des particules -----------------------
write(1984,'(A)') '   <Points>'
write(1984,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
! One can either use min radius or each particle radius
!write(1984,*)     (POSI(IP,1)/RAD(IP), POSI(IP,2)/RAD(IP), POSI(IP,3)/RAD(IP), IP=1,NP)
write(1984,*)     (POSI(IP,1)/RAD, POSI(IP,2)/RAD, POSI(IP,3)/RAD, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
write(1984,'(A)') '   </Points>'
! Sortie des donnees associees aux particules
!-----------------------------------------------
write(1984,'(A)') '   <PointData>'
!1. Numeros identifiants des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="NumId"  format="ascii">'
write(1984,*)     (IP, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!2. rayons des particules 
write(1984,'(A)') '        <DataArray type="Float32" Name="Rayon"  format="ascii">'
write(1984,*)     (1.0, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
write(1984,'(A)') '   </PointData>'
write(1984,'(A)') '   <Cells>'
write(1984,'(A)') '      <DataArray type="Int32" Name="connectivity" format="ascii">'
write(1984,'(A)') '      </DataArray>'
write(1984,'(A)') '       <DataArray type="Int32" Name="offsets" format="ascii">'
write(1984,'(A)') '       </DataArray>'
write(1984,'(A)') '       <DataArray type="UInt8" Name="types" format="ascii">'
write(1984,'(A)') '       </DataArray>'
write(1984,'(A)') '   </Cells>'
write(1984,'(A)') ' </Piece>'  
write(1984,'(A)') '</UnstructuredGrid>'   
write(1984,'(A)') '</VTKFile>'
close(1984)

write(*,*) 'Creation du Fichier ', FILENAME

! ... Format Declarations ...
23 format (1x,'<Piece NumberOfPoints=" ',i12,' " NumberOfCells="0">')
10100 format (A,I8.8,A)
end subroutine PRINT_PARAVIEW_PART_POS_ONLY
!***
