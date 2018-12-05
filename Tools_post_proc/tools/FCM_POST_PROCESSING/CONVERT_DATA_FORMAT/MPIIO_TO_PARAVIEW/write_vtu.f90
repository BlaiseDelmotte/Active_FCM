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
subroutine write_vtu(filevtu) 

  use mod_commun
  implicit none
 
  !.. Scalars ..
  integer :: itraj,nbin1,nbin2,pas
  !
  !.. Arrays ..
  character(len=4) :: nom_cas
  character(len=16), intent(in) :: filevtu
  integer :: np,ig
  !
  !.. Executable Statements ..
  !
  do np=1,npart
    do ig=1,nig
     id_part(np,ig)=1.*npart*(1-ig)+1.*np
    end do
   end do
  !
  open(1984,file=filevtu)
  write(1984,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
  write(1984,'(A)') '<UnstructuredGrid>' 
  write(1984,23) npart
  ! Sortie des positions des particules -----------------------
  write(1984,'(A)') '   <Points>'
  write(1984,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  write(1984,*)     ((xabspart(np,ig), yabspart(np,ig), zabspart(np,ig), np=1,npart),ig=1,nig)
  write(1984,'(A)') '        </DataArray>'
  write(1984,'(A)') '   </Points>'
  ! Sortie des donnees associees aux particules
  !-----------------------------------------------
  write(1984,'(A)') '   <PointData>'
  !1. vitesses des particules
  write(1984,'(A)') '        <DataArray type="Float32" Name="Vitesse" NumberOfComponents="3" format="ascii">'
  write(1984,*)     ((uabspart(np,ig), vabspart(np,ig), wabspart(np,ig), np=1,npart),ig=1,nig)
  write(1984,'(A)') '        </DataArray>'
  !2. rayons des particules 
  write(1984,'(A)') '        <DataArray type="Float32" Name="Rayon"  format="ascii">'
  write(1984,*)     ((rpart, np=1,npart),ig=1,nig)
  write(1984,'(A)') '        </DataArray>'
  !3. Numeros identifiants des particules
  write(1984,'(A)') '        <DataArray type="Float32" Name="NumId"  format="ascii">'
  write(1984,*)     ((id_part(np,ig), np=1,npart),ig=1,nig)
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

  write(*,*) 'Creation du Fichier ', filevtu

  ! ... Format Declarations ...
  23 format (1x,'<Piece NumberOfPoints=" ',i12,' " NumberOfCells="0">')

end subroutine write_vtu
!***
