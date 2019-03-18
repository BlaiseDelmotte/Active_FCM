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
subroutine PRINT_PARAVIEW_PART(TIME,&
                               NP, &
                               POSI, &
                               POSI_NOPER, &
                               VEL, &
                               ROT, &
                               PSWIM, &
                               P_SCAL_PMEAN, &
                               RAD, &
                               FCM_VSW)


implicit none


integer,                     intent(in) :: TIME
integer,                     intent(in) :: NP
real(kind=8), dimension(NP,3), intent(in) :: POSI
real(kind=8), dimension(NP,3), intent(in) :: POSI_NOPER
real(kind=8), dimension(NP,3), intent(in) :: VEL
real(kind=8), dimension(NP,3), intent(in) :: ROT
real(kind=8), dimension(NP,3), intent(in) :: PSWIM
real(kind=8), dimension(NP), intent(in) :: P_SCAL_PMEAN
real(kind=8), dimension(NP), intent(in) :: RAD
real(kind=8), intent(in) :: FCM_VSW
!---------------------------------------------------------------------
character(len=40) :: FILENAME

real(kind=8) :: MINRAD
integer :: IP
!
!.. Executable Statements ..

MINRAD = minval(RAD)

write(FILENAME,10100)'PART_KINEMATICS_t_',TIME,'.vtu'
! POSITION AVEC CONDITONS PERIODIQUES
open(1984,file=trim(FILENAME))
write(1984,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(1984,'(A)') '<UnstructuredGrid>' 
write(1984,23) NP
! Sortie des positions des particules -----------------------
write(1984,'(A)') '   <Points>'
write(1984,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
! One can either use min radius or each particle radius
!write(1984,*)     (POSI(IP,1)/RAD(IP), POSI(IP,2)/RAD(IP), POSI(IP,3)/RAD(IP), IP=1,NP)
write(1984,*)     (POSI(IP,1)/MINRAD, POSI(IP,2)/MINRAD, POSI(IP,3)/MINRAD, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
write(1984,'(A)') '   </Points>'
! Sortie des donnees associees aux particules
!-----------------------------------------------
write(1984,'(A)') '   <PointData>'
!1. vitesses des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Vitesse" NumberOfComponents="3" format="ascii">'
write(1984,*)     (VEL(IP,1)/FCM_VSW, VEL(IP,2)/FCM_VSW, VEL(IP,3)/FCM_VSW, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!2. rotation des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Rotation" NumberOfComponents="3" format="ascii">'
write(1984,*)     (ROT(IP,1), ROT(IP,2), ROT(IP,3), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!3. orientatiion des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Orientation" NumberOfComponents="3" format="ascii">'
write(1984,*)     (PSWIM(IP,1), PSWIM(IP,2), PSWIM(IP,3), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!3b. produit scalaire avec vecteur orientation moyen
write(1984,'(A)') '        <DataArray type="Float32" Name="Scal_Pmean"  format="ascii">'
write(1984,*)     (P_SCAL_PMEAN(IP), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!4. rayons des particules 
write(1984,'(A)') '        <DataArray type="Float32" Name="Rayon"  format="ascii">'
write(1984,*)     (RAD(IP)/MINRAD, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!5. Numeros identifiants des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="NumId"  format="ascii">'
write(1984,*)     (IP, IP=1,NP)
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


write(FILENAME,10100)'PART_KINEMATICS_NOPER_t_',TIME,'.vtu'
! POSITION SANS CONDITONS PERIODIQUES
open(1984,file=trim(FILENAME))
write(1984,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
write(1984,'(A)') '<UnstructuredGrid>' 
write(1984,23) NP
! Sortie des positions des particules sans conditions periodiques-----------------------
write(1984,'(A)') '   <Points>'
write(1984,'(A)') '        <DataArray type="Float32" NumberOfComponents="3" format="ascii">'
! One can either use min radius or each particle radius
!write(1984,*)     (POSI_NOPER(IP,1)/RAD(IP), POSI_NOPER(IP,2)/RAD(IP), POSI_NOPER(IP,3)/RAD(IP), IP=1,NP)
write(1984,*)     (POSI_NOPER(IP,1)/MINRAD, POSI_NOPER(IP,2)/MINRAD, POSI_NOPER(IP,3)/MINRAD, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
write(1984,'(A)') '   </Points>'
! Sortie des donnees associees aux particules
!-----------------------------------------------
write(1984,'(A)') '   <PointData>'
!1. vitesses des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Vitesse" NumberOfComponents="3" format="ascii">'
write(1984,*)     (VEL(IP,1)/FCM_VSW, VEL(IP,2)/FCM_VSW, VEL(IP,3)/FCM_VSW, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!2. rotation des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Rotation" NumberOfComponents="3" format="ascii">'
write(1984,*)     (ROT(IP,1), ROT(IP,2), ROT(IP,3), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!3. orientatiion des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="Orientation" NumberOfComponents="3" format="ascii">'
write(1984,*)     (PSWIM(IP,1), PSWIM(IP,2), PSWIM(IP,3), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!3b. produit scalaire avec vecteur orientation moyen
write(1984,'(A)') '        <DataArray type="Float32" Name="Scal_Pmean"  format="ascii">'
write(1984,*)     (P_SCAL_PMEAN(IP), IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!4. rayons des particules 
write(1984,'(A)') '        <DataArray type="Float32" Name="Rayon"  format="ascii">'
write(1984,*)     (RAD(IP)/MINRAD, IP=1,NP)
write(1984,'(A)') '        </DataArray>'
!5. Numeros identifiants des particules
write(1984,'(A)') '        <DataArray type="Float32" Name="NumId"  format="ascii">'
write(1984,*)     (IP, IP=1,NP)
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
end subroutine PRINT_PARAVIEW_PART
!***
