!***** EnsightCase
! SYNOPSIS
!  To write a Ensight's case file 
! AUTHOR
!  A. Pedrono 
!  Service Codes et Simulations Numeriques : P. Elyakime, H. Neau, A. Pedrono, A. Stoukov
!  Institut de Mecanique des Fluides de Toulouse
!  Adapted from Stephane.Montesino's documents  (hmg.inpg.fr)
! CREATION DATE
!  23/10/2012
!****************************************************************************** 
!
!    VarNameArray.....: Name of the variables  
!    GeoName..........: Name of the geometrie
!    VarTypeArray.....: Variables type 1 => Scalar       3 => Vector
!    ntini.......: filename start number
!    nstop.......: filename end number
!    nstep......: filename increment
!    tab_time....: physical times  
!
!
     subroutine EnsightCase(CaseName,VarNameArray,VarTypeArray,nb_var, &
                            GeoName,ntini,nstop,nstep,tab_time)

     implicit none

     
     character(len=20) :: CaseName
     character(len=80) :: GeoName
     integer, intent(in) :: nb_var
     character(len=8), dimension(nb_var), intent(in) :: VarNameArray
     integer, dimension(nb_var),intent(in) :: VarTypeArray
     integer,intent(in)::nstop,ntini,nstep
     integer::FileUnit,i,nfile,l
     real (kind=8), dimension((nstop-ntini+1)/nstep) :: tab_time

     write(*,'(/2A)') ' Creating case file for Ensight and Paraview: ' &
                       ,VarNameArray


     nfile=(nstop-ntini)/nstep+1 
     print*,'nfile = ', nfile

     FileUnit = 140
     open(FileUnit,file=trim(CaseName)//'.case')

     write(FileUnit,10) trim(GeoName)//'.geo'
  10 format( &
      'FORMAT'            ,/ , &
      'type: ensight gold',//, &
      'GEOMETRY'          ,/ , &
      'model:    ',A         ,//, &
      'VARIABLE')


     do l=1,nb_var
       if(VarTypeArray(l) == 1) &
         write(FileUnit,15)trim(VarNameArray(l)), & 
                           trim(VarNameArray(l))//'********.scl'
       if(VarTypeArray(l) == 3) &
         write(FileUnit,25)trim(VarNameArray(l)), &
                           trim(VarNameArray(l))//'********.vec'
     end do
     if (nfile > 1) then
       write(FileUnit,45) nfile,ntini,nstep
       write(FileUnit,'(f15.3)') (tab_time(i),i=1,nfile)
     end if

     close(FileUnit)

  15 format('scalar per node: ',A,'   ', A)
  25 format('vector per node: ',A,'   ', A)

  45 format( &
     /,'TIME            '      , &
     /,'time set: 1     '      , &
     /,'number of steps: '      ,i4 , &
     /,'filename start number: ',i10 &
     /,'filename increment: '   ,i10 &
     /,'time values: ' &
     )

     end subroutine
