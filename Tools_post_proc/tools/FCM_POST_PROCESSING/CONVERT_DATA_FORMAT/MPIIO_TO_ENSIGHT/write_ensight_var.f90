!***** WriteEnsightVar
! SYNOPSIS
!  Writes result data in Ensight's format  
! AUTHORS
!  A. Pedrono 
!  Service Codes et Simulations Numeriques : P. Elyakime, H. Neau, A. Pedrono, A. Stoukov
!  Institut de Mecanique des Fluides de Toulouse
!  Adapted from Stephane.Montesino's documents  (hmg.inpg.fr)
! CREATION DATE
!  23/10/2012
!
!****************************************************************************** 
!
!
!    ndv..........: number of dimension of the variable (1=>scalar 3=>vector)
!    m1,m2,m3.....: size of the variable in the x1,x2,x3 direction
!    var..........: data to be written
!    Varname......: word used to build filenames
!    WriteBinary..: switch to choose Ensight Format in binary or ascii files
!    imin,imax....: range of writting data in the x1 direction
!    jmin,jmax....: range of writting data in the x2 direction
!    kmin,kmax....: range of writting data in the x3 direction
!    npt..........: time step number    
!
!    
!****************************************************************************** 

     subroutine WriteEnsightVar(ndv,m1,m2,m3,var,VarName,WriteBinary, &
                                imin,imax,jmin,jmax,kmin,kmax,npt)

     implicit none

     integer     ,intent(in)::m1,m2,m3,ndv
     integer     ,intent(in)::imin,imax,jmin,jmax,kmin,kmax
     real*8,dimension(ndv,m1,m2,m3),intent(in)::var
     character*8,intent(in)::Varname
     logical     ,intent(in)::WriteBinary
     integer     ,intent(in)::npt


     character(len=80):: VarFileName,buffer
     character(len=80):: part,block
     character(len=8):: increment
     integer::FileUnit,i,j,k,npart,m,reclength

     FileUnit = 140
     part ='part'
     npart=1
     block='block rectilinear'
     reclength=80*3+4*(1+(imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1)*ndv)
     if (npt < 0) then
       if (ndv.eq.1)VarFileName = trim(Varname)//'.scl'
       if (ndv.eq.3) then
        VarFileName = trim(Varname)//'.vec'
       end if
     else 
       write(increment,'(i8.8)')npt
       if (ndv.eq.1)VarFileName = trim(Varname)//increment//'.scl'
       if (ndv.eq.3) then 
        VarFileName = trim(Varname)//increment//'.vec'
        print*, 'VarFileName = ', VarFileName
       end if
     end if
    !      write(*,'(5x,A)') VarFileName
    
    if(WriteBinary) then
       open (unit=FileUnit,file=VarFileName,status="replace", &
!~              form='UNFORMATTED',access="direct",recl=reclength)
             form='BINARY',access="direct",recl=reclength)
       write(unit=FileUnit,rec=1) VarFileName &
                                 ,part,npart,block &
                                 ,((((SNGl(var(m,i,j,k)) &
                                 ,i=imin,imax) &
                                 ,j=jmin,jmax) &
                                 ,k=kmin,kmax) &
                                 ,m=1,ndv)
     else
       open (unit=FileUnit,file=VarFileName)
       write(buffer,'(a,a)') Varname        
       write(FileUnit,'(A)') buffer         
       write(FileUnit,'(A)') 'part'
       write(FileUnit,'(I10)')npart
       write(FileUnit,'(A)') 'block rectilinear'
       do m=1,ndv
         write(FileUnit,'(e12.5)') &
          (((SNGl(var(m,i,j,k)),i=imin,imax),j=jmin,jmax),k=kmin,kmax)
       enddo
     endif
     close(FileUnit)

     end  subroutine
!
!    
