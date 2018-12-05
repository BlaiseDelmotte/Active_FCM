!***** WriteRectEnsightGeo
! SYNOPSIS
!  Writes result data in Ensight's format for a rectilinear mesh 
! AUTHORS
!  A. Pedrono 
!  Service Codes et Simulations Numeriques : P. Elyakime, H. Neau, A. Pedrono, A. Stoukov
!  Institut de Mecanique des Fluides de Toulouse
!  Adapted from Stephane.Montesino's documents  (hmg.inpg.fr)
! CREATION DATE
!  23/10/2012
! DESCRITION
!    writes mesh data in Ensight's ascii or Binary format for rectilinear geometry
!
!    imax-imin,jmax-jmin,kmax-kmin : number of nodes in the x1,x2,x3 direction
!    x1,x2,x3....: coordinates
!
     subroutine WriteRectEnsightGeo(imin,imax,jmin,jmax,kmin,kmax, &
                                    x1,x2,x3,FileName,WriteBinary)

     implicit none

     integer,intent(in)::imin,imax,jmin,jmax,kmin,kmax
     real (kind=8), dimension(imax),intent(in)::x1
     real (kind=8), dimension(jmax),intent(in)::x2
     real (kind=8), dimension(kmax),intent(in)::x3
     logical               ,intent(in)::WriteBinary
     character(len=80)     ,intent(in)::FileName

     character(len=80)::binary_form
     character(len=80)::file_description1,file_description2
     character(len=80)::node_id,element_id
     character(len=80)::part,description_part,block

     integer::FileUnit,i,j,k,npart,isize,jsize,ksize
     integer::reclength

     FileUnit = 40

     binary_form      ='C Binary'
     file_description1='Ensight Model Geometry File Created by '
     file_description2='WriteRectEnsightGeo Routine'
     node_id          ='node id off'
     element_id       ='element id off'
     part             ='part'
     npart            =1
     description_part ='3D periodic channel'
     block            ='block rectilinear'
     isize=imax-imin+1
     jsize=jmax-jmin+1
     ksize=kmax-kmin+1

     reclength=80*8+4*(4+isize+jsize+ksize)

     if (WriteBinary) then
       open (unit=FileUnit,file=trim(FileName)//'.geo',status="replace", &
            form='UNFORMATTED',access="direct",recl=reclength)
       write(unit=FileUnit,rec=1) binary_form &
                                ,file_description1 &
                                ,file_description2 &
                                ,node_id &
                                ,element_id &
                                ,part,npart &
                                ,description_part &
                                ,block &
                                ,isize,jsize,ksize &
                                ,(sngl(x1(i)),i=imin,imax) &
                                ,(sngl(x2(j)),j=jmin,jmax) &
                                ,(sngl(x3(k)),k=kmin,kmax)
     else
       open (unit=FileUnit,file=trim(FileName)//'.geo')
       write(FileUnit,'(A)') 'Ensight Model Geometry File Created by '
       write(FileUnit,'(A)') 'WriteRectEnsightGeo Routine'
       write(FileUnit,'(A)') 'node id off'
       write(FileUnit,'(A)') 'element id off'
       write(FileUnit,'(A)') 'part'
       write(FileUnit,'(i10)')npart
       write(FileUnit,'(A)') '3D periodic channel'
       write(FileUnit,'(A)') 'block rectilinear'
       write(FileUnit,'(3i10)') isize,jsize,ksize
       write(FileUnit,'(E12.5)') (x1(i),i=imin,imax)
       write(FileUnit,'(E12.5)') (x2(j),j=jmin,jmax)
       write(FileUnit,'(E12.5)') (x3(k),k=kmin,kmax)
     endif

     end subroutine
!**********************************************************************
