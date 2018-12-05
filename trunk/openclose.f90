!!======================================================================
!!
!!
!!======================================================================
subroutine OPENCLOSE(FLAG)

!!======================================================================
!!
!! Units between 300 and 399 are reserved for the fluid
!!               400 and XXX                      particles 
!!
!!======================================================================
use dns_dim            !- Dimension
use param_phys         !- Physical & numerical parameters

implicit none


!-----------------------------------------------------------------------
! ARRAYS STATEMENT
!-----------------------------------------------------------------------
logical, intent(in) :: FLAG

character(len=40) :: FILENAME
integer :: I

!- number of opened files
integer :: NBFILE
!-----------------------------------------------------------------------

NBFILE = 0

if(FLAG) then

!!======================================================================
!! 1. Open files for fluid
!!======================================================================
 if(MYID==0) then
 
  if(LEVEL1_STFLU) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l1.stat'
   open(unit=301, file=trim(FILENAME), status='replace')
   write(301,20001)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL2_STFLU) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l2.stat'
   open(unit=302, file=trim(FILENAME), status='replace')
   write(302,20002)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL3_STFLU) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l3.stat'
   open(unit=303, file=trim(FILENAME), status='replace')
   write(303,20003)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
  if(LEVEL4_STFLU) then
   NBFILE = NBFILE + 1
   FILENAME = 'fluid_l4.stat'
   open(unit=304, file=trim(FILENAME), status='replace')
   write(304,20004)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if


!!======================================================================
!! 2. Open files for scalar
!!======================================================================
if(LEVEL0_STSCL) then 
  if(LEVEL1_STSCL) then
   NBFILE = NBFILE + 1
   FILENAME = 'scalar_l1.stat'
   open(unit=401, file=trim(FILENAME), status='replace')
   write(401,20201)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if

  if(LEVEL2_STSCL) then
   NBFILE = NBFILE + 1
   FILENAME = 'scalar_l2.stat'
   open(unit=402, file=trim(FILENAME), status='replace')
   write(402,20202)
   if(MYID==0) write(*,10700)trim(FILENAME)
  end if
end if




!!======================================================================
!! 3. Open files for particles
!!======================================================================
if(LEVEL0_STPAR) then 

 if(LEVEL1_STPAR) then
  do I =1, NIG

   NBFILE = NBFILE + 1
   write(FILENAME,10600)'part_l1_p',I,'.stat'
   open(unit=500+I, file=trim(FILENAME), status='replace')
   write(500+I,20101)
   if(MYID==0) write(*,10700)trim(FILENAME)

   NBFILE = NBFILE + 1
   write(FILENAME,10600)'part_l1_qp_p',I,'.stat'
   open(unit=520+I, file=trim(FILENAME), status='replace')
   write(520+I,20103)
   if(MYID==0) write(*,10700)trim(FILENAME)

   if(SOLVE_SCALAR) then
    NBFILE = NBFILE + 1
    write(FILENAME,10600)'part_l1_scl_p',I,'.stat'
    open(unit=560+I, file=trim(FILENAME), status='replace')
    write(560+I,20105)
    if(MYID==0) write(*,10700)trim(FILENAME)
   end if

   if(SOLVE_COLLISION > 0) then
    NBFILE = NBFILE + 1
    write(FILENAME,10600)'part_l1_col_p',I,'.stat'
    open(unit=580+I, file=trim(FILENAME), status='replace')
    write(580+I,20107)
    if(MYID==0) write(*,10700)trim(FILENAME)
   end if


  end do
 end if


! if(LEVEL2_STPAR) then
!  do I =1, NIG
!   NBFILE = NBFILE + 1
!   write(FILENAME,10600)'part_l2_Rp_p',I,'.stat'
!   open(unit=600+I, file=trim(FILENAME), status='replace')
!   if(SOLVE_SCALAR) then 
!    write(600+I,20106)
!   else
!    write(600+I,20104)
!   end if 
!   if(MYID==0) write(*,10700)trim(FILENAME)
!  end do
! end if


 end if


end if



if(MYID==0) write(*,10800) 'Files opened --> ',NBFILE

!!====================================================================
!! 4. CLOSE FILES
!!====================================================================
else

if(MYID==0) then

 close(301)
 close(302)
 close(303)
 close(304)


end if


end if


!-----------------------------------------------------------------------
10600 format (A,I2.2,A)
10700 format (1x,'   +   ',A)
10800 format (1x,A,I3)

20001 format('# t, <uf>, <vf>, <wf>, <kf>, <epsf>, <uf.uf>,<vf.vf>,<wf.wf>,<uf.vf>,<uf.wf>,<vf.wf>,eta_K,tau_K,v_K, Taylor g')
20002 format('# t, max(|dui/dxi|), < dui/dxi >, <(dui/dxi)^2>')
20003 format('# r, Rii, Rij, f, g')
20004 format('# t, Sk, Tk, Taylor f, Taylor g, Tay_f/Tay_g/2^0.5')

20101 format('# t, <up>, <vp>, <wp>, <uf@p>, <vf@p>, <wf@p>, Rp,xx, Rp,yy, Rp,zz,&
              Rp,xy, Rp,xz, Rp,yz,Rf@p,xx, Rf@p,yy, Rf@p,zz, Rf@p,xy, Rf@p,xz, Rf@p,yz, &
              Rfp,xx, Rfp,yy, Rfp,zz, Rfp,xy, Rfp,xz, Rfp,yz')
20103 format('# t, qp, qfp, qf@p, <1/tp>, <Rep>, <Cd>, <Vr>')
20104 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p')
20105 format('# t, <Tp>, <Tf@p>, <Tp^2>, <Tf@p^2>, <Tp.up>, <Tp.vp>, <Tp.wp>,&
<Tp.uf@p>,<Tp.vf@p>, <Tp.wf@p>, <Tf@p.up>, <Tf@p.vp>, <Tf@p.wp>,&
<Tf@p.uf@p>,<Tf@p.vf@p>, <Tf@p.wf@p>,<1/tp>, <Nup>, <Tp.Tf@p>')

20106 format('# t, Rpx, Rpy, Rpz, Rp, Rf@px, Rf@py, Rf@pz, Rf@p, RTp, RTf@p, RTfvf, RvfTf')

20107 format('# t, Nclose, Nover, Ncol, 1/tcol, dqp, <overlap>/dp')

20201 format('# t, <theta>, <theta^2>, <uf.theta>, <vf.theta>, <wf.theta>')
20202 format('# t, <theta^2>, diss, prod, Sx, Sy, Sz, Kx, Ky, Kz')

end subroutine OPENCLOSE
