SUBROUTINE CHECK_Q
use, intrinsic :: iso_fortran_env
use m_bzmesh, only: qbmc
use m_functions, only: exit_on_error
use m_params, only: debug_mode
use m_rixs, only: nrixs, rixs
use m_sdt
implicit none
! local vars
character(11), parameter :: srcname=" in CHECK_Q"
integer i1,i2,i3,isp
logical printed
real(REAL64), parameter :: eps=1.e-3_REAL64
real(REAL64) absx,absy,absz,temp(3,3),y(3),det,x(3)

! 'print_ut' is called twice (here and in printdata.f90) since if code
! stops here then information about ut, ut1 and so on would not be printed
call print_ut

  3 format( 4x, 'Diff. for expansion coeffs for isp=', i0, ': ', 3f12.5)
  5 format(/,4x,"Some q vectors do not lie on the k-mesh. Possible q vectors")
 10 format(4x,"ISPEC=",i0,4x,"q_input in r.l.u. is (",3f10.5," )")
 15 format(15x,"q_input in qmbc is (",3f10.5," )")
 20 format(14x,'q (in r.l.u.)',15x,'len(q)',4x,'abs(q-q_input) in (1/a)')
 25 format(4x,3f12.5,4x,f12.5,8x,f12.5)

printed = .false.
do isp = 1, nrixs
  call dinv33( qbmc, 0, temp, det )
  y(:) = matmul( temp, rixs(isp)%q(:) )

  absx = abs( y(1) - nint( y(1) ) )
  absy = abs( y(2) - nint( y(2) ) )
  absz = abs( y(3) - nint( y(3) ) )
  if ( debug_mode ) &
    write( unit = output_unit, fmt = 3 ) isp, absx, absy, absz

! q-vector does not lie on the k-mesh. Here I want to find
! the k-vector closest to q
  if(absx > eps .or. absy > eps .or. absz > eps)then
    if(.not. printed)then
      write(unit=output_unit,fmt=5)
      printed=.true.
    endif

    write(unit=output_unit,fmt=10)isp,matmul(ut,rixs(isp)%q)
    write(unit=output_unit,fmt=15)y(:)
    write(unit=output_unit,fmt=20)
    do i3=floor(y(3)),ceiling(y(3))
      do i2=floor(y(2)),ceiling(y(2))
        do i1=floor(y(1)),ceiling(y(1))
          x(:)=matmul(qbmc,real((/i1,i2,i3/),kind=REAL64))
          write(unit=output_unit,fmt=25)matmul(ut,x)&
               ,sqrt(dot_product(x,x))&
               ,sqrt(dot_product(rixs(isp)%q,x))
        enddo
      enddo
    enddo
  endif
enddo

if(printed) call exit_on_error(srcname)

END SUBROUTINE CHECK_Q

SUBROUTINE GET_IK2KQ
! For each point k of the whole BZ ik2kq(k) returns a pointer to the point k+q

! If ik2ibz(1)>0 then ik2kqibz(k) is also filled which for a point k returns a
! pointer to the point k+q from the irreducible part of BZ

! If ik2ibz(1)=0 then other elements of ik2ibz and the whole array ik2kqibz
! are not accessed
use, intrinsic :: iso_fortran_env
use m_bz
use m_bzmesh
use m_rixs, only: nrixs, rixs
implicit none
! local vars
integer jq(3), isp, ik, i, i1, i2, i3, iq
real(REAL64) v(3), q(3), x

do isp=1,nrixs
  q(:)=rixs(isp)%q(:)
  ik=0
  do i3=1,ndxyz(3)
    do i2=1,ndxyz(2)
      do i1=1,ndxyz(1)
        ik=ik+1
        do i=1,3
          v(i)=(i1-1)*qbmc(i,1)+(i2-1)*qbmc(i,2)+(i3-1)*qbmc(i,3)+q(i) ! v=k+q
        enddo
        do i=1,3
          x=v(1)*rbmc(1,i)+v(2)*rbmc(2,i)+v(3)*rbmc(3,i)
          jq(i)=nint(x)
          jq(i)=mod(jq(i),ndxyz(i))+1
          if(jq(i) <= 0)jq(i)=jq(i)+ndxyz(i)
        enddo                       ! i
        iq=(jq(3)-1)*ndxyz(2)*ndxyz(1)+(jq(2)-1)*ndxyz(1)+jq(1)
        rixs(isp)%ik2kq(ik)=iq
      enddo                         ! i1
    enddo                           ! i2
  enddo                             ! i3
enddo                               ! isp

  5 format(4x,'i=',i0,4x,'x(i)=',f15.5,4x,'ik=',i0)
 10 format(4x,'For ISPEC= ',i0,4x,' the k-mesh node closest to q-vector',&
          ' is (',3f14.8,')')
 15 format(4x,'Maximum value of rixs(',i0,')%ik2ibz= ',i0&
          ,' is larger then nkibz= ',i0)
 20 format(4x,'Maximum value of rixs(',i0,')%ik2kq= ',i0&
          ,' is larger then nkbz= ',i0)

END SUBROUTINE GET_IK2KQ

SUBROUTINE PRINT_UT
use, intrinsic :: iso_fortran_env
use m_bnd, only: rbas, qbas
use m_bzmesh, only: rbmc, qbmc
use m_params, only: iprint
use m_sdt, only: ut, ut1
implicit none

 10 format(/,4x,'Direct and reciprocal space',/,17x,'U matrix',28x&
          ,'Inverse U matrix')
 15 format(1x,3f12.6,3x,':',3f12.6)
 20 format(19x,'RBAS',36x,'QBAS')
 25 format(1x,3f12.6,3x,':',3f12.6)
 30 format(19x,'rbmc',36x,'qbmc')
 35 format(1x,3f12.6,3x,':',3f12.6)

if(iprint > 10)then
  write(unit=output_unit,fmt=10)
  write(unit=output_unit,fmt=15)ut(1,:),ut1(1,:)
  write(unit=output_unit,fmt=15)ut(2,:),ut1(2,:)
  write(unit=output_unit,fmt=15)ut(3,:),ut1(3,:)
  write(unit=output_unit,fmt=20)
  write(unit=output_unit,fmt=25)rbas(1,:),qbas(1,:)
  write(unit=output_unit,fmt=25)rbas(2,:),qbas(2,:)
  write(unit=output_unit,fmt=25)rbas(3,:),qbas(3,:)
  write(unit=output_unit,fmt=30)
  write(unit=output_unit,fmt=35)rbmc(1,:),qbmc(1,:)
  write(unit=output_unit,fmt=35)rbmc(2,:),qbmc(2,:)
  write(unit=output_unit,fmt=35)rbmc(3,:),qbmc(3,:)
endif

END SUBROUTINE PRINT_UT

