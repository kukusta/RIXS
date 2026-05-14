SUBROUTINE CHECK_Q
use, intrinsic :: iso_fortran_env
use m_bzmesh, only: qbmc
use m_functions, only: exit_on_error
use m_rixs, only: nrixs, rixs
use m_sdt
implicit none
! local vars
character(*), parameter :: srcname=" in CHECK_Q"
integer i1, i2, i3, isp
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
use, intrinsic :: iso_fortran_env
use m_bz
use m_bzmesh
use m_rixs, only: nrixs, rixs
implicit none
! local vars
integer jq(3), isp, ik, i, i1, i2, i3, iq
real(REAL64) v(3), q(3), x

do isp = 1, nrixs
  q(:) = rixs(isp)%q(:)
  ik = 0
  do i3 = 1, ndxyz(3)
    do i2 = 1, ndxyz(2)
      do i1 = 1, ndxyz(1)
        ik = ik + 1
        do i = 1, 3
          v(i)=(i1-1)*qbmc(i,1)+(i2-1)*qbmc(i,2)+(i3-1)*qbmc(i,3)+q(i) ! v=k+q
        enddo
        do i = 1, 3
          x = v(1) * rbmc(1, i) + v(2) * rbmc(2, i) + v(3) * rbmc(3, i)
          jq(i) = nint(x)
          jq(i) = mod(jq(i), ndxyz(i)) + 1
          if (jq(i) <= 0) jq(i) = jq(i) + ndxyz(i)
        enddo
        iq = (jq(3) - 1) * ndxyz(2) * ndxyz(1) + (jq(2) - 1) * ndxyz(1) + jq(1)
        rixs(isp)%ik2kq(ik) = iq
      enddo
    enddo
  enddo
enddo

  5 format(4x, 'i=', i0, 4x, 'x(i)=', f15.5, 4x, 'ik=', i0)
 10 format(4x, 'For ISPEC= ', i0, 4x, ' the k-mesh node closest to q-vector', &
           ' is (',3f14.8,')')
 15 format(4x, 'Maximum value of rixs(',i0,')%ik2ibz= ', i0, &
           ' is larger then nkibz= ', i0)
 20 format(4x, 'Maximum value of rixs(', i0, ')%ik2kq= ', i0, &
           ' is larger then nkbz= ', i0)

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

! character(4) rixs(isp)%Reg_BZ
! integer rixs(isp)%Reg_BZ_size(3), rixs(isp)%Reg_BZ_center_qbmc(3)
! real(REAL64) rixs(isp)%Reg_BZ_center(3)

SUBROUTINE CHECK_POINT_IN_BZ
use, intrinsic :: iso_fortran_env
use m_bzmesh, only: qbmc, ndxyz
use m_functions, only: exit_on_error
use m_constants, only: eps3
use m_rixs, only: nrixs, rixs
use m_sdt
implicit none
! local vars
character(*), parameter :: srcname=" in CHECK_POINT_IN_BZ"
integer i, i1, i2, i3, isp
logical does_not_lie_on_kmesh
real(REAL64) absx, absy, absz, temp(3, 3), y(3), det, x(3)

  3 format(   4x, 'Diff. for expansion coeffs for isp=', i0, ': ', 3f12.5)
  5 format(/, 4x, "Some Reg_BZ_center do not lie on the k-mesh. Possible k-points")
 10 format(4x, "ISPEC=", i0, 4x, "k-point is (", 3f8.5, " )")
 15 format(15x, "Reg_BZ_center in qmbc is (", 3f10.5, " )")
 20 format(14x, 'k-point (in r. l. u.)', 10x, 'len(q)', 4x, 'abs(k-k_input) in (1/a)')
 25 format(4x, 3f12.5, 4x, f12.5, 8x, f12.5)

does_not_lie_on_kmesh = .false.
do isp = 1, nrixs
  call dinv33(qbmc, 0, temp, det)
  y(:) = matmul(temp, rixs(isp)%Reg_BZ_center)
  absx = abs( y(1) - nint(y(1)) )
  absy = abs( y(2) - nint(y(2)) )
  absz = abs( y(3) - nint(y(3)) )
! q-vector does not lie on the k-mesh. Here I want to find
! the k-vector closest to q
  if (absx > eps3 .or. absy > eps3 .or. absz > eps3) then
    if (.not. does_not_lie_on_kmesh) then
      write(unit = output_unit, fmt = 5)
      does_not_lie_on_kmesh = .true.
    endif
    write(unit = output_unit, fmt = 10) isp, &
         matmul(ut, rixs(isp)%Reg_BZ_center(:))
    write(unit = output_unit, fmt = 15) y(:)
    write(unit = output_unit, fmt = 20)
    do i3 = floor(y(3)), ceiling(y(3))
      do i2 = floor(y(2)), ceiling(y(2))
        do i1 = floor(y(1)), ceiling(y(1))
          x(:) = matmul(qbmc, real((/i1, i2, i3/), kind = REAL64))
          write(unit = output_unit, fmt = 25) matmul(ut, x) &
               ,sqrt(dot_product(x, x))&
               ,sqrt(dot_product(rixs(isp)%Reg_BZ_center(:), x))
        enddo
      enddo
    enddo
  endif
  if (does_not_lie_on_kmesh) cycle ! no need to save noninteger y(:)
  do i = 1, 3
    rixs(isp)%Reg_BZ_center_qbmc(i) = nint(y(i)) + 1
    if (rixs(isp)%Reg_BZ_center_qbmc(i) <= 0) then
       rixs(isp)%Reg_BZ_center_qbmc(i) = &
       rixs(isp)%Reg_BZ_center_qbmc(i) + ndxyz(i)
    endif
  enddo
enddo

if (does_not_lie_on_kmesh) call exit_on_error(srcname)

END SUBROUTINE CHECK_POINT_IN_BZ

! Reg_BZ == 'p'
! Reg_BZ_center(3) must lie on l-mesh: Reg_BZ_center(:)=0.5 0.5 0.5
! Reg_BZ_size(1:3) - integer number of microcells centered at Reg_BZ_center
! it is better to set Reg_BZ_size(1:3) even but not larger then ndxyz(1:3)
! Reg_BZ_size(1:3) = 4 4 4 (ndiv = 161616) must be even
! then it would be parallelepiped of size +/-2 +/-2 +/-2 centered at 0.5 0.5 0.5
! Reg_BZ_center_qbmc(1:3) - coordinates of region_in_BZ_center in qbmc
! basis surface of parallelepiped is included

SUBROUTINE DATA_REGION_IN_BZ
use, intrinsic :: iso_fortran_env
use m_bz
use m_bzmesh
use m_functions, only: exit_on_error
use m_rixs, only: nrixs, rixs
implicit none
! local vars
character(*), parameter :: srcname=' in DATA_REGION_IN_BZ'
integer c(3), d(3), i1, i2, i3, ik, isp, j1, j2, j3
integer, external :: check_kpoint
logical wrong_parameter

 10 format(/, 4x, 'For ISPEC=', i0, &
          ' wrong number for Reg_BZ_center: ik= ', i0, '. (', 3i6, ')')

wrong_parameter = .false.
do isp = 1, nrixs
  do i3 = 1, 3
    if (rixs(isp)%Reg_BZ_size(i3) < 0) then
      write(unit = output_unit, fmt = '(a, i2, a, i0, a)') 'For SPEC=', isp, &
           ' set for the region_in_BZ_size(', i3, ') < 0 but should be >= 0.'
      wrong_parameter = .true.
    endif
    if (mod(rixs(isp)%Reg_BZ_size(i3), 2) > 0) then
      write(unit = output_unit, fmt = '(4x,5(a,i0))') 'For ISPEC=', isp, &
           ' Reg_BZ_size(', i3, ')=', rixs(isp)%Reg_BZ_size(i3), &
           ' but must be even.'
      wrong_parameter = .true.
    endif
    if (rixs(isp)%Reg_BZ_size(i3) > ndxyz(i3)) then
      write(unit = output_unit, fmt = '(4x, 5(a, i0))') 'For ISPEC=', isp, &
           ' Reg_BZ_size(',i3,')=', rixs(isp)%Reg_BZ_size(i3), &
           ' > ndxyz(',i3,')=', ndxyz(i3)
      wrong_parameter = .true.
    endif
  enddo
enddo
if (wrong_parameter) stop

do isp = 1, nrixs
  rixs(isp)%map_region(:) = 1
  if (rixs(isp)%Reg_BZ(1:1) == 'n') cycle
  if (.not. rixs(isp)%Reg_BZ(1:1) == 'p') cycle
  c(:) = rixs(isp)%Reg_BZ_center_qbmc(:)
  d(:) = rixs(isp)%Reg_BZ_size(:)
!  ik = (i3 - 1) * ndxyz(2) * ndxyz(1) + (i2 - 1) * ndxyz(1) + i1
!  if (ik <= 0 .or. ik > nkbz) then
!    write(unit = output_unit, fmt = 10) isp, ik, i1, i2, i3
!    call exit_on_error(srcname)
!  endif
!  jk = 0
  rixs(isp)%map_region(:) = 0
  do j3 = -d(3) / 2, d(3) / 2
    do j2 = -d(2) / 2, d(2) / 2
      do j1 = -d(1) / 2, d(1) / 2
        i3 = c(3) + j3
        i2 = c(2) + j2
        i1 = c(1) + j1
        if (i3 > ndxyz(3)) i3 = i3 - ndxyz(3)
        if (i3 <= 0) i3 = i3 + ndxyz(3)
        if (i2 > ndxyz(2)) i2 = i2 - ndxyz(2)
        if (i2 <= 0) i2 = i2 + ndxyz(2)
        if (i1 > ndxyz(1)) i1 = i1 - ndxyz(1)
        if (i1 <= 0) i1 = i1 + ndxyz(1)
        ik = (i3 - 1) * ndxyz(2) * ndxyz(1) + (i2 - 1) * ndxyz(1) + i1
        rixs(isp)%map_region(ik) = 1
      enddo
    enddo
  enddo
! write(*,*)'For ISPEC=',isp,' map is '
! write(*,*)'For ISPEC=',isp,' sum is ', sum(rixs(isp)%map_region(:))
! write(*,'(1000(15i4,/))')rixs(isp)%map_region(:)
enddo

END SUBROUTINE DATA_REGION_IN_BZ
