SUBROUTINE SET_ABSORP_RBZ(isp,ia,kpk)
use m_aux, only: lmtindex,i1,i2,i3
use m_bnd
use m_bz
use m_bzmesh, only: nkbz
use m_functions
use m_rixs
implicit none
! input
integer isp,ia
! output
real(REAL64) kpk(nb,rixs(isp)%nkrbz)
! local vars
character(BUFFER_SIZE) msg
complex(REAL64), parameter :: z0=(0._REAL64,0._REAL64)
complex(REAL64) mmek(3),I(3,3)
integer, allocatable :: irbz2bz(:)
integer ib0,kbz,krbz,kibz,igkibz,ib,iat,iat1,iatn,iok

if(ia < 1 .or. ia > n_mme)then
!$omp critical
  write(unit=output_unit,fmt="(4x,a,i0)")&
       'ERROR in SET_ABSORP_RBZ: wrong value ia=',ia
!$omp end critical
  call exit_on_error
endif

iat1 = lmtdata( lmtindex(ia) )%iat1
iatn = lmtdata( lmtindex(ia) )%iatn
ib0 = 0

! nkibz == maxval(ipqibz)
if( rixs(isp)%nkrbz > maxval(ipqibz) .and. rixs(isp)%nkrbz < nkbz ) then
  allocate( irbz2bz( rixs(isp)%nkrbz ), stat = iok, errmsg = msg)
  if( iok /= 0 ) call print_allocation_error( 'rbz2bz', iok, msg )
  irbz2bz(:)=0
  do krbz = 1, rixs(isp)%nkrbz
    do kbz = 1, nkbz
      if( rixs(isp)%ik2rbz(kbz) == krbz )then
        irbz2bz(krbz) = kbz
        exit
      endif
    enddo
  enddo
  if( minval( irbz2bz ) == 0 ) then
    write(unit=output_unit,fmt="(4x,a)") &
      'Some elements of irbz2bz array are not referenced.'
      call exit_on_error
  endif
  if( maxval( irbz2bz ) > nkbz ) then
    write(unit=output_unit,fmt="(4x,a)") &
      'Some elements of irbz2bz array are larger then nkbz.'
      call exit_on_error
  endif
endif

do krbz = 1, rixs(isp)%nkrbz
  if( rixs(isp)%ng == 1 ) then
! no symmetry
    igkibz = 1
    kbz = krbz
    kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
  elseif( rixs(isp)%ng == nopused ) then
! full symmetry
    igkibz=1
    kibz=krbz
  else
! partial symmetry
    kbz = irbz2bz(krbz)
    kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
    igkibz = ig4qibz( i1(kbz), i2(kbz), i3(kbz) )
  endif
  do ib = 1, nb
    do iat = iat1, iatn
      mmek(:)=lmtdata( lmtindex(ia) )%mme(:, iat, ia, ib, kibz)
      if(igkibz /= 1) mmek(:) = mme_trans(igkibz, mmek)
      call tensor_from_mme(mmek, I)
      kpk( ib0+ib, krbz ) = real(I(1,1)+I(2,2),kind=REAL64)
    enddo         ! iat
  enddo           ! ib
enddo
deallocate(irbz2bz,stat=iok)

kpk(:,:)=kpk(:,:)/cmplx( 8*(iatn - iat1 + 1), kind=REAL64)

END SUBROUTINE SET_ABSORP_RBZ
