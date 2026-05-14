SUBROUTINE SET_ABSORP_RBZ( isp, ia, kpk )
use m_aux, only: lmtindex
use m_bnd
use m_bz
use m_bzmesh, only: ndxyz
use m_functions
use m_rixs
implicit none
! input
integer isp, ia
! output
real(REAL64) kpk( nb, rixs(isp)%nkrbz )
! local vars
character(18), parameter :: srcname=' in SET_ABSORP_RBZ'
complex(REAL64), parameter :: z0 = ( 0._REAL64, 0._REAL64 )
complex(REAL64) mmek(3), I( 3, 3 )
integer ib0, kbz, krbz, ib, iat, iat1, iatn , k1, k2, k3
real(REAL64), parameter :: eps = 1.e-9_REAL64
real(REAL64) Ikrbz

if ( ia < 1 .or. ia > n_mme ) then
!$omp critical
  write( unit = output_unit, fmt = "( 4x, a, i0 )" ) &
    'ERROR in SET_ABSORP_RBZ: wrong value ia=', ia
!$omp end critical
  call exit_on_error( srcname )
endif

iat1 = lmtdata( lmtindex(ia) )%iat1
iatn = lmtdata( lmtindex(ia) )%iatn
ib0 = 0
! to check consistensy
kpk( : , : ) = -1._REAL64

kbz = 0
krbz = 0
do k3 = 1, ndxyz(3)
  do k2 = 1, ndxyz(2)
    do k1 = 1, ndxyz(1)
      kbz = kbz + 1
!$$$      if ( nopused == 1 ) then
!$$$        krbz = krbz + 1
!$$$      else
      if ( rixs(isp)%ng == 1 ) then
        krbz = krbz + 1
      else
        if ( rixs(isp)%ig4q( k1, k2, k3 ) /= 1 ) cycle
        krbz = krbz + 1
        if ( krbz /= rixs(isp)%ik2rbz( kbz ) ) then
          write( unit = output_unit, fmt = "( 4x, a )") &
            'krbz /= rixs(isp)%ik2rbz( kbz )'
          call exit_on_error( srcname )
        endif
!$$$        endif
      endif

      do ib = 1, nb
        Ikrbz = 0._REAL64
        do iat = iat1, iatn
!$$$          if ( nopused == 1 ) then
!$$$            mmek(:) = lmtdata( lmtindex(ia) )%mme( : , iat , ia , ib , kbz )
!$$$          else
          mmek(:) = lmtdata( lmtindex(ia) )%mmebz( : , iat, ia, ib, kbz )
!$$$          endif
          call tensor_from_mme( mmek, I )
          Ikrbz = Ikrbz + real( I( 1, 1 ) + I( 2, 2 ), kind = REAL64 )
        enddo         ! iat
        if ( kpk( ib0 + ib, krbz ) < -0.5 ) then
          kpk( ib0 + ib, krbz ) = Ikrbz
        else
          if ( abs( Ikrbz - kpk( ib0 + ib, krbz ) ) > eps ) then
            ! write( unit = output_unit , fmt = "( 4x , a )") &
            !   'ERROR' // srcname // ': Ikrbz /= kpk( ib0 + ib , krbz )'
            write( unit = output_unit, fmt = "( 4x, 'For kbz= ', i0 &
                 , ' krbz= ', i0, ' new mme= ', f0.10, ' old mme= ' &
                 , f0.10 )" ) kbz, krbz, Ikrbz, kpk( ib0 + ib, krbz )
            ! call exit_on_error ( srcname )
          endif
        endif
      enddo           ! ib
    enddo             ! k3
  enddo               ! k2
enddo                 ! k1

if ( minval( kpk ) < -0.5_REAL64 ) then
  write( unit = output_unit, fmt = "( 4x, 3 (a, i0 ) )") &
    'Some elements in kpk array for isp=', isp, &
    ', ia=', ia, ' have not been set.'
  call exit_on_error( srcname )
endif

kpk( : , : ) = kpk( : , : ) / cmplx( 8 * (iatn - iat1 + 1), kind = REAL64 )

END SUBROUTINE SET_ABSORP_RBZ

! nkibz == maxval(ipqibz)
! if( rixs(isp)%nkrbz > maxval(ipqibz) .and. rixs(isp)%nkrbz < nkbz ) then
!   allocate( irbz2bz( rixs(isp)%nkrbz ), stat = iok, errmsg = msg)
!   if( iok /= 0 ) call print_allocation_error( 'rbz2bz', iok, msg )
!   irbz2bz(:)=0
!   do krbz = 1, rixs(isp)%nkrbz
!     do kbz = 1, nkbz
!       if( rixs(isp)%ik2rbz(kbz) == krbz )then
!         irbz2bz(krbz) = kbz
!         exit
!       endif
!     enddo
!   enddo
!   if( minval( irbz2bz ) == 0 ) then
!     write(unit=output_unit,fmt="(4x,a)") &
!       'Some elements of irbz2bz array are not referenced.'
!       call exit_on_error
!   endif
!   if( maxval( irbz2bz ) > nkbz ) then
!     write(unit=output_unit,fmt="(4x,a)") &
!       'Some elements of irbz2bz array are larger then nkbz.'
!       call exit_on_error
!   endif
! endif

!   if( rixs(isp)%ng == 1 ) then
! ! no symmetry
!     igkibz = 1
!     kbz = krbz
!     kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
!   elseif( rixs(isp)%ng == nopused ) then
! ! full symmetry
!     igkibz=1
!     kibz=krbz
!   else
! ! partial symmetry
!     kbz = rixs(isp)%ikbz2rbz(krbz)
!     kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
!     igkibz = ig4qibz( i1(kbz), i2(kbz), i3(kbz) )
!   endif

! !$$$ do krbz = 1 , rixs(isp)%nkrbz
! do kbz = 1 , nkbz

!   ! if ( nopused == 1 ) then
!   !   kbz = krbz
!   !   kibz = kbz
!   !   igkibz = 1
!   ! else
!   !   if ( rixs(isp)%ng == 1 ) then
!   !     kbz = krbz
!   !     kibz = ik2ibz(kbz)
!   !     igkibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!   !   else
!   !     kbz = rixs(isp)%ikbz2rbz(krbz)
!   !     kibz = ipqibz( i1(kbz) , i2(kbz) , i3(kbz) )
!   !     igkibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!   !   endif
!   ! endif

!   if ( nopused == 1 ) then
!     kibz = kbz
!     krbz = kbz
!     igk2ibz = 1
!     igk2rbz = 1
!   else
!     if ( rixs(isp)%ng == 1 ) then
!       kibz = ik2ibz(kbz)
!       krbz = kbz
!       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       igk2rbz = 1
!     else
! ! kbz is in IBZ if ig4qibz(i1(kbz),i2(kbz),i3(kbz)) == 1
! ! kbz is in RBZ if rixs(isp)%ig4q(i1(kbz),i2(kbz),i3(kbz)) == 1
!       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       igk2rbz = rixs(isp)%ig4q( i1(kbz) , i2(kbz) , i3(kbz) )
!       ! if ( igk2ibz == 1 .and. igk2rbz == 1 ) then
!       ! endif
!       ! if ( igk2ibz /= 1 .and. igk2rbz == 1 ) then
!       ! endif
!       ! if ( igk2ibz == 1 .and. igk2rbz /= 1 ) then
!       ! endif
!       ! if ( igk2ibz /= 1 .and. igk2rbz /= 1 ) then
!       ! endif
!       kibz = ipqibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       krbz = rixs(isp)%ik2rbz(kbz)
!     endif
!   endif

!   do ib = 1 , nb
!     do iat = iat1 , iatn
!       mmek(:) = lmtdata( lmtindex(ia) )%mme( : , iat , ia , ib , kibz )

!       if ( igk2ibz /= 1 ) then
!         op( : , : ) = g( : , : , igk2ibz )
!         if ( iopnum(igk2ibz) > 0 ) then
!           mmek(:) = matmul( op , mmek )
!         else
!           mmek(:) = matmul( op , conjg( mmek ) )
!         endif
!       endif

!       if ( igk2rbz /= 1 ) then
!         op( : , : ) = transpose(rixs(isp)%g( : , : , igk2rbz ))
!         if ( iopnum(igk2rbz) > 0 ) then
!           mmek(:) = matmul( op , mmek )
!         else
!           mmek(:) = matmul( op , conjg( mmek ) )
!         endif
!       endif

!       call tensor_from_mme(mmek, I)
!       Ikrbz = real( I(1,1) + I(2,2) , kind = REAL64 )
!       if ( kpk( ib0+ib , krbz ) < -0.5 ) then
!         kpk( ib0+ib , krbz ) = Ikrbz
!       else
!         if ( abs ( Ikrbz - kpk( ib0+ib , krbz ) ) > eps ) then
!           write( unit = output_unit , fmt = "( 4x , a )") &
!             'ERROR' // srcname // ': Ikrbz /= kpk( ib0+ib , krbz )'
!           call exit_on_error ( srcname )
!         endif
!       endif
!     enddo         ! iat
!   enddo           ! ib
! enddo
