SUBROUTINE set_absorp_rbz(isp, ia, kpk)
use m_aux, only: lmtindex
use m_bnd
use m_bz
use m_bzmesh, only: ndxyz, nkbz
use m_constants, only: eps6, z0
use m_functions
use m_rixs
implicit none
! input
integer isp, ia
! output
!vvv real(REAL64) kpk(nb, rixs(isp)%nkrbz)
real(REAL64) kpk(nb, nkbz)
! local vars
character(18), parameter :: srcname=' in SET_ABSORP_RBZ'
complex(REAL64) mmek(3), I(3, 3)
integer ib0, kbz, ib, iat, iat1, iatn , k1, k2, k3
real(REAL64) Ikbz

if (ia < 1 .or. ia > n_mme) then
!$omp critical
  write(unit = output_unit, fmt = "(4x, a, i0)" ) &
    'ERROR in SET_ABSORP_RBZ: wrong value ia=', ia
!$omp end critical
  call exit_on_error(srcname)
endif

iat1 = lmtdata(lmtindex(ia))%iat1
iatn = lmtdata(lmtindex(ia))%iatn
ib0 = 0
! to check consistensy
kpk( : , : ) = -1._REAL64

kbz = 0
!vvv krbz = 0
do k3 = 1, ndxyz(3)
  do k2 = 1, ndxyz(2)
    do k1 = 1, ndxyz(1)
      kbz = kbz + 1
!$$$      if ( nopused == 1 ) then
!$$$        krbz = krbz + 1
!$$$      else
!vvv      if ( rixs(isp)%ng == 1 ) then
!vvv        krbz = krbz + 1
!vvv      else
!vvv        if ( rixs(isp)%ig4q( k1, k2, k3 ) /= 1 ) cycle
!vvv        krbz = krbz + 1
!vvv        if ( krbz /= rixs(isp)%ik2rbz( kbz ) ) then
!vvv          write( unit = output_unit, fmt = "( 4x, a )") &
!vvv            'krbz /= rixs(isp)%ik2rbz( kbz )'
!vvv          call exit_on_error( srcname )
!vvv        endif
!$$$        endif
!vvv      endif
      do ib = 1, nb
!vvv         Ikrbz = 0._REAL64
        Ikbz = 0._REAL64
        do iat = iat1, iatn
!$$$          if ( nopused == 1 ) then
!$$$            mmek(:) = lmtdata( lmtindex(ia) )%mme( : , iat , ia , ib , kbz )
!$$$          else
          mmek(:) = lmtdata(lmtindex(ia))%mmebz(:, iat, ia, ib, kbz)
!$$$          endif
          call tensor_from_mme(mmek, I)
          Ikbz = Ikbz + real(I(1, 1) + I(2, 2), kind = REAL64)
        enddo
        kpk(ib0 + ib, kbz) = Ikbz
!vvv        if (kpk(ib0 + ib, krbz) < -0.5 ) then
!vvv          kpk( ib0 + ib, krbz ) = Ikrbz
!vvv        else
!vvv          if ( abs( Ikrbz - kpk( ib0 + ib, krbz ) ) > eps6 ) then
!vvv            ! write( unit = output_unit , fmt = "( 4x , a )") &
!vvv            !   'ERROR' // srcname // ': Ikrbz /= kpk( ib0 + ib , krbz )'
!vvv            write( unit = output_unit, fmt = "( 4x, 'For kbz= ', i0 &
!vvv                 , ' krbz= ', i0, ' new mme= ', f0.10, ' old mme= ' &
!vvv                 , f0.10 )" ) kbz, krbz, Ikrbz, kpk( ib0 + ib, krbz )
!vvv            ! call exit_on_error ( srcname )
!vvv          endif
!vvv        endif
      enddo
    enddo
  enddo
enddo

if (minval(kpk) < -0.5_REAL64) then
  write( unit = output_unit, fmt = "(4x, 3(a, i0))") &
    'Some elements in kpk array for isp=', isp, &
    ', ia=', ia, ' have not been set.'
  call exit_on_error(srcname)
endif

kpk(:,:) = kpk(:,:) / cmplx(8 * (iatn - iat1 + 1), kind = REAL64)

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
!         if ( abs ( Ikrbz - kpk( ib0+ib , krbz ) ) > eps6 ) then
!           write( unit = output_unit , fmt = "( 4x , a )") &
!             'ERROR' // srcname // ': Ikrbz /= kpk( ib0+ib , krbz )'
!           call exit_on_error ( srcname )
!         endif
!       endif
!     enddo         ! iat
!   enddo           ! ib
! enddo
