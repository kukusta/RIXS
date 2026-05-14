SUBROUTINE KPMESH
use m_bnd
use m_bz
use m_bzmesh
use m_functions
use m_params, only: BUFFER_SIZE, iprint, nbi, nbf
use m_rixs, only: nrixs, rixs
implicit none
character(10), parameter :: srcname=' in KPMESH'
integer ib, iop, is_g4q, ng, k1, k2, k3, iprint_local, isp, kbz, krbz
real(REAL64), parameter :: eps=1.e-5_REAL64
real(REAL64), allocatable :: g_temp( : , : , : )
real(REAL64) op( 3, 3 ), dpi, q2(3), x, norm2, start_ref, finish_ref

if ( iprint > 0 ) then
  write( unit = output_unit , fmt = "( 4x , a )" ) 'Symmetry for RIXS.'
else
  write( unit = output_unit , fmt = "()" )
endif
iprint_local = -1

! get matrices for symmetry operations
ng = nopused
call allocate( g, 'g' // srcname, udim1 = 3, udim2 = 3, udim3 = ng, &
     nullify = .true. )
! allocate(g(3,3,ng),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('g',iok,msg)
do iop = 1, ng
  call g4ig( iopnum(iop), op )
  g( : , : , iop ) = op( : , : )
enddo
if ( iprint > 0 ) then
  write( unit = output_unit, fmt = "( 4x, a, / , 4x, a )" ) &
       'Symmetry for RIXS has been defined.', 'K-mesh for IBZ.'
endif
! ipq should be allocated outside of kmesh
if ( .not. allocated( ipq ) ) then
  call allocate( ipq, 'ipq' // srcname, udim1 = ndxyz(1), &
       udim2 = ndxyz(2), udim3 = ndxyz(3) )
  ! allocate(ipq(ndxyz(1),ndxyz(2),ndxyz(3)),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('ipq in KPMESH',iok,msg)
else
  write( unit = output_unit, fmt = "( 4x, a )" ) &
    'In KPMESH: ipq array is allocated. Why?'
  call exit_on_error( srcname )
endif

! get ipq and ig4q (set is_g4q=1) for IBZ, maximal symmetry from LMTO
is_g4q = 0

!$$$ if ( ng > 1 ) then
  is_g4q = 1
  if ( .not. allocated( ig4q ) ) then
    call allocate( ig4q, 'ig4q' // srcname, udim1 = ndxyz(1), &
            udim2 = ndxyz(2), udim3 = ndxyz(3) )
  ! allocate(ig4q(ndxyz(1),ndxyz(2),ndxyz(3)),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('ig4q in KPMESH',iok,msg)
  else
    write( unit = output_unit, fmt = "( 4x, a )" ) &
         'ig4q array is allocated. Why?'
    call exit_on_error( srcname )
  endif

  if ( .not. allocated( ipqibz ) ) then
    call allocate( ipqibz, 'ipqibz ' // srcname, udim1 = ndxyz(1), &
         udim2 = ndxyz(2), udim3 = ndxyz(3) )
  ! allocate(ipqibz(ndxyz(1),ndxyz(2),ndxyz(3)),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('ipqibz in KPMESH',iok,msg)
  else
    write( unit = output_unit, fmt = "( 4x, a )" ) &
      'ipqibz array is allocated. Why?'
    call exit_on_error( srcname )
  endif

  if ( .not. allocated( ig4qibz ) ) then
    call allocate( ig4qibz, 'ig4qibz' // srcname, udim1 = ndxyz(1), &
         udim2 = ndxyz(2), udim3 = ndxyz(3) )
  ! allocate(ig4qibz(ndxyz(1),ndxyz(2),ndxyz(3)),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('ig4qibz in KPMESH',iok,msg)
  else
    write( unit = output_unit, fmt = "( 4x, a )" ) &
      'ig4qibz array is allocated. Why?'
    call exit_on_error( srcname )
  endif
!$$$ endif  !  ng > 1

call kmesh( rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
       ibzext, is_idold, is_avtet, iprint_local )

! several extra checks
if ( nkbz /= ndxyz(1) * ndxyz(2) * ndxyz(3) ) then
  write( unit = output_unit, fmt = "( 4x , a )" ) &
    'In KPMESH: nkbz != n1 * n2 * n3'
  call exit_on_error( srcname )
endif
if ( nkibz /= npnt ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) 'In KPMESH: nkibz != npnt.'
  call exit_on_error( srcname )
endif

!$$$ if ( ng > 1 ) then
  call allocate( ik2ibz, 'ik2ibz' // srcname, udim1 = nkbz )
! allocate(ik2ibz(nkbz),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('ik2ibz',iok,msg)
  forall( k1 = 1 : ndxyz(1), k2 = 1 : ndxyz(2), k3 = 1 : ndxyz(3) )
    ik2ibz( ( k3 - 1 ) * ndxyz(2) * ndxyz(1) + ( k2 - 1 ) * ndxyz(1) + k1 &
          ) = ipq( k1, k2, k3 )
  end forall
  ipqibz( : , : , : ) = ipq( : , : , : )
  ig4qibz( : , : , : ) = ig4q( : , : , : )
!$$$ endif

if ( iprint > 0 ) &
  write( unit = output_unit, fmt = "( 4x, a )" ) 'K-mesh for IBZ has been set.'

  5 format( 4x, 'ERROR in KPMESH: for iop=', i0, ' op. matrix is degenerate.' )
 10 format( /, 4x, 'For ispec=', i0, ' no symmetry operations found ( ng = 0 ).' )

call cpu_time( start_ref )
if ( iprint > 0 ) &
  write( unit = output_unit, fmt = "( 4x, a )" ) &
       'Start referencing arrays for all RIXS spectra.'

!$$$ if ( ng > 1 )
call allocate( g_temp, 'g_temp' // srcname, udim1 = 3, udim2 = 3, udim3 = nopused )

! allocate(g_temp(3,3,nopused),stat=iok,errmsg=msg)
! if(iok /= 0) &
!   call print_allocation_error('g_temp',iok,msg)

do isp = 1, nrixs
  !$$$if ( ng == 1 )then
!$$$    rixs(isp)%ng = 1
!$$$    is_g4q = 0
!$$$    rixs(isp)%iopnum(1) = 1
!$$$    call kmesh( rbas , ndxyz(1) , ndxyz(2) , ndxyz(3) , g , ng , is_g4q , &
!$$$         is_wgt , ibzext , is_idold , is_avtet , iprint_local )
!$$$    rixs(isp)%nkrbz = nkibz
!$$$    rixs(isp)%ntrbz = ntibz
!$$$    rixs(isp)%fnorm = 2._REAL64 * pcellvol * dpi()
!$$$    call allocate( rixs(isp)%itetr , 'rixs(' // int2string(isp) // ')%itetr' &
!$$$         // srcname , ldim1 = 0 , udim1 = 4 , udim2 = rixs(isp)%ntrbz )
!$$$    call allocate( rixs(isp)%idold , 'rixs(' // int2string(isp) // ')%idold' &
!$$$         //srcname , udim1 = 4 , udim2 = rixs(isp)%ntrbz )
!$$$    rixs(isp)%itetr( : , : ) = itetr( : , : )
!$$$    rixs(isp)%idold( : , : ) = idold( : , : )
!$$$  else
    rixs( isp )%iopnum(:) = 0
  ! allocate(rixs(isp)%ipq(ndxyz(1),ndxyz(2),ndxyz(3))&
  !         ,stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%ipq',iok,msg)
    if ( rixs( isp )%use_symmetry( 1 : 1 ) == 'n' ) then
      rixs( isp )%ng = 1
      g_temp( : , : , 1 ) = g( : , : , 1 )
      rixs( isp )%iopnum( 1 ) = iopnum( 1 )
    else
      rixs(isp)%ng = 0
      do iop = 1, nopused
        norm2 = dot_product( rixs(isp)%q, rixs(isp)%q )
        op( : , : ) = g( : , : , iop )
        if ( maxval( abs(op) ) < eps ) then
          write( unit = output_unit, fmt = 5 )iop
          call exit_on_error( srcname )
        endif
        q2(:) = matmul( op, rixs(isp)%q )
        x = dot_product( rixs(isp)%q, q2 )
        if ( norm2 - abs(x) < eps ) then
          if ( norm2 - x > eps .and. rixs( isp )%use_symmetry( 1:1 ) == 'p' ) cycle
          rixs( isp )%ng = rixs( isp )%ng + 1
          g_temp( : , : , rixs( isp )%ng ) = op( : , : )
          rixs( isp )%iopnum( rixs( isp )%ng ) = iopnum( iop )
        endif
      enddo
      if ( rixs( isp )%ng == 0 ) then
        write( unit = output_unit, fmt = 10 ) isp
        call exit_on_error( srcname )
      endif
    endif   !  rixs(isp)%use_symmetry(1:1) == 'n'

    call allocate( rixs(isp)%g, 'rixs(' // int2string(isp) // ')%g' // &
         srcname, udim1 = 3, udim2 = 3, udim3 = rixs( isp )%ng )
    rixs( isp )%g( : , : , : ) = g_temp( : , : , : rixs( isp )%ng )

  ! allocate(rixs(isp)%g(3,3,rixs(isp)%ng),stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%g',iok,msg)

    is_g4q = 0
    if( rixs(isp)%ng > 1 ) is_g4q = 1
    call kmesh( rbas, ndxyz(1), ndxyz(2), ndxyz(3), rixs( isp )%g, &
         rixs( isp )%ng, is_g4q, is_wgt, ibzext, is_idold, is_avtet, &
         iprint_local )
    rixs( isp )%nkrbz = nkibz
    rixs( isp )%ntrbz = ntibz
!    rixs(isp)%fnorm=omg48*real(rixs(isp)%ng,kind=REAL64)/2._REAL64
    rixs( isp )%fnorm = 2._REAL64 * pcellvol * dpi()
    call allocate( rixs( isp )%itetr, 'rixs(' // trim( int2string(isp) ) // &
         ')%itetr' // srcname, ldim1 = 0, udim1 = 4, udim2 = rixs( isp )%ntrbz )
  ! allocate(rixs(isp)%itetr(0:4,rixs(isp)%ntrbz),stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%itetr',iok,msg)
    call allocate( rixs( isp )%idold, 'rixs(' // trim( int2string(isp) ) // &
         ')%idold' // srcname, udim1 = 4, udim2 = rixs( isp )%ntrbz )
  ! allocate(rixs(isp)%idold(4,rixs(isp)%ntrbz),stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%idold',iok,msg)
    rixs( isp )%itetr( : , : ) = itetr( : , : )
    rixs( isp )%idold( : , : ) = idold( : , : )

    if( rixs( isp )%ng > 1 ) then
      call allocate( rixs( isp )%ig4q, 'rixs(' // trim( int2string(isp) ) // &
           ')%ig4q' // srcname, udim1 = ndxyz(1), udim2 = ndxyz(2), &
           udim3 = ndxyz(3) )
      call allocate( rixs( isp )%ipq, 'rixs(' // trim( int2string(isp) ) // &
           ')%ipq' // srcname, udim1 = ndxyz(1), udim2 = ndxyz(2), &
           udim3 = ndxyz(3) )
      call allocate( rixs( isp )%ik2rbz, 'rixs(' // int2string(isp) // &
           ')%ik2rbz' // srcname, udim1 = nkbz )
  ! allocate(rixs(isp)%ik2rbz(nkbz),stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%ik2rbz',iok,msg)
      rixs( isp )%ipq( : , : , : ) = ipq( : , : , : )
      rixs( isp )%ig4q( : , : , : ) = ig4q( : , : , : )
      kbz = 0
      rixs( isp )%ik2rbz(:) = 0
      do k3 = 1, ndxyz(3)
        do k2 = 1, ndxyz(2)
          do k1 = 1, ndxyz(1)
            kbz = kbz + 1
            rixs( isp )%ik2rbz( kbz ) = ipq( k1, k2, k3 )
          enddo
        enddo
      enddo
! extra check
      if ( maxval( rixs( isp )%ik2rbz ) > rixs( isp )%nkrbz ) then
        write( unit = output_unit, fmt = "( 4x, a, 3( i0, a ) )" ) &
             'ERROR in KPMESH: for isp=', isp, &
             ' maxval(rixs(', isp, ')%ik2rbz) > rixs(', isp, ')%nkrbz.'
        call exit_on_error( srcname )
      endif
      if ( minval( rixs( isp )%ik2rbz ) == 0 ) then
        write( unit = output_unit, fmt = "( 4x, 2( a, i0 ), a )" ) &
             'ERROR in KPMESH: for isp=', isp, &
             ' minval(rixs(', isp, ')%ik2rbz) == 0.'
        call exit_on_error( srcname )
      endif

      call allocate( rixs( isp )%ikbz2rbz, 'rixs(' // trim( int2string(isp) ) &
           // ')%ikbz2rbz' // srcname, udim1 = rixs( isp )%nkrbz, nullify = .true. )

  ! allocate( rixs(isp)%ikbz2rbz( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%ikbz2rbz',iok,msg)
  ! rixs(isp)%ikbz2rbz(:) = 0

      do krbz = 1, rixs(isp)%nkrbz
        do kbz = 1, nkbz
          if ( rixs( isp )%ik2rbz( kbz ) == krbz ) then
            rixs( isp )%ikbz2rbz( krbz ) = kbz
            exit
          endif
        enddo
      enddo
      if ( minval( rixs( isp )%ikbz2rbz ) <= 0 ) then
        write( unit = output_unit, fmt = "( 4x, a )" ) &
          'Some elements of ikbz2rbz array are not referenced.'
        call exit_on_error( srcname )
      endif
      if( maxval( rixs( isp )%ikbz2rbz ) > nkbz ) then
        write( unit = output_unit, fmt = "( 4x, a )" ) &
          'Some elements of ikbz2rbz array are larger then nkbz.'
        call exit_on_error( srcname )
      endif

    endif  !  rixs(isp)%ng > 1

  !$$$ endif  !  ng == 1

  call allocate( rixs(isp)%ik2kq, 'rixs(' // trim (int2string(isp) ) &
       // ')%ik2kq' // srcname, udim1 = nkbz )
  ! allocate(rixs(isp)%ik2kq(nkbz),stat=iok,errmsg=msg)
  ! if(iok /= 0) &
  !   call print_allocation_error(&
  !        'rixs('//trim(int2string(isp))//')%ik2kq',iok,msg)
enddo     ! isp

if ( ng > 1 ) call deallocate( g_temp, 'g_temp' // srcname )
! deallocate(g_temp,stat=iok)

call cpu_time( finish_ref )
if ( iprint > 0 ) &
  write( unit = output_unit, fmt = "( 4x, a, f0.4, ' s' )" ) &
       'Finished referencing arrays for all RIXS spectra. Time required: ', &
       finish_ref - start_ref

! check if all q vectors lie on q mesh
call check_q
! get arrays for referencing band energies. ik2kq( : , : ) and ik2kqibz( : , : )
! contain k-point numbers for all q-vectors and are used in OpenMP loop.
call get_ik2kq

! Restore data for full symmetry from LMTO. It is used in set_ncf_ncl.f90
is_g4q = 0
call kmesh( rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
     ibzext, is_idold, is_avtet, iprint_local )
call set_nbi_nbf

if ( nbi(1) > nbi(2) ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) 'Error: nbi(1) > nbi(2).'
  call exit_on_error( srcname )
endif
if ( nbf(1) > nbf(2) ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) 'Error: nbf(1) > nbf(2).'
  call exit_on_error( srcname )
endif
if ( iprint > 0 ) &
  write( unit = output_unit, fmt = "( 4x, a )" ) 'K-mesh has been defined.'

call allocate( ebz, srcname, udim1 = nb, udim2 = nkbz )
forall( ib = 1 : nb, k1 = 1 : nkbz ) ebz( ib, k1 ) = e( ib, ik2ibz(k1) )

! call check_sym

END SUBROUTINE KPMESH

SUBROUTINE CHECK_SYM
use, intrinsic :: ISO_FORTRAN_ENV
use m_bnd
implicit none

complex(REAL64), parameter :: z0 = ( 0._REAL64, 0._REAL64 )
complex(REAL64), parameter :: zi = ( 0._REAL64, 1._REAL64 )
complex(REAL64) vec(3), vec2(3), rm(3,3), rm2(3,3), rm_iop(3,3)
integer iop, ix, iy
real m( 3, 3 ), vec_re(3), vec_im(3)
real(REAL64) op(3,3)

! symmetry check: initilaze random matrix and apply symmetry operations
! from iopnum. r_matr is random complex matrix with real numbers on its
! main diagonal
call random_seed()
call random_number( vec_re )
call random_number( vec_im )
vec( : ) = cmplx( vec_re, kind = REAL64 ) + zi * cmplx( vec_im, kind = REAL64 )
call tensor_from_mme( vec, rm )

write( unit = output_unit, fmt = "( /, 4x, a )" ) &
  'The first algorithm for symmetrization'
rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
    rm2( : , : ) = rm2( : , : ) + matmul( transpose(op), matmul( rm, op ) )
  else
    rm2( : , : ) = rm2( : , : ) + &
                   matmul( transpose(op), matmul( conjg(rm), op ) )
  endif
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit, &
       fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit , fmt = "( / )")

write( unit = output_unit , fmt = "( 4x, a )") &
  'The second algorithm for symmetrization'
rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
!$$$    vec2( : ) = matmul( transpose(op) , vec )
    vec2( : ) = matmul( op, vec )
  else
!$$$    vec2( : ) = matmul( transpose(op) , conjg( vec ) )
    vec2( : ) = matmul( op, conjg( vec ) )
  endif
  call tensor_from_mme( vec2, rm_iop )
  rm2( : , : ) = rm2( : , : ) + rm_iop( : , : )
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3 
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit &
       , fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit, fmt = "( / )")

write( unit = output_unit, fmt = "( 4x, a )") &
  'The third algorithm for symmetrization'
call random_number( m )
rm( : , : ) = m( : , : )
call random_number( m )
rm( : , : ) = rm( : , : ) + zi * m( : , : )

rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
    rm2( : , : ) = rm2( : , : ) + matmul( transpose(op), matmul( rm, op ) )
  else
    rm2( : , : ) = rm2( : , : ) &
                 + matmul( transpose(op), matmul( conjg(rm), op ) )
  endif
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3 
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit &
       , fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit, fmt = "( / )")

END SUBROUTINE CHECK_SYM
