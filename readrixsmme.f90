SUBROUTINE READRIXSMME
use m_aux ! , only: i1, i2, i3, map
use m_bnd ! , only : npnt , nb , nopused , iopnum
use m_bz
use m_bzmesh, only: nkbz
use m_files, only: rim, mmefile, bns, bnsfile
use m_functions
use m_params ! , only: irixs, iabsorp, iprint, BUFFER_SIZE, check_kstar, debug_mode
use m_rixs
use m_time
use omp_lib
implicit none
! local vars
character(15), parameter :: srcname = ' in READRIXSMME'
character(BUFFER_SIZE) msg
complex(REAL64), parameter :: z0 = ( 0._REAL64, 0._REAL64 )
complex(REAL64), allocatable :: mme( :, :, :, : )
complex(REAL64) mme_kibz(3), mme_kbz(3), tensor(3, 3)
!$$$ integer, allocatable :: isp4k(:), kbz4k(:)
integer iok, isp, k, ia, iat, ib, iop, kbz, ia1, iat1, ian, iatn
integer(INT64) reclength, size
logical exist
real(REAL64) timeini, timefinal, op( 3, 3 )
! data for CHECK_KSTAR
character(8) named
integer idummy,ntxtd,ibz
real ef4

if ( read_rixfile ) return

if ( ( irixs == 0 .or. irixs == 10 ) .and. iabsorp == 0 ) then
  write( unit = output_unit, fmt = 5 )
  return
endif

if ( check_kstar .and. debug_mode ) then
  inquire ( file = bnsfile, exist = exist )
  if ( .not. exist ) then
    write( unit = output_unit, fmt = 9 ) trim( bnsfile )
    call exit_on_error( srcname )
  endif
  open( newunit = bns, file = bnsfile, action = 'read', &
        status = 'unknown', form = 'unformatted' ,iostat = iok )
  if( iok /= 0 ) stop '    Bnsfile cannot be opened, exiting.'

  read( bns, iostat = iok )named(:)
  if( iok /= 0 ) stop '    Cannot read named(:) from bnsfile.'
  if( named(8:8) == '+' ) then
    ntxtd = 4
  else
    ntxtd = 1
  endif
  nopused = -1
  high_sym_pnt_txt_dir(:) = ""
  read( unit = bns, iostat = iok) idummy, ef4, npnt, idummy, nb, &
      idummy, idummy, idummy, ibz, idummy, ( high_sym_pnt(:,k), &
      high_sym_pnt_txt_dir(k)(1:ntxtd), k = 1, idummy ), nopused, &
      iopnum(:)
  if( iok /= 0 ) nopused = -1
  if( ibz /= 0) then
      call exit_on_error( srcname // &
           ': RIXS does not work with extended BZ mesh.' )
  endif
endif

inquire ( file = mmefile, exist = exist )
if ( .not. exist ) then
  write( unit = output_unit, fmt = 9 ) trim(mmefile)
  call exit_on_error( srcname )
endif

  5 format( /, 4x, 'Joint DOS without absorption: RIXS mme are not read.' )
  9 format( 4x, 'File ', a, ' does not exist.' )

timeini = omp_get_wtime()
reclength = int( storage_size(z0), kind = INT64 ) * 3_INT64 * int( nb, &
      kind = INT64 ) * int( n_mme, kind = INT64 ) * int( max_nhsort, &
      kind = INT64 )
open( newunit = rim, file = mmefile, action = 'read', &
      status = 'old', form = 'unformatted', access = 'direct', &
      recl = reclength, iostat = iok )
if ( iok /= 0 ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) &
      ' mmefile cannot be opened, exiting.'
  call exit_on_error( srcname )
endif

inquire ( file = mmefile, size = size )
! check the size of mmefile. In main LMTO program rixemmefile must be
! overwritten so its size can be easy checked
if ( size /= reclength * int( npnt, kind = INT64 ) ) then
  write( unit = output_unit, fmt = 13 ) &
       reclength * int( npnt, kind = INT64 ), size
  call exit_on_error( srcname )
endif

 10 format( 4x, 'Operation #(', i0, '), its matrix' )
 11 format( 4x, 'RIXS mme for k=', i0 )
 12 format( 4x, 'iat=', i2, ', ia=', i2, ', ib=', i2, ', mme(1:3)=', &
            6f10.6 )
 13 format( 4x, "Wrong size of rixsmmefile. Expected size is ", i0, &
            " B but real size is ", i0, " B" )
 15 format( 4x, 'Start reading of matrix elements.' )
 20 format( 4x, 'Size of mmefile is ', i0, ' MB (', a,' B)' )

if ( iprint > 0 ) then
  write( unit = output_unit, fmt = 15 )
  write( unit = output_unit, fmt = 20 ) &
       max( 1, int( real( size, kind = REAL64 ) / &
                    real( 1024**2, kind = REAL64 ) ) ), &
       trim( int2string( size, ' ' ) )
endif

if ( check_kstar .and. debug_mode ) then
  call allocate( mme, 'mme array' // srcname, udim1 = 3, &
       udim2 = max_nhsort, udim3 = n_mme, udim4 = nb )
  do k = 1, npnt
    if ( k == 1 ) then
      read( unit = rim, rec = k, iostat = iok ) mme( : , : , : , : )
      if ( iok /= 0 ) then
        write( unit = output_unit, fmt = 35 ) k
        call exit_on_error( srcname )
      endif
    endif
    if( k > 1 ) then
      call g4ig( iopnum(k), op )
!      op = transpose( op )
      write( unit = output_unit , fmt = 10 ) iopnum(k)
      write( unit = output_unit , fmt = "( 3f15.10 )" ) op(1,:)
      write( unit = output_unit , fmt = "( 3f15.10 )" ) op(2,:)
      write( unit = output_unit , fmt = "( 3f15.10 )" ) op(3,:)
    endif
    write( unit = output_unit , fmt = 11 ) k
    do iat = 1, max_nhsort
      do ia = 1, n_mme
        do ib = 1, nb
          mme_kibz(:) = mme( 1:3, iat, ia, ib )
          mme_kbz(:) = mme_kibz(:)
          if ( k > 1 ) then
            if ( iopnum(k) > 0 ) then
              mme_kbz(:) = matmul( op, mme_kibz )
!              mme_kbz(:) = matmul( mme_kibz, op )
            else
              mme_kbz(:) = matmul( op, conjg( mme_kibz ) )
!              mme_kbz(:) = matmul( conjg( mme_kibz ), op )
!              mme_kbz(:) = matmul( mme_kibz, op )
            endif
          endif
          write( unit = output_unit , fmt = 12 ) iat, ia, ib, mme_kbz
          cycle
          call tensor_from_mme( mme_kbz , tensor )
          do isp = 1 , 3
            write( unit = output_unit , fmt = "(3(' (',2(f15.10),')'))") &
                 tensor(isp,:)
          enddo
       enddo
      enddo
    enddo
  enddo
  call deallocate( mme, 'mme' // srcname )
  call deallocate_global_arrays( success = .true. )
  stop '    STOP in READRIXSMME: CHECK_KSTAR is set'
endif

!$omp parallel default(shared) private( iok, isp, k, mme, iat1, iatn, ia1, ian )
!$omp do schedule( dynamic, 1 )
do  k = 1, npnt
  allocate( mme( 3, max_nhsort, n_mme, nb ), stat = iok, errmsg = msg )
  if ( iok /= 0 ) then
    call print_allocation_error('readrixsmme for k='//trim(int2string(k))&
         ,iok,msg)
    call deallocate_global_arrays
    stop
  endif
!$omp critical
  read( unit = rim, rec = k, iostat = iok ) mme( : , : , : , : )
!$omp end critical
  if ( iok /= 0 ) then
    write( unit = output_unit, fmt = 35 ) k
    call exit_on_error( srcname )
  endif
!$omp critical
  do isp = 1, nlmtdata
    iat1 = lmtdata(isp)%iat1
    iatn = lmtdata(isp)%iatn
    ia1  = lmtdata(isp)%ia1
    ian  = lmtdata(isp)%ian
    lmtdata(isp)%mme( : , iat1 : iatn, ia1 : ian, : , k ) = &
                 mme( : , : iatn - iat1 + 1, ia1 :ian, : )
  enddo
  call deallocate(mme)
!$omp end critical
enddo   ! k from IBZ
!$omp end do
!$omp end parallel

!$$$ if ( nopused > 1 ) then
!$$$  call allocate( isp4k , 'isp4k' // srcname, udim1 = nkbz * nlmtdata )
!$$$  call allocate( kbz4k , 'kbz4k' // srcname, udim1 = nkbz * nlmtdata )  
  do isp = 1, nlmtdata
    call allocate( lmtdata( isp )%mmebz, &
                   'lmtdata(' // trim( int2string(isp) ) // ')%mmebz' //     &
                   srcname, udim1 = 3, ldim2 = lmtdata( isp )%iat1,          &
                   udim2 =  lmtdata( isp )%iatn, ldim3 = lmtdata( isp )%ia1, &
                   udim3 = lmtdata( isp )%ian, udim4 = nb, udim5 = nkbz )
!$$$    isp4k( (isp - 1) * nkbz + 1 : isp * nkbz ) = isp
!$$$    do k = 1, nkbz
!$$$      isp4k( (isp - 1) * nkbz + k ) = isp
!$$$      kbz4k( (isp - 1) * nkbz + k ) = k
!$$$    enddo
  enddo
! $omp parallel default(shared) private( isp, kbz, iop, op, ib, ia, iat, &
! $omp&         mme_kibz, mme_kbz )
! $omp do schedule( dynamic, 1 )
  do kbz = 1, nkbz ! * nlmtdata
  ! do isp = 1 , nlmtdata
  !   do kbz = 1 , nkbz
    ! isp = isp4k(k)
    ! kbz = kbz4k(k)
    iop = ig4qibz( i1(kbz), i2(kbz), i3(kbz) )
    op( : , : ) = g( : , : , iop )
    do ib = 1, nb
!$$$      do ia = lmtdata(isp)%ia1, lmtdata(isp)%ian
      do ia = 1, n_mme
        if ( sum( map( ia, : ) ) == 0)cycle
        isp = lmtindex( ia )
        do iat = lmtdata( isp )%iat1, lmtdata( isp )%iatn
          mme_kibz(:) = lmtdata( isp )%mme( : , iat, ia, ib, ik2ibz(kbz) )
          ! if ( iop == 1 ) then
          !   mme_kbz(:) = mme_kibz(:)
          ! else
          !   if ( iopnum(iop) > 0 ) then
          !     mme_kbz(:) = matmul( op, mme_kibz )
          !   else
          !     mme_kbz(:) = matmul( op, conjg( mme_kibz ) )
          !   endif
          ! endif
          if ( iopnum(iop) > 0 ) then
            mme_kbz(:) = matmul( op, mme_kibz )
          else
            mme_kbz(:) = matmul( op, conjg( mme_kibz ) )
          endif
          lmtdata( isp )%mmebz( : , iat, ia, ib, kbz ) = mme_kbz(:)
        enddo  ! iat
      enddo    ! ia
    enddo      ! ib
  enddo        ! k
! $omp end do
! $omp end parallel
!$$$ endif   ! nopused > 1

 35 format( 4x, 'In READRIXSMME: cannot read from mmefile for k=', i0 )
 40 format( 4x, 'Matrix elements have been read.' )

timefinal = omp_get_wtime()
call set_readmme( timefinal - timeini )

if ( iprint > 0 ) write( unit = output_unit, fmt = 40 )

END SUBROUTINE READRIXSMME
