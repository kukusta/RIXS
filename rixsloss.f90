SUBROUTINE RIXSLOSS
use m_aux
use m_bnd
use m_bz
use m_bzmesh
use m_functions, only: exit_on_error, allocate, deallocate
use m_params
use m_results
use m_rixs
use m_time
use omp_lib
implicit none
! local vars
character(12), parameter :: srcname=' in RIXSLOSS'
character(BUFFER_SIZE) msg, spec_type
integer, parameter :: NTR_STEP = 20
integer iok, itr, isp, jbi, jbf, kibz, kqibz, kbz, kqbz, nCPUitr, hh, mm, krbz
real, allocatable :: a(:)
real(REAL64), allocatable :: ebi(:) , ebf(:), work(:)
real(REAL64), external :: simpson38
real(REAL64), parameter :: eps = 1.e-4_REAL64
real(REAL64) CPUini, CPUfinal, realini, realfinal, CPUitr_begin, CPUitr_end, &
             CPUitr, real_itr, real_itr2, real_itr_end, real_itr_begin, ss

if ( read_rixfile ) return

call cpu_time( CPUini )
realini = omp_get_wtime()

  5 format( /, 4x, a, ' spectra calculations.' )
 10 format( 4x, 'OpenMP do-loop in 1 thread for ', i0, &
            ' k-points in BZ (Joint DOS).' )
 11 format( 4x, 'OpenMP do-loop in ', i0, ' threads for ', i0, &
            ' k-points in BZ (Joint DOS).' )
 15 format( 4x, 'Joint DOS calculations for ispec=', i0, ' thread=', i0 )

! allocate(loss(nw,nrixs),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('loss',iok,msg)
! loss(:,:) = 0._REAL64

call allocate( loss, alias = 'loss' // srcname, udim1 = nw, udim2 = nrixs, &
     nullify = .true. )

if ( absorption_only ) then
  call set_CPU_RIXS( 0._REAL64 )
  call set_real_RIXS( 0._REAL64 )
  return
endif

if ( abs(irixs) /= 0 .and. abs(irixs) /= 10 ) then
  spec_type = 'RIXS loss'
else
  spec_type = 'Joint DOS'
endif
write( unit = output_unit, fmt = 5 ) trim( spec_type )

nCPUitr = 0
CPUitr = 0._REAL64
real_itr = 0._REAL64

if ( irixs == 0 ) then

!$omp parallel default(shared) private( a, jbi, jbf, kbz, krbz, ebi, ebf,&
!$omp&         itr, kibz, kqibz, kqbz, isp, CPUitr_begin, CPUitr_end,    &
!$omp&         real_itr_begin, real_itr_end, hh, mm, ss, work )

!$omp master

  if( omp_get_num_threads() == 1 )then
    write( unit = output_unit, fmt = 10 ) nkbz
  else
    write( unit = output_unit, fmt = 11 ) omp_get_num_threads(), nkbz
  endif

  if ( omp_get_nested() ) then
    write( unit = output_unit, fmt = " ( 4x, a ) ") &
      trim( spec_type ) // ': nested parallelism has been disabled'
    call omp_set_nested( .false. )
  endif

!$omp end master

!$omp do  schedule( dynamic, 1 )
! reduction( + : loss )

  do isp = 1, nrixs

    allocate( ebi( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'ebi (Joint DOS)', iok, msg )

    allocate( ebf( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'ebf (Joint DOS)', iok, msg )

    allocate( a( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'a (Joint DOS)', iok, msg )

    call allocate( work, 'work' // srcname, udim1 = nw, nullify = .true. )

!$omp critical
    write( unit = output_unit, fmt = 15 ) isp, omp_get_thread_num()
!$omp end critical

    do jbi = rixs(isp)%nbi(1), rixs(isp)%nbi(2)
      do jbf = rixs(isp)%nbf(1), rixs(isp)%nbf(2)
        do kbz = 1, nkbz
          kqbz = rixs(isp)%ik2kq(kbz)
!$$$          if ( nopused == 1 ) then ! rixs(isp)%ng
!$$$            krbz = kbz
!$$$            kibz = kbz
!$$$            kqibz = kqbz
!$$$          else
            if ( rixs(isp)%ng == 1 ) then
              krbz = kbz
              kibz = ik2ibz(kbz)
              kqibz = ik2ibz(kqbz)
            else
              kibz = ik2ibz(kbz)
              krbz = rixs(isp)%ik2rbz(kbz)
              kqibz = ik2ibz(kqbz)
            endif
!$$$          endif
          ebi(krbz) = e( jbi, kibz )
          ebf(krbz) = e( jbf, kqibz )
          a(krbz) = 1.
        enddo

!        if ( isp == 6 ) then
!          if ( jbi == rixs(isp)%nbi(1) .and. jbf == rixs(isp)%nbf(1) ) then
!            do kbz = 1, rixs(isp)%nkrbz
!              write( * , * ) kbz, ( ebf(kbz) - ebi(kbz) ) * 13.6
!            enddo
!          endif
!        endif

        if ( debug_mode ) then
!$$$          if ( minval( ebf - ebi ) < 0._REAL64 )then
!$$$            write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$                 "ERROR in RIXSLOSS: for some k points ebf < ebi."
!$$$            call exit_on_error( srcname )
!$$$          endif
          if ( minval( a ) < 0.5 )then
            write( unit = output_unit, fmt = "( 4x, a )" ) &
                 "ERROR in RIXSLOSS: some a(:) elements are not referenced."
            call exit_on_error( srcname )
          endif
        endif
        call bzopt( rixs(isp)%itetr, rixs(isp)%idold, rixs(isp)%ntrbz, &
             qbmc, rixs(isp)%fnorm, ebi, ebf, a, ef, work )
      enddo
    enddo
!$omp critical
    loss( : , isp ) = work( : )
!$omp end critical
    deallocate( ebi, stat = iok, errmsg = msg )
    deallocate( ebf, stat = iok, errmsg = msg )
    deallocate(   a, stat = iok, errmsg = msg )
    call deallocate( work, 'work' // srcname )
  enddo
!$omp end do

!$omp end parallel

else

!$omp parallel default(shared) private( ebi, iok, msg, ebf, a, kbz, kqbz, &
!$omp&         krbz , kibz , kqibz , CPUitr_begin , CPUitr_end ,          &
!$omp&         real_itr_begin , real_itr_end , hh , mm , ss, work )

!$omp master

  if( omp_get_num_threads() == 1 )then
    write( unit = output_unit, fmt = 20 )nkbz, ntr
  else
    write( unit = output_unit, fmt = 21 )omp_get_num_threads(), nkbz, ntr
  endif
  if ( omp_get_nested() ) then
    write( unit = output_unit, fmt = " ( 4x, a ) ") &
      trim( spec_type ) // ': nested parallelism has been disabled'
    call omp_set_nested( .false. )
  endif

!$omp end master

!$omp do schedule( dynamic, 1 )
! reduction( + : loss )

  do itr = 1, ntr

    if( itr > NTR_STEP .and.                                &
        mod( itr - 1, max( ntr / NTR_STEP, 1 ) ) == 0 .and. &
        iprint > 0)then
      real_itr2 = real_itr * real( ntr - itr + 1 ) / real( omp_get_num_threads() )
      hh = floor( real_itr2 / 3600. )
      mm = floor( real_itr2 - real(hh) * 3600. ) / 60.
      ss = real_itr2 - real(hh) * 3600. - real(mm) * 60.
!$omp critical
      write( unit = output_unit, fmt = 25 ) itr, omp_get_thread_num(), &
           100*itr/ntr, hh, mm, ss
!$omp end critical
    endif

! estimate time
    call cpu_time( CPUitr_begin )
    real_itr_begin = omp_get_wtime()

    allocate( ebi( rixs( ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'ebi (RIXS)', iok, msg )

    allocate( ebf( rixs( ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'ebf (RIXS)', iok, msg )

    allocate( a( rixs(ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
    if( iok /= 0 ) call print_allocation_error( 'a (RIXS)', iok, msg )

    if ( abs( irixs ) == 10 ) then
      a( : ) = 1.
    else
      call set_a( itr, a )
    endif

    call allocate( work, 'work' // srcname, udim1 = nw, nullify = .true. )

    ebi(:) = huge( 0._REAL64 ) / 2._REAL64
    ebf(:) = huge( 0._REAL64 ) / 2._REAL64
    do kbz = 1, nkbz
      kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!$$$          if ( nopused == 1 ) then ! rixs(isp)%ng
!$$$            krbz = kbz
!$$$            kibz = kbz
!$$$            kqibz = kqbz
!$$$          else
      if ( rixs(ispec(itr))%ng == 1 ) then
        krbz = kbz
        kibz = ik2ibz(kbz)
        kqibz = ik2ibz(kqbz)
      else
        kibz = ik2ibz(kbz)
        krbz = rixs( ispec(itr) )%ik2rbz(kbz)
        kqibz = ik2ibz(kqbz)
      endif
!$$$          endif
      if ( ebi(krbz) > huge( 0._REAL64 ) / 2.2_REAL64 ) then
        ebi(krbz) = e( ibi(itr), kibz )
      else
        if ( debug_mode ) then
          if( abs( ebi(krbz) - e( ibi(itr), kibz ) ) > eps ) then
            write(*,*) 'abs( ebi(krbz) - e( ibi(itr), kibz ) ) = ', &
                 abs( ebi(krbz) - e( ibi(itr), kibz ) )
            call exit_on_error( srcname // ': ebi')
          endif
        endif
      endif
      if ( ebf(krbz) > huge( 0._REAL64 ) / 2.2_REAL64 ) then
        ebf(krbz) = e( ibf(itr), kqibz )
      else
        if ( debug_mode ) then
          if( abs( ebf(krbz) - e( ibf(itr), kqibz ) ) > eps ) then
            write(*,*) 'abs( ebf(krbz) - e( ibf(itr), kqibz ) ) = ', &
                 abs( ebf(krbz) - e( ibf(itr), kqibz ) )
            call exit_on_error( srcname // ': ebf' )
          endif
        endif
      endif
    enddo

    if ( debug_mode ) then
!$$$      if ( minval( ebf - ebi ) < 0._REAL64 )then
!$$$        write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$             "ERROR in RIXSLOSS: for some k points ebf < ebi."
!$$$        call exit_on_error( srcname )
!$$$      endif
      if ( minval( a ) < 0. ) then
!$omp critical
        write( unit = output_unit, fmt = "( 4x, 2(a) )" ) &
             'ERROR in SET_A: at least for one k-point its value' // &
             ' in mme array is not set'
!$omp end critical
        call exit_on_error( srcname )
      endif
    endif

    call bzopt( rixs( ispec(itr) )%itetr, rixs( ispec(itr) )%idold, &
         rixs( ispec(itr) )%ntrbz, qbmc, rixs( ispec(itr) )%fnorm,  &
         ebi, ebf, a, ef, work )
!$omp critical
    loss( : , ispec(itr) ) = loss( : , ispec(itr) ) + work( : )
!$omp end critical
    deallocate( ebi, stat = iok, errmsg = msg )
    deallocate( ebf, stat = iok, errmsg = msg )
    deallocate(   a, stat = iok, errmsg = msg )
    call deallocate( work, 'work' // srcname )
! estimate time
    call cpu_time( CPUitr_end )
    real_itr_end = omp_get_wtime()
!$omp critical
    CPUitr = ( real( nCPUitr, kind = REAL64 ) * CPUitr + CPUitr_end - &
             CPUitr_begin ) / real( nCPUitr + 1, kind = REAL64 )
    real_itr = ( real( nCPUitr, kind = REAL64 ) * real_itr + real_itr_end - &
               real_itr_begin ) / real( nCPUitr + 1, kind = REAL64 )
    nCPUitr = nCPUitr + 1
!$omp end critical
  enddo
!$omp end do
!$omp end parallel
endif

!$$$ if(abs(irixs) == 10)then
!$$$   do isp=1,nrixs
!$$$     loss(:,isp)=loss(:,isp)/real(2*abs(rixs(isp)%nk),kind=REAL64)
!$$$   enddo
!$$$ endif

write( unit = output_unit, fmt = 30 )
! spectra integration
do isp = 1, nrixs
  write( unit = output_unit, fmt = 35 ) trim(spec_type), isp, &
    simpson38( loss( 1, isp ), nw, wmin, wmax )
enddo

call cpu_time( CPUfinal )
realfinal = omp_get_wtime()
call set_CPU_RIXS( CPUfinal - CPUini )
call set_real_RIXS( realfinal - realini )
isp = 2

 20 format( 4x, 'OpenMP do-loop in 1 thread: ', i0, &
            ' k-points in BZ and ', i0, ' transitions in bzopt.' )
 21 format( 4x, 'OpenMP do-loop in ', i0, ' threads: ', i0, &
            ' k-points in BZ and ', i0, ' transitions in bzopt.' )
 25 format( 4x, 'OpenMP do-loop: itr=', i9, ', thread=', i3, &
            ', progress=', i2, '%, time left=', i4, ' hh ', i2, &
            ' mm ', f4.1, ' ss' )
 30 format( 4x, 'RIXS loss spectra have been calculated.', / )
 35 format( 4x, 'Integral of the ', a, ' spectrum for isp=', i0, &
            ' is ', f0.10 )
 55 format( 4x, "ERROR: mme=", f10.6, " < 0 for krbz=", i0)

END SUBROUTINE RIXSLOSS
