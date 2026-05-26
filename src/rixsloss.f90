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
character(*), parameter :: srcname=' in RIXSLOSS'
character(10) text, spec_type
integer iok, itr, isp, jbi, jbf, kibz, kqibz, kbz, kqbz, nCPUitr, hh, mm, krbz
real, allocatable :: a(:)
real(REAL64), allocatable :: ebi(:) , ebf(:), work(:)
real(REAL64), external :: simpson38
real(REAL64) CPUini, CPUfinal, realini, realfinal, CPUitr_begin, CPUitr_end, &
             CPUitr, real_itr, real_itr2, real_itr_end, real_itr_begin, ss

if (read_rixfile) return

call cpu_time(CPUini)
realini = omp_get_wtime()

  5 format(/, 4x, a, ' spectra calculations.')
 10 format(4x, 'OpenMP do-loop in 1 thread for ', i0, &
           ' k-points in BZ (Joint DOS).')
 11 format(4x, 'OpenMP do-loop in ', i0, ' threads for ', i0, &
           ' k-points in BZ (Joint DOS).')
 15 format(4x, 'Joint DOS calculations for ispec=', i0, ' thread=', i0)

! allocate(loss(nw,nrixs),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('loss',iok,msg)
! loss(:,:) = 0._REAL64

call allocate(loss, alias = 'loss' // srcname, udim1 = nw, udim2 = nrixs, &
     nullify = .true.)

if (absorption_only) then
  call set_CPU_RIXS(0._REAL64)
  call set_real_RIXS(0._REAL64)
  return
endif

if (irixs /= 0 .and. irixs /= 10) then
  spec_type = 'RIXS loss'
else
  spec_type = 'Joint DOS'
endif
write(unit = output_unit, fmt = 5) trim(spec_type)

nCPUitr  = 0
CPUitr   = 0._REAL64
real_itr = 0._REAL64

if (irixs == 0) then
!$omp parallel default(shared) private(a, jbi, jbf, kbz, krbz, ebi, ebf, &
!$omp&         itr, kibz, kqibz, kqbz, isp, CPUitr_begin, CPUitr_end,    &
!$omp&         real_itr_begin, real_itr_end, hh, mm, ss, work)
!$omp master
  if (omp_get_num_threads() == 1) then
    write(unit = output_unit, fmt = 10) nkbz
  else
    write(unit = output_unit, fmt = 11) omp_get_num_threads(), nkbz
  endif
  if (omp_get_max_active_levels() > 1) then
    write(unit = output_unit, fmt = "(4x, a)") &
      trim(spec_type) // ': nested parallelism has been disabled'
    call omp_set_nested(.false.)
  endif
!$omp end master
!$omp do  schedule( dynamic, 1 )
! reduction( + : loss )
  do isp = 1, nrixs
!vvv    allocate( ebi( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'ebi (Joint DOS)', iok, msg )
!vvv    allocate( ebf( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'ebf (Joint DOS)', iok, msg )
!vvv    allocate( a( rixs(isp)%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'a (Joint DOS)', iok, msg )
    call allocate( ebi,  'ebi' // srcname, udim1 = nkbz)
    call allocate( ebf,  'ebi' // srcname, udim1 = nkbz)
    call allocate(   a,    'a' // srcname, udim1 = nkbz)
    call allocate(work, 'work' // srcname, udim1 = nw, nullify = .true.)
!$omp critical
    write(unit = output_unit, fmt = 15) isp, omp_get_thread_num()
!$omp end critical
    do jbi = rixs(isp)%nbi(1), rixs(isp)%nbi(2)
      do jbf = rixs(isp)%nbf(1), rixs(isp)%nbf(2)
        do kbz = 1, nkbz
!$$$          if ( nopused == 1 ) then ! rixs(isp)%ng
!$$$            krbz = kbz
!$$$            kibz = kbz
!$$$            kqibz = kqbz
!$$$          else
!vvv            if ( rixs(isp)%ng == 1 ) then
!vvv              krbz = kbz
!vvv              kibz = ik2ibz(kbz)
!vvv              kqibz = ik2ibz(kqbz)
!vvv            else
!vvv              kibz = ik2ibz(kbz)
!vvv              krbz = rixs(isp)%ik2rbz(kbz)
!vvv              kqibz = ik2ibz(kqbz)
!vvv            endif
!$$$          endif
!vvv          ebi(krbz) = e( jbi, kibz )
!vvv          ebf(krbz) = e( jbf, kqibz )
!vvv          a(krbz) = 1.
          kibz  = ik2ibz(kbz)
          kqbz  = rixs(isp)%ik2kq(kbz)
          kqibz = ik2ibz(kqbz)
          ebi(kbz) = e(jbi, kibz)
          ebf(kbz) = e(jbf, kqibz)
          a(kbz) = 1.
        enddo
!vvv        call bzopt( rixs(isp)%itetr, rixs(isp)%idold, rixs(isp)%ntrbz, &
!vvv             qbmc, rixs(isp)%fnorm, ebi, ebf, a, ef, work )
        call bzopt(itetr, idold, ntibz, qbmc, rixs(isp)%fnorm, ebi, ebf, a, &
             ef, work)
      enddo
    enddo
!$omp critical
    write(*,*)'nrixs=', nrixs, nkibz, nkbz, ntibz, ntbz
    loss(:, isp) = work(:)
!$omp end critical
!vvv    deallocate( ebi, stat = iok, errmsg = msg )
!vvv    deallocate( ebf, stat = iok, errmsg = msg )
!vvv    deallocate(   a, stat = iok, errmsg = msg )
    call deallocate( ebi,  'ebi' // srcname)
    call deallocate( ebf,  'ebf' // srcname)
    call deallocate(   a,    'a' // srcname)
    call deallocate(work, 'work' // srcname)
  enddo
!$omp end do
!$omp end parallel
else
!$omp parallel default(shared) private( ebi, iok, ebf, a, kbz, kqbz, &
!$omp&         krbz , kibz , kqibz , CPUitr_begin , CPUitr_end ,          &
!$omp&         real_itr_begin , real_itr_end , hh , mm , ss, work )
!$omp master
  text(:) = ' threads: '
  if (omp_get_num_threads() == 1) text(:) = ' thread: '
  write(unit = output_unit, fmt = 21) omp_get_num_threads(), text, nkbz, ntr
  if (omp_get_max_active_levels() > 1) then
    write(unit = output_unit, fmt = "(4x, a)") &
      trim(spec_type) // ': nested parallelism has been disabled'
    call omp_set_nested(.false.)
  endif
!$omp end master
!$omp do schedule(static, 1) ! dynamic
! reduction( + : loss )
  do itr = 1, ntr
    if (itr > NTR_STEP .and. mod(itr - 1, max(ntr / NTR_STEP, 1)) == 0 .and. &
        iprint > 0) then
      real_itr2 = real_itr * real(ntr - itr + 1) / real(omp_get_num_threads())
      hh = floor(real_itr2 / 3600.)
      mm = floor(real_itr2 - real(hh) * 3600.) / 60.
      ss = real_itr2 - real(hh) * 3600. - real(mm) * 60.
!$omp critical
      write(unit = output_unit, fmt = 25) itr, omp_get_thread_num(), &
           100 * itr / ntr, hh, mm, ss
!$omp end critical
    endif
    call cpu_time(CPUitr_begin)
    real_itr_begin = omp_get_wtime()

!vvv    allocate( ebi( rixs( ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'ebi (RIXS)', iok, msg )
!vvv    allocate( ebf( rixs( ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'ebf (RIXS)', iok, msg )
!vvv    allocate( a( rixs(ispec(itr) )%nkrbz ), stat = iok, errmsg = msg )
!vvv    if( iok /= 0 ) call print_allocation_error( 'a (RIXS)', iok, msg )

    call allocate( ebi,  'ebi (RIXS)' // srcname, udim1 = nkbz)
    call allocate( ebf,  'ebf (RIXS)' // srcname, udim1 = nkbz)
    call allocate(   a,     'a(RIXS)' // srcname, udim1 = nkbz)
    call allocate(work, 'work (RIXS)' // srcname, udim1 = nw, nullify = .true.)
    if (irixs == 10) then
      a(:) = 1.
    else
      call set_a(itr, a)
    endif
!vvv    ebi(:) = huge(0._REAL64) / 2._REAL64
!vvv    ebf(:) = huge(0._REAL64) / 2._REAL64
    do kbz = 1, nkbz
!$$$          if ( nopused == 1 ) then ! rixs(isp)%ng
!$$$            krbz = kbz
!$$$            kibz = kbz
!$$$            kqibz = kqbz
!$$$          else
!vvv      if ( rixs(ispec(itr))%ng == 1 ) then
!vvv        krbz = kbz
!vvv        kibz = ik2ibz(kbz)
!vvv        kqibz = ik2ibz(kqbz)
!vvv      else
!vvv        kibz = ik2ibz(kbz)
!vvv        krbz = rixs( ispec(itr) )%ik2rbz(kbz)
!vvv        kqibz = ik2ibz(kqbz)
!vvv      endif

        kibz  = ik2ibz(kbz)
        kqbz  = rixs(ispec(itr))%ik2kq(kbz)
        kqibz = ik2ibz(kqbz)
        ebi(kbz) = e(ibi(itr), kibz)
        ebf(kbz) = e(ibf(itr), kqibz)
        ! if (1.01*1.57 / 13.605698 >= ebf(kbz) - ebi(kbz) .and. &
        !      0.99*1.57 / 13.605698 <= ebf(kbz) - ebi(kbz)) then
        !   if (a(kbz) > 0.1) write(*,*) ibi(itr), ibf(itr), i1(kbz), i2(kbz), i3(kbz), a(kbz)
        ! endif

!$$$          endif
!vvv      if ( ebi(krbz) > huge( 0._REAL64 ) / 2.2_REAL64 ) then
!vvv        ebi(krbz) = e( ibi(itr), kibz )
!vvv      endif
!vvv      if ( ebf(krbz) > huge( 0._REAL64 ) / 2.2_REAL64 ) then
!vvv        ebf(krbz) = e( ibf(itr), kqibz )
!vvv      endif
    enddo
!vvv    call bzopt( rixs( ispec(itr) )%itetr, rixs( ispec(itr) )%idold, &
!vvv         rixs( ispec(itr) )%ntrbz, qbmc, rixs( ispec(itr) )%fnorm,  &
!vvv         ebi, ebf, a, ef, work )
    call bzopt(itetr, idold, ntibz, qbmc, rixs(ispec(itr))%fnorm,  &
         ebi, ebf, a, ef, work)
!$omp critical
    loss(:, ispec(itr)) = loss(:, ispec(itr)) + work(:)
!$omp end critical
!vvv    deallocate( ebi, stat = iok, errmsg = msg )
!vvv    deallocate( ebf, stat = iok, errmsg = msg )
!vvv    deallocate(   a, stat = iok, errmsg = msg )
    call deallocate( ebi,  'ebi (RIXS)' // srcname)
    call deallocate( ebf,  'ebf (RIXS)' // srcname)
    call deallocate(   a,    'a (RIXS)' // srcname)
    call deallocate(work, 'work (RIXS)' // srcname)
! estimate time
    call cpu_time(CPUitr_end)
    real_itr_end = omp_get_wtime()
!$omp critical
    CPUitr = (real(nCPUitr, kind = REAL64) * CPUitr + CPUitr_end - &
              CPUitr_begin) / real(nCPUitr + 1, kind = REAL64)
    real_itr = (real(nCPUitr, kind = REAL64) * real_itr + real_itr_end - &
                real_itr_begin) / real(nCPUitr + 1, kind = REAL64)
    nCPUitr = nCPUitr + 1
!$omp end critical
  enddo
!$omp end do
!$omp end parallel
endif

write(unit = output_unit, fmt = 30)
! spectra integration
do isp = 1, nrixs
  write(unit = output_unit, fmt = 35) trim(spec_type), isp, &
       simpson38(loss(1, isp), nw, wmin, wmax)
enddo

call cpu_time(CPUfinal)
realfinal = omp_get_wtime()
call set_CPU_RIXS(CPUfinal - CPUini)
call set_real_RIXS(realfinal - realini)

 21 format(4x, 'OpenMP do-loop in ', i0, a, i0, ' k-points in BZ and ', i0, &
           ' bzopt transitions')
 25 format(4x, 'OpenMP do-loop: itr=', i9, ', thread=', i3, ', progress=', i2, &
           '%, time left=', i4, ' hh ', i2, ' mm ', f4.1, ' ss')
 30 format(4x, 'RIXS loss spectra have been calculated.', /)
 35 format(4x, 'Integral of the ', a, ' spectrum for isp=', i0, ' is ', f12.4)

END SUBROUTINE RIXSLOSS
