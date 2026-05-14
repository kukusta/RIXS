SUBROUTINE ABSORPTION
use m_aux, only: map, lmtindex
use m_bnd
use m_bz
use m_bzmesh
use m_functions
use m_params
use m_results
use m_rixs
use m_time
use omp_lib
implicit none
! local vars
character(BUFFER_SIZE) msg
character(9) text
integer ia, iok, nf, npol, iabem, is_g4q, ng, isp, ib, k, iprint_local
real(REAL64), allocatable :: kpk(:,:), dos(:), bz(:), work(:,:), ebzk(:,:)
real(REAL64), external :: simpson38
real(REAL64) CPUini, CPUfinal, dec, realini, realfinal

if (read_rixfile) return
if (iabsorp < 1 .or. iabsorp > 2) return

write(unit = output_unit, fmt = "()")
iprint_local = -1
npol = 1
iabem = 0
ne = int((emax - emin) / de) + 1
! nf=nint((ef-emin)/de)+1
nf = nint((ef - emin + eps5) / de) + 1
call allocate(absorp, 'absorp in ABSORPTION', udim1 = ne, udim2 = nrixs, &
     udim3 = iabsorp, nullify = .true.)
call allocate(work, 'work in ABSORPTION', udim1 = ne, udim2 = n_mme, &
     nullify = .true.)

! k-point integration over IBZ
call cpu_time(CPUini)
realini = omp_get_wtime()
write(unit = output_unit, fmt = 5)
is_g4q = 1              ! reference ig4q array
ng = nopused
call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
     ibzext, is_idold, is_avtet, iprint_local)

! omg48s=real(ng,kind=REAL64)/real(2*nfu,kind=REAL64)
! vol1t=1._REAL64/vol1t/real(6*nkbz*ng,kind=REAL64)

omg48s = 0.5_REAL64
vol1t = volomg

! $omp parallel default(shared) private(ia, dec, bz, kpk, dos, iok, msg)
! $omp master
text(:) = ' threads'
if (omp_get_num_threads() == 1) text(:) = ' thread'
write(unit = output_unit, fmt = 11) omp_get_num_threads(), &
     text( : len_trim(text)), n_mme
if (omp_get_max_active_levels() > 1) then
  msg(:) = 'ABSORPTION over IBZ: nested parallelism has been disabled'
  write(unit = output_unit, fmt = 10) trim(msg)       
  call omp_set_nested(.false.)
endif
! $omp end master
! $omp barrier
! $omp do schedule(dynamic,1)
do ia = 1, n_mme
! $omp critical
  if (iprint > 0) then
    if (sum(map(ia,:)) == 0) then
      if (iprint > 10) then
        write(unit = output_unit, fmt = 15) ia, omp_get_thread_num()
      endif
    else
      write(unit = output_unit, fmt = 20) ia, omp_get_thread_num()
    endif
  endif
! $omp end critical
  if (sum(map(ia,:)) == 0) cycle
  call allocate(kpk, 'kpk in ABSORPTION', udim1 = nb, udim2 = nkibz)
  call set_absorp_ibz(ia, kpk)
! shift the "photon energy" scale in order to account for the splitting
! dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%ec
! trying to use enlmin instead of ec
  dec = lmtdata(lmtindex(ia))%esp(ia) - lmtdata(lmtindex(ia))%enlmin
  dec = lmtdata(lmtindex(ia))%esp(ia) - lmtdata(lmtindex(ia))%ec
! integrate and sum contributions from spin split core levels.
! kpk is input. bz is output. dos is not used now.
! limzn = limallzn = nbmax = nb and nbmin = 1.
  call allocate( bz,  'bz in ABSORPTION', udim1 = ne)
  call allocate(dos, 'dos in ABSORPTION', udim1 = ne)
  call tetdos(ne, npol - 1, npol, nb, nb, emin, de, bz, dos, e, kpk, 1, nb)
  call x_intp(iabem, dec, nf, de, ne, bz, work(1, ia))
  call deallocate(kpk, 'kpk in ABSORPTION')
  call deallocate( bz, ' bz in ABSORPTION')
  call deallocate(dos, 'dos in ABSORPTION')
enddo
! $omp enddo
! $omp end parallel

do ia = 1, n_mme
  do isp = 1, nrixs
    if (map(ia, isp) > 0) then
      absorp(:, isp, 1) = absorp(:, isp, 1) + work(:, ia)
    endif
  enddo
enddo

  5 format(4x, 'ABSORPTION: Integration over the IBZ.')
 10 format(4x, a)
 11 format(4x, 'OpenMP do-loop in ', i0, a, ' for ', i0, ' core levels.')
 12 format(4x, 'Calculations in 1 thread for ', i0, ' core levels.')
 15 format(4x, 'Cycling  do-loop  for ia=', i0, ', thread=', i0, '.')
 20 format(4x, 'Doing integration for ia=', i0, ', thread=', i0, '.')
 25 format(/, 4x, 'ABSORPTION: Integration over the BZ.')
 35 format(4x, 'Integral of the absoprtion spectrum for isp=', i0, ' is ', &
           f0.10)

! time needed to integrate over IBZ
call cpu_time(CPUfinal)
realfinal = omp_get_wtime()
call set_CPU_ABSORP_IBZ(CPUfinal - CPUini)
call set_real_ABSORP_IBZ(realfinal - realini)
! spectra integration
do isp = 1, nrixs
  write(unit = output_unit, fmt = 35) isp, &
       simpson38(absorp(1, isp, 1), ne, emin, emax)
enddo
if (iabsorp /= 2) return

! k-point integration over BZ
call cpu_time(CPUini)
realini = omp_get_wtime()
work(:,:) = 0._REAL64
write(unit = output_unit, fmt = 25)
write(unit = output_unit, fmt = 12) n_mme
call allocate( bz, alias =  'bz in ABSORPTION BZ', udim1 = ne)
call allocate(dos, alias = 'dos in ABSORPTION BZ', udim1 = ne)
do ia = 1, n_mme
  do isp = 1, nrixs
    if (map(ia, isp) == 0) cycle
    write(unit = output_unit, fmt = "(4x, 2(a, i3))") &
         'Calculation for ia=', ia, ', ispec=', isp
! local arrays
!vvv    call allocate(ebzk, alias = 'ebzk in ABSORPTION', udim1 = nb, &
!vvv         udim2 = rixs(isp)%nkrbz)
!vvv    call allocate(kpk, alias = 'kpk in ABSORPTION', udim1 = nb, &
!vvv         udim2 = rixs(isp)%nkrbz)
    call allocate(ebzk, 'ebzk in ABSORPTION', udim1 = nb, udim2 = nkbz)
    call allocate( kpk, ' kpk in ABSORPTION', udim1 = nb, udim2 = nkbz)
! I want to be on the safe side here
    is_g4q = 0              ! do not reference ig4q array
!vvv    if (nopused == 1) then
      call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, &
           is_wgt, ibzext, is_idold, is_avtet, iprint_local)
!vvv    else
!vvv      call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), rixs(isp)%g, &
!vvv           rixs(isp)%ng, is_g4q, is_wgt, ibzext, is_idold, is_avtet, &
!vvv           iprint_local)
!vvv    endif
    omg48s = 4._REAL64
    vol1t = volomg
!vvv    if (nopused == 1) then
!vvv        ebzk(:,:) = e(:,:)
!vvv    else
!vvv      if (rixs(isp)%ng == 1) then
    forall(ib = 1:nb, k = 1:nkbz) ebzk(ib, k) = e(ib, ik2ibz(k))
!vvv      else
!vvv        forall(ib = 1:nb , k = 1:nkbz) &
!vvv              ebzk(ib, rixs(isp)%ik2rbz(k)) = e(ib, ik2ibz(k))
!vvv      endif
!vvv    endif
    call set_absorp_rbz(isp, ia, kpk)
! shift the "photon energy" scale in order to account for the splitting
    dec = lmtdata(lmtindex(ia))%esp(ia) - lmtdata(lmtindex(ia))%ec
! integrate and sum contributions from spin split core levels.
! kpk is input. bz is output. dos is not used now.
! limzn = limallzn = nbmax = nb and nbmin = 1.
    bz(:) = 0._REAL64
    call tetdos(ne, npol - 1, npol, nb, nb, emin, de, bz, dos, ebzk, kpk, 1, nb)
    call x_intp(iabem, dec, nf, de, ne, bz, work(1, ia))
    call deallocate(ebzk)
    call deallocate(kpk)
  enddo
enddo
do ia = 1, n_mme
  do isp = 1, nrixs
    if (map(ia,isp) > 0) then
      absorp(:, isp, 2) = absorp(:, isp, 2) + work(:, ia)
    endif
  enddo
enddo
call deallocate(bz   , alias = 'bz in ABSORPTION'  )
call deallocate(dos  , alias = 'dos in ABSORPTION' )
call deallocate(work , alias = 'work in ABSORPTION')
! time needed to integrate over RBZ
call cpu_time(CPUfinal)
realfinal = omp_get_wtime()
call set_CPU_ABSORP_RBZ(CPUfinal - CPUini)
call set_real_ABSORP_RBZ(realfinal - realini)
! spectra integration
do isp = 1, nrixs
  write(unit = output_unit, fmt = 35)isp, &
       simpson38(absorp(1, isp, 2), ne, emin, emax) / 8.
enddo
return

!$$$  AHTUNG: THIS CODE IS NOT USED NOW. IT IS SAVED FOR FUTURE

! local arrays
call allocate(ebzk, 'ebzk (ABSORPTION)', udim1 = nb, udim2 = nkbz)
work(:,:)=0._REAL64

! k-point integration over BZ
call cpu_time(CPUini)
realini=omp_get_wtime()
write(unit=output_unit,fmt=25)
is_g4q=0              ! not to reference ig4q array
ng=1
call kmesh(rbas,ndxyz(1),ndxyz(2),ndxyz(3),g,ng,is_g4q,is_wgt&
     ,ibzext,is_idold,is_avtet,iprint_local)

! omg48s=real(ng,kind=REAL64)/real(2*nfu,kind=REAL64)
! vol1t=1._REAL64/vol1t/real(6*nkbz*ng,kind=REAL64)

omg48s=4._REAL64
vol1t=volomg

forall(ib=1:nb,k=1:nkbz) ebzk(ib,k)=e(ib,ik2ibz(k))

!$omp parallel default(shared) private(ia,dec,bz,kpk,dos,iok,msg)
!$omp master
if(omp_get_num_threads() == 1)then
  write(unit=output_unit,fmt=10)n_mme
else
  write(unit=output_unit,fmt=11)omp_get_num_threads(),n_mme
endif
!$omp end master
!$omp do schedule(dynamic,1)
do ia=1,n_mme
!$omp critical
  if(iprint > 0)then
    if(sum(map(ia,:)) == 0)then
      if(iprint > 10)write(unit=output_unit,fmt=15)ia,omp_get_thread_num()
    else
      write(unit=output_unit,fmt=20)ia,omp_get_thread_num()
    endif
  endif
!$omp end critical
  if(sum(map(ia,:)) == 0)cycle
  allocate(kpk(nb,nkbz),stat=iok,errmsg=msg)
  if(iok /= 0)call print_allocation_error('kpk for BZ thread #'//&
       trim(int2string(omp_get_thread_num())),iok,msg)
  call set_absorp_bz(ia,kpk)
  ! shift the "photon energy" scale in order to account for the splitting
  dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%ec
  ! integrate and sum contributions from spin split core levels.
  ! kpk is input. bz is output. dos is not used now.
  ! limzn=limallzn=nbmax=nb and nbmin=1.
  allocate(bz(ne),stat=iok,errmsg=msg)
  if(iok /= 0)call print_allocation_error('bz for BZ thread #'//&
       trim(int2string(omp_get_thread_num())),iok,msg)
  allocate(dos(ne),stat=iok,errmsg=msg)
  if(iok /= 0)call print_allocation_error('dos for BZ thread #'//&
       trim(int2string(omp_get_thread_num())),iok,msg)
  bz(:)=0._REAL64
  call tetdos(ne,npol-1,npol,nb,nb,emin,de,bz,dos,ebzk,kpk,1,nb)
  call x_intp(iabem,dec,nf,de,ne,bz,work(1,ia))
  deallocate(kpk,stat=iok)
  deallocate(bz,stat=iok)
  deallocate(dos,stat=iok)
enddo                             ! ia
!$omp end do
!$omp end parallel

do ia=1,n_mme
  do isp=1,nrixs
    if(map(ia,isp) > 0)then
      absorp(:,isp,2)=absorp(:,isp,2)+work(:,ia)
    endif
  enddo
enddo

deallocate(work,stat=iok)
deallocate(ebzk,stat=iok)

! time needed to integrate over BZ
call cpu_time(CPUfinal)
realfinal=omp_get_wtime()
call set_CPU_ABSORP_BZ(CPUfinal-CPUini)
call set_real_ABSORP_BZ(realfinal-realini)

END SUBROUTINE ABSORPTION

