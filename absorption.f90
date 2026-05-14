SUBROUTINE ABSORPTION
use m_aux, only: map, lmtindex
use m_bnd
use m_bz
use m_bzmesh
use m_functions
use m_params ! , only: iabsorp, iprint, volomg
use m_results
use m_rixs
use m_time
use omp_lib
implicit none
! local vars
character(BUFFER_SIZE) msg
integer ia,iok,nf,npol,iabem,is_g4q,ng,isp,ib,k,iprint_local
real(REAL64), allocatable :: kpk(:,:),dos(:),bz(:),work(:,:),ebzk(:,:)
real(REAL64), external :: simpson38
real(REAL64) CPUini,CPUfinal,dec,realini,realfinal

if ( read_rixfile ) return
if(iabsorp < 1 .or. iabsorp > 2) return

write(unit=output_unit,fmt="()")
iprint_local=-1
npol=1
iabem=0
ne=int((emax-emin)/de)+1
! nf=nint((ef-emin)/de)+1
nf=nint((ef-emin+1.e-5_REAL64)/de)+1

call allocate( absorp , alias = 'absorp in ABSORPTION' , udim1 = ne &
     , udim2 = nrixs , udim3 = iabsorp , nullify = .true. )

! allocate(absorp(ne,nrixs,iabsorp),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('absorp',iok,msg)
! absorp(:,:,:)=0._REAL64

call allocate( work , alias = 'work in ABSORPTION' , udim1 = ne &
     , udim2 = n_mme , nullify = .true. )

! allocate(work(ne,n_mme),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('work',iok,msg)
! work(:,:)=0._REAL64

! k-point integration over IBZ
call cpu_time(CPUini)
realini=omp_get_wtime()
write(unit=output_unit,fmt=5)
is_g4q=1              ! reference ig4q array
ng=nopused
call kmesh(rbas,ndxyz(1),ndxyz(2),ndxyz(3),g,ng,is_g4q,is_wgt,&
     ibzext,is_idold,is_avtet,iprint_local)

! omg48s=real(ng,kind=REAL64)/real(2*nfu,kind=REAL64)
! vol1t=1._REAL64/vol1t/real(6*nkbz*ng,kind=REAL64)

omg48s=0.5_REAL64
vol1t=volomg

!$omp parallel default(shared) private(ia,dec,bz,kpk,dos,iok,msg)
!$omp master
if(omp_get_num_threads() == 1)then
  write(unit=output_unit,fmt=10)n_mme
else
  write(unit=output_unit,fmt=11)omp_get_num_threads(),n_mme
endif
if ( omp_get_nested() ) then
  write( unit = output_unit , fmt = " ( 4x , a ) ") &
    'ABSORPTION over IBZ: nested parallelism has been disabled'
  call omp_set_nested( .false. )
endif
!$omp end master
!$omp barrier
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

  call allocate( kpk , alias = 'kpk in ABSORPTION' , udim1 = nb &
       , udim2 = nkibz)

  ! allocate(kpk(nb,nkibz),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('kpk for IBZ thread #'//&
  !      trim(int2string(omp_get_thread_num())),iok,msg)

  call set_absorp_ibz(ia,kpk)
  ! shift the "photon energy" scale in order to account for the splitting

!  dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%ec
! trying to use enlmin instead of ec

  dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%enlmin

  dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%ec

  ! integrate and sum contributions from spin split core levels.
  ! kpk is input. bz is output. dos is not used now.
  ! limzn=limallzn=nbmax=nb and nbmin=1.

  call allocate( bz , alias = 'bz in ABSORPTION' , udim1 = ne &
       , nullify = .true. )

  ! allocate(bz(ne),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('bz for IBZ thread #'//&
  !      trim(int2string(omp_get_thread_num())),iok,msg)
  ! bz(:)=0._REAL64

  call allocate( dos , alias = 'dos in ABSORPTION' , udim1 = ne )

  ! allocate(dos(ne),stat=iok,errmsg=msg)
  ! if(iok /= 0)call print_allocation_error('dos for IBZ thread #'//&
  !      trim(int2string(omp_get_thread_num())),iok,msg)

  call tetdos(ne,npol-1,npol,nb,nb,emin,de,bz,dos,e,kpk,1,nb)
  call x_intp(iabem,dec,nf,de,ne,bz,work(1,ia))

  call deallocate( kpk , alias = 'kpk in ABSORPTION' )
  call deallocate( bz , alias = 'bz in ABSORPTION' )
  call deallocate( dos , alias = 'dos in ABSORPTION' )

  ! deallocate(kpk,stat=iok)
  ! deallocate(bz,stat=iok)
  ! deallocate(dos,stat=iok)
enddo                             ! ia
!$omp enddo
!$omp end parallel

do ia=1,n_mme
  do isp=1,nrixs
    if(map(ia,isp) > 0)then
      absorp(:,isp,1)=absorp(:,isp,1)+work(:,ia)
    endif
  enddo
enddo

  5 format(4x,'ABSORPTION: Integration over the IBZ.')
 10 format(4x,'OpenMP do-loop in 1 thread for ',i0,' core levels.')
 11 format(4x,'OpenMP do-loop in ',i0,' threads for ',i0,' core levels.')
 12 format(4x,'Calculations in 1 thread for ',i0,' core levels.')
 15 format(4x,'Cycling  do-loop  for ia=',i0,', thread=',i0,'.')
 20 format(4x,'Doing integration for ia=',i0,', thread=',i0,'.')
 25 format(/,4x,'ABSORPTION: Integration over the RBZ.')
 35 format(4x,'Integral of the absoprtion spectrum for isp=',i0,' is ',f0.10)

! time needed to integrate over IBZ
call cpu_time(CPUfinal)
realfinal=omp_get_wtime()
call set_CPU_ABSORP_IBZ(CPUfinal-CPUini)
call set_real_ABSORP_IBZ(realfinal-realini)
! spectra integration
do isp=1,nrixs
  write(unit=output_unit,fmt=35)isp,simpson38(absorp(1,isp,1),ne,emin,emax)
enddo
if(iabsorp /= 2)return

! k-point integration over RBZ
call cpu_time(CPUini)
realini=omp_get_wtime()
work(:,:)=0._REAL64
write(unit=output_unit,fmt=25)
write(unit=output_unit,fmt=12)n_mme

call allocate( bz , alias = 'bz in ABSORPTION RBZ' , udim1 = ne )

! allocate(bz(ne),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('bz for RBZ',iok,msg)

call allocate( dos , alias = 'dos in ABSORPTION RBZ' , udim1 = ne )

! allocate(dos(ne),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('dos for RBZ',iok,msg)

do ia = 1 , n_mme
  do isp = 1 , nrixs
    if ( map( ia , isp ) == 0 ) cycle
    write( unit = output_unit , fmt ="( 4x , 2 ( a , i3 ) ) " ) &
         'Calculation for ia=' , ia , ', ispec=' , isp

! local arrays
    call allocate( ebzk , alias = 'ebzk in ABSORPTION' , udim1 = nb &
         , udim2 = rixs(isp)%nkrbz )
    ! allocate(ebzk(nb,rixs(isp)%nkrbz),stat=iok,errmsg=msg)
    ! if(iok /= 0)call print_allocation_error('ebzk',iok,msg)

    call allocate( kpk , alias = 'kpk in ABSORPTION' , udim1 = nb &
         , udim2 = rixs(isp)%nkrbz )
    ! allocate(kpk(nb,rixs(isp)%nkrbz),stat=iok,errmsg=msg)
    ! if(iok /= 0)call print_allocation_error('kpk for RBZ',iok,msg)

! I want to be on the safe side here
    is_g4q=0              ! not to reference ig4q array
    if ( nopused == 1 ) then
      call kmesh( rbas , ndxyz(1) , ndxyz(2) , ndxyz(3) , g , ng , is_g4q , &
           is_wgt , ibzext , is_idold , is_avtet , iprint_local )
    else
      call kmesh( rbas , ndxyz(1) , ndxyz(2) , ndxyz(3) , rixs(isp)%g &
           , rixs(isp)%ng , is_g4q , is_wgt , ibzext , is_idold , is_avtet &
           , iprint_local )
    endif
    omg48s=4._REAL64
    vol1t=volomg
    if ( nopused == 1 ) then
        ebzk( : , : ) = e( : , : )
    else
      if ( rixs(isp)%ng == 1 ) then
        forall( ib = 1 : nb , k = 1 : nkbz ) &
          ebzk( ib , k ) = e( ib , ik2ibz(k) )
      else
        forall( ib = 1 : nb , k = 1 : nkbz ) &
          ebzk( ib , rixs(isp)%ik2rbz(k) ) = e( ib , ik2ibz(k) )
      endif
    endif
    call set_absorp_rbz( isp , ia , kpk )
  ! shift the "photon energy" scale in order to account for the splitting
    dec=lmtdata(lmtindex(ia))%esp(ia)-lmtdata(lmtindex(ia))%ec
  ! integrate and sum contributions from spin split core levels.
  ! kpk is input. bz is output. dos is not used now.
  ! limzn=limallzn=nbmax=nb and nbmin=1.
    bz(:)=0._REAL64
    call tetdos(ne,npol-1,npol,nb,nb,emin,de,bz,dos,ebzk,kpk,1,nb)
    call x_intp(iabem,dec,nf,de,ne,bz,work(1,ia))
    call deallocate(ebzk)
    call deallocate(kpk)
    ! deallocate(ebzk,stat=iok)
    ! deallocate(kpk,stat=iok)
  enddo
enddo                             ! ia
do ia=1,n_mme
  do isp=1,nrixs
    if(map(ia,isp) > 0)then
      absorp(:,isp,2)=absorp(:,isp,2)+work(:,ia)
    endif
  enddo
enddo
call deallocate( bz , alias = 'bz in ABSORPTION' )
call deallocate( dos , alias = 'dos in ABSORPTION' )
call deallocate( work , alias = 'work in ABSORPTION' )

! deallocate(bz,stat=iok)
! deallocate(dos,stat=iok)
! deallocate(work,stat=iok)

! time needed to integrate over RBZ
call cpu_time(CPUfinal)
realfinal=omp_get_wtime()
call set_CPU_ABSORP_RBZ(CPUfinal-CPUini)
call set_real_ABSORP_RBZ(realfinal-realini)
! spectra integration
do isp=1,nrixs
  write(unit=output_unit,fmt=35)isp,simpson38(absorp(1,isp,2),ne,emin,emax)/8.
enddo
return

!$$$  AHTUNG: THIS CODE IS NOT USED NOW. IT IS SAVED FOR FUTURE

! local arrays
allocate(ebzk(nb,nkbz),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('ebzk',iok,msg)
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

