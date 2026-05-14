SUBROUTINE read_RIXS_mme
use m_params
use m_time
use omp_lib
implicit none
real(REAL64) timeini, timefinal

if (read_rixfile) return

if ((irixs == 0 .or. irixs == 10) .and. iabsorp == 0) then
  write(unit = output_unit, fmt = 10) 'Joint DOS without absorption:' // &
       ' RIXS mme are not read'
  return
endif
if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Ready to read matrix elements'
endif
timeini = omp_get_wtime()
call open_mmefile
call perform_check_kstar  ! if check_kstar is set RIXS would stop here
call allocate_lmtdata_mme
call read_mme_for_IBZ
call set_mme_for_BZ
timefinal = omp_get_wtime()
call set_readmme(timefinal - timeini)
if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Matrix elements have been read'
endif
 10 format(4x, a)

END SUBROUTINE read_RIXS_mme

SUBROUTINE open_mmefile
use, intrinsic :: iso_fortran_env
use m_bnd
use m_constants, only: z0
use m_files, only: rim, mmefile
use m_functions, only: int2string, exit_on_error
use m_params
use m_rixs, only: max_nhsort, n_mme
implicit none
character(*), parameter :: srcname = ' in READ_RIXS_MME'
integer iok
integer(INT64) reclength, size
logical exist

if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Ready to open mmefile'
endif
inquire (file = mmefile, exist = exist)
if (.not. exist) then
  write(unit = output_unit, fmt = 10) 'File ' // trim(mmefile) // &
       ' does not exist'
  call exit_on_error(srcname)
endif
reclength = int(storage_size(z0), kind = INT64) * 3_INT64 * &
     int(nb, kind = INT64) * int(n_mme, kind = INT64) *     &
     int(max_nhsort, kind = INT64)
open(newunit = rim, file = mmefile, action = 'read', status = 'old', &
     form = 'unformatted', access = 'direct', recl = reclength, iostat = iok)
if (iok /= 0) then
  write(unit = output_unit, fmt = 10) 'mmefile cannot be opened, exiting'
  call exit_on_error(srcname)
endif
inquire (file = mmefile, size = size)
! check the size of mmefile. In main LMTO program rixsmmefile must be
! overwritten so its size can be easy checked
if (size /= reclength * int(npnt, kind = INT64)) then
  write(unit = output_unit, fmt = 13) reclength * int(npnt, kind = INT64), size
  call exit_on_error(srcname)
endif
if (iprint > 0) then
  write(unit = output_unit, fmt = 20) max(1, int(real(size, kind = REAL64) / &
       real(1024**2, kind = REAL64))), trim(int2string(size, ' '))
endif

10 format(4x, a)
13 format(4x, "Wrong size of rixsmmefile. Expected size is ", i0, &
          " B but real size is ", i0, " B")
20 format(4x, 'Size of mmefile is ', i0, ' MB (', a,' B)')

END SUBROUTINE open_mmefile

SUBROUTINE perform_check_kstar
use, intrinsic :: iso_fortran_env
use m_bnd, only: iopnum, nb, npnt
use m_files, only: rim
use m_functions, only: allocate, exit_on_error, deallocate
use m_params, only: check_kstar
use m_rixs, only: max_nhsort, n_mme
implicit none
character(*), parameter :: srcname = ' in DO_CHECK_KSTAR'
complex(REAL64), allocatable :: mme(:,:,:,:)
complex(REAL64) mme_kibz(3), mme_kbz(3), tensor(3, 3)
integer ia, iat, ib, iok, isp, k
real(REAL64) op(3, 3)

  9 format(4x, 'File ', a, ' does not exist.')
 10 format(4x, 'Operation #(', i0, '), its matrix')
 15 format(3f15.10)
 11 format(4x, 'RIXS mme for k=', i0)
 12 format(4x, 'iat=', i2, ', ia=', i2, ', ib=', i2, ', mme(1:3)=', 6f10.6)
 35 format(4x, 'In READRIXSMME: cannot read from mmefile for k=', i0)
  ! data for CHECK_KSTAR
  ! character(8) named
  ! integer idummy,ntxtd,ibz
  ! real ef4
  ! inquire (file = bnsfile, exist = exist)
  ! if (.not. exist) then
  !   write(unit = output_unit, fmt = 9) trim(bnsfile)
  !   call exit_on_error(srcname)
  ! endif
  ! open(newunit = bns, file = bnsfile, action = 'read', status = 'unknown', &
  !      form = 'unformatted' ,iostat = iok)
  ! if (iok /= 0) stop '    Bnsfile cannot be opened, exiting.'
  ! read(bns, iostat = iok) named(:)
  ! if(iok /= 0) stop '    Cannot read named(:) from bnsfile.'
  ! if (named(8:8) == '+') then
  !   ntxtd = 4
  ! else
  !   ntxtd = 1
  ! endif
  ! nopused = -1
  ! high_sym_pnt_txt_dir(:) = ""
  ! read( unit = bns, iostat = iok) idummy, ef4, npnt, idummy, nb, &
  !     idummy, idummy, idummy, ibz, idummy, (high_sym_pnt(:, k), &
  !     high_sym_pnt_txt_dir(k)(1:ntxtd), k = 1, idummy), nopused, &
  !     iopnum(:)
  ! if (iok /= 0) nopused = -1
  ! if (ibz /= 0) then
  !     call exit_on_error(srcname // &
  !                        ': RIXS does not work with extended BZ mesh.')
  ! endif

if (.not. check_kstar) return

call allocate(mme, 'mme array' // srcname, udim1 = 3, udim2 = max_nhsort, &
     udim3 = n_mme, udim4 = nb)
do k = 1, npnt
  if (k == 1) then
    read(unit = rim, rec = k, iostat = iok) mme(:,:,:,:)
    if (iok /= 0) then
      write(unit = output_unit, fmt = 35) k
      call exit_on_error(srcname)
    endif
  endif
  if (k > 1) then
    call g4ig(iopnum(k), op)
!      op = transpose( op )
    write(unit = output_unit, fmt = 10) iopnum(k)
    write(unit = output_unit, fmt = 15) op(1,:)
    write(unit = output_unit, fmt = 15) op(2,:)
    write(unit = output_unit, fmt = 15) op(3,:)
  endif
  write(unit = output_unit, fmt = 11) k
  do iat = 1, max_nhsort
    do ia = 1, n_mme
      do ib = 1, nb
        mme_kibz(:) = mme(:, iat, ia, ib)
        mme_kbz(:)  = mme_kibz(:)
        if (k > 1) then
          if (iopnum(k) > 0) then
            mme_kbz(:) = matmul(op, mme_kibz)
!              mme_kbz(:) = matmul( mme_kibz, op )
          else
            mme_kbz(:) = matmul(op, conjg(mme_kibz))
!              mme_kbz(:) = matmul( conjg( mme_kibz ), op )
!              mme_kbz(:) = matmul( mme_kibz, op )
          endif
        endif
        write(unit = output_unit, fmt = 12) iat, ia, ib, mme_kbz
        cycle
        call tensor_from_mme(mme_kbz, tensor)
        do isp = 1, 3
          write(unit = output_unit, fmt = 20) tensor(isp,:)
        enddo
     enddo
    enddo
  enddo
enddo
call deallocate(mme, 'mme' // srcname)
call exit_on_error

 20 format(3(' (',2(f15.10),')'))

END SUBROUTINE perform_check_kstar

SUBROUTINE allocate_lmtdata_mme
use m_aux, only: map
use m_bnd, only: nb, npnt
use m_bzmesh, only: nkbz
use m_functions, only: allocate, int2string
use m_rixs, only: lmtdata, nlmtdata
implicit none
character(*), parameter :: srcname = ' in READ_RIXS_DATA'
integer isp, ia1, ian, iat1, iatn

do isp = 1, nlmtdata
!  if (.not. lmtdata(isp)%needed_in_RIXS) cycle
  if (sum(map(lmtdata(isp)%ia1, : )) == 0) cycle
  iat1 = lmtdata(isp)%iat1
  iatn = lmtdata(isp)%iatn
  ia1  = lmtdata(isp)%ia1
  ian  = lmtdata(isp)%ian
  call allocate(lmtdata(isp)%mme, 'lmtdata(' // trim(int2string(isp)) // &
       ')%mme' // srcname, udim1 = 3, ldim2 = iat1, udim2 = iatn, &
       ldim3 = ia1, udim3 = ian, udim4 = nb, udim5 = npnt)
  call allocate(lmtdata(isp)%mmebz, 'lmtdata(' // trim(int2string(isp)) // &
       ')%mmebz' // srcname, udim1 = 3, ldim2 = iat1, udim2 = iatn, &
       ldim3 = ia1, udim3 = ian, udim4 = nb, udim5 = nkbz)
enddo

END SUBROUTINE allocate_lmtdata_mme

SUBROUTINE read_mme_for_IBZ
use, intrinsic :: iso_fortran_env
use m_aux, only: map
use m_bnd
use m_files, only: rim
use m_functions
use m_params, only: iprint
use m_rixs
implicit none
character(*), parameter :: srcname = ' in READ_RIXS_MME'
complex(REAL64), allocatable :: mme(:, :, :, :)
integer iok, ia1, ian, iat1, iatn, isp, k

if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Ready to read matrix elements for IBZ'
endif
do k = 1, npnt
  call allocate(mme, 'mme for k= ' // trim(int2string(k)) // srcname, &
       udim1 = 3, udim2 = max_nhsort, udim3 = n_mme, udim4 = nb)
  read(unit = rim, rec = k, iostat = iok) mme(:, :, :, :)
  if (iok /= 0) then
    write(unit = output_unit, fmt = 35) k
    call exit_on_error(srcname)
  endif
  do isp = 1, nlmtdata
!    if (.not. lmtdata(isp)%needed_in_RIXS) cycle
    if (sum(map(lmtdata(isp)%ia1, : )) == 0) cycle
    iat1 = lmtdata(isp)%iat1
    iatn = lmtdata(isp)%iatn
    ia1  = lmtdata(isp)%ia1
    ian  = lmtdata(isp)%ian
    lmtdata(isp)%mme(:, iat1 : iatn,            ia1 : ian, :, k) = &
                 mme(:,      : iatn - iat1 + 1, ia1 : ian, :)
  enddo
!  call analyze_mme(k, mme)
  call deallocate(mme, 'mme for k= ' // trim(int2string(k)) // srcname)
enddo

 10 format(4x, a)
 35 format(4x, 'In READRIXSMME: cannot read from mmefile for k= ', i0)

END SUBROUTINE read_mme_for_IBZ

SUBROUTINE set_mme_for_BZ
use, intrinsic :: iso_fortran_env
use m_aux
use m_bnd, only: nb, g, iopnum
use m_bz
use m_bzmesh, only: nkbz
use m_functions, only: allocate, int2string
use m_params, only: iprint
use m_rixs
use omp_lib
implicit none
character(*), parameter :: srcname = ' in READ_RIXS_MME'
complex(REAL64) mme_kibz(3), mme_kbz(3)
integer ja, iat, ib, isp, iop, kbz
real(REAL64) op(3, 3)

if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Ready to set matrix elements for BZ'
endif

!$$$ if ( nopused > 1 ) then
!$$$  call allocate( isp4k , 'isp4k' // srcname, udim1 = nkbz * nlmtdata )
!$$$  call allocate( kbz4k , 'kbz4k' // srcname, udim1 = nkbz * nlmtdata )  

  do isp = 1, nlmtdata
!    if (.not. lmtdata(isp)%needed_in_RIXS) cycle
    if (sum( map(lmtdata(isp)%ia1, : ) ) == 0) cycle
!$$$    isp4k( (isp - 1) * nkbz + 1 : isp * nkbz ) = isp
!$$$    do k = 1, nkbz
!$$$      isp4k( (isp - 1) * nkbz + k ) = isp
!$$$      kbz4k( (isp - 1) * nkbz + k ) = k
!$$$    enddo
  enddo
! $omp parallel default(shared) private(isp, kbz, iop, op, ib, ja, iat, &
! $omp&         mme_kibz, mme_kbz)
! $omp do schedule(dynamic, 1)
  do kbz = 1, nkbz ! * nlmtdata
  ! do isp = 1 , nlmtdata
  !   do kbz = 1 , nkbz
    ! isp = isp4k(k)
    ! kbz = kbz4k(k)
    iop = ig4qibz(i1(kbz), i2(kbz), i3(kbz))
    op(:,:) = g(:,:, iop)
    do ib = 1, nb
!$$$      do ja = lmtdata(isp)%ia1, lmtdata(isp)%ian
      do ja = 1, n_mme
        if (sum(map(ja, :)) == 0)cycle
        isp = lmtindex(ja)
        do iat = lmtdata(isp)%iat1, lmtdata(isp)%iatn
          mme_kibz(:) = lmtdata(isp)%mme(:, iat, ja, ib, ik2ibz(kbz))
          ! if ( iop == 1 ) then
          !   mme_kbz(:) = mme_kibz(:)
          ! else
          !   if ( iopnum(iop) > 0 ) then
          !     mme_kbz(:) = matmul( op, mme_kibz )
          !   else
          !     mme_kbz(:) = matmul( op, conjg( mme_kibz ) )
          !   endif
          ! endif
          if (iopnum(iop) > 0) then
            mme_kbz(:) = matmul(op, mme_kibz)
          else
            mme_kbz(:) = matmul(op, conjg(mme_kibz))
          endif
          lmtdata(isp)%mmebz(:, iat, ja, ib, kbz) = mme_kbz(:)
        enddo  ! iat
      enddo    ! ja
    enddo      ! ib
  enddo        ! k
! $omp end do
! $omp end parallel
!$$$ endif   ! nopused > 1

 10 format(4x, a)

END SUBROUTINE set_mme_for_BZ
