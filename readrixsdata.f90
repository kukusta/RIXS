SUBROUTINE read_RIXS_data
use m_params !, only: irixs, iabsorp, read_rixfile
implicit none
character(*), parameter :: srcname = ' in READRIXSDATA'

if (read_rixfile) return

call open_datafile
call read_COMMON_DATA
! Joint DOS without absorption => no mme needed
if (irixs == 0 .and. iabsorp == 0 .and. .not. new_inr) return
call set_nlmtdata
call read_data_to_lmtdata
call search_for_duplicates
call read_MISC
call set_map_array
! call set_needed_in_RIXS
! min=huge(0._REAL64)
! do isp=1,nrixs
!   call set_ncf_ncl(ncf,ncl)
!   do ib=ncf,nb
!     do k=1,npnt
!       diff=e(ib,k)-rixs(isp)%ec
!       if(diff < min)then
!         ib_min=ib
!         k_min=k
!       endif
!     enddo
!   enddo
! enddo
! if(iprint > 0) write(unit=output_unit,fmt=45)
!  45 format(4x,'Resonant energies have been calculated.')

END SUBROUTINE read_RIXS_data

SUBROUTINE open_datafile
use, intrinsic :: iso_fortran_env
use m_files, only: rid, datafile
use m_functions, only: exit_on_error
implicit none
character(18), parameter :: srcname = ' in READ_RIXS_DATA'
integer iok
logical exist

inquire(file = datafile, exist = exist)
if (.not. exist) then
  write(unit = output_unit, fmt = 5) trim(datafile)
  call exit_on_error(srcname)
endif
open(newunit = rid, file = datafile, action = 'read', status = 'unknown', &
     position = 'rewind', iostat = iok)
if(iok /= 0) then
  write(unit = output_unit, fmt = 10) 'Can not open datafile ' // trim(datafile)
  call exit_on_error(srcname)
endif

 5 format(/, 4x, 'File ', a, ' does not exist.')
10 format(4x, a)

END SUBROUTINE open_datafile

SUBROUTINE read_COMMON_DATA
use m_constants, only: BUFFER_SIZE, eps5
use m_files, only: rid
use m_functions, only: exit_on_error, allocate
use m_params
use m_rixs
implicit none
character(*), parameter :: srcname = ' in READ_RIXS_DATA'
character(BUFFER_SIZE) msg
integer iok

namelist /common_data/ nrixslmt, ispc, max_nhsort, nlabel, n_mme, &
     emin, emax, de, nfu, volomg

read(unit = rid, nml = common_data, iomsg = msg, iostat = iok)
if (iok /= 0) then
  write(unit = output_unit, fmt = 5) trim(msg)
  call exit_on_error(srcname)
endif

 5 format(4x, 'Error in datafile in COMMON_DATA section: ', a)

if (abs(dw) < eps5) dw = de
if (abs(wmin) < eps5 .and. abs(wmax) < eps5) then
  wmin = 0._REAL64
  wmax = emax - emin
  if (int((wmax - wmin) / dw) + 1 >= huge(0)) then
    msg(:) = 'ERROR: number of points on omega mesh is out of range for INT32'
    write(unit = output_unit, fmt = 10) trim(msg)
  endif
  nw = int((wmax - wmin) / dw) + 1
endif

10 format(4x, a)

END SUBROUTINE read_COMMON_DATA

SUBROUTINE set_nlmtdata
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
use m_files, only: rid
use m_functions
use m_rixs, only: nlmtdata, nrixslmt
implicit none
character(4) sname
character(6) txtel
character(*), parameter :: srcname = ' in READ_RIXS_DATA'
character(BUFFER_SIZE) msg
integer iok, isp, nn, nk, writemme, iat1, iatn, label, ia1, ian, isort
logical lstop, skip
real(REAL64) ec, enlmin, enlmax, dnhsort

namelist /static/ nn, nk, writemme, iat1, iatn, label, txtel, sname, ec, &
     enlmin, enlmax, dnhsort, ia1, ian, isort, skip

nlmtdata = 0
lstop = .false.
do isp = 1, nrixslmt
  skip = .false.
  read(unit = rid, nml = static, iomsg = msg, iostat = iok)
  if (iok == IOSTAT_END) exit
  if (iok /= 0) then
    write(unit = output_unit, fmt = 15) isp, trim(msg)
    lstop = .true.
  endif
  if (skip) write(unit = output_unit, fmt = 20) isp
  if (writemme == 0 .or. skip .or. lstop) cycle
  nlmtdata = nlmtdata + 1
enddo
if (nlmtdata == 0) then
  write(unit = output_unit, fmt = 10) 'ERROR: No valid data found in datafile'
  lstop = .true.
endif
if (lstop) call exit_on_error(srcname)

 10 format(4x, a)
 15 format(4x, 'Error in datafile in STATIC section for ISPEC=', i0, ': ', a)
 20 format(4x, 'WARNING: you are going to skip LMTO matrix elements' // &
          ' for ISPEC=', i0,' (skip= T in RIXS datafile)')

END SUBROUTINE set_nlmtdata

SUBROUTINE read_data_to_lmtdata
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
use m_files, only: rid
use m_functions
use m_rixs, only: nlmtdata, lmtdata, nrixslmt
implicit none
character(4) sname
character(6) txtel
character(16), parameter :: srcname = ' in READRIXSDATA'
character(BUFFER_SIZE) msg
integer iok, isp, nn, nk, writemme, iat1, iatn, label, ia1, ian, isort, ilmtdata
logical lstop, skip
real(REAL64), allocatable :: esp(:), rat(:,:)
real(REAL64) ec, enlmin, enlmax, dnhsort

namelist /static/ nn, nk, writemme, iat1, iatn, label, txtel, sname, ec, &
     enlmin, enlmax, dnhsort, ia1, ian, isort, skip
namelist /esp_rat/ esp, rat

rewind(rid)
call allocate(lmtdata, 'lmttemp' // srcname, udim1 = nlmtdata)
lstop = .false.
ilmtdata = 0
do isp = 1, nrixslmt
  skip = .false.
  read(unit = rid, nml = static, iomsg = msg, iostat = iok)
  if (iok == IOSTAT_END) exit
  if (writemme == 0 .or. skip) cycle
  ilmtdata = ilmtdata + 1
  lmtdata(ilmtdata)%nn             = nn
  lmtdata(ilmtdata)%nk             = nk
  lmtdata(ilmtdata)%iwritemme      = writemme
  lmtdata(ilmtdata)%ia1            = ia1
  lmtdata(ilmtdata)%ian            = ian
  lmtdata(ilmtdata)%iat1           = iat1
  lmtdata(ilmtdata)%iatn           = iatn
  lmtdata(ilmtdata)%label          = label
  lmtdata(ilmtdata)%txtel(:)       = txtel(:)
  lmtdata(ilmtdata)%sname(:)       = sname(:)
  lmtdata(ilmtdata)%ec             = ec
  lmtdata(ilmtdata)%enlmin         = enlmin
  lmtdata(ilmtdata)%enlmax         = enlmax
  lmtdata(ilmtdata)%dnhsort        = dnhsort
  lmtdata(ilmtdata)%isort          = isort
!  lmtdata(ilmtdata)%needed_in_RIXS = .false.
  call allocate(esp, 'esp' // srcname, ldim1 = ia1, udim1 = ian)
  call allocate(rat, 'rat' // srcname, udim1 = 3, ldim2 = iat1, udim2 = iatn)
  read(unit = rid, nml = esp_rat, iomsg = msg, iostat = iok)
  if (iok /= 0) then
    write(unit = output_unit, fmt = 30) isp, trim(msg)
    lstop = .true.
  endif
  if (lstop) then
    call deallocate(esp, 'esp' // srcname)
    call deallocate(rat, 'rat' // srcname)
    cycle
  endif
  call allocate(lmtdata(ilmtdata)%esp, 'lmtdata(' // trim(int2string(ilmtdata))&
       // ')%esp' //srcname, ldim1 = ia1, udim1 = ian)
  call allocate(lmtdata(ilmtdata)%rat, 'lmtdata(' // trim(int2string(ilmtdata))&
       // ')%rat' //srcname, udim1 = 3, ldim2 = iat1, udim2 = iatn)
  lmtdata(ilmtdata)%esp(:)   = esp(:)
  lmtdata(ilmtdata)%rat(:,:) = rat(:,:)
  call deallocate(esp, 'esp' // srcname)
  call deallocate(rat, 'rat' // srcname)
enddo
if (lstop) call exit_on_error(srcname)

 10 format(4x, a)
 15 format(4x, 'Error in datafile in STATIC section for ISPEC=', i0, ': ', a)
 20 format(4x, 'WARNING: you are going to skip LMTO matrix elements' // &
          ' for ISPEC=', i0,' (skip= T in RIXS datafile)')
 30 format(4x, 'Error in datafile in ESP_RAT section for ISPEC=', i0, ': ', a)

END SUBROUTINE read_data_to_lmtdata

SUBROUTINE search_for_duplicates
use, intrinsic :: iso_fortran_env
use m_functions, only: exit_on_error
use m_rixs, only: nlmtdata, lmtdata
implicit none
character(18), parameter :: srcname = ' in READ_RIXS_DATA'
integer isp, jsp
logical lstop

do isp = 1, nlmtdata
  do jsp = isp + 1, nlmtdata
    if (lmtdata(jsp)%nn    == lmtdata(isp)%nn              .and. &
        lmtdata(jsp)%nk    == lmtdata(isp)%nk              .and. &
        lmtdata(jsp)%isort == lmtdata(isp)%isort           .and. &
        lmtdata(jsp)%iwritemme*lmtdata(isp)%iwritemme /= 0       &
      ) then
      write(unit = output_unit, fmt = 5) isp, jsp
      lstop = .true.
    endif
  enddo
enddo
if (lstop) then
  write(unit = output_unit, fmt = 10) 'To continue remove duplicates' // &
     ' in LMTO input and recalculate matrix elements'
  call exit_on_error
endif

 5 format(   4x, 'Duplicated mme found in datafile: ISPEC=', i2, ' and ', i2)
10 format(/, 4x, a)

END SUBROUTINE search_for_duplicates

SUBROUTINE read_MISC
use, intrinsic :: iso_fortran_env
use m_bnd, only: natom, nopused
use m_constants, only: BUFFER_SIZE
use m_files, only: rid, mmefile
use m_functions, only: allocate, exit_on_error
use m_sdt, only: igtta
implicit none
character(BUFFER_SIZE) msg, mmefile_save
character(18), parameter :: srcname = ' in READ_RIXS_DATA'
integer iok

namelist /misc/ igtta, mmefile

call allocate(igtta, 'igtta' // srcname, udim1 = nopused, udim2 = natom)
! mmefile from datafile is used only if it has not been set in inptfile
mmefile_save(:)=""
if (len_trim(mmefile) /= 0) mmefile_save(:) = mmefile(:)
read(unit = rid, nml = misc, iomsg = msg, iostat = iok)
if (iok /= 0) then
  write(unit = output_unit, fmt = 35) trim(msg)
  call exit_on_error(srcname)
endif
if (len_trim(mmefile_save) /= 0) mmefile(:) = mmefile_save(:)

 35 format(4x, 'Error in datafile in MISC section: ', a)

END SUBROUTINE read_MISC

! SUBROUTINE set_needed_in_RIXS
! use m_rixs
! implicit none
! character(4) sname_RIXS, sname_LMT
! character(2) txtel_RIXS, txtel_LMT
! integer isp, jsp

! do isp = 1, nrixs
!    txtel_RIXS(:) = rixs(isp)%txtel(:2)
!    sname_RIXS(:) = rixs(isp)%sname(:)
!    do jsp = 1, nlmtdata
!       txtel_LMT(:) = lmtdata(jsp)%txtel(:2)
!       sname_LMT(:) = lmtdata(jsp)%sname(:)
!       if (sname_RIXS  == sname_LMT .and. txtel_RIXS == txtel_LMT) then
!         lmtdata(jsp)%needed_in_RIXS = .true.
!       endif
!   enddo
! enddo

! END SUBROUTINE set_needed_in_RIXS

LOGICAL FUNCTION SAVED(sname,txtel,isort)
use m_rixs, only: lmtdata, nlmtdata
implicit none
!input vars
character(4) sname
character(6) txtel
integer isort
! local vars
integer isp

saved=.false.
do isp=1,nlmtdata
  if(lmtdata(isp)%sname(:) == sname(:) .and. &
     lmtdata(isp)%txtel(:2) == txtel(:2) .and. &
     lmtdata(isp)%isort == isort) saved=.true.
enddo

END FUNCTION SAVED

!$$$ close(unit=rid,iostat=iok)
! this check is performed in READRIXSMME subroutine
! if(.not. newini)then
!   inquire(file=mmefile,exist=exist)
!   if(.not. exist)then
!     write(unit=output_unit,fmt=5)trim(mmefile)
!     call deallocate_global_arrays
!     stop
!   endif
! endif
! call allocate(lmtdata, 'lmtdata' // srcname, udim1 = )
! allocate(lmtdata(nrixslmt2), stat = iok, errmsg = msg)
! if (iok /= 0) call print_allocation_error('lmtdata', iok, msg)

! ! used to analyze absorption edges below
! rixs(1:nrixs)%ec=1._REAL64

! nlmtdata=0
! do isp = 1, nrixs
!   first = .true.
!   do jsp = 1, nrixslmt2
!     if (rixs(isp)%sname(:) == lmttemp(jsp)%sname(:)) then
!       if (rixs(isp)%txtel(:2) == lmttemp(jsp)%txtel(:2)) then
!         rixs(isp)%nn    = lmttemp(jsp)%nn
!         rixs(isp)%nk    = lmttemp(jsp)%nk
!         rixs(isp)%label = lmttemp(jsp)%label
! ! these four lines are not clear. When RIXS is calculated for several atom
! ! sorts with different ec, dnhsort, enlmax values these values do not have
! ! unambigous values. They can be different for different sorts.
!         rixs(isp)%dnhsort = lmttemp(jsp)%dnhsort
!         rixs(isp)%enlmax  = lmttemp(jsp)%enlmax
!         rixs(isp)%isort   = lmttemp(jsp)%isort
! ! sometimes there are several nonequivalent atoms with the same Z
! ! here warning is printed if this takes place
!         if (rixs(isp)%ec > 0._REAL64) then
! ! always lmt(jsp)%ec < 0
!           rixs(isp)%ec = lmttemp(jsp)%ec
!           txtel(:)     = lmttemp(jsp)%txtel(:)
!         else
!           if (first) then
!             first = .false.
!             write(unit = output_unit, fmt = 50) isp, txtel(:)
!             write(unit = output_unit, fmt = 55, advance = 'no') txtel(:)
!           endif
!           write(unit = output_unit, fmt = 60, advance = 'no') lmttemp(jsp)%txtel(:)
!         endif

!  50 format(/, 4x, 'For SPEC=', i0, " the absorption edge is defined for '", a, "'")
!  55 format(4x, "Inequivalent species with the same Z are '", a, "'")
!  60 format(", '", a, "'")

! ! if for some spectrum input energy is less than 0 it is given in eV
! ! relative to E_F. If it is positive it is given in absolute units.
! ! It's value in absolute units is usually large and difficult to define
! ! If input energy is exactly zero, then it means some mistake (error
! ! message is printed).
! !$$$ ! very small negative value, error for -0.001 eV < input energy < 0 eV
! !$$$ ! if input energy is not set in inrfile than by default input energy is
! !$$$ ! set equal to -0.0001 eV and causes an error here. input energy must
! !$$$         if(rixs(isp)%en > -7.e-5_REAL64 .and. rixs(isp)%en < -1.e-6_REAL64)then
! !$$$           write(unit=output_unit,fmt=37)isp
! !$$$           call deallocate_global_arrays
! !$$$ 37 format(/, 4x, 'For ISPEC=', i0, ' the energy of incoming x-ray is not set.')
! !$$$           stop
! !$$$         endif
! ! negative or zero input energy, relative to E_F
!         if(rixs(isp)%en < 1.e-6_REAL64)then
!           rixs(isp)%en=ef-rixs(isp)%ec-rixs(isp)%en
!         endif
!         if(.not. saved(lmttemp(jsp)%sname,lmttemp(jsp)%txtel,lmttemp(jsp)%isort))then
!           nlmtdata=nlmtdata+1
!           ia1=lmttemp(jsp)%ia1
!           ian=lmttemp(jsp)%ian
!           iat1=lmttemp(jsp)%iat1
!           iatn=lmttemp(jsp)%iatn
!           allocate(lmtdata(nlmtdata)%esp(ia1:ian)&
!                   ,lmtdata(nlmtdata)%rat(3,iat1:iatn)&
!                   ,lmtdata(nlmtdata)%mme(3,iat1:iatn,ia1:ian,nb,npnt)&
!                   ,stat=iok,errmsg=msg)
!           if(iok /= 0)call print_allocation_error('lmtdata%esp,rat',iok,msg)
!           lmtdata(nlmtdata)%sname(:)=lmttemp(jsp)%sname(:)
!           lmtdata(nlmtdata)%txtel(:)=lmttemp(jsp)%txtel(:)
!           lmtdata(nlmtdata)%nn=lmttemp(jsp)%nn
!           lmtdata(nlmtdata)%nk=lmttemp(jsp)%nk
!           lmtdata(nlmtdata)%isort=lmttemp(jsp)%isort
!           lmtdata(nlmtdata)%ia1=lmttemp(jsp)%ia1
!           lmtdata(nlmtdata)%ian=lmttemp(jsp)%ian
!           lmtdata(nlmtdata)%iat1=lmttemp(jsp)%iat1
!           lmtdata(nlmtdata)%iatn=lmttemp(jsp)%iatn
!           lmtdata(nlmtdata)%label=lmttemp(jsp)%label
!           lmtdata(nlmtdata)%ec=lmttemp(jsp)%ec
!           lmtdata(nlmtdata)%enlmin=lmttemp(jsp)%enlmin
!           lmtdata(nlmtdata)%enlmax=lmttemp(jsp)%enlmax
!           lmtdata(nlmtdata)%dnhsort=lmttemp(jsp)%dnhsort
!           lmtdata(nlmtdata)%esp(:)=lmttemp(jsp)%esp(:)
!           lmtdata(nlmtdata)%rat(:,:)=lmttemp(jsp)%rat(:,:)
!         endif
!       endif
!     endif
!   enddo
!   if (.not. first) write(unit = output_unit, fmt =" ()")
! enddo

! stop = .false.
! do isp = 1, nrixs
!   if (rixs(isp)%nn == 0) then
!     write(unit = output_unit, fmt = 40) rixs(isp)%txtel(:2), trim(rixs(isp)%sname)
!     stop = .true.
!   endif
! enddo
! if(stop) stop

!  40 format(4x,'No mme have been found for ',a,' at ',a,'-edge')

! if (.not. newini) then
!   do isp = 1, nlmttemp
!     deallocate(lmttemp(isp)%esp,lmttemp(isp)%rat,stat=iok)
!   enddo
!   deallocate(lmttemp,stat=iok)
! endif

! ! do isp=1,nrixs
! !   ! rixs(isp)%q(:)=matmul(ut1,rixs(isp)%q)
! !   ! rixs(isp)%q(:)=matmul(ut,rixs(isp)%q)
! !   rixs(isp)%q(:)=matmul(qbas,rixs(isp)%q)/2._REAL64/dpi()
! ! enddo

