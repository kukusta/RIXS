SUBROUTINE readrix
use m_bnd, only: ef
use m_files
use m_functions
use m_params
use m_results
use m_rixs, only: nrixs, rixs, nlabel, de, emin, emax
implicit none
character(*), parameter :: srcname = ' in READ_RIXFILE'
integer i, iok, isp, nerixs
logical negative_energy_loss_in_rixfile
real wmin_4, wmax_4, dw_4, emin_4, emax_4, de_4, ef_4, ec_4, dnhsort_4, enlmax_4
real, allocatable :: rixs_array(:), absorp_array(:)

if (.not. read_rixfile) return

call open_rixsfile
read(rix) nrixs, nlabel, iabsorp
iabsorp = iabsorp - 1
read(rix) nw, wmin_4, wmax_4, dw_4
if (real(wmin_4, kind = REAL64) /= wmin .or. real(wmax_4, kind = REAL64) &
     /= wmax .or. real(dw_4, kind = REAL64) /= dw) then
  write(unit = output_unit, fmt = 10) 'wmin or wmax, or dw in inputfile' // &
       ' and rixfile differ. Their value from rixfile will be used'
endif
call allocate(loss, 'loss' // srcname, udim1 = nw, udim2 = nrixs)
negative_energy_loss_in_rixfile = .false.
if (wmin_4 < 0.) negative_energy_loss_in_rixfile = .true.
if (negative_energy_loss_in_rixfile) then
  wmin = -real(wmax_4, kind = REAL64)
  wmax = -real(wmin_4, kind = REAL64)
else
  wmin = real(wmin_4, kind = REAL64)
  wmax = real(wmax_4, kind = REAL64)
endif
dw = real(dw_4, kind = REAL64)
call allocate(rixs_array, 'rixs_array' // srcname, udim1 = nw)
read(rix) nerixs, emin_4, emax_4, de_4, ef_4
emin = real(emin_4, kind = REAL64)
emax = real(emax_4, kind = REAL64)
  de = real(de_4,   kind = REAL64)
  ef = real(ef_4,   kind = REAL64)
if (iabsorp >= 1) then
  call allocate(absorp, 'absorp' // srcname, udim1 = nerixs, udim2 = nrixs, &
       udim3 = iabsorp)
  call allocate(absorp_array, 'absorp_array' // srcname, udim1 = nerixs)
endif
if (negative_energy_loss_in_rixfile) then
  call allocate(loss_neg_w, 'loss_neg_w' // srcname, udim1 = nw, udim2 = nrixs)
endif
do isp = 1, nrixs
  read(rix) rixs(isp)%label, rixs(isp)%isort, ec_4, dnhsort_4, enlmax_4, &
            rixs(isp)%nn, rixs(isp)%nk, rixs(isp)%txtel(:)
  rixs(isp)%ec = real(ec_4, kind = REAL64)
  rixs(isp)%dnhsort = real(dnhsort_4, kind = REAL64)
  rixs(isp)%enlmax = real(enlmax_4, kind = REAL64)
  read(rix, iostat = iok) rixs_array(:)
  if (iok /= 0) then
    write(unit = output_unit, fmt = 10) 'ERROR: I can not read rixs loss' // &
         ' for ISPEC=' // trim(int2string(isp))
    call exit_on_error(srcname)
  endif
  if (negative_energy_loss_in_rixfile) then
    loss_neg_w(:, isp) = real(rixs_array(:), kind = REAL64)
    do i = 1, nw
      loss(i, isp) = loss_neg_w(nw - (i - 1), isp)
    enddo
  else
    loss(:, isp) = real(rixs_array(:), kind = REAL64)
  endif
  do i = 1, iabsorp
    read(rix) absorp_array(:)
    absorp(:, isp, i) = real(absorp_array(:), kind = REAL64)
  enddo
enddo
call deallocate(rixs_array, 'rixs_array' // srcname)
if (iabsorp >= 1) call deallocate(absorp_array, 'absorp_array' // srcname)

10 format(4x, a)

END SUBROUTINE readrix

SUBROUTINE open_rixsfile
use, intrinsic :: iso_fortran_env
use m_files, only: rix, rixsfile
use m_functions, only: exit_on_error
implicit none
character(*), parameter :: srcname=' in READRIX'
integer iok
logical exist
  
inquire (file = rixsfile, exist = exist)
if (.not. exist) then
  write(unit = output_unit, fmt = 10) 'File ' // trim(rixsfile) // &
       ' does not exist' 
  call exit_on_error(srcname)
endif
open(newunit = rix, file = rixsfile, action = 'read', status = 'unknown', &
     form = 'unformatted', iostat = iok)
if (iok /= 0) call exit_on_error(srcname)

10 format(4x, a)

END SUBROUTINE open_rixsfile
