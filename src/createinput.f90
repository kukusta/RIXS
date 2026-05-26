SUBROUTINE create_input
implicit none
! local vars
integer ncl, ncf
logical exist

call inquire_or_create_inputfile(exist)
call read_bnd_and_rixs_data(ncl, ncf)
call create_header
call write_PARAMETERS(ncl, ncf)
call write_SPECS(ncl, ncf)
call write_hsym_points_to_inputfile
call print_epilog_and_stop(exist)

END SUBROUTINE create_input

SUBROUTINE inquire_or_create_inputfile(exist)
use, intrinsic :: iso_fortran_env
use m_files, only: inp, inptfile
use m_params, only: iprint
implicit none
character key, key2
integer iok
logical exist

inquire (file = inptfile, exist = exist)
iprint = 1
key = 'N'
key2 = ' '
if (exist) then
  do
    write(unit = output_unit, fmt = 1, advance = 'no') trim(inptfile)
    read(unit = input_unit, fmt = '(a1)') key2
    if (key2 /= ' ') key = key2
    if (key == 'y' .or. key == 'Y') exit
    if (key == 'n' .or. key == 'N') then
      write(unit = output_unit, fmt = 10) trim(inptfile)
      stop
    endif
  enddo
endif
open(newunit = inp, file = inptfile, action = 'write', status = 'unknown', &
     iostat = iok)
if (iok /= 0) stop '    Input file cannot be opened, exiting.'

  1 format(4x, 'File ', a, ' exists. Do you want to overwrite it? [y/N] ')
 10 format(4x, 'File ', a, ' has not been changed.')

END SUBROUTINE inquire_or_create_inputfile

SUBROUTINE read_bnd_and_rixs_data(ncl, ncf)
use m_bnd, only: nb
use m_params, only: iprint, check_kstar
implicit none
integer ncl, ncf

! get data from main lmto program and set some defaults
if (.not. check_kstar) then
  call readbnd
  iprint = -1
  call set_ncf_ncl(ncf, ncl)
  iprint = 1
else
  ncl = 0
  ncf = 0
  nb  = 0
endif
call read_RIXS_data

END SUBROUTINE read_bnd_and_rixs_data

SUBROUTINE create_header
use m_files, only: inp
implicit none
character(10) date, time
character(80) hostname
integer iok

! I switched this header off
if (.false.) then
  call date_and_time(date = date, time = time)
  write(unit = inp, fmt = 30)
  write(unit = inp, fmt = 35, advance = 'no') date(1:4), date(5:6), date(7:8), &
       time(1:2), time(3:4), time(5:6)
  call get_environment_variable('HOSTNAME', value = hostname, status = iok)
  if (iok == 0) then
    write(unit = inp, fmt = "(' on ',a,/)") trim(hostname)
  else
    write(unit = inp, fmt = "(2(/))")
  endif
endif

 30 format(4x, 'This input file for RIXS calculations was generated')
 35 format(4x, 'automatically on ', a, '-', a, '-', a, 1x, a, ':', a, ':', a)

END SUBROUTINE create_header

SUBROUTINE write_PARAMETERS(ncl, ncf)
use m_bnd, only: nb
use m_files
use m_params
use units, only: sry2ev
implicit none
integer ncl, ncf

irixs = 1
write(unit = inp, fmt = "(4x, '&PARAMETERS')")
write(unit = inp, fmt = 2) irixs, iabsorp, iprint
if (check_kstar) then
  write(unit = inp, fmt = 3) 0, ncl, ncf, nb
else
  write(unit = inp, fmt = 3) 1, ncl, ncf, nb
endif
write(unit = inp, fmt = 5) wmin * sry2ev, wmax * sry2ev, -dw * sry2ev
write(unit = inp, fmt = 6) trim(rixsfile), trim(map2deinfile), &
     trim(map2dqfile), trim(datafile), trim(bndfile), trim(sdtfile), &
     trim(mmefile)
write(unit = inp, fmt = 7) ! check_kstar
write(unit = inp, fmt = 120)

 2 format(4x, 'irixs= ', i0, 4x, 'iabsorp= ', i0, 4x, 'iprint= ', i0)
 3 format(4x, '! nbi(1:2)=', 2(i5), 4x, 'nbf(1:2)=', 2(i5))
 5 format(4x, '! wmin=', f8.3, 4x, 'wmax=', f10.6, 4x, 'dw=', f10.5)
 6 format(4x, '! rixsfile= "', a, '"', 4x, 'map2deinfile= "', a, '"', &
          4x, 'map2dqfile= "', a, '"', /, 4x, '! datafile= "', a, '"', &
          4x, 'bndfile= "', a, '"', 4x, 'sdtfile= "', a, '"', /, &
          4x, '! mmefile= "', a, '"')
 7 format(4x, 'k_in(:)= 0. -1. 1.', 4x, 'q_axis(:)= 0. 0. 1.')
!         , /, 4x, 'check_kstar= ', L )

 120 format(4x, '/', /)

END SUBROUTINE write_PARAMETERS

SUBROUTINE write_SPECS(ncl, ncf)
use m_bnd, only: nb
use m_files, only: inp
use m_params, only: check_kstar, wmin, wmax
use m_rixs, only: nlmtdata, lmtdata
use units, only: sry2ev
implicit none
integer isp, ncl, ncf
real Gamma, energy

do isp = 1, nlmtdata
  energy = ((wmax + wmin) / 2. - lmtdata(isp)%ec) * sry2ev
  call spwidth_cp01(lmtdata(isp)%txtel(:2), lmtdata(isp)%nn, &
       lmtdata(isp)%nk, Gamma)
  write(unit = inp, fmt = 10) '&SPEC'
  write(unit = inp, fmt = 12) lmtdata(isp)%txtel(:2), lmtdata(isp)%sname(:)
  if (check_kstar) then
    write(unit = inp, fmt = 3) 0, ncl, ncf, nb
  else
    write(unit = inp, fmt = 3) 1, ncl, ncf, nb
  endif
  write(unit = inp, fmt = 15) 0., 0., 0.
  write(unit = inp, fmt = 18) energy, Gamma
  write(unit = inp, fmt = 21)
  write(unit = inp, fmt = 120)
enddo

  3 format(4x, '! nbi(1:2)=', 2(i5), 4x, 'nbf(1:2)=', 2(i5))
 10 format(4x, a)
 12 format(4x, "atom= '", a, "'", 4x, "specname= '", a, "'", 4x, "skip= F", 4x,&
           "! use_symmetry= 'n'")
 15 format(4x, 'q(:)= ', 3(f7.3), 4x, '! in G1, G2, G3')
 18 format(4x, 'input_energy= ', f8.2, 4x, 'Gamma= ', f6.3, 4x, '! both in eV')
 21 format(4x, "channel= 'p-s '")
!  21 format(4x,'save2dat= F',4x,'map2d_ein= F',4x,'map2d_q= F')
 22 format(/                                                          &
    ,4x,'IRIXS:',/                                                    &
    ,4x,'0 - Joint DOS (BZ integration with constant mme)',/          &
    ,4x,'1 - RIXS loss spectra',/                                     &
    ,4x,'2 - RIXS with matrix elements averaged over polarizations',/,/ &
    ! ,4x,'3 - dispersion of RIXS loss spectra',/,/                     &
    ,4x,'IABSORP:',/                                                  &
    ,4x,'0 - absorption spectra are not calculated',/                 &
    ,4x,'1 - IBZ integration',/                                       &
    ,4x,'2 - IBZ and BZ integration')
 120 format(4x, '/', /)

END SUBROUTINE write_SPECS

SUBROUTINE write_hsym_points_to_inputfile
use m_bnd, only: high_sym_pnt, high_sym_pnt_txt_dir
use m_files, only: inp
use m_params, only:  write_hsym_points
implicit none
integer ipnt, ipnt2
logical cycle

write_hsym_points = .false.
if (write_hsym_points) then
  write(unit = inp, fmt = 25)
  do ipnt = 1, 48
    if (len_trim(high_sym_pnt_txt_dir(ipnt)) == 0) exit
    cycle = .false.
    do ipnt2 = 1, ipnt - 1
      if (high_sym_pnt_txt_dir(ipnt2) == high_sym_pnt_txt_dir(ipnt)) then
        cycle = .true.
        exit
      endif
    enddo
    if (cycle) cycle
    write(unit = inp, fmt = 45) high_sym_pnt(:,ipnt), high_sym_pnt_txt_dir(ipnt)
  enddo
endif

 25 format(/, 4x, 'HIGH SYMMETRY POINTS:')
 45 format(4x, 3(f14.10), 4x, "'", a, "'")

END SUBROUTINE write_hsym_points_to_inputfile

SUBROUTINE print_epilog_and_stop(exist)
use, intrinsic :: iso_fortran_env
use m_files, only: inptfile
use m_functions, only: exit_on_error
implicit none
logical exist

if (exist) then
  write(unit = output_unit, fmt = 9) trim(inptfile)
else
  write(unit = output_unit, fmt = 11) trim(inptfile)
endif
call exit_on_error

  9 format(4x, 'File ', a, ' has been overwritten.')
 11 format(4x, 'File ', a, ' has been created.')

END SUBROUTINE print_epilog_and_stop

