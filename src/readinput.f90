SUBROUTINE readinput
use, intrinsic :: iso_fortran_env
use m_params, only: new_inr
use m_rixs, only: rixs
implicit none
! integer i

! character(4) rixs(isp)%Reg_BZ
! character(6) rixs(isp)%Reg_BZ_states
! integer rixs(isp)%Reg_BZ_size(3), rixs(isp)%Reg_BZ_center_qbmc(3)
! real(REAL64) rixs(isp)%Reg_BZ_center(3)

call get_cmd_arguments_and_basename
call set_defaults_for_FILES
if (new_inr) call create_input  ! create new input file and stop
call open_inputfile
call read_PARAMETERS
call normalize_vectors
call check_data_in_PARAMETERS
call print_calculation_mode
call get_number_of_spectra
call get_data_for_RIXS_spectra
! Later I use 'nn' field to check if mme for ispec where written to mmefile
! rixs(:)%nn=0
rixs(:)%ntr     = 0
rixs(:)%min_mme = 0._REAL64
rixs(:)%max_mme = 0._REAL64
rixs(:)%av_mme  = 0._REAL64

! if (.not. (irixs == 0 .or. abs(irixs) == 10)) then
!   do i = 1, nrixs
!     if (rixs(i)%Gamma < eps2) &
!       write(unit = output_unit, fmt = 85) i, rixs(i)%Gamma * dry2ev
!   enddo
!   endif

!  85 format(4x, 'For ISPEC=', i0, ' Gamma value of', f6.4, ' eV is too small')

END SUBROUTINE READINPUT

SUBROUTINE get_cmd_arguments_and_basename
use, intrinsic :: iso_fortran_env
use m_constants
use m_files, only: basename
use m_functions, only: exit_on_error
use m_params, only: new_inr
implicit none
character(BUFFER_SIZE) buf
character(15), parameter :: srcname = ' in READINPUT'
integer ilmt, iinr, iarg

call get_command_argument(1, buf)
ilmt = index(buf, '.lmt', back = .true.)
iinr = index(buf, '.inr', back = .true.)
basename(:) = ""
if (ilmt /= 0) then
  basename(: ilmt) = buf(: ilmt)
  new_inr = .true.
else
  if (iinr /= 0) then
    basename(: iinr) = buf(: iinr)
    new_inr = .false.
  else
    write(unit = output_unit, fmt = 10) &
         'RIXS or LMTO inputfile must always be the first argument.'
    call exit_on_error(srcname)
  endif
endif
do iarg = 2, command_argument_count()
  call get_command_argument(iarg, buf)
enddo

10 format(4x, a)

END SUBROUTINE get_cmd_arguments_and_basename

SUBROUTINE set_defaults_for_FILES
use m_files
implicit none
! Initially all files are closed
bnd = 0
bns = 0
inp = 0
rid = 0
sdt = 0
rix = 0
rim = 0
m2d_q = 0
m2d_ein = 0
! inptfile or lmtfile is always the first argument
inptfile(:)     = trim(basename) // 'inr'
bndfile(:)      = trim(basename) // 'bnd'
bnsfile(:)      = trim(basename) // 'bns'
datafile(:)     = trim(basename) // 'rid'
sdtfile(:)      = trim(basename) // 'sdt'
rixsfile(:)     = trim(basename) // 'rix'
map2deinfile(:) = trim(basename) // 'e2d'
map2dqfile(:)   = trim(basename) // 'q2d'
! Empty filename means that mmefile would be read from datafile.
! To be reread it should be empty after reading inptfile as well,
! i.e. 'mmefile' variable in &PARAMETERS section should be commented.
mmefile(:)      = ""

END SUBROUTINE set_defaults_for_FILES

SUBROUTINE open_inputfile
use, intrinsic :: iso_fortran_env
use m_files, only: inptfile, inp
implicit none
integer iok
logical exist

inquire(file = inptfile, exist = exist)
if (.not. exist) then
  write(unit = output_unit , fmt = 12) trim(inptfile)
  stop
endif
open(newunit = inp, file = inptfile, action = 'read', position = 'rewind', &
     status = 'unknown', iostat = iok)
if (iok /= 0) stop '    Input file cannot be opened, exiting.'

 12 format(4x, 'In READINPUT: file ', a, ' does not exist.')

END SUBROUTINE open_inputfile

SUBROUTINE read_PARAMETERS
use, intrinsic :: iso_fortran_env
use m_files
use m_functions, only: exit_on_error
use m_params
implicit none
character(BUFFER_SIZE) buf
character(BUFFER_SIZE), external :: to_lowercase
character(*), parameter :: srcname = ' in read_PARAMETERS'
character(4) Reg_BZ
character(6) Reg_BZ_states
integer iok, nbi(2), nbf(2), Reg_BZ_size(3)
logical map2d_ein, map2d_q, save2dat
real(REAL64) Reg_BZ_center(3)

namelist /parameters/ irixs, iabsorp, iprint, nbi, nbf, wmin, wmax, dw,    &
     rixsfile, bndfile, bnsfile, datafile, sdtfile, mmefile, map2deinfile, &
     map2dqfile, q_axis, k_in, absorption_only, save2dat, map2d_ein,& !, A_bas
     map2d_q, check_kstar, read_rixfile, Reg_BZ, Reg_BZ_size, Reg_BZ_center, &
     ntr_step, negative_energy_loss, Reg_BZ_states, &
     allow_negative_RIXS_values, switch_ini_final

iabsorp    = 0
iprint     = 1
irixs      = 1
ntr_step   = 10
wmin       = 0.
wmax       = 0.
dw         = 0.
q_axis(:)  = 0.
! A_bas(:,:) = 0.
k_in(:)    = 0.
nbi(:)     = 0
nbf(:)     = 0
! forall(i = 1:3) A_bas(i,i) = 1.
Reg_BZ(:)            = 'none'
Reg_BZ_size(:)       = 0
Reg_BZ_center(:)     = 0.
Reg_BZ_states(:)     = 'below '
absorption_only            = .false.
save2dat                   = .false.
map2d_ein                  = .false.
map2d_q                    = .false.
check_kstar                = .false.
negative_energy_loss       = .false.
allow_negative_RIXS_values = .false.
read_rixfile               = .false.
switch_ini_final           = .false.
read(unit = inp, nml = parameters, iomsg = buf, iostat = iok)
if (iok /= 0) then
  write(unit = output_unit, fmt = 10) 'Error in input file ' // &
       'in &PARAMETERS: ' // trim(buf)
  call exit_on_error(srcname)
endif
nbi_global(:)                 = nbi(:)
nbf_global(:)                 = nbf(:)
save2dat_global               = save2dat
map2d_ein_global              = map2d_ein
map2d_q_global                = map2d_q
Reg_BZ_global(:)        = to_lowercase(Reg_BZ)
Reg_BZ_size_global(:)   = Reg_BZ_size(:)
Reg_BZ_center_global(:) = Reg_BZ_center(:)
Reg_BZ_states_global(:) = to_lowercase(Reg_BZ_states)

 10 format(4x, a)

END SUBROUTINE read_PARAMETERS

SUBROUTINE normalize_vectors
use m_constants, only: eps3
use m_params, only: k_in, q_axis! , A_bas
implicit none

if (norm2(k_in     )  > eps3) k_in(:)    = k_in(:)    / norm2(k_in     )
if (norm2(q_axis   )  > eps3) q_axis(:)  = q_axis(:)  / norm2(q_axis   )
! if (norm2(A_bas(:,1)) > eps3) A_bas(:,1) = A_bas(:,1) / norm2(A_bas(:,1))
! if (norm2(A_bas(:,2)) > eps3) A_bas(:,2) = A_bas(:,2) / norm2(A_bas(:,2))
! if (norm2(A_bas(:,3)) > eps3) A_bas(:,3) = A_bas(:,3) / norm2(A_bas(:,3))

END SUBROUTINE normalize_vectors

SUBROUTINE print_calculation_mode
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
use m_params
implicit none
character(BUFFER_SIZE) buf

if (read_rixfile) then
  write(unit = output_unit, fmt = 10) 'Read rixsfile and write it again only'
  return
endif

buf(:) = 'disabled'
if (irixs == 0 .or. irixs == 10) buf(:) = 'enabled'
write(unit = output_unit, fmt = 40) trim(buf)

buf(:) = 'disabled'
if (irixs /= 0 .and. irixs /= 10) buf(:) = 'enabled'
write(unit = output_unit, fmt = 45) trim(buf)

buf(:) = 'disabled'
if (irixs == 2) buf(:) = 'enabled'
write(unit = output_unit, fmt = 50) trim(buf)

! buf(:) = 'disabled'
! if (irixs > 0) buf(:) = 'enabled'
! write(unit = output_unit, fmt = 55) trim(buf)

buf(:) = 'disabled'
if (irixs == 3) buf(:) = 'enabled'
write(unit = output_unit, fmt = 57) trim(buf)

buf(:) = 'disabled'
if (iabsorp == 1 .or. iabsorp == 2) buf(:) = 'enabled'
write(unit = output_unit, fmt = 60) trim(buf)
buf(:) = 'disabled'
if (iabsorp == 2) buf(:) = 'enabled'
write(unit = output_unit,fmt = 65) trim(buf)

 10 format(4x, a)
 40 format(4x, 'Joint DOS spectra',      15x, a)
 45 format(4x, 'RIXS loss spectra',      15x, a)
 50 format(4x, 'RIXS with averaged mme', 10x, a)
! 55 format(4x, 'RIXS fractional translations', 4x, a)
 57 format(4x, 'RIXS Lorentzian',        17x, a)
 60 format(4x, 'Absorption over IBZ',    13x, a)
 65 format(4x, 'Absorption over BZ ',    13x, a, /)

END SUBROUTINE print_calculation_mode

SUBROUTINE get_number_of_spectra
use, intrinsic :: iso_fortran_env
use m_constants, only: MAX_NRIXS
! use m_params, only: read_rixfile
use m_rixs, only: rixs_spectrum, nrixs
implicit none
integer nspecread
logical eof, skip, stop
type(rixs_spectrum) DTO

! if (read_rixfile) return

nrixs = 0
nspecread = 0
stop = .false.
eof = .false.
do while (.not. eof)
  skip = .false.
  nspecread = nspecread + 1
  call read_SPEC(nspecread, DTO, eof, skip, stop)
  if (eof) exit
  if (.not. skip) nrixs = nrixs + 1
enddo
if (stop) stop
if (nrixs == 0) then
  write(unit = output_unit, fmt = 10) 'No valid data in input file found.'
  stop
endif
if (nrixs > MAX_NRIXS) then
  write(unit = output_unit, fmt = 90) nrixs, MAX_NRIXS
  stop
endif

 10 format(4x, a)
 90 format(4x, 'Number of RIXS spectra is too large: nrixs=', i0, &
          ' MAX_NRIXS=', i0)

END SUBROUTINE get_number_of_spectra

SUBROUTINE get_data_for_RIXS_spectra
use, intrinsic :: iso_fortran_env
use m_files, only: inp
use m_functions, only: allocate, exit_on_error
use m_rixs, only: rixs, nrixs
implicit none
character(*), parameter :: srcname = ' in READ_INPUT'
integer isp, nspecread
logical eof, skip, stop

! if (read_rixfile) return

call allocate(rixs, 'rixs' // srcname, udim1 = nrixs)
rewind(inp)
isp  = 1
eof  = .false.
skip = .false.
do
  call read_SPEC(nspecread, rixs(isp), eof, skip, stop)
  if (skip) cycle
  if (eof) return
  call set_nn_nk(isp)
  if (isp == nrixs) return
  isp = isp + 1
enddo

END SUBROUTINE get_data_for_RIXS_spectra

SUBROUTINE read_SPEC(nspecread, DTO, eof, skip, stop)
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
use m_files, only: inp
use m_functions, only: int2string
use m_params
use m_rixs, only: rixs_spectrum
use units, only: dry2ev
implicit none
character(BUFFER_SIZE) buf
character(BUFFER_SIZE), external :: to_lowercase
character(6) atom, Reg_BZ_states
character(4) specname, channel, use_symmetry, Reg_BZ
character(*), parameter :: srcname = ' in read_SPEC'
integer iok, nspecread, nbi(2), nbf(2), Reg_BZ_size(3)
logical eof, skip, stop, map2d_ein, map2d_q, save2dat
real(REAL64) peak_near_E(16), q(3), input_energy, Gamma, e_in(3), e_out(3)
real(REAL64) Reg_BZ_center(3)
type(rixs_spectrum) DTO

namelist /spec/ atom, specname, q, input_energy, Gamma, channel, nbi, &
     nbf, skip, e_in, e_out, save2dat, map2d_ein, map2d_q, use_symmetry, &
     Reg_BZ, Reg_BZ_size, Reg_BZ_center, Reg_BZ_states, peak_near_E

atom           = ""
specname       = ""
e_in(:)        = 0._REAL64
e_out(:)       = 0._REAL64
q(:)           = 0._REAL64
input_energy   = 0._REAL64
Gamma          = 0._REAL64
channel        = ""
nbi(:)         = nbi_global(:)
nbf(:)         = nbf_global(:)
use_symmetry   = 'n'
skip           = .false.
save2dat       = save2dat_global
map2d_ein      = map2d_ein_global
map2d_q        = map2d_q_global
peak_near_E(:) = -1._REAL64
Reg_BZ(:)        = Reg_BZ_global(:)
Reg_BZ_size(:)   = Reg_BZ_size_global(:)
Reg_BZ_center(:) = Reg_BZ_center_global(:)
Reg_BZ_states(:) = Reg_BZ_states_global(:)
read(unit = inp, nml = spec, iomsg = buf, iostat = iok)
if (iok == 0 .and. skip) return
if (iok == IOSTAT_END) then
  eof = .true.
  return
endif
if (iok /= 0) then
  write(unit = output_unit, fmt = 5) nspecread, trim(buf)
  stop = .true.
  return
endif
if (input_energy < 0.) then
  buf(:) = 'For ISPEC= ' // trim(int2string(nspecread)) // ' input_energy < 0' &
        // ' (is set relative to E_Fermi)'
  write(unit = output_unit, fmt = 10) trim(buf)
  buf(:) = 'Right now I do not know how to do this unambiguously'
  write(unit = output_unit, fmt = 10) trim(buf)
  buf(:) = 'Please use absolute value for input energy'
  write(unit = output_unit, fmt = 10) trim(buf)
  stop = .true.
  return
endif
DTO%txtel(:)               = atom(:)
DTO%sname(:)               = specname(:)
DTO%q(:)                   = q(:)
DTO%en                     = input_energy / dry2ev
DTO%Gamma                  = Gamma / dry2ev
DTO%channel(:)             = channel(:)
DTO%nbi(:)                 = nbi(:)
DTO%nbf(:)                 = nbf(:)
DTO%e_in(:)                = cmplx(e_in, kind = REAL64)
DTO%e_out(:)               = cmplx(e_out, kind = REAL64)
DTO%use_symmetry(:)        = use_symmetry(:)
DTO%save2dat               = save2dat
DTO%map2d_ein              = map2d_ein
DTO%map2d_q                = map2d_q
DTO%Reg_BZ(:)        = to_lowercase(Reg_BZ)
DTO%Reg_BZ_size(:)   = Reg_BZ_size(:)
DTO%Reg_BZ_states(:) = to_lowercase(Reg_BZ_states)
DTO%Reg_BZ_center(:) = Reg_BZ_center(:)
DTO%peak_near_E(:)   = peak_near_E(:) / dry2ev

  5 format(4x, 'Error in input file in &SPEC for ISPEC= ', i0, ': ', a)
 10 format(4x, a)

END SUBROUTINE read_SPEC

SUBROUTINE check_data_in_PARAMETERS
use, intrinsic :: iso_fortran_env
use m_constants, only: eps3, eps5
use m_params
use units, only: dry2ev
implicit none
real(REAL64), external :: ddet33

! if (read_rixfile) return

if (nbi_global(1) < 0) stop '    Error in &PARAMETERS: nbi(1) < 0.'
if (nbi_global(2) < 0) stop '    Error in &PARAMETERS: nbi(2) < 0.'
if (nbi_global(1) > nbi_global(2)) then
  stop '    Error in &PARAMETERS: nbi(1) > nbi(2).'
endif
if (nbf_global(1) < 0) stop '    Error in &PARAMETERS: nbf(1) < 0.'
if (nbf_global(2) < 0) stop '    Error in &PARAMETERS: nbf(2) < 0.'
if (nbf_global(1) > nbf_global(2)) then
  stop '    Error in &PARAMETERS: nbf(1) > nbf(2).'
endif
if (irixs < 0) then
  write(unit = output_unit, fmt = 15) 'In &PARAMETERS: irixs < 0. ' // &
       'Set negative_energy_loss= T instead.', &
       ' irixs=', abs(irixs), ' would be used in further calculations.'
  irixs = abs(irixs)
  negative_energy_loss = .true.
endif
if (irixs > 3 .and. irixs /= 10) then
  write(unit = output_unit, fmt = 20) irixs
  stop
endif
if (.not. (irixs == 0 .or. irixs == 10)) then
  if (norm2(q_axis) < eps3) then
    write(unit = output_unit, fmt = 10) &
         'ERROR in READINPUT: length of q_axis is too small'
    stop
  endif
  ! if (ddet33(A_bas) < 1._REAL64 - eps5) then
  !   write(unit = output_unit, fmt = 10) 'ERROR in READINPUT:' // &
  !        ' matrix A does not define a valid Cartesian coordinate system.'
  !   stop
  ! endif
endif
if (iabsorp < 0 .or. iabsorp > 2) then
  write(unit = output_unit, fmt = 10) &
       'ERROR in READINPUT: iabsorp < 0 or iabsorp > 2 found.'
  stop
endif
if (wmin < 0._REAL64) stop '    ERROR in READINPUT: wmin < 0.'
if (wmin > wmax) stop '    ERROR in READINPUT: wmin > wmax.'
! if energies are given in eV
if (dw < 0._REAL64) then
  dw   =  -dw / dry2ev
  wmin = wmin / dry2ev
  wmax = wmax / dry2ev
endif
if (wmin > 0._REAL64 .or. wmax > 0._REAL64) then
  if (abs(dw) < eps5) then
    write(unit = output_unit, fmt = 10) &
         'ERROR in READINPUT: wmin or wmax > 0 but dw = 0.'
    stop
  endif
  if ((wmax - wmin) / dw >= real(huge(0), kind = REAL64)) then
    write(unit = output_unit, fmt = 29) (wmax - wmin) / dw
    stop
  endif
  nw = int((wmax - wmin) / dw) + 1
endif
if (norm2(k_in) < eps3 .and. .not. (irixs == 0 .or. irixs == 10)) then
  write(unit = output_unit, fmt = 10) 'Input vector k_in = 0.'
  stop
endif

 10 format(4x, a)
 15 format(4x, a, /, a, i2, a)
 20 format(4x, 'ERROR in READINPUT: wrong irixs vaule found: irixs=', i2)
 29 format(4x, 'ERROR in READINPUT: number of points on omega mesh '&
          ,e11.4,' is out of range for INT32.')

END SUBROUTINE check_data_in_PARAMETERS

SUBROUTINE set_nn_nk(isp)
use m_rixs, only: rixs
implicit none
integer isp, nn_index, nk_index
character    :: tn(1:5)  = [ 'K', 'L', 'M', 'N', 'O' ]
character(3) :: tk(-4:3) = ['VII','V  ','III','I  ','   ','II ','IV ','VI ']

do nn_index = 1, 5
  if (rixs(isp)%sname(1:1) == tn(nn_index)) rixs(isp)%nn = nn_index
enddo
do nk_index = -4, 3
  if (rixs(isp)%sname(2:4) == tk(nk_index)) rixs(isp)%nk = nk_index
enddo
if (rixs(isp)%nk == 0) rixs(isp)%nk = -1

END SUBROUTINE set_nn_nk
