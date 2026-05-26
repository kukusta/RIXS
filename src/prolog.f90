SUBROUTINE prolog
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
implicit none
character(BUFFER_SIZE) buf
character(6) version

call get_command_argument(1, buf)
if (command_argument_count() < 1 .or. &
    trim(buf) == '-h' .or. trim(buf) == '--help') then
  call usage
endif

call s4ver(version)
write(unit = output_unit, fmt = 35) version(:)

 35 format(  4x, 'RESONANT INELASTIC X-RAY SCATTERING CALCULATIONS', &
          /, 4x, 'based on PY-LMTO band structure code Ver. ', a, /)

END SUBROUTINE prolog

SUBROUTINE usage
use, intrinsic :: iso_fortran_env

if (command_argument_count() < 1) then
  write(unit = output_unit, fmt = 10) 'ERROR: No input file found.'
endif
write(unit = output_unit, fmt = 10) 'USAGE: rixs <inputfile.inr>'
write(unit = output_unit, fmt = 12) 'rixs -h or --help prints this message'
stop

10 format(4x, a)
12 format(11x, a)

END SUBROUTINE usage
