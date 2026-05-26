SUBROUTINE PRINT_TIME_STATISTICS
use m_constants, only: BUFFER_SIZE
use m_params, only: irixs, iabsorp, absorption_only, read_rixfile
use m_time
implicit none
! local vars
character(BUFFER_SIZE) fmt
integer i, hh, mm, lower, upper
integer, parameter :: NMAX=7
integer, external :: memory_usage
real(REAL64) time(NMAX), ss
character(33) text(NMAX)

if (read_rixfile) return

write(unit = output_unit, fmt = 5) max(memory_usage(), 1)
text(:)=[ 'time to read mmefile             ', &
          'real time for RIXS               ', &
          'CPU time for RIXS                ', &
          'real time for absorption over IBZ', &
          'CPU time for absorption over IBZ ', &
          'real time for absorption over RBZ', &
          'CPU time for absorption over RBZ '  ]
if (irixs == 0 .or. irixs == 10) then
  text(2) = 'real time for Joint DOS          '
  text(3) = 'CPU time for Joint DOS           '
endif
time(:)=[ get_readmme(),         &
          get_real_RIXS(),       &
          get_CPU_RIXS(),        &
          get_real_ABSORP_IBZ(), &
          get_CPU_ABSORP_IBZ(),  &
          get_real_ABSORP_RBZ(), &
          get_CPU_ABSORP_RBZ()   ]

write(unit = output_unit, fmt = 10)

  5 format(/, 4x, 'MEMORY USAGE: max_used_memory= ', i0, ' MB.')
 10 format(/, 4x, 'TIME STATISTICS')

lower=1
if((irixs == 0 .or. irixs == 10) .and. iabsorp == 0)lower=2
if (iabsorp == 1) then
  upper = 5
else if (iabsorp == 2) then
  upper = 7
else
  upper = 3
endif

do i = lower, upper
  if (i == 2 .and. absorption_only) cycle
  if (i == 3 .and. absorption_only) cycle
  hh = floor(time(i) / 3600.)
  mm = floor(time(i) - real(hh) * 3600.) / 60.
  ss = time(i) - real(hh) * 3600. - real(mm) * 60.
  fmt(:) = "(4x, '" // text(i) // "', 2x, i3, ' hh ', i2, ' mm ', f5.2, ' ss')"
  write(unit = output_unit, fmt = fmt) hh, mm, ss
enddo

END SUBROUTINE PRINT_TIME_STATISTICS
