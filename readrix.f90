SUBROUTINE READRIX
use m_files
use m_functions
use m_params
use m_results, only: loss
use m_rixs, only: nrixs, nlabel
implicit none
character(11), parameter :: srcname=' in READRIX'
integer i, iok, isp
logical exist
real wmin4, wmax4, dw4
real, allocatable :: loss2(:)

if( .not. read_rixfile ) return

inquire( file = rixsfile, exist = exist )
if( .not. exist ) then
  write( unit = output_unit, fmt = &
       "( 4x, 'File ', a, ' does not exist' )" ) trim( rixsfile )
  return
endif

open( newunit = rix, file = rixsfile, action = 'read', status = 'unknown', &
      form = 'unformatted', iostat = iok )
if( iok /= 0 ) stop '    Error while opening rixfile.'

read(rix) nrixs, nlabel, iabsorp
iabsorp = iabsorp - 1
read(rix) nw, wmin4, wmax4, dw4
wmin = real( wmin4, kind = REAL64)
wmax = real( wmax4, kind = REAL64)
dw = real(dw4, kind = REAL64)
read(rix) ! nerixs,real(emin),real(emax),real(de),real(ef)

call allocate( loss, 'loss' // srcname, udim1 = nw, udim2 = nrixs )
call allocate( loss2, 'loss2' // srcname, udim1 = nw )

do isp = 1, nrixs
  read(rix) ! rixs(isp)%label, rixs(isp)%isort, real( rixs(isp)%ec ), &
            ! real( rixs(isp)%dnhsort ), real( rixs(isp)%enlmax ), &
            ! rixs(isp)%nn, rixs(isp)%nk, rixs(isp)%txtel(:)
  read( rix, iostat = iok ) loss2(:)
  if( iok /= 0 ) then
    write( unit = output_unit, fmt = &
         "( 4x, 'Can not read rixs loss for ISPEC #', i0 )" ) isp
    exit
  endif
  loss( : , isp ) = real( loss2(:), kind = REAL64 )
  do i = 1, iabsorp
    read(rix) ! real( absorp(:,isp,i) )
  enddo
enddo

call deallocate( loss2, 'loss2' // srcname )

END SUBROUTINE READRIX
