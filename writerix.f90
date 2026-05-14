SUBROUTINE write_results

call write_rixfile
call find_local_maxima
call write_datafiles
call write_2d_maps

END SUBROUTINE write_results

SUBROUTINE write_rixfile
use m_bnd, only: ef
use m_files, only: rix, rixsfile
use m_functions
use m_params
use m_results
use m_rixs
use units
implicit none
! local vars
character(BUFFER_SIZE) spec_type
integer channel_label, iok, i, isp, nerixs
logical exist

if (iprint > 0) write(unit = output_unit, fmt = 4) 'Starting WRITE_RIXFILE'
if (irixs /= 0 .and. irixs /= 10) then
  spec_type = 'RIXS loss'
else
  spec_type = 'Joint DOS'
endif

inquire (file = rixsfile, exist = exist)
if (exist) then
  write(unit = output_unit, fmt = 10) 'File ' // trim(rixsfile) // &
       ' exists and will be overwritten.'
endif
open(newunit = rix, file = rixsfile, status = 'replace', &
     form = 'unformatted', iostat = iok)
if (iok /= 0) then
  write(unit = output_unit, fmt = 10) 'ERROR: I can not open rixfile'
  call exit_on_error
endif
! if (iabsorp == 1 .or. iabsorp == 2) then
!   write(rix) nrixs, nlabel, 1 + iabsorp
! else
!   write(rix) nrixs, nlabel,1+0
! endif
write(rix) nrixs, nlabel, 1 + iabsorp
nerixs = int((emax - emin) / de + 1)
if (negative_energy_loss) then
  write(rix) nw, -real(wmax), -real(wmin), real(dw)
else
  write(rix) nw, real(wmin), real(wmax), real(dw)
endif
write(rix) nerixs, real(emin), real(emax), real(de), real(ef)
! check if there are negative values of RIXS loss intensity
do isp = 1, nrixs
  do i = 1, nw
    if (loss(i, isp) < 0._REAL64 .and. .not. allow_negative_RIXS_values) then
      if (abs(loss(i, isp) ) > eps5) then
        write(unit = output_unit, fmt = 5) isp, trim(spec_type), i, loss(i, isp)
      endif
      loss(i, isp) = 0._REAL64
    endif
  enddo
enddo
if (negative_energy_loss .and. .not. allocated(loss_neg_w)) then
  call allocate(loss_neg_w, 'loss_neg_w', udim1 = nw, udim2 = nrixs)
endif
!  write(unit = output_unit, fmt = "(/)")
do isp = 1, nrixs
  channel_label = 1
  if (rixs(isp)%channel == 's-p ') channel_label = 2
  if (rixs(isp)%channel == 'p-s ') channel_label = 3
  if (rixs(isp)%channel == 'p-p ') channel_label = 4
  if (rixs(isp)%channel == 'l-l ') channel_label = 5
  if (rixs(isp)%channel == 'l-r ') channel_label = 6
  if (rixs(isp)%channel == 'r-l ') channel_label = 7
  if (rixs(isp)%channel == 'r-r ') channel_label = 8
!  rixs(isp)%isort = 0
!  write(rix) rixs(isp)%label, rixs(isp)%isort, real(rixs(isp)%ec),    &
  write(rix) rixs(isp)%label, channel_label, real(rixs(isp)%ec),    &
       real(rixs(isp)%dnhsort), real(rixs(isp)%enlmax), rixs(isp)%nn, &
       rixs(isp)%nk, rixs(isp)%txtel(:)
  if (negative_energy_loss) then
    do i = 1, nw
      loss_neg_w(i, isp) = loss(nw - (i - 1), isp)
    enddo
    write(rix) real(loss_neg_w(:, isp))
  else
    write(rix) real(loss(:, isp))
  endif
  do i = 1, iabsorp
    write(rix) real(absorp(:, isp, i))
  enddo
enddo
if (iprint > 0) write(unit = output_unit, fmt = 10) 'Rixsfile has been written'

  4 format(/, 4x, a)
  5 format(4x, 'For ISPEC=', i2, 1x, a, '(', i0, ')= ', f9.5, &
          ' and would be set to 0.' )
  9 format(4(4x, f11.6), 4x, f20.6)
 10 format(4x, a)

! if ( .not. ( irixs == 0 .or. abs(irixs) == 10 ) ) then
!   write( unit = output_unit , fmt = 30 )
!   do isp = 1 , nrixs
!     write( unit = output_unit , fmt = 40 ) isp , rixs(isp)%min_mme , &
!       rixs(isp)%max_mme , rixs(isp)%av_mme / real( rixs(isp)%ntr , &
!       kind = REAL64 ) , rixs(isp)%ntr
!   enddo
! endif
! 30 format(/, 4x, 'Maximal values of matrix elements:', /, &
!           4x, 'ISPEC',4x,'min(abs(mme))',4x,'max(abs(mme))', &
!           4x, 'avr(abs(mme))', 4x, 'rixs%ntr')
! 40 format(4x, i3, 4x, f0.10, 4x, f0.10, 4x, f0.10, 4x, i13 )

END SUBROUTINE write_rixfile

SUBROUTINE find_local_maxima
use, intrinsic :: iso_fortran_env
use m_constants, only: BUFFER_SIZE
use m_functions
use m_params
use m_results, only: loss
use m_rixs
use units, only: dry2ev
implicit none
character(BUFFER_SIZE) msg
integer imax, isp, ipnt, icenter, ipos
real(REAL64) v1, v2, v3, v4, v5

do isp = 1, nrixs
  do imax = 1, 16 ! size(rixs(isp)%peak_near_E)
    if (rixs(isp)%peak_near_E(imax) < 0.) cycle
    icenter = nint((rixs(isp)%peak_near_E(imax) - wmin) / dw) + 1
    if (icenter < 0 .or. icenter > nw) cycle
    do ipnt = -15, 15
      ipos = icenter + ipnt
      if (ipos < 5 .or. ipos > nw - 5) cycle
      v1 = loss(ipos - 2, isp)
      v2 = loss(ipos - 1, isp)
      v3 = loss(ipos - 0, isp)
      v4 = loss(ipos + 1, isp)
      v5 = loss(ipos + 2, isp)
      if (v1 < v2 .and. v2 < v3 .and. v3 > v4 .and. v4 > v5) then
        write(msg, fmt = 5) (wmin + dw * (ipos - 1)) * dry2ev
        msg(:) = 'Local maximum for ISPEC= ' // trim(int2string(isp)) // &
             ' at E=' // trim(msg) // ' eV'
        write(unit = output_unit, fmt = 10) trim(msg)
      endif
    enddo
  enddo
enddo
 5 format(f6.3)
10 format(4x, a)
END SUBROUTINE find_local_maxima

SUBROUTINE write_datafiles
use m_constants, only: BUFFER_SIZE
use m_files, only: rixsfile
use m_functions
use m_params
use m_results
use m_rixs, only: nrixs, rixs
use units
implicit none
character(BUFFER_SIZE) datfile
integer i, iok, isp, idat
logical exist

do isp = 1, nrixs
  if (.not. rixs(isp)%save2dat) cycle
  i = index(rixsfile, '.rix')
  datfile(:) = rixsfile(: i - 1) // "_ispec_" // trim(int2string(isp)) // ".dat"
  inquire(file = datfile, exist = exist)
  if (exist) write(unit = output_unit, fmt = 10) trim(datfile)
  open(newunit = idat, file = datfile, status = 'replace', iostat = iok)
  if (iok /= 0) then
    write(unit = output_unit, fmt = 10) &
         'ERROR: I can not open datfile ' // trim(datfile)
    call exit_on_error
  endif
! (rixs(isp)%en - ef + rixs(isp)%ec) * dry2ev,
  do i = 1, nw
    if (negative_energy_loss) then
      write(unit = idat, fmt = 8) -(wmin + (i - 1) * dw) * dry2ev, &
           real(loss_neg_w(i, isp))
    else
      write(unit = idat, fmt = 8) (wmin + (i - 1) * dw) * dry2ev, &
           real(loss(i, isp))
    endif
  enddo
  call close(unit = idat)
enddo

 8 format(4x, f11.6, 4x, f20.6)
10 format(4x, a)

END SUBROUTINE write_datafiles

SUBROUTINE write_2d_maps
! use m_bnd, only: ef
use m_files
use m_functions
use m_params
use m_results
use m_rixs
use units
implicit none
! local vars
! character(BUFFER_SIZE) spec_type
integer iok, i, isp ! , nerixs
logical exist
real maxRIXS, maxRIXSispec, maxRIXSenergy

if (m2d_ein ==0 .and. m2d_q == 0) return
if (irixs == 1 .or. irixs == 2) then
! open m2d files only if some spectra should be saved
  m2d_ein = 0
  m2d_q   = 0
  do isp = 1, nrixs
    if (rixs(isp)%map2d_ein) m2d_ein = m2d_ein + 1
    if (rixs(isp)%map2d_q  ) m2d_q   = m2d_q   + 1
  enddo
  if (m2d_ein > 0) then
    inquire(file = map2deinfile, exist = exist)
    if (exist) then
      write(unit = output_unit, fmt= 10) 'File ' // trim(map2deinfile) // &
           ' exists and will be overwritten'
    endif
    open(newunit = m2d_ein, file = map2deinfile, status = 'replace', iostat=iok)
    if(iok /= 0)then
      write(unit = output_unit, fmt = 10) &
           'ERROR: I can not open map2deinfile=' // trim(map2deinfile)
    endif
    write(unit = m2d_ein, fmt = 2)
    write(unit = output_unit, fmt = 10) &
         "Data for 2D map will be written to file " // trim(map2deinfile)
  endif
  if (m2d_q > 0) then
    inquire(file = map2dqfile, exist = exist)
    if (exist) then
      write(unit = output_unit, fmt = 10) &
           'File ' // trim(map2dqfile) // ' exists and will be overwritten'
    endif
    open(newunit = m2d_q, file = map2dqfile, status = 'replace', iostat = iok)
    if (iok /= 0) then
      write(unit = output_unit, fmt = 10) &
           'ERROR: I can not open map2dqfile=' // trim(map2dqfile)
    endif
        write(unit = m2d_q, fmt = 3)
    write(unit = output_unit, fmt = 10) &
         "Data for 2D map will be written to file " // trim(map2dqfile)
  endif
endif

2 format(  "# 1st column is energy loss", &
        /, "# 2nd column is normalized RIXS (max. RIXS is 1.0)")
3 format(  "# 1st column is energy loss", &
        /, "# 2nd column is qx, qy, qz", &
        /, "# 3rd column is normalized RIXS (max. RIXS is 1.0)")

! find max RIXS value for m2d file and check if there are
! negative values of RIXS loss intensity
maxRIXS = -1.e10
do isp = 1, nrixs
  maxRIXSispec = -1.e10
  maxRIXSenergy = 0.
  if (negative_energy_loss) then
    do i = 1, nw
      if (loss_neg_w(i, isp) > maxRIXS) maxRIXS = loss_neg_w(i, isp)
      if (loss_neg_w(i, isp) > maxRIXSispec) then
        maxRIXSispec  = loss_neg_w(i, isp)
        maxRIXSenergy = -real(wmin + (i - 1) * dw)
      endif
    enddo
  else
    do i = 1, nw
      if (loss(i, isp) > maxRIXS) then
        maxRIXS = loss(i, isp)
      endif
      if (loss(i, isp) > maxRIXSispec) then
        maxRIXSispec = loss(i, isp)
        maxRIXSenergy = real(wmin + (i - 1) * dw)
      endif
    enddo
  endif
  write(unit = output_unit, fmt = 7) isp, maxRIXSispec, maxRIXSenergy * sry2ev
enddo

7 format(4x, 'Maximum value for ISPEC=', i0, ' is ', f15.5, &
         ' at E= ', f10.5, ' eV' )
8 format(4x, f11.6, 4x, f20.6)
9 format(4(4x, f11.6), 4x, f20.6)
! 10 format(4x,'Datefile ',a,' exists and will be overwritten')

if (irixs == 1 .or. irixs == 2 .and. (m2d_ein /= 0 .or. m2d_q /= 0)) then
  write(unit = output_unit, fmt = 10) 'Data for 2D plot would be normalized.'
endif
do isp = 1, nrixs
! save data to map2d files
  if (irixs == 1 .or. irixs == 2) then
    if (rixs(isp)%map2d_ein) then
      do i = 1, nw
        if (negative_energy_loss) then
          write(unit = m2d_ein, fmt = 8) -(wmin + (i - 1) * dw) * dry2ev, &
               real(loss_neg_w(i, isp)) / maxRIXS
        else
          write(unit = m2d_ein, fmt = 8) (wmin + (i - 1) * dw) * dry2ev, &
               real(loss(i, isp)) / maxRIXS
        endif
      enddo
    endif
    if (rixs(isp)%map2d_q) then
      do i = 1, nw
        if (negative_energy_loss) then
          write(unit = m2d_q, fmt = 9) -(wmin + (i - 1) * dw) * dry2ev, &
               rixs(isp)%q(:), real(loss_neg_w(i, isp)) / maxRIXS
        else
          write(unit = m2d_q, fmt = 9) (wmin + (i - 1) * dw) * dry2ev, &
               rixs(isp)%q(:), real(loss(i,isp)) / maxRIXS
        endif
      enddo
    endif
  endif
enddo

 10 format(4x, a)
 15 format( 4x, 'Can not deallocate loss2(:), error code is ', i0)
 30 format( /, 4x, 'Maximal values of matrix elements:', /, &
           4x, 'ISPEC',4x,'min(abs(mme))',4x,'max(abs(mme))', &
           4x, 'avr(abs(mme))', 4x, 'rixs%ntr')
 40 format( 4x, i3, 4x, f0.10, 4x, f0.10, 4x, f0.10, 4x, i13 )

END SUBROUTINE write_2d_maps
