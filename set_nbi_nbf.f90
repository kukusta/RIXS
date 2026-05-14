SUBROUTINE SET_NBI_NBF
use, intrinsic :: iso_fortran_env
use m_bzmesh
use m_functions
use m_params, only: iprint, nbi_global, nbf_global
use m_bnd, only: nb
use m_rixs, only: nrixs, rixs
implicit none
! local vars
character(*), parameter :: srcname = ' in SET_NBI_NBF'
integer ncf, ncl, isp, i
logical lstop

! nspin == 1 and ndimspin == 1 is assumed
! ncf is the lowest NOT fully occupied band
! ncl the last NOT completely empty band

lstop = .false.
call set_ncf_ncl(ncf, ncl)

  5 format(4x, 'Correct nbi(1) from input data.')
 10 format(4x, 'Error in input data: nbi(1)=', i0, ' is too large.')
 15 format(4x, 'Correct nbi(2) from input data.')
 20 format(4x, 'Error in input data: nbi(2)=', i0, ' is too large.')
 25 format(4x, 'Error in input data: nbf(1)=', i0, ' is too small.')
 30 format(4x, 'Correct nbf(1) from input data.')
 35 format(4x, 'Error in input data: nbf(1)=', i0, ' is too large.')

! check if input data are consistent with ncl, ncf.
! check nbi(1)
if (nbi_global(1) == 0) then
  nbi_global(1) = 1
else
  if (nbi_global(1) >= 1 .and. nbi_global(1) <= ncl) then
    if (iprint > 0) write(unit = output_unit, fmt = 5)
  endif
  if (nbi_global(1) >= ncl + 1) then
    write(unit = output_unit, fmt = 10) ncl
    lstop = .true.
  endif
endif
! check nbi(2)
if (nbi_global(2) == 0) then
  nbi_global(2) = ncl
else
  if (nbi_global(2) >= 1 .and. nbi_global(2) <= ncl) then
    if (iprint > 0) write(unit = output_unit, fmt = 15)
  endif
  if (nbi_global(2) >= ncl+1) then
    write(unit = output_unit, fmt = 20) nbi_global(2)
    lstop = .true.
  endif
endif
! check nbf(1)
if (nbf_global(1) == 0) then
  nbf_global(1) = ncf
else
  if (nbf_global(1) >= 1 .and. nbf_global(1) <= ncf-1) then
    write(unit = output_unit, fmt = 25) nbf_global(1)
    lstop = .true.
  endif
  if (nbf_global(1) >= ncf .and. nbf_global(1) <= nb) then
    if (iprint > 0) write(unit = output_unit, fmt = 30)
  endif
  if (nbf_global(1) >= nb + 1) then
    write(unit = output_unit, fmt = 35) nbf_global(1)
    lstop = .true.
  endif
endif
! check nbf(2)
if (nbf_global(2) == 0) then
  nbf_global(2) = nb
else
  if (nbf_global(2) >= 1 .and. nbf_global(2) <= ncf-1) then
    write(unit = output_unit, fmt = 40) nbf_global(2)
    lstop = .true.
  endif
  if (nbf_global(2) >= ncf .and. nbf_global(2) <= nb) then
    if (iprint > 0) write(unit = output_unit, fmt = 45)
  endif
  if (nbf_global(2) >= nb + 1) then
    write(unit = output_unit, fmt = 50) nbf_global(2)
    lstop = .true.
  endif
endif

 40 format(4x, 'Error in input data: nbf(2)=', i0, ' is too small.')
 45 format(4x, 'Correct nbf(2) from input data.')
 50 format(4x, 'Error in input data: nbf(2)=', i0, ' is too large.')
 55 format(4x, 'Recommended band range from bndfile: nbi(1:2) = 1..', &
           i0, ' nbf(1:2) = ', i0, '..', i0, '.')

if (lstop) then
  write(unit = output_unit, fmt = 55) ncl, ncf, nb
  call exit_on_error(srcname)
endif

 60 format(4x, 'Error in input data for SPEC=', i0, ': nbi(', i0, ') < 0.')
 65 format(4x, 'Error in input data for SPEC=', i0, ': nbi(', i0,   &
          ') is too large. Maximum nbi(', i0, ')=', i0)
 70 format(4x, 'Error in input data for SPEC=', i0, ': nbf(', i0, ') < 0.')

do isp = 1, nrixs
  do i = 1, 2
    if (rixs(isp)%nbi(i) < 0) then
      write(unit = output_unit, fmt = 60) isp, i
      lstop = .true.
    endif
    if (rixs(isp)%nbi(i) == 0) rixs(isp)%nbi(i) = nbi_global(i)
    if (rixs(isp)%nbi(i) > 0) then
      if (rixs(isp)%nbi(i) > nbi_global(2)) then
        write(unit = output_unit, fmt = 65) isp, i, i, nbi_global(2)
        lstop = .true.
      endif
    endif
    if (rixs(isp)%nbf(i) < 0) then
      write(unit = output_unit, fmt = 70) isp, i
      lstop = .true.
    endif
    if (rixs(isp)%nbf(i) == 0) rixs(isp)%nbf(i)=nbf_global(i)
    if (rixs(isp)%nbf(i) > 0) then
      if (rixs(isp)%nbf(i) < nbf_global(1)) then
        write(unit = output_unit, fmt = 75) isp, i, i, nbf_global(1)
        lstop = .true.
      endif
      if (rixs(isp)%nbf(i) > nbf_global(2)) then
        write(unit = output_unit, fmt = 80) isp, i, i, nbf_global(2)
        lstop = .true.
      endif
    endif
  enddo
  if (rixs(isp)%nbi(1) > rixs(isp)%nbi(2)) then
    write(unit = output_unit, fmt = 85) isp
    lstop = .true.
  endif
  if (rixs(isp)%nbf(1) > rixs(isp)%nbf(2)) then
    write(unit = output_unit, fmt = 90) isp
    lstop = .true.
  endif
enddo
if (lstop) call exit_on_error(srcname)

if (iprint > 10) then
  do isp = 1, nrixs
    write(unit = output_unit, fmt = 95) isp, rixs(isp)%nbi(:), rixs(isp)%nbf(:)
  enddo
endif

 75 format(4x, 'Error in input data for SPEC=', i0, ': nbf(', i0,  &
          ') is too small. Minimal nbf(', i0, ')=', i0)
 80 format(4x, 'Error in input data for SPEC=', i0, ': nbf(', i0,  &
          ') is too large. Maximum nbf(', i0, ')=', i0)
 85 format(4x, 'Error in input data for SPEC=', i0, ': nbi(1) > nbi(2).')
 90 format(4x, 'Error in input data for SPEC=', i0, ': nbf(1) > nbf(2).')
 95 format(4x, 'For ISPEC=', i3, '  nbi(1:2)= ', 2(i4, 2x), '  nbf(1:2)= ', &
          2(i4, 2x))

 ! 95 format(4x,'For ISPEC=',i0,' initial states: ',i4,' -',i4&
 !          ,', final states: ',i4,' -',i4)

END SUBROUTINE SET_NBI_NBF
