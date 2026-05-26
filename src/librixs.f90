SUBROUTINE X_INTP(iabem,dec,nf,de,ne,input,output)
! interpolates spectrum to a common omega mesh. For absorption
! (iabem == -1) and emission (iabem == 1) only transitions
! to unoccupied or occupied states are considered.
! de can be dexmcd .or. derexs. or. derixs
! ne can be nexmcd .or. nerexs. or. nerixs
use, intrinsic :: iso_fortran_env
implicit none
! input and output
integer iabem, nf, ne
real(REAL64) de, dec, input(ne), output(ne)
! local vars
integer i, i1, in, j
real(REAL64) x, x1

j = dec / de
if (dec < 0._REAL64) j = j - 1    ! x is always positive
x1 = dec / de - j
x  = 1._REAL64 - x1
i1 = max( 1 , 1 - j )
in = min( ne , ne - j + 1 )
if ( iabem == -1 ) then          ! to unoccupied states only
  i1 = max( i1 , nf - j )
else if ( iabem == 1 ) then      ! to occupied states
  in = min( in , nf - j )
endif
do i = i1 , in - 1
  output(i) = output(i) + input(i + j) * x + input(i + j + 1) * x1
enddo
output(in) = output(in) + input(in + j) * x

END SUBROUTINE X_INTP

SUBROUTINE TENSOR_FROM_MME(mme, tensor)
use, intrinsic :: iso_fortran_env
implicit none
complex(REAL64) :: mme(3)      ! ,intent(in)
complex(REAL64) :: tensor(3, 3)! ,intent(out)
! local vars
complex(REAL64) mxt, myt, mzt

mxt = mme(1)
myt = mme(2)
mzt = mme(3)
tensor(1,1) = mxt * conjg(mxt)
tensor(1,2) = mxt * conjg(myt)
tensor(1,3) = mxt * conjg(mzt)
tensor(2,1) = conjg(tensor(1,2))
tensor(2,2) = myt * conjg(myt)
tensor(2,3) = myt * conjg(mzt)
tensor(3,1) = conjg(tensor(1,3))
tensor(3,2) = conjg(tensor(2,3))
tensor(3,3) = mzt * conjg(mzt)

END SUBROUTINE TENSOR_FROM_MME

CHARACTER(BUFFER_SIZE) FUNCTION to_lowercase(input)
use m_constants, only: BUFFER_SIZE
implicit none
character(*) input
character(26) :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
character(26) :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz'
integer i, pos

to_lowercase = trim(input)
do i = 1, len_trim(input)
   pos = index(UPPER_CASE, input(i:i))
   if (pos > 0) to_lowercase(i:i) = LOWER_CASE(pos:pos)
enddo

END FUNCTION to_lowercase

SUBROUTINE set_aux_arrays
use m_params
implicit none

if (read_rixfile) return
if (check_kstar ) return
if (irixs == 0 .and. iabsorp == 0) return

call set_lmtindex
! call set_map_array <-- now it is called from read_RIXS_data
call set_i1_i2_i3
if (irixs == 0) return
call set_ntr
call set_ispec_ibi_ibf

END SUBROUTINE set_aux_arrays

SUBROUTINE set_lmtindex
use m_aux, only: lmtindex
use m_functions, only: allocate
use m_rixs, only: nlmtdata, lmtdata, n_mme
implicit none
character(10), parameter :: srcname = ' in SETAUX'
integer jsp

call allocate(lmtindex, 'lmtindex' // srcname, udim1 = n_mme, nullify = .true.)
do jsp = 1, nlmtdata
  lmtindex( lmtdata(jsp)%ia1 : lmtdata(jsp)%ian ) = jsp
enddo

END SUBROUTINE set_lmtindex

SUBROUTINE set_map_array
use, intrinsic :: iso_fortran_env
use m_aux, only: map
use m_functions, only: allocate, exit_on_error
use m_rixs, only: nrixs, rixs, n_mme, lmtdata, nlmtdata
implicit none
character(*), parameter :: srcname = ' in SET_MAP'
integer jsp, jk
logical lstop

call allocate(map, 'map' // srcname, udim1 = n_mme, udim2 = nrixs, &
     nullify = .true.)
do jsp = 1, nrixs
  do jk = 1, nlmtdata
    if (lmtdata(jk)%txtel(:2) == rixs(jsp)%txtel(:2)) then
      if (lmtdata(jk)%sname(:) == rixs(jsp)%sname(:)) then
        map( lmtdata(jk)%ia1 : lmtdata(jk)%ian, jsp ) = 1
      endif
    endif
  enddo
enddo
lstop = .false.
do jsp = 1, nrixs
  if (sum(map( :, jsp)) == 0) then
    write(unit = output_unit, fmt = 5) jsp
    lstop = .true.
  endif
enddo
if (lstop) then
  call print_map_array
  call exit_on_error
endif

 5 format(4x, 'ERROR in SET_MAP: no core levels have been', &
         ' found for ispec=', i0)

END SUBROUTINE set_map_array

SUBROUTINE print_map_array
use, intrinsic :: iso_fortran_env
use m_rixs, only: nrixs, n_mme
use m_aux, only: map
implicit none
integer ja, jsp

write(unit = output_unit, fmt = 10) 'MAP ARRAY (ia, ispec_RIXS)'
write(unit = output_unit, fmt = 12) (jsp, jsp = 1, nrixs)
do ja = 1, n_mme
  write(unit = output_unit, fmt = 14) ja, map(ja, : )
enddo
 10 format(4x, a)
 12 format(9x, 100i4)
 14 format(4x, i4, ':', 100(i4))

END SUBROUTINE print_map_array

SUBROUTINE set_i1_i2_i3
use m_aux, only: i1, i2, i3
use m_bzmesh, only: nkbz, ndxyz
use m_functions, only: allocate
implicit none
character(*), parameter :: srcname = ' in SET_I1_I2_I3'
integer j1, j2, j3, jk

call allocate(i1, 'i1' // srcname, udim1 = nkbz)
call allocate(i2, 'i2' // srcname, udim1 = nkbz)
call allocate(i3, 'i3' // srcname, udim1 = nkbz)
jk = 1
do j3 = 1, ndxyz(3)
  do j2 = 1, ndxyz(2)
    do j1 = 1, ndxyz(1)
      i1(jk) = j1
      i2(jk) = j2
      i3(jk) = j3
      jk     = jk + 1
    enddo
  enddo
enddo

END SUBROUTINE set_i1_i2_i3

SUBROUTINE set_ntr
use, intrinsic :: iso_fortran_env
use m_aux
use m_functions, only: exit_on_error
use m_rixs, only: nrixs, rixs, n_mme
implicit none
character(*), parameter :: srcname = ' in SET_NTR'
integer ja_counter, jbi, jbf, jsp, nbi1, nbi2, nbf1, nbf2

! call print_map_array
ntr = 0
do jsp = 1, nrixs
  ja_counter = 0
  nbi1 = rixs(jsp)%nbi(1)
  nbi2 = rixs(jsp)%nbi(2)
  nbf1 = rixs(jsp)%nbf(1)
  nbf2 = rixs(jsp)%nbf(2)
!$$$  do ja = 1, n_mme
!$$$    if (map(ja, jsp) == 0) cycle
!$$$    ja_counter = ja_counter + 1
    do jbi = nbi1, nbi2
      do jbf = nbf1, nbf2
        if (ntr == huge(0)) then
          write(unit = output_unit, fmt = 10) 'Number of transactions ' // &
               'is out of range for INT32'
          call exit_on_error(srcname)
        endif
        ntr = ntr + 1
        rixs(jsp)%ntr = rixs(jsp)%ntr + 1
      enddo
    enddo
!$$$  enddo
!$$$  write(unit = output_unit, fmt = 20) 'For ISPEC = ', jsp, ' ntr = ', &
!$$$       rixs(jsp)%ntr, ' = ', ja_counter, ' * ', (nbi2 - nbi1 + 1) ,' * ', &
!$$$       (nbf2 - nbf1 + 1) 
enddo
if (ntr == 0) then
  write(unit = output_unit, fmt = 25) ntr, nrixs, n_mme
  call print_map_array
  call exit_on_error
endif
if (ntr /= sum(rixs(:)%ntr)) then
  write(unit = output_unit, fmt = 30) ntr, sum(rixs(:)%ntr)
  call exit_on_error
endif

 10 format(4x, a)
!$$$ 20 format(4x, 5(a, i0))
 25 format(4x, 'ERROR in SET_MAP_ARRAY: ntr=', i0, 4x, 'nrixs=', i0, &
           4x, 'n_mme=', i0)
 30 format(4x, 'ERROR in SET_MAP_ARRAY: ntr /= sum(rixs(:)%ntr): ', i0, &
           ' /= ', i0)

END SUBROUTINE set_ntr

SUBROUTINE set_ispec_ibi_ibf
use m_aux
use m_functions, only: allocate
use m_rixs, only: nrixs, rixs! , n_mme
implicit none
character(*), parameter :: srcname = ' in SET_ISPEC_IBI_IBF'
integer itr, jsp, jbi, jbf

!$$$ call allocate(ia   , 'ia'    // srcname, udim1 = ntr)
call allocate(ispec, 'ispec' // srcname, udim1 = ntr)
call allocate(ibi  , 'ibi'   // srcname, udim1 = ntr)
call allocate(ibf  , 'ibf'   // srcname, udim1 = ntr)

itr = 1
do jsp = 1, nrixs
!$$$  do ja = 1, n_mme
!$$$    if (map(ja, jsp) == 0) cycle
    do jbi = rixs(jsp)%nbi(1), rixs(jsp)%nbi(2)
      do jbf = rixs(jsp)%nbf(1), rixs(jsp)%nbf(2)
!$$$        ia(itr)=ja
        ispec(itr) = jsp
        ibi(itr)   = jbi
        ibf(itr)   = jbf
        itr        = itr + 1
      enddo
    enddo
!$$$  enddo
enddo

END SUBROUTINE set_ispec_ibi_ibf

SUBROUTINE PRINT_ALLOCATION_ERROR(name, code, msg)
use, intrinsic :: iso_fortran_env
implicit none
character(*) name, msg
integer code

write(unit = output_unit, fmt = 5 )
write(unit = output_unit, fmt = 10) trim(name)
write(unit = output_unit, fmt = 15) code
write(unit = output_unit, fmt = 20) trim(msg )

  5 format(/, 4x, 'An ALLOCATION ERROR OCCURIED')
 10 format(   4x, 'Array: ',         a )
 15 format(   4x, 'Error code: ',    i0)
 20 format(   4x, 'Error message: ', a )

END SUBROUTINE PRINT_ALLOCATION_ERROR

SUBROUTINE PRINT_DEALLOCATION_ERROR(name, code, msg)
use , intrinsic :: iso_fortran_env
implicit none
character(*) name , msg
integer code

write(unit = output_unit, fmt = 5 )
write(unit = output_unit, fmt = 10) trim(name)
write(unit = output_unit, fmt = 15) code
write(unit = output_unit, fmt = 20) trim(msg )

  5 format(/ , 4x, 'THE DEALLOCATION ERROR OCCURIED')
 10 format(    4x, 'Array: ',         a )
 15 format(    4x, 'Error code: ',    i0)
 20 format(    4x, 'Error message: ', a )

END SUBROUTINE PRINT_DEALLOCATION_ERROR
