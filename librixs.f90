SUBROUTINE X_INTP(iabem,dec,nf,de,ne,input,output)
! interpolates spectrum to a common omega mesh. For absorption
! (iabem == -1) and emission (iabem == 1) only transitions
! to unoccupied or occupied states are considered.
! de can be dexmcd .or. derexs. or. derixs
! ne can be nexmcd .or. nerexs. or. nerixs
use, intrinsic :: iso_fortran_env
implicit none
! input and output
integer iabem , nf , ne
real(REAL64) de , dec , input(ne) , output(ne)
! local vars
integer i , i1 , in , j
real(REAL64) x , x1

j = dec / de
if ( dec < 0._REAL64 ) j = j - 1    ! x is always positive
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
  output(i) = output(i) + input( i + j ) * x + input( i + j + 1 ) * x1
enddo
output(in) = output(in) + input( in + j ) * x

END SUBROUTINE X_INTP

SUBROUTINE TENSOR_FROM_MME(mme,tensor)
use, intrinsic :: iso_fortran_env
implicit none
complex(REAL64) :: mme(3)      ! ,intent(in)
complex(REAL64) :: tensor(3,3) ! ,intent(out)
! local vars
complex(REAL64) mxt,myt,mzt

mxt=mme(1)
myt=mme(2)
mzt=mme(3)
tensor(1,1)=mxt*conjg(mxt)
tensor(1,2)=mxt*conjg(myt)
tensor(1,3)=mxt*conjg(mzt)
tensor(2,1)=conjg(tensor(1,2))
tensor(2,2)=myt*conjg(myt)
tensor(2,3)=myt*conjg(mzt)
tensor(3,1)=conjg(tensor(1,3))
tensor(3,2)=conjg(tensor(2,3))
tensor(3,3)=mzt*conjg(mzt)

END SUBROUTINE TENSOR_FROM_MME

SUBROUTINE SETAUX
use m_aux
use m_functions
use m_params !, only: iabsorp, irixs, iprint
use m_bzmesh, only: nkbz, ndxyz
use m_rixs, only: nrixs, rixs, n_mme, lmtdata, nlmtdata! ,ispc
implicit none
! local vars
character(*), parameter :: srcname = ' in SETAUX'
character(BUFFER_SIZE) msg
integer iok, itr, jsp, j1, j2, j3, ja, jbi, jbf, jk
logical lstop

if ( read_rixfile ) return
if ( check_kstar .and. debug_mode ) return

if ( irixs == 0 .and. iabsorp == 0 ) return

allocate(lmtindex(n_mme),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('lmtindex',iok,msg)
lmtindex(:) = 0
do jsp = 1, nlmtdata
  lmtindex( lmtdata(jsp)%ia1 : lmtdata(jsp)%ian ) = jsp
enddo

!$$$ allocate(map(n_mme,nrixs),stat=iok,errmsg=msg)
!$$$ if(iok /= 0)call print_allocation_error('map',iok,msg)
!$$$ map( : , : ) = 0

call allocate( map, 'map' // srcname, udim1 = n_mme, udim2 = nrixs, &
     nullify = .true. )

do jsp = 1, nrixs
  do jk = 1, nlmtdata
    if( lmtdata(jk)%txtel(:2) == rixs(jsp)%txtel(:2) ) then
      if( lmtdata(jk)%sname(:) == rixs(jsp)%sname(:) ) then
        map( lmtdata(jk)%ia1 : lmtdata(jk)%ian, jsp ) = 1
      endif
    endif
  enddo
enddo

lstop = .false.
do jsp = 1, nrixs
  if ( sum( map( :, jsp ) ) == 0 ) then
    write( unit = output_unit, fmt = 5 ) jsp
    lstop = .true.
  endif
enddo
if ( lstop ) then
  write( unit = output_unit, fmt = 10 )
  do ja = 1, n_mme
    write( unit = output_unit, fmt = 11) map( ja, : )
  enddo
  call exit_on_error
endif

  5 format( 4x, 'ERROR in SETAUX: no core levels have been' &
          , ' found for ispec=', i0 )
 10 format( 4x, 'MAP ARRAY ( ia, isp )' )
 11 format( 4x, 100(i4) )

if ( debug_mode ) then
  write( unit = output_unit, fmt = 10 )
  do ja = 1, n_mme
    write( unit = output_unit, fmt = 11 ) map( ja, : )
  enddo
endif

allocate(i1(nkbz),i2(nkbz),i3(nkbz),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('i1,i2,i3',iok,msg)

jk = 0
do j3 = 1, ndxyz(3)
  do j2 = 1, ndxyz(2)
    do j1 = 1, ndxyz(1)
      jk = jk + 1
      i1(jk) = j1
      i2(jk) = j2
      i3(jk) = j3
    enddo                        ! j1
  enddo                          ! j2
enddo                            ! j3
if ( jk /= nkbz ) then
  write( unit = output_unit, fmt = 15 )
  call deallocate_global_arrays
  stop
endif

 15 format( 4x, 'In SETAUX: wrong nkbz value.' )

if ( irixs == 0 ) return

ntr = 0
do jsp = 1, nrixs
!$$$  do ja=1,n_mme
!$$$    if(map(ja,jsp) == 0)cycle
    do jbi = rixs(jsp)%nbi(1), rixs(jsp)%nbi(2)
      do jbf = rixs(jsp)%nbf(1), rixs(jsp)%nbf(2)
        if ( ntr == huge(0) ) then
          write( unit = output_unit, fmt = 20 )
          call deallocate_global_arrays
          stop
        endif
        ntr = ntr + 1
        rixs(jsp)%ntr = rixs(jsp)%ntr + 1
      enddo
    enddo
!$$$  enddo
enddo
if ( ntr == 0 ) then
  write( unit = output_unit, fmt = 25 ) ntr, nrixs, n_mme
  write( unit = output_unit, fmt = 10 )
  do ja = 1, n_mme
    write( unit = output_unit, fmt = 11 ) map( ja, : )
  enddo
  call deallocate_global_arrays
  stop
endif
if ( ntr /= sum( rixs(:)%ntr ) ) then
  write ( unit = output_unit, fmt = 30 ) ntr, sum( rixs(:)%ntr )
  call deallocate_global_arrays
  stop
endif

 20 format( 4x, 'Number of transactions is out of range for INT32.' )
 25 format( 4x, 'ERROR in SETAUX: ntr=', i0, 4x, 'nrixs=', i0, 4x &
          , 'n_mme=', i0 )
 30 format( 4x, 'ERROR in SETAUX: ntr /= sum(rixs(:)%ntr: ', i0, ' /= ' &
          , i0 )

!$$$ allocate(ispec(ntr),ia(ntr),ibi(ntr),ibf(ntr),stat=iok,errmsg=msg)
!$$$ if(iok /= 0)call print_allocation_error('ispec,ia,ibi,ibf',iok,msg)

allocate(ispec(ntr),ibi(ntr),ibf(ntr),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('ispec,ibi,ibf',iok,msg)

itr = 0
do jsp = 1, nrixs
!$$$  do ja=1,n_mme
!$$$    if(map(ja,jsp) == 0)cycle
    do jbi = rixs(jsp)%nbi(1), rixs(jsp)%nbi(2)
      do jbf = rixs(jsp)%nbf(1), rixs(jsp)%nbf(2)
        itr = itr + 1
        ispec(itr) = jsp
!$$$         ia(itr)=ja
        ibi(itr) = jbi
        ibf(itr) = jbf
      enddo
    enddo
!$$$  enddo
enddo

END SUBROUTINE SETAUX

SUBROUTINE PRINT_ALLOCATION_ERROR(name,code,msg)
use, intrinsic :: iso_fortran_env
implicit none
character(*) name, msg
integer code

write( unit = output_unit, fmt = 5 )
write( unit = output_unit, fmt = 10 ) trim( name )
write( unit = output_unit, fmt = 15 ) code
write( unit = output_unit, fmt = 20 ) trim( msg )

  5 format( /, 4x, 'An ALLOCATION ERROR OCCURIED' )
 10 format( 4x, 'Array: ', a )
 15 format( 4x, 'Error code: ', i0 )
 20 format( 4x, 'Error message: ', a )

 !  5 format(4x,'An ALLOCATION ERROR OCCURIED')
 ! 10 format(12x,'Array: ',a)
 ! 15 format(7x,'Error code: ',i0)
 ! 20 format(4x,'Error message: ',a)

END SUBROUTINE PRINT_ALLOCATION_ERROR

SUBROUTINE PRINT_DEALLOCATION_ERROR(name,code,msg)
use , intrinsic :: iso_fortran_env
implicit none
character(*) name , msg
integer code

write( unit = output_unit , fmt = 5 )
write( unit = output_unit , fmt = 10 ) trim( name )
write( unit = output_unit , fmt = 15 ) code
write( unit = output_unit , fmt = 20 ) trim( msg )

  5 format( / , 4x , 'THE DEALLOCATION ERROR OCCURIED' )
 10 format( 4x , 'Array: ' , a )
 15 format( 4x , 'Error code: ' , i0 )
 20 format( 4x , 'Error message: ' , a )

END SUBROUTINE PRINT_DEALLOCATION_ERROR
