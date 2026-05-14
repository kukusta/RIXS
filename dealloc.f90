SUBROUTINE EXIT_ON_ERROR( from , stop )
  use, intrinsic :: ISO_FORTRAN_ENV
  use m_functions , only : deallocate_global_arrays
  character(*), optional :: from
  logical, optional :: stop
! local vars
  logical stop_local
  if( present(from) ) &
    write( unit = output_unit , fmt = "( / , 4x , a )") &
      'EXIT_ON_ERROR has been called'//trim(from)
  if ( .not. present( stop ) ) then
    stop_local = .true.
  else
    stop_local = stop
  endif
  call deallocate_global_arrays
  if ( stop_local ) stop
END SUBROUTINE EXIT_ON_ERROR

SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS( success )
use m_bzmesh
use m_aux
use m_bnd
use m_bz
use m_files
use m_functions , only : deallocate , close , int2string
use m_results
use m_rixs
use m_sdt
implicit none
logical , optional :: success
! local vars
character(28), parameter :: srcname=" in DEALLOCATE_GLOBAL_ARRAYS"
integer isp

! these arrays may be deallocated in READRIXSDATA
if ( allocated(lmttemp) ) then
  do isp = 1, nrixslmt
    if ( allocated( lmttemp(isp)%esp ) ) &
      call deallocate( lmttemp(isp)%esp, 'lmt(' // trim( int2string(isp) ) &
           // ')%esp' // srcname )
    if ( allocated( lmttemp(isp)%rat ) ) &
      call deallocate( lmttemp(isp)%rat, 'lmt(' // trim( int2string(isp) ) &
           // ')%rat' // srcname )
    if ( allocated( lmttemp(isp)%mme ) ) &
      call deallocate( lmttemp(isp)%mme, 'lmt(' // trim( int2string(isp) ) &
           // ')%mme' // srcname )
    if ( allocated( lmttemp(isp)%mmebz ) ) &
      call deallocate( lmttemp(isp)%mmebz, 'lmt(' // trim( int2string(isp) ) &
           // ')%mmebz' // srcname )
    ! deallocate(lmt(isp)%esp,stat=iok)
    ! deallocate(lmt(isp)%rat,stat=iok)
    ! deallocate(lmt(isp)%mme,stat=iok)
  enddo
  ! deallocate(lmt,stat=iok)
  call deallocate( lmttemp, 'lmt' // srcname )
endif

if ( allocated(lmtdata) ) then
  do isp = 1, nlmtdata
    if ( allocated( lmtdata(isp)%mme ) ) &
      call deallocate( lmtdata(isp)%mme, 'lmtdata(' &
           // trim( int2string(isp) ) // ')%mme' // srcname )
    if ( allocated( lmtdata(isp)%mmebz ) ) &
      call deallocate( lmtdata(isp)%mmebz, 'lmtdata(' &
           // trim( int2string(isp) ) // ')%mmebz' // srcname )
    if ( allocated( lmtdata(isp)%esp ) ) &
      call deallocate( lmtdata(isp)%esp, 'lmtdata(' &
           // trim( int2string(isp) ) // ')%esp' // srcname )
    if ( allocated( lmtdata(isp)%rat ) ) &
      call deallocate( lmtdata(isp)%rat, 'lmtdata(' &
           // trim( int2string(isp) ) // ')%rat' // srcname )
    ! deallocate(lmtdata(isp)%mme,stat=iok)
    ! deallocate(lmtdata(isp)%esp,stat=iok)
    ! deallocate(lmtdata(isp)%rat,stat=iok)
  enddo
  call deallocate( lmtdata , 'lmtdata' // srcname )
  ! deallocate(lmtdata,stat=iok)
endif

if ( allocated(rixs) ) then
  do isp = 1 , nrixs
    if ( allocated( rixs(isp)%g ) ) &
      call deallocate( rixs(isp)%g , 'rixs(' // trim( int2string(isp) ) &
           // ')%g' // srcname )
    if ( allocated( rixs(isp)%itetr ) ) &
      call deallocate( rixs(isp)%itetr , 'rixs(' // trim( int2string(isp) ) &
           // ')%itetr' // srcname )
    if ( allocated( rixs(isp)%idold ) ) &
      call deallocate( rixs(isp)%idold , 'rixs(' // trim( int2string(isp) ) &
           // ')%idold' // srcname )
    if ( allocated( rixs(isp)%ipq ) ) &
      call deallocate( rixs(isp)%ipq , 'rixs(' // trim( int2string(isp) ) &
           // ')%ipq' // srcname )
    if ( allocated( rixs(isp)%ik2rbz ) ) &
      call deallocate( rixs(isp)%ik2rbz , 'rixs(' // trim( int2string(isp) ) &
           // ')%ik2rbz' // srcname )
    if ( allocated( rixs(isp)%ikbz2rbz ) ) &
      call deallocate( rixs(isp)%ikbz2rbz , 'rixs(' // trim( int2string(isp) )&
           // ')%ikbz2rbz' // srcname )
    if ( allocated( rixs(isp)%ik2kq ) ) &
      call deallocate( rixs(isp)%ik2kq , 'rixs(' // trim( int2string(isp) ) &
           // ')%ik2kq' // srcname )
    ! deallocate(rixs(isp)%g,stat=iok)
    ! deallocate(rixs(isp)%itetr,stat=iok)
    ! deallocate(rixs(isp)%idold,stat=iok)
    ! deallocate(rixs(isp)%ipq,stat=iok)
    ! deallocate(rixs(isp)%ik2rbz,stat=iok)
    ! deallocate(rixs(isp)%ik2kq,stat=iok)
  enddo
  call deallocate( rixs , 'rixs' // srcname )
  ! deallocate(rixs,stat=iok)
endif

! deallocate(vg,stat=iok)
! deallocate(igtta,stat=iok)
! deallocate(e,stat=iok)
! deallocate(g,stat=iok)
! deallocate(ipqibz,stat=iok)
! deallocate(ig4qibz,stat=iok)
! deallocate(ik2ibz,stat=iok)
! deallocate(ispec,stat=iok)
! deallocate(ia,stat=iok)
! deallocate(ibi,stat=iok)
! deallocate(ibf,stat=iok)
! deallocate(i1,stat=iok)
! deallocate(i2,stat=iok)
! deallocate(i3,stat=iok)
! deallocate(lmtindex,stat=iok)
! deallocate(map,stat=iok)
! deallocate(loss,stat=iok)
! deallocate(absorp,stat=iok)
! deallocate(ipq,stat=iok)
! deallocate(pnt,stat=iok)
! deallocate(wgt,stat=iok)
! deallocate(itetr,stat=iok)
! deallocate(idold,stat=iok)

if( allocated(vg)       ) call deallocate( vg       , 'vg'       // srcname )
if( allocated(igtta)    ) call deallocate( igtta    , 'igtta'    // srcname )
if( allocated(e)        ) call deallocate( e        , 'e'        // srcname )
if( allocated(g)        ) call deallocate( g        , 'g'        // srcname )
if( allocated(ipqibz)   ) call deallocate( ipqibz   , 'ipqibz'   // srcname )
if( allocated(ig4qibz)  ) call deallocate( ig4qibz  , 'ig4qibz'  // srcname )
if( allocated(ik2ibz)   ) call deallocate( ik2ibz   , 'ik2ibz'   // srcname )
if( allocated(ispec)    ) call deallocate( ispec    , 'ispec'    // srcname )
!$$$ if( allocated(ia)       ) call deallocate( ia       , 'ia'       // srcname )
if( allocated(ibi)      ) call deallocate( ibi      , 'ibi'      // srcname )
if( allocated(ibf)      ) call deallocate( ibf      , 'ibf'      // srcname )
if( allocated(i1)       ) call deallocate( i1       , 'i1'       // srcname )
if( allocated(i2)       ) call deallocate( i2       , 'i2'       // srcname )
if( allocated(i3)       ) call deallocate( i3       , 'i3'       // srcname )
if( allocated(lmtindex) ) call deallocate( lmtindex , 'lmtindex' // srcname )
if( allocated(map)      ) call deallocate( map      , 'map'      // srcname )
if( allocated(loss)     ) call deallocate( loss     , 'loss'     // srcname )
if( allocated(absorp)   ) call deallocate( absorp   , 'absorp'   // srcname )
if( allocated(ipq)      ) call deallocate( ipq      , 'ipq'      // srcname )
if( allocated(pnt)      ) call deallocate( pnt      , 'pnt'      // srcname )
if( allocated(wgt)      ) call deallocate( wgt      , 'wgt'      // srcname )
if( allocated(itetr)    ) call deallocate( itetr    , 'itetr'    // srcname )
if( allocated(idold)    ) call deallocate( idold    , 'idold'    // srcname )

!$$$ write( unit = output_unit, fmt = "()" )
if ( bnd /= 0 ) then
  call close( unit = bnd )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'BNDFILE has not been opened'
endif
if ( bns /= 0 ) call close( unit = bns )
if ( inp /= 0 ) then
  call close( unit = inp )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'INPUTFILE has not been opened'
endif
if ( rid /= 0 ) then
  call close( unit = rid )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'RIXSDATAFILE has not been opened'
endif
if ( sdt /= 0 ) then
  call close( unit = sdt )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'SDTFILE has not been opened'
endif
if ( rix /= 0 ) then
  call close( unit = rix )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'RIXSFILE has not been opened'
endif
if ( rim /= 0 ) then
  call close( unit = rim )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'RIXSMMEFILE has not been opened'
endif
if ( m2d_q /= 0 ) then
  call close( unit = m2d_q )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'MAP2DQFILE has not been opened'
endif
if ( m2d_ein /= 0 ) then
  call close( unit = m2d_ein )
else
!$$$   write( unit = output_unit, fmt = "( 4x, a )" ) &
!$$$     'MAP2DEINFILE has not been opened'
endif

if ( run /= 0 ) then
  if ( present( success ) ) then
    if ( success )then
      call close( unit = run , status = 'delete' )
    else
      call close( unit = run )
    endif
  else
    call close( unit = run )
  endif
endif

END SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS
