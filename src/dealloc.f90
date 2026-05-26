SUBROUTINE EXIT_ON_ERROR(from)
use, intrinsic :: ISO_FORTRAN_ENV
! use m_functions, only: deallocate_global_arrays
character(*), optional :: from

if (present(from)) then
  write(unit = output_unit, fmt = "(/ , 4x , a)") &
       'EXIT_ON_ERROR has been called' // trim(from)
endif
call deallocate_global_arrays
call close_all_files
stop

END SUBROUTINE EXIT_ON_ERROR

SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS

call deallocate_lmtdata
call deallocate_rixs
call deallocate_everything_else

END SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS

SUBROUTINE DEALLOCATE_LMTDATA
use m_constants, only: BUFFER_SIZE
use m_functions, only: deallocate, int2string
use m_rixs
implicit none
character(BUFFER_SIZE) str
character(28), parameter :: srcname=" in DEALLOCATE_GLOBAL_ARRAYS"
integer isp
 
if (allocated(lmtdata)) then
  do isp = 1, nlmtdata
     if (allocated(lmtdata(isp)%mme)) then
      str(:) = 'lmtdata(' // trim(int2string(isp)) // ')%mme' // srcname
      call deallocate(lmtdata(isp)%mme, str)
    endif
    if (allocated(lmtdata(isp)%mmebz)) then
      str(:) = 'lmtdata(' // trim(int2string(isp)) // ')%mmebz' // srcname
      call deallocate(lmtdata(isp)%mmebz, str)
    endif
    if (allocated(lmtdata(isp)%esp)) then
      str(:) = 'lmtdata(' // trim(int2string(isp)) // ')%esp' // srcname
      call deallocate(lmtdata(isp)%esp, str)
    endif
    if (allocated(lmtdata(isp)%rat)) then
      str(:) = 'lmtdata(' // trim(int2string(isp)) // ')%rat' // srcname
      call deallocate(lmtdata(isp)%rat, str)
    endif
  enddo
  call deallocate(lmtdata, 'lmtdata' // srcname)
endif

END SUBROUTINE DEALLOCATE_LMTDATA

SUBROUTINE DEALLOCATE_RIXS
use m_constants, only: BUFFER_SIZE
use m_functions, only: deallocate, int2string
use m_rixs
implicit none
character(BUFFER_SIZE) message
character(28), parameter :: srcname=" in DEALLOCATE_GLOBAL_ARRAYS"
integer isp

if (allocated(rixs)) then
  do isp = 1, nrixs
    ! if (allocated(rixs(isp)%g)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%g' // srcname
    !   call deallocate(rixs(isp)%g, str)
    ! endif
    ! if (allocated(rixs(isp)%itetr)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%itetr' // srcname
    !   call deallocate(rixs(isp)%itetr, str)
    ! endif
    ! if (allocated(rixs(isp)%idold)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%idold' // srcname
    !   call deallocate(rixs(isp)%idold, str)
    ! endif
    ! if (allocated(rixs(isp)%ipq)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%ipq' // srcname
    !   call deallocate(rixs(isp)%ipq, str)
    ! endif
    ! if (allocated(rixs(isp)%ik2rbz)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%ik2rbz' // srcname
    !   call deallocate(rixs(isp)%ik2rbz, str)
    ! endif
    ! if (allocated(rixs(isp)%ikbz2rbz)) then
    !   str(:) = 'rixs(' // trim(int2string(isp)) // ')%ikbz2rbz' // srcname
    !   call deallocate(rixs(isp)%ikbz2rbz, str)
    ! endif
    if (allocated(rixs(isp)%ik2kq)) then
      message = 'rixs(' // trim(int2string(isp)) // ')%ik2kq' // srcname
      call deallocate(rixs(isp)%ik2kq, message)
    endif
  enddo
  call deallocate(rixs, 'rixs' // srcname)
endif

END SUBROUTINE DEALLOCATE_RIXS

SUBROUTINE DEALLOCATE_EVERYTHING_ELSE
use m_aux
use m_bnd
use m_bz
use m_bzmesh
use m_functions, only: deallocate
use m_results
use m_sdt
implicit none
character(*), parameter :: srcname=" in DEALLOCATE_GLOBAL_ARRAYS"

if (allocated(vg)        ) call deallocate(vg        , 'vg'         // srcname)
if (allocated(igtta)     ) call deallocate(igtta     , 'igtta'      // srcname)
if (allocated(e)         ) call deallocate(e         , 'e'          // srcname)
if (allocated(g)         ) call deallocate(g         , 'g'          // srcname)
if (allocated(ipqibz)    ) call deallocate(ipqibz    , 'ipqibz'     // srcname)
if (allocated(ig4qibz)   ) call deallocate(ig4qibz   , 'ig4qibz'    // srcname)
if (allocated(ik2ibz)    ) call deallocate(ik2ibz    , 'ik2ibz'     // srcname)
if (allocated(ispec)     ) call deallocate(ispec     , 'ispec'      // srcname)
if (allocated(ibi)       ) call deallocate(ibi       , 'ibi'        // srcname)
if (allocated(ibf)       ) call deallocate(ibf       , 'ibf'        // srcname)
if (allocated(i1)        ) call deallocate(i1        , 'i1'         // srcname)
if (allocated(i2)        ) call deallocate(i2        , 'i2'         // srcname)
if (allocated(i3)        ) call deallocate(i3        , 'i3'         // srcname)
if (allocated(lmtindex)  ) call deallocate(lmtindex  , 'lmtindex'   // srcname)
if (allocated(map)       ) call deallocate(map       , 'map'        // srcname)
if (allocated(loss)      ) call deallocate(loss      , 'loss'       // srcname)
if (allocated(loss_neg_w)) call deallocate(loss_neg_w, 'loss_neg_w' // srcname)
if (allocated(absorp)    ) call deallocate(absorp    , 'absorp'     // srcname)
if (allocated(ipq)       ) call deallocate(ipq       , 'ipq'        // srcname)
if (allocated(pnt)       ) call deallocate(pnt       , 'pnt'        // srcname)
if (allocated(wgt)       ) call deallocate(wgt       , 'wgt'        // srcname)
if (allocated(itetr)     ) call deallocate(itetr     , 'itetr'      // srcname)
if (allocated(idold)     ) call deallocate(idold     , 'idold'      // srcname)

END SUBROUTINE DEALLOCATE_EVERYTHING_ELSE

SUBROUTINE CLOSE_ALL_FILES
use m_files
use m_functions, only: close
implicit none

if (bnd /= 0)     call close(unit = bnd)
if (bns /= 0)     call close(unit = bns)
if (inp /= 0)     call close(unit = inp)
if (rid /= 0)     call close(unit = rid)
if (sdt /= 0)     call close(unit = sdt)
if (rix /= 0)     call close(unit = rix)
if (rim /= 0)     call close(unit = rim)
if (m2d_q /= 0)   call close(unit = m2d_q)
if (m2d_ein /= 0) call close(unit = m2d_ein)

END SUBROUTINE CLOSE_ALL_FILES
