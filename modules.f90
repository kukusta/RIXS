! MODULE M_RIXS2
!   use, intrinsic :: iso_fortran_env
!   implicit none

!   type, public :: lmto_spectrum
!     character(4) sname
!     character(6) txtel
!     integer nn,nk,iwritemme
!     integer isort
!     integer ia1,ian
!     integer iat1,iatn
!     integer label
!     real(REAL64) ec,enlmin,enlmax,dnhsort
!     complex(REAL64), allocatable :: mme(:,:,:,:,:) !3,iat1:iatn,ia1:ian,nb,npnt
!     ! mmebz(1:3,iat1:iatn,ia1:ian,nb,1:nkbz)
!     complex(REAL64), allocatable :: mmebz(:,:,:,:,:)
!     real(REAL64), allocatable :: esp(:)
!     real(REAL64), allocatable :: rat(:,:)
!   end type lmto_spectrum

!   type, public :: rixs_spectrum
!     complex(REAL64) e_in(3),e_out(3)
!     character(4) sname,channel,use_symmetry
!     character(6) txtel
!     integer isort,nn,nk,iwritemme,label
!     integer nbi(2),nbf(2),ntr
!     integer ng,ntrbz,nkrbz,iopnum(64)
!     logical save2dat,map2d_ein,map2d_q
!     real(REAL64) ec,enlmin,enlmax,dnhsort,Gamma,q(3)
!     real(REAL64) en,fnorm,k_in(3),k_out(3)
! ! symmetry operations for this spectrum
!     real(REAL64), allocatable :: g(:,:,:)  ! (1:3,1:3,nopused)
! ! list of original k-points forming tetrahedra. For integration in optics.
!     integer, allocatable :: idold(:,:)     ! (4,rixs(isp)%ntrbz)
! ! ipq for rixs(isp)%ng
!     integer, allocatable :: ipq(:,:,:)     ! (ndxyz(1),ndxyz(2),ndxyz(3))
! ! list of irreducible k-points forming tetrahedra
!     integer, allocatable :: itetr(:,:)     ! (0:4,rixs(isp)%ntrbz)
! ! these arrays allow to use symmetry in bzopt
!     integer, allocatable :: ik2rbz(:)      ! 1:nkbz
!     integer, allocatable :: ik2kq(:)       ! 1:nkbz
!     integer, allocatable :: ikbz2rbz(:)    ! 1:rixs(isp)%nkrbz
!     integer, allocatable :: ig4q(:,:,:)    ! 1:nkbz,1:nkbz,1:nkbz
! ! auxiliary statistics
!     real(REAL64) max_mme , min_mme , av_mme
!   end type rixs_spectrum

!   integer nrixs,nrixslmt,ispc,nlabel,ne,n_mme,nfu,nlmtdata,max_nhsort
!   real(REAL64) de,emin,emax,k_in(3),Qaxis(3),Abas(3,3)
!   type(rixs_spectrum), allocatable :: rixs(:)    ! 1:nrixs
!   type(lmto_spectrum), allocatable :: lmttemp(:) ! 1:nrixslmt
!   type(lmto_spectrum), allocatable :: lmtdata(:) ! 1:nlmtdata
! END MODULE M_RIXS2

MODULE M_FILES
  use, intrinsic :: iso_fortran_env
  implicit none

  integer, parameter :: filename_length=256
  integer, parameter :: aliasname_length=16

  TYPE FILE
    character( filename_length ) fname
    character( aliasname_length ) alias
    integer unit
! pointers for bidirectional list
    type(file), pointer :: prev => null()
    type(file), pointer :: next => null()
  END TYPE FILE

  type(file), target, allocatable :: headfile
  type(file), pointer :: tail => null()

  character(80)  basename
  character(80)  bndfile, bnsfile, inptfile, datafile, sdtfile, &
                 rixsfile, runfile, map2deinfile, map2dqfile
  integer        bnd, bns, inp, rid, sdt, rix, run, m2d_ein, m2d_q
! mmefile can be written to scratch
  character(256) mmefile
  integer        rim

!   INTERFACE
!     SUBROUTINE ADD_FILE(fname,unit,alias)
! !      use m_files
!       implicit none
!       character(filename_length) :: fname
!       character(aliasname_length), optional :: alias
!       integer unit
! ! local vars
!       character(256) msg
!       integer iok
!       type(file), pointer :: this
!       type(file), allocatable :: newnode
!     END SUBROUTINE ADD_FILE
!   END INTERFACE

!   INTERFACE
!     INTEGER FUNCTION GET_UNIT(alias)
! !      use m_files
!       implicit none
!       character(aliasname_length) alias
!       type(file), pointer :: this
!     END FUNCTION GET_UNIT
!   END INTERFACE

!   INTERFACE
!     SUBROUTINE PRINT_FILES
! !      use m_files
!       implicit none
!       type(file), pointer :: this
!     END SUBROUTINE PRINT_FILES
!   END INTERFACE

!   INTERFACE
!     SUBROUTINE REMOVE_FILES
! !      use m_files
!     END SUBROUTINE REMOVE_FILES
!   END INTERFACE

END MODULE M_FILES

MODULE M_PARAMS
use, intrinsic :: iso_fortran_env
implicit none
integer, parameter :: MAX_NRIXS = 99
integer, parameter :: BUFFER_SIZE = 256
integer iprint, iabsorp, irixs, nbi(2), nbf(2), nw
logical newini, debug_mode, absorption_only, write_hsym_points, &
     jointbnd_only, check_kstar, read_rixfile
real(REAL64) wmin, wmax, dw, volomg
common /cor/ wmin, wmax, dw, nw
END MODULE M_PARAMS

MODULE M_SDT
use, intrinsic :: iso_fortran_env
implicit none
real(REAL64) ut( 3, 3 ), ut1( 3, 3 )
real(REAL64), allocatable :: vg( : , : ) ! 1:3,1:nopused
integer, allocatable :: igtta( : , : ) ! 1:nopused,1:natom
END MODULE M_SDT

MODULE M_BND
use, intrinsic :: iso_fortran_env
implicit none
character(4) high_sym_pnt_txt_dir(48)
integer natom, npnt, nb, nopused, iopnum(64)
real high_sym_pnt( 3, 48 )
real(REAL64) ef, rbas( 3, 3 ), qbas( 3, 3 )
real(REAL64), allocatable :: e( : , : )      ! 1:nb,1:nkibz
real(REAL64), allocatable :: ebz( : , : )    ! 1:nb,1:nkibz
real(REAL64), allocatable :: g( : , : , : )  ! 1:3,1:3,1:nopused
END MODULE M_BND

MODULE M_BZ
implicit none
integer, parameter :: ibzext = 0
integer, parameter :: is_avtet = 0
integer, parameter :: is_idold = 1
integer, parameter :: is_wgt = 0
! save data for IBZ
integer, allocatable :: ipqibz( : , : , : )  ! 1:n1,1:n2,1:n3
integer, allocatable :: ig4qibz( : , : , : )  ! 1:n1,1:n2,1:n3
integer, allocatable :: ik2ibz(:)      ! 1:nkbz
END MODULE M_BZ

MODULE M_AUX
! Auxiliary arrays for OpenMP calculations
use, intrinsic :: iso_fortran_env
implicit none
integer nthreads, ntr
integer, allocatable :: ispec(:)    ! 1:ntr
!$$$ integer, allocatable :: ia(:)       ! 1:ntr
integer, allocatable :: ibi(:)      ! 1:ntr
integer, allocatable :: ibf(:)      ! 1:ntr
integer, allocatable :: i1(:)       ! 1:nkbz
integer, allocatable :: i2(:)       ! 1:nkbz
integer, allocatable :: i3(:)       ! 1:nkbz

integer, allocatable :: lmtindex(:) ! 1:n_mme
integer, allocatable :: map( : , : )    ! 1:n_mme,1:nrixs
END MODULE M_AUX

MODULE M_RESULTS
use, intrinsic :: iso_fortran_env
implicit none
real(REAL64), allocatable :: loss( : , : )     ! 1:nw,1:nrixs
real(REAL64), allocatable :: absorp( : , : , : ) ! 1:nw,1:nrixs,1:iabsorp
END MODULE M_RESULTS

MODULE M_FUNCTIONS
  use m_params, only: BUFFER_SIZE
  implicit none

  INTERFACE ALLOCATE
    module procedure alloc1dint32, alloc1dreal32, alloc1dreal64, &
           alloc2dint32, alloc2dreal32, alloc2dreal64, alloc3dcomplex64, &
           alloc3dint32, alloc3dreal64, alloc4dcomplex64, alloc5dcomplex64
  END INTERFACE ALLOCATE

  INTERFACE DEALLOCATE
    module procedure dealloc1dint32, dealloc1dreal32, dealloc1dreal64, &
           dealloc1dlmtspec, dealloc1drixsspec, dealloc2dint32, &
           dealloc2dreal32, dealloc2dreal64, dealloc3dint32, &
           dealloc3dreal64, dealloc4dcomplex64, dealloc5dcomplex64
  END INTERFACE DEALLOCATE

  INTERFACE INT2STRING
    module procedure int32_2string,int64_2string
  END INTERFACE INT2STRING

  INTERFACE MME_TRANS
    module procedure mme_trans_real64,mme_trans_cmplx64
  END INTERFACE MME_TRANS

  INTERFACE MME_TRANS_NEW
    module procedure mme_trans_new_cmplx64
  END INTERFACE MME_TRANS_NEW

  INTERFACE MME_TRANS_INV
    module procedure mme_trans_inv_real64,mme_trans_inv_cmplx64
  END INTERFACE MME_TRANS_INV

  INTERFACE DEALLOCATE_GLOBAL_ARRAYS
    SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS( success )
      logical , optional :: success
      character(28) , parameter :: srcname = " in DEALLOCATE_GLOBAL_ARRAYS"
      integer isp
    END SUBROUTINE DEALLOCATE_GLOBAL_ARRAYS
  END INTERFACE DEALLOCATE_GLOBAL_ARRAYS

  ! INTERFACE
  !   subroutine close(unit,status)
  !     use, intrinsic :: ISO_FORTRAN_ENV
  !     use m_params, only: BUFFER_SIZE
  !     character(BUFFER_SIZE), optional :: status
  !     integer unit
  !     character(BUFFER_SIZE) msg
  !     integer iok
  !   end subroutine close
  ! END INTERFACE

  INTERFACE EXIT_ON_ERROR
    subroutine exit_on_error( from , stop )
      character(*), optional :: from
      logical, optional :: stop
    end subroutine exit_on_error
  END INTERFACE EXIT_ON_ERROR

CONTAINS

  CHARACTER(BUFFER_SIZE) FUNCTION INT32_2STRING(int,separator)
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character, optional :: separator
    integer(INT32) int
! local vars
    integer i

    write(int32_2string,*)int
    int32_2string=adjustl(int32_2string)
    if(int >= 0 .and. int < 10) &
      int32_2string(:)='0'//int32_2string(:len(int32_2string)-1)

    if(present(separator))then
      do i=len(trim(int32_2string))-2,2,-3
        int32_2string(:)=int32_2string(:i-1)//separator//&
                      int32_2string(i:len(int32_2string)-1)
      enddo
    endif
  END FUNCTION INT32_2STRING

  CHARACTER(BUFFER_SIZE) FUNCTION INT64_2STRING(int,separator)
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character, optional :: separator
    integer(INT64) int
! local vars
    integer i

    write(int64_2string,*)int
    int64_2string=adjustl(int64_2string)
    if(int >= 0 .and. int < 10) &
      int64_2string(:)='0'//int64_2string(:len(int64_2string)-1)

    if(present(separator))then
      do i=len(trim(int64_2string))-2,2,-3
        int64_2string(:)=int64_2string(:i-1)//separator//&
                      int64_2string(i:len(int64_2string)-1)
      enddo
    endif
  END FUNCTION INT64_2STRING

  SUBROUTINE ALLOC1DINT32( array , alias , ldim1 , udim1 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1
    logical, optional :: nullify
    integer, allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
   if( .not. present(udim1) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 ) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias_local , iok , msg )
    if( present(nullify) .and. nullify ) array(:) = 0
  END SUBROUTINE ALLOC1DINT32

  SUBROUTINE ALLOC1DREAL32( array, alias, ldim1, udim1, nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1
    logical, optional :: nullify
    real, allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg, alias_local
    integer iok, ldim1_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
   if( .not. present(udim1) ) then
!$omp critical
      write( unit = output_unit, fmt = "( 4x, a )" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 ), stat = iok, errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias_local, iok, msg )
    if( present(nullify) .and. nullify ) array(:) = 0.
  END SUBROUTINE ALLOC1DREAL32

  SUBROUTINE ALLOC1DREAL64( array , alias , ldim1 , udim1 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1
    logical, optional :: nullify
    real(REAL64), allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
   if( .not. present(udim1) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 ) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias_local , iok , msg )
    if( present(nullify) .and. nullify ) array(:) = 0._REAL64
  END SUBROUTINE ALLOC1DREAL64

  SUBROUTINE ALLOC2DINT32( array , alias , ldim1 , udim1 &
             , ldim2 , udim2 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1, ldim2, udim2
    logical, optional :: nullify
    integer, allocatable :: array(:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local, ldim2_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if( .not. present(ldim2) ) then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
   if( .not. present(udim1) .or. .not. present(udim2) ) then
!$omp critical
      msg(:)='All array upper indices should be set to allocate array '
      if( present(alias) ) msg(:)=trim(msg)//trim(alias)
      write( unit = output_unit , fmt = "(4x,a)" ) trim(msg)        
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 ) &
            , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_allocation_error( alias_local , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:) = 0
  END SUBROUTINE ALLOC2DINT32

  SUBROUTINE ALLOC2DREAL32( array, alias, ldim1, udim1, ldim2, udim2, nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1, ldim2, udim2
    logical, optional :: nullify
    real, allocatable :: array( : , : )
! local vars
    character(BUFFER_SIZE) msg, alias_local
    integer iok, ldim1_local, ldim2_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if( .not. present(ldim2) ) then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
   if( .not. present(udim1) .or. .not. present(udim2) ) then
!$omp critical
      msg(:) = 'All array upper indices should be set to allocate array '
      if( present(alias) ) msg(:) = trim(msg) // trim(alias)
      write( unit = output_unit, fmt = "( 4x, a )" ) trim(msg)        
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1, ldim2_local : udim2 ), &
         stat = iok, errmsg = msg)
    if( iok /= 0 ) &
      call print_allocation_error( alias_local, iok, msg )
    if( present(nullify) .and. nullify ) array(:,:) = 0.
  END SUBROUTINE ALLOC2DREAL32

  SUBROUTINE ALLOC2DREAL64( array , alias , ldim1 , udim1 &
             , ldim2 , udim2 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1, ldim2, udim2
    logical, optional :: nullify
    real(REAL64), allocatable :: array(:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local, ldim2_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if( .not. present(ldim1) ) then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if( .not. present(ldim2) ) then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
   if( .not. present(udim1) .or. .not. present(udim2) ) then
!$omp critical
      msg(:)='All array upper indices should be set to allocate array '
      if( present(alias) ) msg(:)=trim(msg)//trim(alias)
      write( unit = output_unit , fmt = "(4x,a)" ) trim(msg)        
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 ) &
            , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_allocation_error( alias_local , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:) = 0._REAL64
  END SUBROUTINE ALLOC2DREAL64

  SUBROUTINE ALLOC3DCOMPLEX64( array , alias , ldim1 , udim1 &
             , ldim2 , udim2 , ldim3 , udim3 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(24) , parameter :: srcname="ALLOCATE_ALLOC3DCOMPLEX64"
    character(*) , optional :: alias
    complex(REAL64) , allocatable :: array(:,:,:)
    complex(REAL64) , parameter :: z0 = ( 0._REAL64 , 0._REAL64 )
    integer(INT32) , optional :: ldim1, udim1, ldim2, udim2, ldim3, udim3
    logical , optional :: nullify
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local, ldim2_local, ldim3_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if(.not.present(ldim1))then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if(.not.present(ldim2))then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
    if(.not.present(ldim3))then
      ldim3_local=1
    else
      ldim3_local=ldim3
    endif
    if( .not. present(udim1) .or. .not. present(udim2) .or. &
        .not. present(udim3) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error(srcname)
    endif
    if ( udim1 < ldim1_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim1 < ldim1: ldim1=' &
        ,ldim1_local,' udim1=',udim1
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim2 < ldim2_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim2 < ldim2: ldim2=' &
        ,ldim2_local,' udim2=',udim2
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim3 < ldim3_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim3 < ldim3: ldim3=' &
        ,ldim3_local,' udim1=',udim3
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 &
            , ldim3_local : udim3 ) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:,:) = z0
  END SUBROUTINE ALLOC3DCOMPLEX64

  SUBROUTINE ALLOC3DINT32( array , alias , ldim1 , udim1 &
             , ldim2 , udim2 , ldim3 , udim3 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1, ldim2, udim2, ldim3, udim3
    logical, optional :: nullify
    integer, allocatable :: array(:,:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local, ldim2_local, ldim3_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if(.not.present(ldim1))then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if(.not.present(ldim2))then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
    if(.not.present(ldim3))then
      ldim3_local=1
    else
      ldim3_local=ldim3
    endif
    if( .not. present(udim1) .or. .not. present(udim2) .or. &
        .not. present(udim3) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 &
            , ldim3_local : udim3 ) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:,:) = 0
  END SUBROUTINE ALLOC3DINT32

  SUBROUTINE ALLOC3DREAL64( array , alias , ldim1 , udim1 &
             , ldim2 , udim2 , ldim3 , udim3 , nullify )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(22), parameter :: srcname="ALLOCATE_ALLOC3DREAL64"
    character(*), optional :: alias
    integer(INT32), optional :: ldim1, udim1, ldim2, udim2, ldim3, udim3
    logical, optional :: nullify
    real(REAL64), allocatable :: array(:,:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok, ldim1_local, ldim2_local, ldim3_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if(.not.present(ldim1))then
      ldim1_local=1
    else
      ldim1_local=ldim1
    endif
    if(.not.present(ldim2))then
      ldim2_local=1
    else
      ldim2_local=ldim2
    endif
    if(.not.present(ldim3))then
      ldim3_local=1
    else
      ldim3_local=ldim3
    endif
    if( .not. present(udim1) .or. .not. present(udim2) .or. &
        .not. present(udim3) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error(srcname)
    endif
    if ( udim1 < ldim1_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim1 < ldim1: ldim1=' &
        ,ldim1_local,' udim1=',udim1
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim2 < ldim2_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim2 < ldim2: ldim2=' &
        ,ldim2_local,' udim2=',udim2
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim3 < ldim3_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim3 < ldim3: ldim3=' &
        ,ldim3_local,' udim1=',udim3
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 &
            , ldim3_local : udim3 ) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:,:) = 0._REAL64
  END SUBROUTINE ALLOC3DREAL64

  SUBROUTINE ALLOC4DCOMPLEX64( array , alias , ldim1 , udim1 , ldim2 , &
             udim2 , ldim3 , udim3 , ldim4 , udim4 , nullify )
    use , intrinsic :: ISO_FORTRAN_ENV
    use m_params , only : BUFFER_SIZE
    implicit none
    character(22) , parameter :: srcname = "ALLOCATE_ALLOC4DREAL64"
    character(*) , optional :: alias
    complex(REAL64) , parameter :: z0 = ( 0._REAL64 , 0._REAL64)
    integer(INT32) , optional :: ldim1 , udim1 , ldim2 , udim2 , ldim3 , &
        udim3 , ldim4 , udim4
    logical, optional :: nullify
    complex(REAL64) , allocatable :: array( : , : , : , : )
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok , ldim1_local , ldim2_local , ldim3_local , ldim4_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if ( .not. present(ldim1) ) then
      ldim1_local = 1
    else
      ldim1_local = ldim1
    endif
    if ( .not. present(ldim2) ) then
      ldim2_local = 1
    else
      ldim2_local = ldim2
    endif
    if ( .not. present(ldim3) ) then
      ldim3_local = 1
    else
      ldim3_local = ldim3
    endif
    if ( .not. present(ldim4) ) then
      ldim4_local = 1
    else
      ldim4_local = ldim4
    endif
    if( .not. present(udim1) .or. .not. present(udim2) .or. &
        .not. present(udim3) .or.  .not. present(udim4) ) then
!$omp critical
      write( unit = output_unit , fmt = "( 4x, a )" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error(srcname)
    endif
    if ( udim1 < ldim1_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim1 < ldim1: ldim1=' &
        ,ldim1_local,' udim1=',udim1
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim2 < ldim2_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim2 < ldim2: ldim2=' &
        ,ldim2_local,' udim2=',udim2
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim3 < ldim3_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim3 < ldim3: ldim3=' &
        ,ldim3_local,' udim3=',udim3
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim4 < ldim4_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim4 < ldim4: ldim4=' &
        ,ldim4_local,' udim4=',udim4
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2, &
              ldim3_local : udim3 , ldim4_local : udim4 ) , &
              stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias , iok , msg )
    if( present(nullify) .and. nullify ) array( : , : , : , : ) = z0
  END SUBROUTINE ALLOC4DCOMPLEX64

  SUBROUTINE ALLOC5DCOMPLEX64( array , alias , ldim1 , udim1 , ldim2 , udim2 , &
             ldim3 , udim3 , ldim4 , udim4 , ldim5 , udim5 , nullify )
    use , intrinsic :: ISO_FORTRAN_ENV
    use m_params , only : BUFFER_SIZE
    implicit none
    character(22) , parameter :: srcname = "ALLOCATE_ALLOC5DREAL64"
    character(*) , optional :: alias
    complex(REAL64) , parameter :: z0 = ( 0._REAL64 , 0._REAL64)
    integer(INT32) , optional :: ldim1 , udim1 , ldim2 , udim2 , ldim3 , &
        udim3 , ldim4 , udim4 , ldim5 , udim5
    logical, optional :: nullify
    complex(REAL64) , allocatable :: array( : , : , : , : , : )
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok , ldim1_local , ldim2_local , ldim3_local , ldim4_local , &
            ldim5_local
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    if ( .not. present(ldim1) ) then
      ldim1_local = 1
    else
      ldim1_local = ldim1
    endif
    if ( .not. present(ldim2) ) then
      ldim2_local = 1
    else
      ldim2_local = ldim2
    endif
    if ( .not. present(ldim3) ) then
      ldim3_local = 1
    else
      ldim3_local = ldim3
    endif
    if ( .not. present(ldim4) ) then
      ldim4_local = 1
    else
      ldim4_local = ldim4
    endif
    if ( .not. present(ldim5) ) then
      ldim5_local = 1
    else
      ldim5_local = ldim5
    endif
    if( .not. present(udim1) .or. .not. present(udim2) .or. &
        .not. present(udim3) .or.  .not. present(udim4) .or. &
        .not. present(udim5) ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,a)" ) &
        'All array upper indices should be set to allocate array.'
!$omp end critical
      call exit_on_error(srcname)
    endif
    if ( udim1 < ldim1_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim1 < ldim1: ldim1=' &
        ,ldim1_local,' udim1=',udim1
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim2 < ldim2_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim2 < ldim2: ldim2=' &
        ,ldim2_local,' udim2=',udim2
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim3 < ldim3_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim3 < ldim3: ldim3=' &
        ,ldim3_local,' udim3=',udim3
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim4 < ldim4_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim4 < ldim4: ldim4=' &
        ,ldim4_local,' udim4=',udim4
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    if ( udim5 < ldim5_local ) then
!$omp critical
      write( unit = output_unit , fmt = "(4x,2(a,i0))" ) &
        'Cannot allocate '//trim(alias)//': udim5 < ldim5: ldim5=' &
        ,ldim5_local,' udim5=',udim5
      call exit_on_error( srcname , .true. )
!$omp end critical
    endif
    allocate( array( ldim1_local : udim1 , ldim2_local : udim2 &
            , ldim3_local : udim3 , ldim4_local : udim4 &
            , ldim5_local : udim5) , stat = iok , errmsg = msg )
    if( iok /= 0 ) &
      call print_allocation_error( alias , iok , msg )
    if( present(nullify) .and. nullify ) array(:,:,:,:,:) = z0
  END SUBROUTINE ALLOC5DCOMPLEX64

  SUBROUTINE DEALLOC1DINT32( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer, allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC1DINT32

  SUBROUTINE DEALLOC1DREAL32( array, alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    real, allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array, stat = iok, errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local, iok, msg )
  END SUBROUTINE DEALLOC1DREAL32

  SUBROUTINE DEALLOC1DREAL64( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    real(REAL64), allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC1DREAL64

  SUBROUTINE DEALLOC1DLMTSPEC( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    use m_rixs, only: lmto_spectrum
    implicit none
    character(*), optional :: alias
    type(lmto_spectrum), allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC1DLMTSPEC

  SUBROUTINE DEALLOC1DRIXSSPEC( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    use m_rixs, only: rixs_spectrum
    implicit none
    character(*), optional :: alias
    type(rixs_spectrum), allocatable :: array(:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC1DRIXSSPEC

  SUBROUTINE DEALLOC2DINT32( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer, allocatable :: array(:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC2DINT32

  SUBROUTINE DEALLOC2DREAL32( array, alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    real, allocatable :: array( : , : )
! local vars
    character(BUFFER_SIZE) msg, alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array, stat = iok, errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local, iok, msg )
  END SUBROUTINE DEALLOC2DREAL32

  SUBROUTINE DEALLOC2DREAL64( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    real(REAL64), allocatable :: array(:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC2DREAL64

  SUBROUTINE DEALLOC3DINT32( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    integer, allocatable :: array(:,:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC3DINT32

  SUBROUTINE DEALLOC3DREAL64( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    real(REAL64), allocatable :: array(:,:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC3DREAL64

  SUBROUTINE DEALLOC4DCOMPLEX64( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    complex(REAL64), allocatable :: array( : , : , : , : )
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC4DCOMPLEX64

  SUBROUTINE DEALLOC5DCOMPLEX64( array , alias )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only: BUFFER_SIZE
    implicit none
    character(*), optional :: alias
    complex(REAL64), allocatable :: array(:,:,:,:,:)
! local vars
    character(BUFFER_SIZE) msg , alias_local
    integer iok
    if( .not. present(alias) ) then
      alias_local(:) = 'unknown array'
    else
      alias_local(:) = alias(:)
    endif
    deallocate( array , stat = iok , errmsg = msg)
    if( iok /= 0 ) &
      call print_deallocation_error( alias_local , iok , msg )
  END SUBROUTINE DEALLOC5DCOMPLEX64

  SUBROUTINE CLOSE( unit , status )
    use, intrinsic :: ISO_FORTRAN_ENV
    use m_params, only : BUFFER_SIZE
    character(*), optional :: status
    integer unit
! local vars
    character(BUFFER_SIZE) msg
    integer iok
    logical opened

    if ( unit == 0 ) then
      write( unit = output_unit, fmt = 1)
      return
    endif

   1 format(4x,'WARNING: you try to close unit=0.' &
           ,'This is system unit and will not be closed now')

    inquire( unit = unit, opened = opened )
    if ( .not. opened ) then
      write( unit = output_unit, fmt = 3 ) unit
      return
    endif

   3  format( 4x , 'Unit=' , i0 , &
            ' is not connected to any file so it can not be closed' )

    if ( present( status ) ) then
      close( unit = unit, status = status, iomsg = msg, iostat = iok )
    else
      close( unit = unit, iomsg = msg, iostat = iok )
    endif
    if ( iok /= 0 ) then
      write( unit = output_unit, fmt = 5  )
      write( unit = output_unit, fmt = 10 ) unit
      write( unit = output_unit, fmt = 15 ) iok
      write( unit = output_unit, fmt = 20 ) trim(msg)
    endif

   5 format( /, 4x, 'THE ERROR OCCURIED WHILE TRYING TO CLOSE FILE' )
  10 format( 4x, 'Unit: ', i0 )
  15 format( 4x, 'Error code: ', i0 )
  20 format( 4x, 'Error message: ', a )

  END SUBROUTINE CLOSE

  FUNCTION MME_TRANS_REAL64(igk,mmein) RESULT(mmeout)
    use, intrinsic :: iso_fortran_env
!    use m_trans_matrix
    integer, intent(in) :: igk
    real(REAL64), intent(in) :: mmein(3)
    real(REAL64) mmeout(3)
! local vars
    real(REAL64) op(3,3)

    if(igk /= 1) op(:,:)=trans_matrix(igk)
    if(igk == 1)then
      mmeout(:)=mmein(:)
    else
      mmeout(:)=matmul(op,mmein)
    endif
  END FUNCTION MME_TRANS_REAL64

  FUNCTION MME_TRANS_CMPLX64(igk,mmein) RESULT(mmeout)
    use, intrinsic :: iso_fortran_env
    use m_bnd, only: iopnum
!    use m_trans_matrix
    integer, intent(in) :: igk
    complex(REAL64), intent(in) :: mmein(3)
    complex(REAL64) mmeout(3)
! local vars
    real(REAL64) op(3,3)

    if(igk /= 1) op(:,:)=trans_matrix(igk)
    if(iopnum(igk) > 0)then
      if(igk == 1)then
        mmeout(:)=mmein(:)
      else
        mmeout(:)=matmul(op,mmein)
      endif
    else
      mmeout(:)=matmul(op,conjg(mmein))
    endif
  END FUNCTION MME_TRANS_CMPLX64

  FUNCTION MME_TRANS_NEW_CMPLX64(op,mmein) RESULT(mmeout)
    use, intrinsic :: iso_fortran_env
    real(REAL64) op(3,3)
    complex(REAL64), intent(in) :: mmein(3)
    complex(REAL64) mmeout(3)
    mmeout(:)=matmul(op,conjg(mmein))
    mmeout(:)=matmul(op,mmein)
  END FUNCTION MME_TRANS_NEW_CMPLX64

  FUNCTION MME_TRANS_INV_REAL64(igk,mmein) RESULT(mmeout)
    use, intrinsic :: iso_fortran_env
!    use m_trans_matrix
    integer, intent(in) :: igk
    real(REAL64), intent(in) :: mmein(3)
    real(REAL64) mmeout(3)
! local vars
    real(REAL64) op(3,3)

    mmeout(:)=mmein(:)
    if(igk /= 1)then
      op(:,:)=transpose(trans_matrix(igk))
      mmeout(:)=matmul(op,mmein)
    endif
  END FUNCTION MME_TRANS_INV_REAL64

  FUNCTION MME_TRANS_INV_CMPLX64(igk,mmein) RESULT(mmeout)
    use, intrinsic :: iso_fortran_env
    use m_bnd, only: iopnum
!    use m_trans_matrix
    integer, intent(in) :: igk
    complex(REAL64), intent(in) :: mmein(3)
    complex(REAL64) mmeout(3)
! local vars
    real(REAL64) op(3,3)

    mmeout(:)=mmein(:)
    if(igk /= 1)then
      op(:,:)=transpose(trans_matrix(igk))
      if(iopnum(igk) > 0)then
        mmeout(:)=matmul(op,mmein)
      else
        mmeout(:)=matmul(op,conjg(mmein))
      endif
    endif
  END FUNCTION MME_TRANS_INV_CMPLX64

  FUNCTION TRANS_MATRIX(igk) RESULT(OP)
    use, intrinsic :: iso_fortran_env
    use m_bnd, only: g
    implicit none
    integer igk
    real(REAL64) op(3,3)

    if(igk <= 0) stop '    ERROR in TRANS_MATRIX:  igk <= 0.'
    op(:,:)=g(1:3,1:3,igk)
    op(:,:)=transpose(op)
  END FUNCTION TRANS_MATRIX

END MODULE M_FUNCTIONS

MODULE M_TIME
  use, intrinsic :: iso_fortran_env
  implicit none
  real(REAL64), private :: readmme
  real(REAL64), private :: real_RIXS,real_ABSORP_IBZ,real_ABSORP_RBZ &
      ,real_ABSORP_BZ
  real(REAL64), private :: CPU_RIXS,CPU_ABSORP_IBZ,CPU_ABSORP_RBZ &
      ,CPU_ABSORP_BZ
! not used now
  integer(INT64) nb_ticks_initial,nb_ticks_final
  integer(INT64), private :: nb_ticks_sec,nb_ticks_max

CONTAINS

  REAL(REAL64) FUNCTION GET_REAL_RIXS()
    get_real_RIXS=real_RIXS
  END FUNCTION GET_REAL_RIXS

  SUBROUTINE SET_REAL_RIXS(input)
    real(REAL64) input
    real_RIXS=input
  END SUBROUTINE SET_REAL_RIXS

  REAL(REAL64) FUNCTION GET_REAL_ABSORP_IBZ()
    get_real_ABSORP_IBZ=real_ABSORP_IBZ
  END FUNCTION GET_REAL_ABSORP_IBZ

  SUBROUTINE SET_REAL_ABSORP_IBZ(input)
    real(REAL64) input
    real_ABSORP_IBZ=input
  END SUBROUTINE SET_REAL_ABSORP_IBZ

  REAL(REAL64) FUNCTION GET_REAL_ABSORP_BZ()
    get_real_ABSORP_BZ=real_ABSORP_BZ
  END FUNCTION GET_REAL_ABSORP_BZ

  SUBROUTINE SET_REAL_ABSORP_BZ(input)
    real(REAL64) input
    real_ABSORP_BZ=input
  END SUBROUTINE SET_REAL_ABSORP_BZ

  REAL(REAL64) FUNCTION GET_REAL_ABSORP_RBZ()
    get_real_ABSORP_RBZ=real_ABSORP_RBZ
  END FUNCTION GET_REAL_ABSORP_RBZ

  SUBROUTINE SET_REAL_ABSORP_RBZ(input)
    real(REAL64) input
    real_ABSORP_RBZ=input
  END SUBROUTINE SET_REAL_ABSORP_RBZ

  REAL(REAL64) FUNCTION GET_CPU_RIXS()
    get_CPU_RIXS=CPU_RIXS
  END FUNCTION GET_CPU_RIXS

  SUBROUTINE SET_CPU_RIXS(input)
    real(REAL64) input
    CPU_RIXS=input
  END SUBROUTINE SET_CPU_RIXS

  REAL(REAL64) FUNCTION GET_CPU_ABSORP_IBZ()
    get_CPU_ABSORP_IBZ=CPU_ABSORP_IBZ
  END FUNCTION GET_CPU_ABSORP_IBZ

  SUBROUTINE SET_CPU_ABSORP_IBZ(input)
    real(REAL64) input
    CPU_ABSORP_IBZ=input
  END SUBROUTINE SET_CPU_ABSORP_IBZ

  REAL(REAL64) FUNCTION GET_CPU_ABSORP_BZ()
    get_CPU_ABSORP_BZ=CPU_ABSORP_BZ
  END FUNCTION GET_CPU_ABSORP_BZ

  SUBROUTINE SET_CPU_ABSORP_BZ(input)
    real(REAL64) input
    CPU_ABSORP_BZ=input
  END SUBROUTINE SET_CPU_ABSORP_BZ

  REAL(REAL64) FUNCTION GET_CPU_ABSORP_RBZ()
    get_CPU_ABSORP_RBZ=CPU_ABSORP_RBZ
  END FUNCTION GET_CPU_ABSORP_RBZ

  SUBROUTINE SET_CPU_ABSORP_RBZ(input)
    real(REAL64) input
    CPU_ABSORP_RBZ=input
  END SUBROUTINE SET_CPU_ABSORP_RBZ

  REAL(REAL64) FUNCTION GET_READMME()
    get_readmme=readmme
  END FUNCTION GET_READMME

  SUBROUTINE SET_READMME(input)
    real(REAL64) input
    readmme=input
  END SUBROUTINE SET_READMME

! not used now
  INTEGER(INT64) FUNCTION GET_NB_TICKS_SEC()
    get_nb_ticks_sec=nb_ticks_sec
  END FUNCTION GET_NB_TICKS_SEC

  SUBROUTINE SET_NB_TICKS_SEC(input)
    integer(INT64) input
    nb_ticks_sec=input
  END SUBROUTINE SET_NB_TICKS_SEC

  INTEGER(INT64) FUNCTION GET_NB_TICKS_MAX()
    get_nb_ticks_max=nb_ticks_max
  END FUNCTION GET_NB_TICKS_MAX

  SUBROUTINE SET_NB_TICKS_MAX(input)
    integer(INT64) input
    nb_ticks_max=input
  END SUBROUTINE SET_NB_TICKS_MAX

END MODULE M_TIME

