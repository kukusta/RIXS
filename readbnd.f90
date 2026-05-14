SUBROUTINE READBND
use m_files, only: bnd,bndfile
use m_bzmesh, only: nkibz,ndxyz
use m_functions
use m_bnd
use m_params, only: check_kstar, debug_mode, read_rixfile
implicit none
! local vars
character(11), parameter :: srcname=' in READBND'
character(8) named
character(256) msg
integer idummy,k,iok,ntxtd,ndiv,ibz
logical exist
real, allocatable :: e4(:)
real ef4,rbas4(3,3),qbas4(3,3)    

if ( read_rixfile ) return
if ( check_kstar .and. debug_mode ) return

inquire(file=bndfile,exist=exist)
if(.not. exist)then
  write(unit=output_unit,fmt=5)trim(bndfile)
  call exit_on_error( srcname )
endif

 5 format(4x,'File ',a,' does not exist.')

open(newunit=bnd,file=bndfile,action='read',status='unknown',&
    form='unformatted',iostat=iok)
if(iok /= 0) stop '    Bndfile cannot be opened, exiting.'

read(bnd,iostat=iok)named(:)
if(iok /= 0) stop '    Cannot read named(:) from bndfile.'
if(named(8:8).eq.'+')then
  ntxtd=4
else
  ntxtd=1
endif
nopused=-1
high_sym_pnt_txt_dir(:)=""
read(unit=bnd,iostat=iok)natom,ef4,npnt,ndiv,nb,idummy,idummy,idummy,ibz&
    ,idummy,(high_sym_pnt(:,k),high_sym_pnt_txt_dir(k)(1:ntxtd),k=1,idummy)&
    ,nopused,iopnum(:),rbas4(:,:),qbas4(:,:)
if(iok /= 0)nopused=-1
if(ibz /= 0) stop '    RIXS does not work with extended BZ mesh.'
nkibz=npnt
ef=real(ef4,kind=REAL64)
rbas(:,:)=real(rbas4(:,:),kind=REAL64)
qbas(:,:)=real(qbas4(:,:),kind=REAL64)
allocate(e4(nb),e(nb,npnt),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('e,e4',iok,msg)

do k=1,npnt
  read(bnd,iostat=iok)
  read(bnd,iostat=iok)e4(:)
  if(iok /= 0) stop '    READBND: can not read e(nb,k).'
  e(:,k)=real(e4(:),kind=REAL64)
enddo

if(ndiv < 10000) stop '    RIXS works with Bloechl BZ division only.'
call ndiv2nnn(ndiv,ndxyz(1),ndxyz(2),ndxyz(3))

END SUBROUTINE READBND

