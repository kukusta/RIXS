SUBROUTINE READSDT
use m_bnd, only: nopused,iopnum
use m_files, only: sdt,sdtfile
use m_params
use m_sdt
implicit none
! local vars
integer iok, it(64), inver, iaddparam, icub, iop, j, ngsdt, ngsdt2
logical exist
real(REAL64) rbassdt( 3, 3 ), vg48( 3, 48 )

if ( read_rixfile ) return
if ( check_kstar .and. debug_mode ) return

inquire(file=sdtfile,exist=exist)
if(.not. exist)then
  write(unit=output_unit,fmt=5)trim(sdtfile)
  stop
endif

  5 format(4x,'File ',a,' does not exist.')

open(newunit=sdt,file=sdtfile,action='read',status='unknown'&
    ,form='unformatted',position='rewind',iostat=iok)
if(iok /= 0) stop '    Sdtfile cannot be opened, exiting.'

! rewind(unit=sdt)
read(unit=sdt,iostat=iok)
if(iok /= 0)then
  close(unit=sdt,iostat=iok)
  stop '    General error in sdtfile.'
endif
!read(unit=sdt,iostat=iok)
read(unit=sdt,iostat=iok)rbassdt(:,:),ut(:,:),ut1(:,:)

if(iok /= 0)then
  close(unit=sdt,iostat=iok)
  stop '    General error in sdtfile.'
endif
read(unit=sdt,iostat=iok)inver,iaddparam,ngsdt,icub,it(:),ngsdt2&
    ,(vg48(:,j),j=1,ngsdt)
select case(iok)
  case(IOSTAT_EOR)
    stop '    In READSDT: end-of-record condition occured.'
  case(IOSTAT_END)
    stop '    In READSDT: end-of-file condition occured.'
  case default
end select
if(iok /= 0 .or. ngsdt /= ngsdt2)then
  close(unit=sdt,iostat=iok)
  stop '    General error in sdtfile.'
endif
!$$$ close(unit=sdt,iostat=iok)
allocate(vg(3,nopused),stat=iok)
if(iok /= 0)stop '    In READSDT: cannot allocate memory for vg array.'
do iop=1,nopused
  do j=1,ngsdt
    if(abs(iopnum(iop)) == it(j)) vg(:,iop)=vg48(:,j)
  enddo
enddo

! temporarily
! write(*,*)'    Possible fractional translations are:'
! do iop=1,nopused
!   write(*,"('    iop=',i3,' vg(:)= ',3f15.10)")iop,vg(:,iop)
! enddo
vg(:,:)=0._REAL64

END SUBROUTINE READSDT

