SUBROUTINE READRIXSDATA
use m_bnd, only: natom,nopused,nb,npnt,ef
use m_files
use m_functions
use m_params
use m_rixs
use m_sdt, only: igtta
use units, only: dry2ev
implicit none
! local vars
character(4) sname
character(6) txtel
character(256) msg,mmefile_save
complex(REAL64), parameter :: z0=(0._REAL64,0._REAL64)
integer iok,isp,jsp,nn,nk,writemme,iat1,iatn,label,ia1,ian,isort&
     ,nrixslmt2
logical exist,first,stop,skip
logical, external :: incompatible,saved
real(REAL64), allocatable :: esp(:),rat(:,:)
real(REAL64) ec,enlmin,enlmax,dnhsort

namelist /common_data/ nrixslmt,ispc,max_nhsort,nlabel,n_mme,emin,emax&
     ,de,nfu,volomg
namelist /static/ nn,nk,writemme,iat1,iatn,label,txtel,sname,ec,enlmin&
     ,enlmax,dnhsort,ia1,ian,isort,skip
namelist /esp_rat/ esp,rat
namelist /misc/ igtta,mmefile

if ( read_rixfile ) return

inquire(file=datafile,exist=exist)
if(.not. exist)then
  write(unit=output_unit,fmt=5)trim(datafile)
  stop
endif

  5 format(/,4x,'File ',a,' does not exist.')

open(newunit=rid,file=datafile,action='read',status='unknown'&
    ,position='rewind',iostat=iok)
if(iok /= 0) stop '    datafile cannot be opened, exiting.'

read(unit=rid,nml=common_data,iomsg=msg,iostat=iok)
if(iok /= 0)then
  write(unit=output_unit,fmt=10)trim(msg)
  stop
endif

 10 format(4x,'Error in datafile in &COMMON_DATA section: ',a)

if(dw == 0._REAL64)dw=de
if(wmin == 0._REAL64 .and. wmax == 0._REAL64)then
  wmin=0._REAL64
  wmax=emax-emin
  if((wmax-wmin)/dw >= real(huge(0),kind=REAL64)) &
    stop '    ERROR: number of points on omega mesh is out of range for INT32.'
  nw=int((wmax-wmin)/dw)+1
endif

! Joint DOS without absorption => no mme needed
if(irixs == 0 .and. iabsorp == 0)return

allocate(lmttemp(nrixslmt),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('lmttemp',iok,msg)

 15 format(4x,'Error in datafile in &STATIC section for ISPEC #',i0,': ',a)
 20 format(4x,'Matrix elements for SPEC=',i0,' are supposed to be skipped'&
          ," ('skip' flag is set to TRUE).",/,4x,'If this is done '&
          ,'intentionally then debug mode should be enabled'&
          ,' (set debug_mode= T in PARAMETERS).')
 25 format(4x,'For SPEC #',i0,' nn= ',i0,' nk= ',i0,' and sname= ',a&
          ,' are incompatible.')
 30 format(4x,'Error in datafile in &ESP_RAT section for SPEC #',i0,': ',a)

! if(iprint > 0)write(unit=output_unit,fmt="()")
nrixslmt2=0
stop=.false.
do isp=1,nrixslmt
  skip=.false.
  read(unit=rid,nml=static,iomsg=msg,iostat=iok)
  if(iok /= 0)then
    write(unit=output_unit,fmt=15)isp,trim(msg)
!ccc    stop=.false. why???
    stop=.true.
    cycle
  endif
  if(skip .and. .not. debug_mode)then
    write(unit=output_unit,fmt=20)isp
    call deallocate_global_arrays
    stop
  endif
  if(writemme == 0 .or. skip)cycle

  if(incompatible(nn,nk,sname))then
    write(unit=output_unit,fmt=25)isp,nn,nk,sname
    stop
  endif

  nrixslmt2=nrixslmt2+1
  lmttemp(nrixslmt2)%nn=nn
  lmttemp(nrixslmt2)%nk=nk
  lmttemp(nrixslmt2)%iwritemme=writemme
  lmttemp(nrixslmt2)%ia1=ia1
  lmttemp(nrixslmt2)%ian=ian
  lmttemp(nrixslmt2)%iat1=iat1
  lmttemp(nrixslmt2)%iatn=iatn
  lmttemp(nrixslmt2)%label=label
  lmttemp(nrixslmt2)%txtel(:)=txtel(:)
  lmttemp(nrixslmt2)%sname(:)=sname(:)
  lmttemp(nrixslmt2)%ec=ec
  lmttemp(nrixslmt2)%enlmin=enlmin
  lmttemp(nrixslmt2)%enlmax=enlmax
  lmttemp(nrixslmt2)%dnhsort=dnhsort
  lmttemp(nrixslmt2)%isort=isort

  allocate(esp(ia1:ian),rat(3,iat1:iatn),lmttemp(nrixslmt2)%esp(ia1:ian) &
          ,lmttemp(nrixslmt2)%rat(3,iat1:iatn),stat=iok,errmsg=msg)
  if(iok /= 0) call print_allocation_error('lmttemp%esp,lmttemp%rat',iok,msg)

  read(unit=rid,nml=esp_rat,iomsg=msg,iostat=iok)
  if(iok /= 0)then
    write(unit=output_unit,fmt=30)isp,trim(msg)
    call deallocate_global_arrays
    stop
  endif
  lmttemp(nrixslmt2)%esp(:)=esp(:)
  lmttemp(nrixslmt2)%rat(:,:)=rat(:,:)
  deallocate(esp,rat,stat=iok)
enddo
! search for duplicates
do isp=1,nrixslmt
  do jsp=isp+1,nrixslmt
    if(lmttemp(jsp)%nn    == lmttemp(isp)%nn    .and.     &
       lmttemp(jsp)%nk    == lmttemp(isp)%nk    .and.     &
       lmttemp(jsp)%isort == lmttemp(isp)%isort .and.     &
       lmttemp(jsp)%iwritemme*lmttemp(isp)%iwritemme /= 0 &
      )then
      write(unit=output_unit,fmt=33)isp,jsp
      stop=.true.
    endif
  enddo
enddo
if(stop)then
  call deallocate_global_arrays
  stop
endif

 33 format(4x,"Duplicated mme found in lmtdatafile: ISPEC=",i2," and ",i2)

allocate(igtta(nopused,natom),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('igtta',iok,msg)
! mmefile from datafile is used only if it has not been set in inptfile
mmefile_save(:)=""
if(len_trim(mmefile) /= 0)mmefile_save(:)=mmefile(:)
read(unit=rid,nml=misc,iomsg=msg,iostat=iok)
if(iok /= 0)then
  write(unit=output_unit,fmt=35)trim(msg)
  call deallocate_global_arrays
  stop
endif
if(len_trim(mmefile_save) /= 0)mmefile(:)=mmefile_save(:)
!$$$ close(unit=rid,iostat=iok)
! this check is performed in READRIXSMME subroutine
! if(.not. newini)then
!   inquire(file=mmefile,exist=exist)
!   if(.not. exist)then
!     write(unit=output_unit,fmt=5)trim(mmefile)
!     call deallocate_global_arrays
!     stop
!   endif
! endif

 35 format(4x,'Error in datafile in &MISC section: ',a)
 37 format(/,4x,'For SPEC=',i0,' the energy of incoming x-ray is not set.')

allocate(lmtdata(nrixslmt2),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('lmtdata',iok,msg)

! used to analyze absorption edges below
rixs(1:nrixs)%ec=1._REAL64

nlmtdata=0
do isp=1,nrixs
  first=.true.
  do jsp=1,nrixslmt2
    if(rixs(isp)%sname(:) == lmttemp(jsp)%sname(:))then
      if(rixs(isp)%txtel(:2) == lmttemp(jsp)%txtel(:2))then
        rixs(isp)%nn=lmttemp(jsp)%nn
        rixs(isp)%nk=lmttemp(jsp)%nk
        rixs(isp)%label=lmttemp(jsp)%label
! these four lines are not clear. When RIXS is calculated for several atom
! sorts with different ec, dnhsort, enlmax values these values do not have
! unambigous values. They can be different for different sorts.
        rixs(isp)%dnhsort=lmttemp(jsp)%dnhsort
        rixs(isp)%enlmax=lmttemp(jsp)%enlmax
        rixs(isp)%isort=lmttemp(jsp)%isort
! sometimes there are several nonequivalent atoms with the same Z
! here warning is printed if this takes place
        if(rixs(isp)%ec > 0._REAL64)then
! always lmt(jsp)%ec < 0
          rixs(isp)%ec=lmttemp(jsp)%ec
          txtel(:)=lmttemp(jsp)%txtel(:)
        else
          if(first)then
            first=.false.
            write(unit=output_unit,fmt=50)isp,txtel(:)
            write(unit=output_unit,fmt=55,advance='no')txtel(:)
          endif
          write(unit=output_unit,fmt=60,advance='no')lmttemp(jsp)%txtel(:)
        endif

 50 format(/,4x,'For SPEC=',i0," the absorption edge is defined for '",a,"'")
 55 format(4x,"Inequivalent species with the same Z are '",a,"'")
 60 format(", '",a,"'")

! if for some spectrum input energy is less than 0 it is given in eV
! relative to E_F. If it is positive it is given in absolute units.
! It's value in absolute units is usually large and difficult to define
! If input energy is exactly zero, then it means some mistake (error
! message is printed).
!$$$ ! very small negative value, error for -0.001 eV < input energy < 0 eV
!$$$ ! if input energy is not set in inrfile than by default input energy is
!$$$ ! set equal to -0.0001 eV and causes an error here. input energy must
!$$$         if(rixs(isp)%en > -7.e-5_REAL64 .and. rixs(isp)%en < -1.e-6_REAL64)then
!$$$           write(unit=output_unit,fmt=37)isp
!$$$           call deallocate_global_arrays
!$$$           stop
!$$$         endif
! negative or zero input energy, relative to E_F
        if(rixs(isp)%en < 1.e-6_REAL64)then
          rixs(isp)%en=ef-rixs(isp)%ec-rixs(isp)%en
        endif
        if(.not. saved(lmttemp(jsp)%sname,lmttemp(jsp)%txtel,lmttemp(jsp)%isort))then
          nlmtdata=nlmtdata+1
          ia1=lmttemp(jsp)%ia1
          ian=lmttemp(jsp)%ian
          iat1=lmttemp(jsp)%iat1
          iatn=lmttemp(jsp)%iatn
          allocate(lmtdata(nlmtdata)%esp(ia1:ian)&
                  ,lmtdata(nlmtdata)%rat(3,iat1:iatn)&
                  ,lmtdata(nlmtdata)%mme(3,iat1:iatn,ia1:ian,nb,npnt)&
                  ,stat=iok,errmsg=msg)
          if(iok /= 0)call print_allocation_error('lmtdata%esp,rat',iok,msg)
          lmtdata(nlmtdata)%sname(:)=lmttemp(jsp)%sname(:)
          lmtdata(nlmtdata)%txtel(:)=lmttemp(jsp)%txtel(:)
          lmtdata(nlmtdata)%nn=lmttemp(jsp)%nn
          lmtdata(nlmtdata)%nk=lmttemp(jsp)%nk
          lmtdata(nlmtdata)%isort=lmttemp(jsp)%isort
          lmtdata(nlmtdata)%ia1=lmttemp(jsp)%ia1
          lmtdata(nlmtdata)%ian=lmttemp(jsp)%ian
          lmtdata(nlmtdata)%iat1=lmttemp(jsp)%iat1
          lmtdata(nlmtdata)%iatn=lmttemp(jsp)%iatn
          lmtdata(nlmtdata)%label=lmttemp(jsp)%label
          lmtdata(nlmtdata)%ec=lmttemp(jsp)%ec
          lmtdata(nlmtdata)%enlmin=lmttemp(jsp)%enlmin
          lmtdata(nlmtdata)%enlmax=lmttemp(jsp)%enlmax
          lmtdata(nlmtdata)%dnhsort=lmttemp(jsp)%dnhsort
          lmtdata(nlmtdata)%esp(:)=lmttemp(jsp)%esp(:)
          lmtdata(nlmtdata)%rat(:,:)=lmttemp(jsp)%rat(:,:)
        endif
      endif
    endif
  enddo
  if(.not.first)write(unit=output_unit,fmt="()")
enddo

stop=.false.
do isp=1,nrixs
  if(rixs(isp)%nn == 0)then
    write(unit=output_unit&
         ,fmt=40)rixs(isp)%txtel(:2),trim(rixs(isp)%sname)
    stop=.true.
  endif
enddo
if(stop) stop

 40 format(4x,'No mme have been found for ',a,' at ',a,'-edge')

if(.not. newini)then
  do isp=1,nrixslmt
    deallocate(lmttemp(isp)%esp,lmttemp(isp)%rat,stat=iok)
  enddo
  deallocate(lmttemp,stat=iok)
endif

! do isp=1,nrixs
!   ! rixs(isp)%q(:)=matmul(ut1,rixs(isp)%q)
!   ! rixs(isp)%q(:)=matmul(ut,rixs(isp)%q)
!   rixs(isp)%q(:)=matmul(qbas,rixs(isp)%q)/2._REAL64/dpi()
! enddo

! min=huge(0._REAL64)
! do isp=1,nrixs
!   call set_ncf_ncl(ncf,ncl)
!   do ib=ncf,nb
!     do k=1,npnt
!       diff=e(ib,k)-rixs(isp)%ec
!       if(diff < min)then
!         ib_min=ib
!         k_min=k
!       endif
!     enddo
!   enddo
! enddo
! if(iprint > 0) write(unit=output_unit,fmt=45)
!  45 format(4x,'Resonant energies have been calculated.')

END SUBROUTINE READRIXSDATA

LOGICAL FUNCTION INCOMPATIBLE(nn,nk,sname)
implicit none
integer, intent(in) :: nn,nk
character(4), intent(in) ::  sname
character :: tn(1:5)=['K','L','M','N','O']
character(3) :: tk(-4:3)=['VII','V  ','III','I  ','   ','II ','IV ','VI ']

incompatible=.false.
if(nn == 1)then
  if(tn(nn)//'   ' /= sname)incompatible=.true.
else
  if(tn(nn)//tk(nk) /= sname)incompatible=.true.
endif

END FUNCTION INCOMPATIBLE

LOGICAL FUNCTION SAVED(sname,txtel,isort)
use m_rixs, only: lmtdata, nlmtdata
implicit none
!input vars
character(4) sname
character(6) txtel
integer isort
! local vars
integer isp

saved=.false.
do isp=1,nlmtdata
  if(lmtdata(isp)%sname(:) == sname(:) .and. &
     lmtdata(isp)%txtel(:2) == txtel(:2) .and. &
     lmtdata(isp)%isort == isort) saved=.true.
enddo

END FUNCTION SAVED

