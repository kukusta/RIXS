SUBROUTINE CREATEINPUT
use m_bnd, only: nb, high_sym_pnt, high_sym_pnt_txt_dir
use m_files
use m_functions
use m_params ! , only: iabsorp, irixs, iprint, wmin, wmax, dw, write_hsym_points
use m_rixs, only: nrixslmt, lmttemp
use units, only: sry2ev
implicit none
! local vars
character key, key2
character(8) date
character(10) time
character(80) hostname
integer iok, isp, ncl, ncf, ipnt, ipnt2
real Gamma, energy
logical exist, cycle

inquire( file = inptfile, exist = exist )
iprint = 1
key = 'N'
key2 = ' '
if( exist )then
  do
    write( unit = output_unit, fmt = 1, advance = 'no' ) trim( inptfile )
    read(unit=input_unit,fmt='(a1)')key2
    if(key2 /= ' ')key=key2
    if(key == 'y' .or. key == 'Y') exit
    if(key == 'n' .or. key == 'N')then
      write(unit=output_unit,fmt=10)trim(inptfile)
      stop
    endif
  enddo
endif

  1 format(4x,'File ',a,' exists. Do you want to overwrite it? [y/N] ')
 10 format(4x,'File ',a,' has not been changed.')

! get data from main lmto program and set some defaults
if ( .not. check_kstar ) then
  call readbnd
  iprint=-1
!  call readrixsdata
  call set_ncf_ncl(ncf,ncl)
  iprint=1
else
  ncl = 0
  ncf = 0
  nb = 0
!  call readrixsdata
endif
call readrixsdata
write_hsym_points=.false.

open(newunit=inp,file=inptfile,action='write',status='unknown',iostat=iok)
if(iok /= 0) stop '    Input file cannot be opened, exiting.'

! I switched this header off
if ( .false. ) then
  call date_and_time(date=date,time=time)
  write(unit=inp,fmt=30)
  write(unit=inp,fmt=35,advance='no')date(1:4),date(5:6),date(7:8)&
       ,time(1:2),time(3:4),time(5:6)
  call get_environment_variable('HOSTNAME',value=hostname,status=iok)
  if(iok == 0)then
    write(unit=inp,fmt="(' on ',a,/)")trim(hostname)
  else
    write(unit=inp,fmt="(2(/))")
  endif
endif

 30 format(4x,'This input file for RIXS calculations was generated')
 35 format(4x,'automatically on ',a,'-',a,'-',a,1x,a,':',a,':',a)

write(unit=inp,fmt="(4x,'&PARAMETERS')")
write(unit=inp,fmt=2)irixs,iabsorp,iprint
if ( check_kstar ) then
  write(unit=inp,fmt=3)0,ncl,ncf,nb
else
  write(unit=inp,fmt=3)1,ncl,ncf,nb
endif
write(unit=inp,fmt=5)wmin,wmax,dw
write(unit=inp,fmt=6)trim(rixsfile),trim(map2deinfile),trim(map2dqfile)&
     ,trim(datafile),trim(bndfile),trim(sdtfile),trim(mmefile)
write(unit=inp,fmt=7) ! check_kstar
write(unit=inp,fmt=120)

 2 format(4x,'irixs= ',i0,4x,'iabsorp= ',i0,4x,'iprint= ',i0)
 3 format(4x,'! nbi(1:2)=',2(i5),4x,'nbf(1:2)=',2(i5))
 5 format(4x,'! wmin=',4x,g0.3,4x,'wmax=',4x,g0.12,4x,'dw=',4x,g0.5)
 6 format(4x,'! rixsfile= "',a,'"',4x,'map2deinfile= "',a,'"'&
         ,4x,'map2dqfile= "',a,'"',/,4x,'! datafile= "',a,'"'&
         ,4x,'bndfile= "',a,'"',4x,'sdtfile= "',a,'"',/&
         ,4x,'! mmefile= "',a,'"')
 7 format(4x,'k_in(1:3)= 1. 0. -1.',6x,'! input k-vector in lab. frame'&
         ,/,4x,'Qaxis(1:3)= 0. 0. 1.',6x,'! in lab. frame'&
         ,/,4x,'Abas(1:3,1:3)= 1. 0. 0.   ! x-vector of local frame'&
         ,/,19x,'0. 1. 0.   ! y-vector of local frame'&
         ,/,19x,'0. 0. 1.   ! z-vector of local frame' )
!         , /, 4x, 'check_kstar= ', L )

do isp=1,nrixslmt
  if(lmttemp(isp)%iwritemme == 0)cycle
  energy=((wmax+wmin)/2.-lmttemp(isp)%ec)*sry2ev
  call spwidth_cp01(lmttemp(isp)%txtel(:2),lmttemp(isp)%nn, &
       lmttemp(isp)%nk,Gamma)
  write(unit=inp,fmt="(4x,'&SPEC')")
  write(unit=inp,fmt=12)lmttemp(isp)%txtel(:2),lmttemp(isp)%sname(:)
  if ( check_kstar ) then
    write( unit = inp, fmt = 3 ) 0, ncl, ncf, nb
  else
    write(unit=inp,fmt=3)1,ncl,ncf,nb
  endif
  write(unit=inp,fmt=15)0.,0.,0.
  write(unit=inp,fmt=18)energy,Gamma
  write(unit=inp,fmt=21)
  write(unit=inp,fmt=120)
enddo
! write(unit=inp,fmt=100)lmttemp(nrixslmt)%txtel(:2),lmttemp(nrixslmt)%sname(:)
! write(unit=inp,fmt=3)1,ncl,ncf,nb
! write(unit=inp,fmt=105)10,'s-p '
! write(unit=inp,fmt=15)0.,0.,0.
! write(unit=inp,fmt=110)0.1,0.2,0.3
! write(unit=inp,fmt=18)energy,Gamma
! write(unit=inp,fmt=120)
! write(unit=inp,fmt=22)

 12 format(4x,'atom= "',a,'"',4x,'specname= "',a,'"',4x,"skip= F",4x &
          ,"use_symmetry= 'n'")
 15 format(4x,'q(1:3)= ',3(f7.3),4x,'! in r.l.u. in local frame')
 18 format(4x,'input_energy= ',f8.2,4x,'Gamma= ',f6.3,4x &
          ,'! both in eV')
 21 format(4x,'channel= "p-s "')
!  21 format(4x,'channel= "p-s "',4x,'save2dat= F',4x,'map2d_ein= F'&
!           ,4x,'map2d_q= F')
 22 format(/                                                          &
    ,4x,'IRIXS:',/                                                    &
    ,4x,'0 - Joint DOS (BZ integration with constant mme)',/          &
    ,4x,'1 - RIXS loss spectra',/                                     &
    ,4x,'2 - RIXS with matrix elements averaged over polarizations',/,/ &
    ! ,4x,'3 - dispersion of RIXS loss spectra',/,/                     &
    ,4x,'IABSORP:',/                                                  &
    ,4x,'0 - absorption spectra are not calculated',/                 &
    ,4x,'1 - IBZ integration',/                                       &
    ,4x,'2 - IBZ and BZ integration')
! 100 format(4x,'&Q_DISPERSION' &
!           ,/,4x,'atom= "',a,'"',4x,'specname= "',a,'"',4x,'skip= F',4x&
!           ,"use_symmetry= 'no'")
! 105 format(4x,'nq= ',i0,4x,'channel= "',a,'"')
! 110 format(4x,'dq(1:3)= ',3(f7.3),4x,'! in units of 2*pi/[a,b,c]')
 120 format(4x,'/',/)

if(write_hsym_points)then
  write(unit=inp,fmt=25)
  do ipnt=1,48
    if(len_trim(high_sym_pnt_txt_dir(ipnt)) == 0)exit
    cycle=.false.
    do ipnt2=1,ipnt-1
      if(high_sym_pnt_txt_dir(ipnt2) == high_sym_pnt_txt_dir(ipnt))then
        cycle=.true.
        exit
      endif
    enddo
    if(cycle)cycle
    write(unit=inp,fmt=45)high_sym_pnt(:,ipnt),high_sym_pnt_txt_dir(ipnt)
  enddo
endif

 25 format(/,4x,'HIGH SYMMETRY POINTS:')
 45 format(4x,3(f14.10),4x,"'",a,"'")

if(exist)then
  write(unit=output_unit,fmt=9)trim(inptfile)
else
  write(unit=output_unit,fmt=11)trim(inptfile)
endif
! close(unit=inp,iostat=iok)

  9 format(4x,'File ',a,' has been overwritten.')
 11 format(4x,'File ',a,' has been created.')

call deallocate_global_arrays
stop

END SUBROUTINE CREATEINPUT

