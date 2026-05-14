SUBROUTINE READINPUT
use m_files
use m_functions
use m_params
use m_rixs, only: k_in, nrixs, rixs, Abas, Qaxis
use units, only: dry2ev
implicit none
! local vars
character(6) atom, version
character(4) specname, channel, use_symmetry
character(256) buf
integer i, iarg, iinr, ilmt, iok, nspecread, nbi_save(2), nbf_save(2)
logical exist, skip, save2dat, map2d_ein, map2d_q
logical save2dat_from_parameters, map2d_ein_from_parameters, &
        map2d_q_from_parameters
real input_energy, Gamma, q(3), e_in(3), e_out(3)
real(REAL64), external :: ddet33

namelist /parameters/ irixs, iabsorp, iprint, nbi, nbf, wmin, wmax, dw, &
     k_in, rixsfile, bndfile, bnsfile, datafile, sdtfile, mmefile, &
     map2deinfile, map2dqfile, debug_mode, Abas, Qaxis, absorption_only, &
     save2dat, map2d_ein, map2d_q, check_kstar, read_rixfile
namelist /spec/ atom, specname, q, input_energy, Gamma, channel, nbi, &
     nbf, skip, e_in, e_out, save2dat, map2d_ein, map2d_q, use_symmetry

if(command_argument_count() < 1)then
  stop ' No input found, exiting.'
endif

call s4ver(version)
write( unit = output_unit, fmt = 35 ) version(:)

 35 format( /, 4x, 'RESONANT INELASTIC X-RAY SCATTERING CALCULATIONS', &
            /, 4x, 'based on PY-LMTO band structure code Ver. ', a, / )

! Defaults, PART I
nbi(:) = 0
nbf(:) = 0
iabsorp = 0
iprint = 1
irixs = 1
wmin = 0._REAL64
wmax = 0._REAL64
dw = 0._REAL64
Abas(:,:) = 0._REAL64
qaxis(:) = 0._REAL64
k_in(:) = real( (/ 0., 0., 0. /) , kind = REAL64 )
debug_mode = .false.
absorption_only = .false.
save2dat = .false.
map2d_ein = .false.
map2d_q = .false.
check_kstar = .false.
! Initially all the files are closed
bnd = 0
bns = 0
inp = 0
rid = 0
sdt = 0
rix = 0
rim = 0
m2d_q = 0
m2d_ein = 0

! inptfile or lmtfile is always the first argument
basename(:)     = ""
inptfile(:)     = ""
bndfile(:)      = ""
bnsfile(:)      = ""
datafile(:)     = ""
sdtfile(:)      = ""
rixsfile(:)     = ""
runfile(:)      = ""
map2deinfile(:) = ""
map2dqfile(:)   = ""
mmefile(:)      = ""
call get_command_argument( 1, buf )
ilmt = index( buf, '.lmt', back = .true. )
iinr = index( buf, '.inr', back = .true. )
if ( ilmt /= 0 ) then
  basename( : ilmt ) = buf( : ilmt )
else
  if ( iinr /= 0 ) then
    basename( : iinr ) = buf( : iinr )
  else
    basename(:) = buf(:)
  endif
endif
do iarg = 2, command_argument_count()
  call get_command_argument( iarg, buf )
  if ( trim(buf) == '-k' ) check_kstar = .true.
enddo
! lmtfile as input
if(ilmt > 0)then
  newini = .true.
  inptfile(:)     = trim( basename ) // 'inr'
  bndfile(:)      = trim( basename ) // 'bnd'
  bnsfile(:)      = trim( basename ) // 'bns'
  datafile(:)     = trim( basename ) // 'rid'
  sdtfile(:)      = trim( basename ) // 'sdt'
  rixsfile(:)     = trim( basename ) // 'rix'
  runfile(:)      = trim( basename ) // 'run'
  map2deinfile(:) = trim( basename ) // 'e2d'
  map2dqfile(:)   = trim( basename ) // 'q2d'
  mmefile(:)      = ""
  call createinput
else
  newini = .false.
  if( iinr == 0 ) then
    write( unit = output_unit, fmt = 5 ) trim( basename )
    stop
  endif
  inptfile(:)     = trim( basename ) // 'inr'
! call add_file(buf(:i)//'inr',0,'inptfile')
! call add_file(buf(:i)//'inr',0)
  bndfile(:)      = trim( basename ) // 'bnd'
  datafile(:)     = trim( basename ) // 'rid'
  sdtfile(:)      = trim( basename ) // 'sdt'
  rixsfile(:)     = trim( basename ) // 'rix'
  runfile(:)      = trim( basename ) // 'run'
  map2deinfile(:) = trim( basename ) // 'e2d'
  map2dqfile(:)   = trim( basename ) // 'q2d'
! Empty filename means that mmefile would be read from datafile.
! To be reread it should be empty after reading inptfile as well,
! i.e. 'mmefile' variable in &PARAMETERS section should be commented.
  mmefile(:) = "" ! buf(:i)//'rim'
endif

 5 format( 4x, 'Wrong input file: ', a )

do i = 2, command_argument_count()
  call get_command_argument(i,buf)
  if ( trim(buf) == '-kstar' ) check_kstar = .true.
enddo

call print_files

inquire( file = inptfile, exist = exist )
if( .not. exist ) then
  write( unit = output_unit , fmt = 10 ) trim( inptfile )
  stop
endif

 10 format(4x,'In READINPT: file ',a,' does not exist.')
 15 format(4x,'Error in input file in &PARAMETERS: ',a)

! if(.not. allocated(headfile))&
!   write(*,fmt="(4x,'TEST MESSAGE: headfile is not allocated',/)")

open(newunit=run,file=runfile,status='unknown',iostat=iok)
if(iok /= 0) stop '    Runfile cannot be created, exiting.'

open(newunit=inp,file=inptfile,action='read',position='rewind'&
    ,status='unknown',iostat=iok)
if(iok /= 0) stop '    Input file cannot be opened, exiting.'

! read data from inptfile and overwrite defaults if neccessary
read( unit = inp, nml = parameters, iomsg = buf, iostat = iok )
if(iok /= 0)then
  write(unit=output_unit,fmt=15)trim(buf)
  stop
endif
save2dat_from_parameters=save2dat
map2d_ein_from_parameters=map2d_ein
map2d_q_from_parameters=map2d_q
! try to check for some errors
if(nbi(1) < 0) stop '    Error in input data: nbi(1) < 0.'
if(nbi(2) < 0) stop '    Error in input data: nbi(2) < 0.'
if(nbi(1) > nbi(2)) stop '    Error in input data: nbi(1) > nbi(2).'
if(nbf(1) < 0) stop '    Error in input data: nbf(1) < 0.'
if(nbf(2) < 0) stop '    Error in input data: nbf(2) < 0.'
if(nbf(1) > nbf(2)) stop '    Error in input data: nbf(1) > nbf(2).'
nbi_save(:)=nbi(:)
nbf_save(:)=nbf(:)
if(abs(irixs) > 3 .and. abs(irixs) /= 10)then
  write(unit=output_unit,fmt=20)
  call deallocate_global_arrays
  stop
endif
! if(irixs /= 0 .and. abs(irixs) /= 2)then
if( .not. (irixs /= 0 .or. abs(irixs) /= 10))then
  if(sqrt(dot_product(Qaxis,Qaxis)) < 0.01_REAL64)then
    write(unit=output_unit,fmt=22)
    call deallocate_global_arrays
    stop
  endif
  if(ddet33(Abas) < 0.999_REAL64)then
    write(unit=output_unit,fmt=24)
    call deallocate_global_arrays
    stop
  endif
endif
! name of bnsfile is defined using bndfile
i = index( bndfile, '.', back = .true. )
bnsfile(:) = bndfile( : i ) // 'bns'

 20 format(4x,'ERROR in READINPUT: abs(irixs) > 3 and abs(irixs) /= 10 found.')
 22 format(4x,'ERROR in READINPUT: Qaxis length equals to 0.')
 24 format(4x,'ERROR in READINPUT: matrix A does not define a valid '&
          ,'Cartesian coordinate system.')
 25 format(4x,'WARNING: iabsorp > 2 found. iabsorp=2 is set.')
 27 format(4x,'WARNING: iabsorp < 0 found. iabsorp=0 is set.')

if(iabsorp > 2)then
  iabsorp=2
  write(unit=output_unit,fmt=25)
endif

if(iabsorp < 0)then
  iabsorp=0
  write(unit=output_unit,fmt=27)
endif

if(wmin < 0._REAL64) stop '    ERROR  in READINPUT: wmin < 0.'
if(wmin > wmax) stop '    ERROR  in READINPUT: wmin > wmax.'
! if energies are given in eV
if(dw < 0._REAL64)then
  dw=-dw/dry2ev
  wmin=wmin/dry2ev
  wmax=wmax/dry2ev
endif
if(wmin > 0._REAL64 .or. wmax > 0._REAL64)then
  if(dw == 0._REAL64) &
    stop '    ERROR  in READINPUT: wmin or wmax > 0 but dw = 0.'
  if((wmax-wmin)/dw >= real(huge(0),kind=REAL64))then
    write(unit=output_unit,fmt=29)(wmax-wmin)/dw
    call deallocate_global_arrays
    stop
  endif
  nw=int((wmax-wmin)/dw)+1
endif
if(sqrt(dot_product(k_in,k_in)) < 1.e-2 .and. &
       .not. (irixs /= 0 .or. abs(irixs) /= 10))then
  write(unit=output_unit,fmt=30)
  call deallocate_global_arrays
  stop
endif

 29 format(4x,'ERROR in READINPUT: number of points on omega mesh '&
          ,e11.4,' is out of range for INT32.')
 30 format(4x,'ERROR in READINPUT: norm(k_in) < 1.e-2.')
 40 format(4x,'Joint DOS spectra',15x,a)
 45 format(4x,'RIXS loss spectra',15x,a)
 50 format(4x,'RIXS with averaged mme',10x,a)
 55 format(4x,'RIXS fractional translations',4x,a)
 57 format(4x,'RIXS Lorentzian',17x,a)
 60 format(4x,'Absorption over IBZ',13x,a)
 65 format(4x,'Absorption over RBZ',13x,a)

buf(:)='disabled'
if(irixs == 0 .or. irixs == 10)buf(:)='enabled'
write(unit=output_unit,fmt=40)trim(buf)
buf(:)='disabled'
if(abs(irixs) > 0 .and. abs(irixs) /= 10)buf(:)='enabled'
write(unit=output_unit,fmt=45)trim(buf)
buf(:)='disabled'
if(abs(irixs) == 2)buf(:)='enabled'
write(unit=output_unit,fmt=50)trim(buf)
! buf(:)='disabled'
! if(irixs > 0)buf(:)='enabled'
! write(unit=output_unit,fmt=55)trim(buf)
buf(:)='disabled'
if(irixs == 3)buf(:)='enabled'
write(unit=output_unit,fmt=57)trim(buf)
buf(:)='disabled'
if(iabsorp == 1 .or. iabsorp == 2)buf(:)='enabled'
write(unit=output_unit,fmt=60)trim(buf)
buf(:)='disabled'
if(iabsorp == 2)buf(:)='enabled'
write(unit=output_unit,fmt=65)trim(buf)

nrixs=0
nspecread=0
do
  atom=""
  specname=""
  e_in(:) = 0.
  e_out(:) = 0.
  q(:)=0.
! default input energy in eV is set to E_F level (absorption edge)
  input_energy=0._REAL64
  Gamma=0.
  channel=""
  nbi(:)=0
  nbf(:)=0
  use_symmetry='n'
  save2dat=save2dat_from_parameters
  map2d_ein=map2d_ein_from_parameters
  map2d_q=map2d_q_from_parameters
  skip = .false.
  nspecread=nspecread+1
  read(unit=inp,nml=spec,iomsg=buf,iostat=iok)
  if(iok == IOSTAT_END)exit
  if(iok /= 0)then
    write(unit=output_unit,fmt=75)nspecread,trim(buf)
    cycle
  endif
  if(skip)cycle
  call add_spec( atom, specname, q, input_energy, Gamma, channel, nbi, &
       nbf, e_in, e_out, use_symmetry, save2dat, map2d_ein, map2d_q )
enddo

 69 format(4x,'In input file for &SPEC #',i0,": value use_symmetry= 'm'"&
          ,' is set. This case is not implemented yet.')
 70 format(4x,'Error in input file in &SPEC #',i0,": wrong value "&
          ,"use_symmetry= '",a,"'. use_symmetry= 'none' will be used.")
 75 format(4x,'Error in input file in &SPEC #',i0,': ',a)
 80 format(4x,'Error in input file in &Q_DISPERTION: ',a)
 85 format(4x,'WARNING: for SPEC=',i0,' Gamma value ',f7.4,' is too small')

! nrixs is set in subroutine 'add_spec'
if(nrixs == 0) stop '    No valid data in input file found.'

! Later I use 'nn' field to check if mme for ispec where written to mmefile
rixs(:)%nn=0

rixs(:)%ntr = 0
rixs(:)%min_mme = 0._REAL64
rixs(:)%max_mme = 0._REAL64
rixs(:)%av_mme = 0._REAL64

if ( .not. ( irixs == 0 .or. abs(irixs) == 10 ) ) then
  do i = 1 , nrixs
    if ( rixs(i)%Gamma < 0.1_REAL64 ) &
      write( unit = output_unit , fmt = 85 ) i , rixs(i)%Gamma
  enddo
endif

nbi(:) = nbi_save(:)
nbf(:) = nbf_save(:)

END SUBROUTINE READINPUT

SUBROUTINE ADD_SPEC( atom, specname, q, input_energy, Gamma, channel, &
     nbi, nbf, e_in, e_out, use_symmetry, save2dat, map2d_ein, map2d_q )
use m_functions
use m_params, only: MAX_NRIXS
use m_rixs
use units, only: sry2ev
implicit none
character(4) specname,channel,use_symmetry
character(6) atom
! character(*) atom,specname,channel
character(256) msg
complex(REAL64), parameter :: z0=(0._REAL64,0._REAL64)
integer iok,nbi(2),nbf(2)
logical save2dat,map2d_ein,map2d_q
real q(3),input_energy,Gamma,e_in(3),e_out(3)
type(rixs_spectrum) this,rixs2(nrixs+1)

this%txtel(:)=atom(:)
this%sname(:)=specname(:)
this%q(:)=real(q(:),kind=REAL64)
this%en=real(input_energy/sry2ev,kind=REAL64)
! this%en=real(input_energy,kind=REAL64)
this%Gamma=real(Gamma/sry2ev,kind=REAL64)
this%channel(:)=channel(:)
this%nbi(:)=nbi(:)
this%nbf(:)=nbf(:)
this%e_in(:)=cmplx(e_in,kind=REAL64)
this%e_out(:)=cmplx(e_out,kind=REAL64)
this%use_symmetry(:)=use_symmetry(:)
this%save2dat=save2dat
this%map2d_ein=map2d_ein
this%map2d_q=map2d_q
! this%denom_min=z0
rixs2(nrixs+1)=this

if(nrixs > 0)rixs2(:nrixs)=rixs(:)
if(allocated(rixs))then
  deallocate(rixs,stat=iok)
  if(iok /= 0) stop '    Can not deallocate rixs array.'
endif
nrixs=nrixs+1
allocate(rixs(nrixs),stat=iok,errmsg=msg)
if(iok /= 0)call print_allocation_error('rixs',iok,msg)
rixs(:)=rixs2(:)

if(nrixs > MAX_NRIXS)then
  write(unit=output_unit,fmt=90)nrixs,MAX_NRIXS
  call exit_on_error
endif

 90 format( 4x, 'Number of RIXS spectra is too large: nrixs=', i0, &
            ' MAX_NRIXS=', i0 )

END SUBROUTINE ADD_SPEC
