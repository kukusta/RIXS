SUBROUTINE WRITERIX
use m_bnd , only : ef
use m_files
use m_functions
use m_params
use m_results
use m_rixs
use units
implicit none
! local vars
character(BUFFER_SIZE) msg , datfile , spec_type
integer iok , i , idat , isp , nerixs
logical :: exist
real , allocatable :: loss2(:)
real maxRIXS, maxRIXSispec, maxRIXSenergy
real(REAL64), parameter :: r0 = 1.e-5

if(iprint > 0) &
  write(unit=output_unit,fmt='(4x,"Starting WRITERIX")')

if ( abs(irixs) /= 0 .and. abs(irixs) /= 10 ) then
  spec_type = 'RIXS loss'
else
  spec_type = 'Joint DOS'
endif

if( .not. read_rixfile ) then
  inquire(file=rixsfile,exist=exist)
! old code
!$$$ i=1
!$$$ do while(exist)
!$$$   new_rixsfile(:)=trim(rixsfile)//'.~'//trim(int2string(i))//'~'
!$$$   inquire(file=new_rixsfile,exist=exist)
!$$$   if(.not. exist)then
!$$$     write(unit=output_unit,fmt="(4x,'File ',a,' exists, RIXS spectra '&
!$$$          &,'will be written to file ',a)")trim(rixsfile),trim(new_rixsfile)
!$$$     rixsfile(:)=new_rixsfile(:)
!$$$     exit
!$$$   endif
!$$$   i=i+1
!$$$ enddo
! I replaced it with new code. Reason: there is no 'rename' function in
! Fortran standart, so old code did not really rename rixsfile but stored new
! spectra to new_rixsfile
  if(exist) &
    write(unit=output_unit &
         &,fmt="(4x,'File ',a,' exists and will be overwritten')" &
         )trim(rixsfile)
  open(newunit=rix,file=rixsfile,status='replace',form='unformatted' &
      ,iostat=iok)
  if(iok /= 0) stop '    Error while opening rixfile.'

  if(abs(irixs) == 1 .or. abs(irixs) == 2)then
! open m2d files only if some spectra should be saved
    m2d_ein=0
    m2d_q=0
    do isp=1,nrixs
      if(rixs(isp)%map2d_ein) m2d_ein=m2d_ein+1
      if(rixs(isp)%map2d_q) m2d_q=m2d_q+1
    enddo

    if(m2d_ein > 0)then
      inquire(file=map2deinfile,exist=exist)
      if(exist) &
        write(unit=output_unit &
             &,fmt="(4x,'File ',a,' exists and will be overwritten')" &
             )trim(map2deinfile)
      open(newunit=m2d_ein,file=map2deinfile,status='replace',iostat=iok)
      if(iok /= 0) stop '    Error while opening map2deinfile.'
      write(unit=m2d_ein,fmt=2)"input energy relative to Fermi energy"
      write(unit=output_unit,fmt=3)trim(map2deinfile)
    endif

    if(m2d_q > 0)then
      inquire(file=map2dqfile,exist=exist)
      if(exist) &
        write(unit=output_unit &
             &,fmt="(4x,'File ',a,' exists and will be overwritten')" &
             )trim(map2dqfile)
      open( newunit = m2d_q , file = map2dqfile , status = 'replace' , &
           iostat = iok )
      if ( iok /= 0 ) stop '    Error while opening map2dqfile.'
      write( unit = m2d_q , fmt = 2 ) "qx, qy, qz"
      write( unit = output_unit , fmt = 3 ) trim(map2dqfile)
    endif
  endif

 2 format("# 1st column is energy loss",/&
         ,"# 2nd column is ",a,/&
         ,"# 3rd column is normalized RIXS (max. RIXS is 1.0)")
 3 format(4x,"Data for 2D map will be written to file ",a)

  if(iabsorp == 1 .or. iabsorp == 2)then
    write(rix)nrixs,nlabel,1+iabsorp
  else
    write(rix)nrixs,nlabel,1+0
  endif
  nerixs=int((emax-emin)/de+1)
  if(irixs < 0)then
    write(rix)nw,-real(wmax),-real(wmin),real(dw)
  else
    write(rix)nw,real(wmin),real(wmax),real(dw)
  endif
  write(rix)nerixs,real(emin),real(emax),real(de),real(ef)

  if(irixs < 0)then
    allocate(loss2(nw),stat=iok,errmsg=msg)
    if(iok /= 0)call print_allocation_error('loss2',iok,msg)
  endif

! check if there are negative values of RIXS loss intensity
  do isp = 1, nrixs
    do i = 1, nw
      if ( loss(i, isp) < 0._REAL64 ) then
        if( abs( loss(i, isp) ) > r0 ) &
          write( unit = output_unit, fmt = 5 ) isp, trim(spec_type), i, &
  &            loss(i,isp)
        loss(i, isp) = 0._REAL64
      endif
    enddo
  enddo
  write( unit = output_unit, fmt = "(/)" )

 5 format( 4x, 'ERROR: for isp=', i0, ' ', a, '(', i0, ')= ', g0.4, &
         ' and would be set to 0' )

! find max RIXS value for m2d file and check if there are
! negative values of RIXS loss intensity
  maxRIXS = -1.e10
  do isp = 1, nrixs
    maxRIXSispec = -1.e10
    maxRIXSenergy = 0.
    if( irixs < 0 ) then
      do i = 1, nw
        loss2(i) = real( loss(nw - (i-1) , isp) )
        if( loss2(i) > maxRIXS ) then
          maxRIXS = loss2(i)
        endif
        if( loss2(i) > maxRIXSispec ) then
          maxRIXSispec = loss2(i)
          maxRIXSenergy = -real( wmin + (i - 1) * dw )
        endif
      enddo
    else
      do i = 1, nw
        if( loss(i, isp) > maxRIXS ) then
          maxRIXS = loss(i, isp)
        endif
        if( loss(i, isp) > maxRIXSispec ) then
          maxRIXSispec = loss(i, isp)
          maxRIXSenergy = real( wmin + (i - 1) * dw )
        endif
      enddo
    endif
    write( unit = output_unit, fmt = 7 ) isp, maxRIXSispec, &
         maxRIXSenergy * sry2ev
  enddo

  7 format( 4x, 'Maximum value for ISPEC=', i0, ' is ', f15.5, &
          ' at E= ', f10.5, ' eV' )
  8 format(2(4x,f11.6),4x,f20.6)
  9 format(4(4x,f11.6),4x,f20.6)
 10 format(4x,'Datefile ',a,' exists and will be overwritten')

  if( abs(irixs) == 1 .or. abs(irixs) == 2 .and.  &
                (m2d_ein /= 0 .or. m2d_q /= 0) )  &
    write( unit = output_unit, fmt = "( 4x, a )") &
         'Data for 2D plot would be normalized.'

  do isp = 1, nrixs
    rixs(isp)%isort = 0
    write(rix) rixs(isp)%label, rixs(isp)%isort, real( rixs(isp)%ec ), &
         real( rixs(isp)%dnhsort ), real( rixs(isp)%enlmax ), &
         rixs(isp)%nn, rixs(isp)%nk, rixs(isp)%txtel(:)
    if ( irixs < 0 ) then
      do i = 1, nw
        loss2(i) = real( loss(nw-(i-1), isp) )
      enddo
      write(rix) loss2(:)
    else
      write(rix) real( loss(:,isp) )
    endif
    do i = 1, iabsorp
      write(rix) real( absorp(:,isp,i) )
    enddo

! save data to map2d files
    if(abs(irixs) == 1 .or. abs(irixs) == 2)then
      if(rixs(isp)%map2d_ein)then
        do i=1,nw
          if(irixs < 0)then
            write(unit=m2d_ein,fmt=8)-(wmin+(i-1)*dw)*dry2ev&
                 ,(rixs(isp)%en-ef+rixs(isp)%ec)*dry2ev&
                 ,loss2(i)/maxRIXS
          else
            write(unit=m2d_ein,fmt=8)(wmin+(i-1)*dw)*dry2ev&
                 ,(rixs(isp)%en-ef+rixs(isp)%ec)*dry2ev&
                 ,real(loss(i,isp))/maxRIXS
          endif
        enddo
      endif
      if(rixs(isp)%map2d_q)then
        do i=1,nw
          if(irixs < 0)then
            write(unit=m2d_q,fmt=9)-(wmin+(i-1)*dw)*dry2ev&
                 ,rixs(isp)%q(:),loss2(i)/maxRIXS
          else
            write(unit=m2d_q,fmt=9)(wmin+(i-1)*dw)*dry2ev&
                 ,rixs(isp)%q(:),real(loss(i,isp))/maxRIXS
          endif
        enddo
      endif
    endif
  enddo
  if ( irixs < 0 ) then
    deallocate( loss2 , stat = iok )
    if ( iok /= 0 ) write( unit = output_unit , fmt = 15 ) iok
  endif
endif
do isp = 1, nrixs
! save data to datfile if required
  if(rixs(isp)%save2dat)then
    i=index(rixsfile,'.rix')
    datfile(:)=rixsfile(:i-1)//"_ispec_"//trim(int2string(isp))//".dat"
    inquire(file=datfile,exist=exist)
    if(exist) write(unit=output_unit,fmt=10)trim(datfile)
    open(newunit=idat,file=datfile,status='replace',iostat=iok)
    if(iok /= 0) stop '    Error while opening datfile.'
!    write( unit = idat , fmt = "( 4x, 3f8.3, 4a )" ) 0.07, 0.3, 1.0, &
!         "    'L'", "    'L'", '    F', &
!         ' / offset, sigma_0, omega_0 in eV; type_sigma, type_smearing;' &
!         // ' ask to overwrite outfile?'
    do i=1,nw
      if(irixs < 0)then
        write(unit=idat,fmt=8)-(wmin+(i-1)*dw)*dry2ev&
             ,(rixs(isp)%en-ef+rixs(isp)%ec)*dry2ev,loss2(i)
      else
        write( unit = idat, fmt = 8 ) (wmin + (i - 1) * dw ) * dry2ev, &
!$$$               ,(rixs(isp)%en-ef+rixs(isp)%ec)*dry2ev,
             (wmin + (i - 1) * dw ) * dry2ev / 1.25, real( loss( i, isp ) )
      endif
    enddo
    close(unit=idat,iostat=iok)
  endif
enddo

! if ( .not. ( irixs == 0 .or. abs(irixs) == 10 ) ) then
!   write( unit = output_unit , fmt = 30 )
!   do isp = 1 , nrixs
!     write( unit = output_unit , fmt = 40 ) isp , rixs(isp)%min_mme , &
!       rixs(isp)%max_mme , rixs(isp)%av_mme / real( rixs(isp)%ntr , &
!       kind = REAL64 ) , rixs(isp)%ntr
!   enddo
! endif

if ( iprint > 0 ) write( unit = output_unit , fmt = 50 )

 15 format( 4x, 'Can not deallocate loss2(:), error code is ', i0)
 30 format( /, 4x, 'Maximal values of matrix elements:', /, &
           4x, 'ISPEC',4x,'min(abs(mme))',4x,'max(abs(mme))', &
           4x, 'avr(abs(mme))', 4x, 'rixs%ntr')
 40 format( 4x, i3, 4x, f0.10, 4x, f0.10, 4x, f0.10, 4x, i13 )
 50 format( /, 4x, 'Rixsfile is written, rixs is done' )

END SUBROUTINE WRITERIX

