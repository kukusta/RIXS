      program JOINT_BNS
      use, intrinsic :: ISO_FORTRAN_ENV
      implicit none

      character named1*8,txtdir1(48)*4
      character named2*8,txtdir2(48)*4
      character(256) bnsfile, bnsfile1, bnsfile2
      character(256) fbfile, fbfile1, fbfile2
      integer natom1, npnt1, ndiv1, n61, nspin1, ndimspin1, icss1, ibz1
      integer nbndr1, nopused1
      integer natom2, npnt2, ndiv2, n62, nspin2, ndimspin2, icss2, ibz2
      integer nbndr2, nopused2
      integer bns, bns1, bns2
      integer i, ien, ind, ind_max, ind_min, iok, ispin, j, k, nspin
      integer n6, npnt, ntxtd
      real(REAL64) epsk
      real ef, ef1, ef2, par11, boa1, coa1, par12, boa2, coa2
      real pnt01(3,64), pnt02(3,64)
      real, allocatable :: e(:,:,:),e1(:,:,:),e2(:,:,:)
      real, allocatable :: q1(:,:),q2(:,:)

      if( command_argument_count() < 2 )then
        call usage
        stop
      endif
! input file for energies below ef
      call get_command_argument(1, bnsfile1)
! input file for energies above ef
      call get_command_argument(2, bnsfile2)
      open(newunit = bns1, file = bnsfile1, action = 'read',            &
     &    status = 'unknown', form = 'unformatted', iostat = iok)
      if(iok /= 0) stop ' Bnsfile1 cannot be opened, exiting.'
      open(newunit = bns2, file = bnsfile2, action = 'read',            &
     &    status = 'unknown', form = 'unformatted', iostat = iok)
      if(iok /= 0) stop ' Bnsfile2 cannot be opened, exiting.'
! output band file
      ind = index(bnsfile2, '.', back = .true.)
      if(ind <= 0) then
        write(*,*) 'wrong filename: ', trim(bnsfile2)
        stop
      endif
      bnsfile(:) = bnsfile2( : ind - 1 )//'_JBNS.bns'
      open(newunit = bns, file = bnsfile, action = 'write',             &
     &    status = 'unknown', form = 'unformatted', iostat = iok)
      if(iok /= 0) stop ' Bnsfile cannot be opened, exiting.'
! fat bands for energies below ef
      call get_command_argument(3, fbfile1)
! fat bands for energies below ef
      call get_command_argument(4, fbfile2)

      epsk=10._REAL64*tiny(0._REAL64)
! read bnsfile1 and bnsfile2 and check their headers
      read(bns1, iostat = iok) named1, par11, boa1, coa1
      if(named1(8:8) == '+')then
        ntxtd = 4
      else
        ntxtd = 1
      endif
      nopused1 = -1
      read(bns1, iostat = iok) natom1, ef1, npnt1, ndiv1, n61, nspin1,  &
     &    ndimspin1, icss1, ibz1, nbndr1, (pnt01(:, j), txtdir1(j)(1 :  &
     &    ntxtd), j = 1, nbndr1), nopused1
      if( nopused1 /= -1 ) then
        write(*,*)'nopused /= -1 in bnsfile: ',trim(bnsfile1)
        stop 'is it the valid bnsfile? Exiting on error.'
      endif

      read(bns2, iostat = iok) named2, par12, boa2, coa2
      if(named2(8:8) == '+')then
        ntxtd=4
      else
        ntxtd=1
      endif
      nopused2 = -1
      read(bns2, iostat = iok) natom2, ef2, npnt2, ndiv2, n62, nspin2,  &
     &    ndimspin2, icss2, ibz2, nbndr2, (pnt02(:, j), txtdir2(j)(1 :  &
     &    ntxtd), j = 1, nbndr2), nopused2
      if( nopused2 /= -1 ) then
        write(*,*)'nopused /= -1 in bnsfile: ',trim(bnsfile2)
        stop 'is it the valid bnsfile? Exiting on error.'
      endif

! check parameters
      if(natom1 /= natom2) stop ' natom1 /= natom2'
      if(abs (ef1 - ef2 ) > epsk ) stop ' ef1 /= ef2'
      if(npnt1 /= npnt2) stop ' npnt1 /= npnt2'
      if(ndiv1 /= ndiv2) stop ' ndiv1 /= ndiv2'
      if(n61 /= n62) stop ' n61 /= n62'
      if(nspin1 /= nspin2) stop ' nspin1 /= nspin2'
      if(ndimspin1 /= ndimspin2) stop ' ndimspin1 /= ndimspin2'
      if(icss1 /= icss2) stop ' icss1 /= icss2'
      if(ibz1 /= ibz2) stop ' ibz1 /= ibz2'
      if(nbndr1 /= nbndr2) stop ' nbndr1 /= nbndr2'

      allocate(e1(nspin1, npnt1, n61), stat = iok)
      if(iok /= 0)then
        write(*,*)'can not allocate e1 array'
        stop
      endif
      allocate(q1(npnt1, 3), stat = iok)
      if(iok /= 0)then
        write(*,*)'can not allocate q1 array'
        stop
      endif
      do k = 1, npnt1
        read(bns1, iostat = iok) q1(k, :)
        if(iok /= 0)then
          write(*,*)'can not read q1(:) for k=', k
          stop
        endif
        do i = 1, 3
          if( abs( q1(k, i) ) < epsk ) q1(k, i) = 0.
        enddo
        do ispin = 1, nspin1
          read(bns1, iostat = iok)e1(ispin, k, :)
          if(iok /= 0)then
            write(*,*)'cannot read energies for ispin=', ispin,', k=', k
            stop
          endif
        enddo
      enddo

      allocate(e2(nspin2, npnt2, n62), stat = iok)
      if(iok /= 0)then
        write(*,*)'can not allocate e2 array'
        stop
      endif
      allocate(q2(npnt2, 3), stat = iok)
      if(iok /= 0)then
        write(*,*)'can not allocate q2 array'
        stop
      endif
      do k = 1, npnt2
        read(bns2, iostat = iok) q2(k, :)
        if(iok /= 0)then
          write(*,*)'can not read q2(:) for k=', k
          stop
        endif
        do i = 1, 3
          if( abs( q2(k, i) ) < epsk ) q2(k, i) = 0.
        enddo
        do ispin = 1, nspin2
          read(bns2, iostat = iok)e2(ispin, k, :)
          if(iok /= 0)then
            write(*,*)'cannot read energies for ispin=', ispin,', k=', k
            stop
          endif
        enddo
      enddo

      nspin = nspin1
      n6 = n61
      npnt = npnt1
      ef = ef1
      allocate(e(nspin, npnt, n6), stat = iok)
      if(iok /= 0)then
        write(*,*)'can not allocate e array'
        stop
      endif
      do ispin = 1, nspin
        ind_max = 0
        do k = 1, npnt
          do ien = 1, n6
            if(e1(ispin, k, ien) > ef)then
              if(ind_max < ien - 1) ind_max = ien - 1
              exit
            endif
          enddo       
        enddo
        do k = 1, npnt
          do ien = 1, n6
            if(e2(ispin, k, ien) > ef)then
              ind = ien - 1
              exit
            endif
          enddo
          if(ind > ind_max) stop ' ERROR: ind > ind_max'
          do ien = 1, ind_max
            e(ispin, k, ien) = min(e1(ispin, k, ien),ef)
          enddo
        enddo
        ind_min = n6
        do k = 1, npnt
          do ien = 1, n6
            if(e2(ispin, k, ien) > ef)then
              if(ind_min > ien) ind_min = ien
              exit
            endif
          enddo       
        enddo
        if(ind_max - ind_min < 0) stop ' ind_max - ind_min < 0'
        do k = 1, npnt
          do ien = ind_min, n6 - (ind_max - ind_min + 1)
              e(ispin, k, ind_max - ind_min + 1 + ien ) =               &
     &        max(e2(ispin, k, ien),ef)
          enddo
        enddo
      enddo

! write data to final bnsfile
      write(bns, iostat = iok) named1, par11, boa1, coa1
      if(named1(8:8) == '+')then
        ntxtd = 4
      else
        ntxtd = 1
      endif
      write(bns, iostat = iok) natom1, ef1, npnt1, ndiv1, n61, nspin1,  &
     &    ndimspin1, icss1, ibz1, nbndr1, (pnt01(:, j), txtdir1(j)(1 :  &
     &    ntxtd), j = 1, nbndr1)
      do k = 1, npnt
        write(bns, iostat = iok) q1(k, :)
        if(iok /= 0)then
          write(*,*)'can not write q(:) for k=', k
          stop
        endif
        do ispin = 1, nspin
          write(bns, iostat = iok)e(ispin, k, :)
          if(iok /= 0)then
            write(*,*)'cannot write energies for ispin=',ispin,', k=',k
            stop
          endif
        enddo
      enddo

      close(bns)
      close(bns1)
      close(bns2)
      deallocate(e)
      deallocate(q1)
      deallocate(e1)
      deallocate(q2)
      deallocate(e2)

      write(*,*)'joint band structure has been written to '             &
     &     ,trim(bnsfile)

      end program JOINT_BNS

      subroutine usage
      implicit none
      write(*,*)"check input files"
      end subroutine usage

!-----------lmt.fbs-----------
      subroutine fbs_u2f(lunin,lunout,ierr)
      dimension q(3),irl(0:8)
      allocatable ql(:),ils(:)

      ierr=0
      read(lunin,err=10,end=10)npntfb,nbndfb,nqlfb,nsort,nspeig,irel
      allocate(ql(nqlfb),ils(nsort),stat=iok)
      if(iok.ne.0)goto 10
      read(lunin,err=10,end=10)ils(1:nsort),nspb
      write(lunout,*)npntfb,nbndfb,nqlfb,nsort,nspeig,irel
      write(lunout,*)ils(1:nsort),nspb
      if(nspb.eq.-1)then
        do isrt=1,nsort
          read(lunin,err=10,end=10)irl(0:ils(isrt))
          write(lunout,'(14i4)')irl(0:ils(isrt))
        enddo
        nspb=1
      endif
      do isp=1,nspeig
        do k=1,npntfb
          read(lunin,err=10,end=10)q(1:3)
          write(lunout,"(3f10.6,' /k',i5)")q(1:3),k
          do ib=1,nbndfb
            read(lunin,err=10,end=10)n,e
            write(lunout,*)n,e
            do ispb=1,nspb
              read(lunin,err=10,end=10)ql(1:nqlfb)
              write(lunout,*)ql(1:nqlfb)
            enddo
          enddo                         ! ib
        enddo                           ! k
      enddo                             ! isp
      deallocate(ql,ils,stat=iok)
      return
 10   continue
      if(allocated(ql))deallocate(ql,ils,stat=iok)
      ierr=1
      end
