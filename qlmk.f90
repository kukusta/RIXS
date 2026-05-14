      function lreclcp()

      use m_param
      use m_atsrtl, only: ils

      lrecl=0
      if(ibnd_only.eq.0.and.iwrscf.ne.-2)then
        nls=nsort+sum(ils(1:nsort))
        lrecl=mxqlm*nls*limzn           ! real(8)
      endif
      lreclcp=lrecl*lr8                 ! record length in bytes
      end

      subroutine cptest(lrecl)

      use m_param
      include 'PREC.FI'	
      include 'FILES.FI'

      logical opened
     
      if(ibnd_only.gt.0.or.iwrscf.eq.-2)return
      inquire(unit=14,opened=opened)    !,recl=itrecl)
      if(opened)then                    ! lrecl may change
        close(14,status='delete',iostat=iostat)
      endif
      open(unit=14,file=fileq,access='direct',recl=lrecl
     $     ,status='unknown',form='unformatted') !,buffered='yes')
      end

      subroutine qlread(k,qlmk,ispin)

! The routine reads QLM file when IWRBN>0 or restores qlmk from alm(k)

      use m_param
      use m_atsrtl, only: ils
      use m_alm

      INCLUDE 'PREC.FI'

      common /reclen/ lrecl_cp

      dimension qlmk(0:lmax,nsort,limzn,*)

      nqlm=mxqlm
      nsp=max(ipertr,1)
      if(iwrbn.eq.0)then                ! copy data from arrays
        do i=1,nqlm
          iq=(ispin-1)*nqlm+i
          do ib=1,nband
            j=0
            do isort=1,nsort
              do l=0,ils(isort)
                j=j+1
                qlmk(l,isort,ib,i)=alm(j,ib,k,iq)
              enddo
            enddo
          enddo
        enddo     
      else                              ! Read QLM from file
! may be it is safer to check if the fileq is opened
        ir=(k-1)*nspeig+ispin           !+1
        read(14,rec=ir)((((qlmk(l,is,i,j),l=0,ils(is)),is=1,nsort)
     $       ,i=1,nband),j=1,nqlm)
      endif
      end

      subroutine qlwrt(k,qlmk,ils,ispin)
!*******************************************************
!  This program writes QLM to control point file       *
!  or copies qlmk to alm(k)                            *
!*******************************************************

      use m_param
      use m_alm

      INCLUDE 'PREC.FI'
      INCLUDE 'FILES.FI'

      integer ils(*)
      dimension qlmk(0:lmax,nsort,limzn,*)

      if(iwrbn.lt.0)return
      if(ibnd_only.gt.0.or.iwrscf.eq.-2)return
      nqlm=mxqlm
      if(iwrbn.eq.0)then                ! copy data from arrays
        do i=1,nqlm
          iq=(ispin-1)*nqlm+i
          do ib=1,nband
            j=0
            do isort=1,nsort
              do l=0,ils(isort)
                j=j+1
                alm(j,ib,k,iq)=qlmk(l,isort,ib,i)
              enddo
            enddo
          enddo
        enddo     
      else
!$OMP CRITICAL(QLWRT_W)
        ir=(k-1)*nspeig+ispin           !+1
        write(14,rec=ir)((((qlmk(l,is,i,j),l=0,ils(is)),is=1,nsort)
     $       ,i=1,nband),j=1,nqlm) 
!$OMP END CRITICAL(QLWRT_W)
      endif
      end
