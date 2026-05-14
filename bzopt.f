      subroutine bzopt(itetr,idold,ntet,qb,fnorm,ebi,ebf,a,ef,xi)
c
      implicit double precision (a-h,o-z)
c
      dimension ebi(*),ebf(*),itetr(0:4,ntet),idold(4,ntet),qb(*) ! input
      real a(*)
      dimension xi(*)                   ! output
      dimension eti(4),etf(4),etf1(4),at(4),at1(4),qt(3,4),qt1(3,4) ! local
     $     ,ip(4),del(4)
c
      common /cor/ emin,emax,de,npe
c
      efi=ef
      eff=ef
c
      do 500 i=1,ntet
c assignment of values for one tetrahedron
        nile=0
        nige=0
        nfle=0
        nwge=0
        nwle=0
        blwnorm=itetr(0,i)*fnorm
        do k=1,4
          i1=itetr(k,i)
          eti(k)=ebi(i1)
          etf1(k)=ebf(i1)
          del(k)=etf1(k)-eti(k)
          at1(k)=a(i1)*blwnorm
          call qqblohl(qt1(1,k),qb,idold(k,i))
          if(eti(k).le.ef)nile=nile+1
          if(eti(k).ge.ef)nige=nige+1
          if(etf1(k).le.ef)nfle=nfle+1
          if(del(k).ge.emax)nwge=nwge+1
          if(del(k).le.emin)nwle=nwle+1
        enddo
        if(nige.eq.4.or.nfle.eq.4)goto 500
        if(nwge.eq.4.or.nwle.eq.4)goto 500
c
        if(nile.eq.4)then               ! occupied tetrahedron
          call corect(eti,etf1,at1,qt1,eff,xi)
          goto 500
        endif
c ordering of et(i) coupled with at(i),delt(i) and qt(i,k)
        do k=1,4
          ip(k)=k
        enddo
        do k=1,3
          do l=1,4-k
            if(eti(l+1).lt.eti(l))then
              x=eti(l)
              eti(l)=eti(l+1)
              eti(l+1)=x
              ii=ip(l)
              ip(l)=ip(l+1)
              ip(l+1)=ii
            endif
          enddo
        enddo
        do k=1,4
          l=ip(k)
          at(k)=at1(l)
          etf(k)=etf1(l)
          qt(1,k)=qt1(1,l)
          qt(2,k)=qt1(2,l)
          qt(3,k)=qt1(3,l)
        enddo
c
        if(efi.le.eti(2)) then          ! full corner of the tetrahedron
c only initial state 1 is occupied
          do k=2,4
            x=(efi-eti(1))/(eti(k)-eti(1))
            eti(k)=efi
            etf(k)=etf(1)+x*(etf(k)-etf(1))
            at(k)=at(1)+x*(at(k)-at(1))
            qt(1,k)=qt(1,1)+x*(qt(1,k)-qt(1,1))
            qt(2,k)=qt(2,1)+x*(qt(2,k)-qt(2,1))
            qt(3,k)=qt(3,1)+x*(qt(3,k)-qt(3,1))
          enddo
          call corect(eti,etf,at,qt,eff,xi)
        elseif(efi.ge.eti(3)) then       ! empty corner of the tetrahedron
c only initial state 4 is unoccupied
          call corect(eti,etf,at,qt,eff,xi)
          do k=1,3
            x=(efi-eti(4))/(eti(k)-eti(4))
            eti(k)=efi
            etf(k)=etf(4)+x*(etf(k)-etf(4))
            at(k)=-(at(4)+x*(at(k)-at(4)))
            qt(1,k)=qt(1,4)+x*(qt(1,k)-qt(1,4))
            qt(2,k)=qt(2,4)+x*(qt(2,k)-qt(2,4))
            qt(3,k)=qt(3,4)+x*(qt(3,k)-qt(3,4))
          enddo
          at(4)=-at(4)
          call corect(eti,etf,at,qt,eff,xi)
        else
c initial states 1 and 2 are occupied
c empty edge of the tetrahedron
          do k=3,4
            x=(efi-eti(1))/(eti(k)-eti(1))
            etf1(k)=etf(1)+x*(etf(k)-etf(1))
            at1(k)=at(1)+x*(at(k)-at(1))
            qt1(1,k)=qt(1,1)+x*(qt(1,k)-qt(1,1))
            qt1(2,k)=qt(2,1)+x*(qt(2,k)-qt(2,1))
            qt1(3,k)=qt(3,1)+x*(qt(3,k)-qt(3,1))
            x=(efi-eti(2))/(eti(k)-eti(2))
            eti(k)=efi
            etf(k)=etf(2)+x*(etf(k)-etf(2))
            at(k)=at(2)+x*(at(k)-at(2))
            qt(1,k)=qt(1,2)+x*(qt(1,k)-qt(1,2))
            qt(2,k)=qt(2,2)+x*(qt(2,k)-qt(2,2))
            qt(3,k)=qt(3,2)+x*(qt(3,k)-qt(3,2))
          enddo
          call corect(eti,etf,at,qt,eff,xi)
c
          eti(2)=efi
          etf(2)=etf1(3)
          at(2)=at1(3)
          qt(1,2)=qt1(1,3)
          qt(2,2)=qt1(2,3)
          qt(3,2)=qt1(3,3)
          call corect(eti,etf,at,qt,eff,xi)
c
          eti(3)=efi
          etf(3)=etf1(4)
          at(3)=at1(4)
          qt(1,3)=qt1(1,4)
          qt(2,3)=qt1(2,4)
          qt(3,3)=qt1(3,4)
          call corect(eti,etf,at,qt,eff,xi)
        endif
c
 500  continue
      end
c
c---------------------------------------------------------------------
      subroutine corect(eti0,etf0,at0,qt0,ef,xi)
c
      implicit double precision (a-h,o-z)
c
      dimension xi(*),eti0(4),etf0(4),at0(4),qt0(3,4)
      dimension etf(4),del(4),at(4),qt(3,4),del1(4),at1(4),qt1(3,4)
     $     ,ip(4)
c
      nge=0
      nle=0
      delmax=-1.d10
      do k=1,4
        del1(k)=etf0(k)-eti0(k)
        if(del1(k).gt.delmax)delmax=del1(k)
        if(etf0(k).le.ef)nle=nle+1
        if(etf0(k).ge.ef)nge=nge+1
      enddo
      if(nle.eq.4)return                ! all final states are occupied
      if(delmax.lt.0.d0)return          ! negativ omega
      if(nge.eq.4)then                  ! all final states are emty
        call onetet(xi,del1,at0,qt0)
        return
      endif
c assignment of values for one tetrahedron
      do k=1,4
        etf(k)=etf0(k)
        ip(k)=k
      enddo
      do k=1,3
        do l=1,4-k
          if(etf(l+1).gt.etf(l))then
            x=etf(l)
            etf(l)=etf(l+1)
            etf(l+1)=x
            ii=ip(l)
            ip(l)=ip(l+1)
            ip(l+1)=ii
          endif
        enddo
      enddo
      do k=1,4
        l=ip(k)
        at(k)=at0(l)
        del(k)=del1(l)
        qt(1,k)=qt0(1,l)
        qt(2,k)=qt0(2,l)
        qt(3,k)=qt0(3,l)
      enddo
      if(ef.gt.etf(2))then          ! one corner is occupied
        do k=2,4
          x=(ef-etf(1))/(etf(k)-etf(1))
          del(k)=del(1)+x*(del(k)-del(1))
          at(k)=at(1)+x*(at(k)-at(1))
          qt(1,k)=qt(1,1)+x*(qt(1,k)-qt(1,1))
          qt(2,k)=qt(2,1)+x*(qt(2,k)-qt(2,1))
          qt(3,k)=qt(3,1)+x*(qt(3,k)-qt(3,1))
        enddo
        call onetet(xi,del,at,qt)
      elseif(ef.le.etf(3))then          ! one corner is empty
        call onetet(xi,del,at,qt)
        do k=1,3
          x=(ef-etf(4))/(etf(k)-etf(4))
          del(k)=del(4)+x*(del(k)-del(4))
          at(k)=-(at(4)+x*(at(k)-at(4)))
          qt(1,k)=qt(1,4)+x*(qt(1,k)-qt(1,4))
          qt(2,k)=qt(2,4)+x*(qt(2,k)-qt(2,4))
          qt(3,k)=qt(3,4)+x*(qt(3,k)-qt(3,4))
        enddo
        at(4)=-at(4)
        call onetet(xi,del,at,qt)
      else
        do k=3,4
          x=(ef-etf(1))/(etf(k)-etf(1))
          del1(k)=del(1)+x*(del(k)-del(1))
          at1(k)=at(1)+x*(at(k)-at(1))
          qt1(1,k)=qt(1,1)+x*(qt(1,k)-qt(1,1))
          qt1(2,k)=qt(2,1)+x*(qt(2,k)-qt(2,1))
          qt1(3,k)=qt(3,1)+x*(qt(3,k)-qt(3,1))
          x=(ef-etf(2))/(etf(k)-etf(2))
          del(k)=del(2)+x*(del(k)-del(2))
          at(k)=at(2)+x*(at(k)-at(2))
          qt(1,k)=qt(1,2)+x*(qt(1,k)-qt(1,2))
          qt(2,k)=qt(2,2)+x*(qt(2,k)-qt(2,2))
          qt(3,k)=qt(3,2)+x*(qt(3,k)-qt(3,2))
        enddo
        call onetet(xi,del,at,qt)
c
        del(2)=del1(3)
        at(2)=at1(3)
        qt(1,2)=qt1(1,3)
        qt(2,2)=qt1(2,3)
        qt(3,2)=qt1(3,3)
        call onetet(xi,del,at,qt)
c
        del(3)=del1(4)
        at(3)=at1(4)
        qt(1,3)=qt1(1,4)
        qt(2,3)=qt1(2,4)
        qt(3,3)=qt1(3,4)
        call onetet(xi,del,at,qt)
      endif
      end
c---------------------------------------------------------------------
      subroutine onetet(xi,del0,at,qt)
c
      implicit double precision (a-h,o-z)
      parameter(split=1.d-10)
c
      dimension xi(*),del0(4),at(4),qt(3,4)
      dimension del(4),a(4),q(9)
      common /cor/ emin,emax,de,npe
c
      eint(e,delx,ax,c1,s1)=c1*(ax+s1*(e-delx))*(e-delx)*(e-delx)
c
      do i=1,4
        del(i)=del0(i)
        a(i)=at(i)
      enddo
c$$$      write(20,'('' de'',4g18.11)')(del(k),k=1,4)
c$$$      write(20,'(''  a'',4g18.11)')(a(k),k=1,4)
c ordering  del(i) coupled with a(i)
      do k=1,3
        do i=1,4-k
          if(del(i+1).lt.del(i))then
            x=del(i)
            del(i)=del(i+1)
            del(i+1)=x
            x=a(i)
            a(i)=a(i+1)
            a(i+1)=x
          endif
        enddo
      enddo
      if(del(4).le.0.d0)return
      nmin=(del(1)-emin)/de+2.000001d0
      if(nmin.gt.npe)return
c volume of the tetrahedron
      k=0
      do i=1,3
        q(k+1)=qt(1,i)-qt(1,4)
        q(k+2)=qt(2,i)-qt(2,4)
        q(k+3)=qt(3,i)-qt(3,4)
        k=k+3
      enddo
      vol=q(1)*q(5)*q(9)+q(4)*q(8)*q(3)+q(7)*q(2)*q(6)
     $     -q(3)*q(5)*q(7)-q(6)*q(8)*q(1)-q(9)*q(2)*q(4)
      vol=abs(vol)/2.d0                 ! 3*volume of the tetrahedron
c
      if(nmin.lt.1)nmin=2 !1
      nmax=(del(2)-emin)/de+1.000001d0
      if(nmax.gt.npe)nmax=npe
      if(nmin.le.nmax)then
        c1=vol/((del(2)-del(1))*(del(3)-del(1))*(del(4)-del(1)))
        s1=((a(2)-a(1))/(del(2)-del(1))
     $     +(a(3)-a(1))/(del(3)-del(1))
     $     +(a(4)-a(1))/(del(4)-del(1)))/3.d0
        do k=nmin,nmax
          e=emin+(k-1)*de
          xi(k)=xi(k)+eint(e,del(1),a(1),c1,s1)
        enddo
        nmin=nmax+1
      endif
c
      nmax=(del(3)-emin)/de+1.000001d0
      if(nmax.gt.npe)nmax=npe
      if(nmin.le.nmax)then
        if(del(2)-del(1).gt.split) then
          c1= vol/((del(2)-del(1))*(del(3)-del(1))*(del(4)-del(1)))
          s1=((a(2)-a(1))/(del(2)-del(1))
     $       +(a(3)-a(1))/(del(3)-del(1))
     $       +(a(4)-a(1))/(del(4)-del(1)))/3.d0
          c2=-vol/((del(1)-del(2))*(del(3)-del(2))*(del(4)-del(2)))
          s2=((a(1)-a(2))/(del(1)-del(2))
     $       +(a(3)-a(2))/(del(3)-del(2))
     $       +(a(4)-a(2))/(del(4)-del(2)))/3.d0
          do k=nmin,nmax
            e=emin+(k-1)*de
            xi(k)=xi(k)+eint(e,del(1),a(1),c1,s1)
     $               -eint(e,del(2),a(2),c2,s2)
          enddo
        elseif(del(4)-del(3).gt.split) then
          c4=-vol/((del(1)-del(4))*(del(2)-del(4))*(del(3)-del(4)))
          s4=((a(1)-a(4))/(del(1)-del(4))
     $       +(a(2)-a(4))/(del(2)-del(4))
     $       +(a(3)-a(4))/(del(3)-del(4)))/3.d0
          c3= vol/((del(1)-del(3))*(del(2)-del(3))*(del(4)-del(3)))
          s3=((a(1)-a(3))/(del(1)-del(3))
     $       +(a(2)-a(3))/(del(2)-del(3))
     $       +(a(4)-a(3))/(del(4)-del(3)))/3.d0
          do k=nmin,nmax
            e=emin+(k-1)*de
            xi(k)=xi(k)+eint(e,del(4),a(4),c4,s4)
     $               -eint(e,del(3),a(3),c3,s3)
          enddo
        else                          ! two pairs of equal del
c splitting  del(i)
          del(1)=del(2)-split
          c1= vol/((del(2)-del(1))*(del(3)-del(1))*(del(4)-del(1)))
          s1=((a(2)-a(1))/(del(2)-del(1))
     $       +(a(3)-a(1))/(del(3)-del(1))
     $       +(a(4)-a(1))/(del(4)-del(1)))/3.d0
          c2=-vol/((del(1)-del(2))*(del(3)-del(2))*(del(4)-del(2)))
          s2=((a(1)-a(2))/(del(1)-del(2))
     $       +(a(3)-a(2))/(del(3)-del(2))
     $       +(a(4)-a(2))/(del(4)-del(2)))/3.d0
          do k=nmin,nmax
            e=emin+(k-1)*de
            xi(k)=xi(k)+eint(e,del(1),a(1),c1,s1)
     $               -eint(e,del(2),a(2),c2,s2)
          enddo
        endif
        nmin=nmax+1
      endif
c
      nmax=(del(4)-emin)/de+1.000001d0
      if(nmax.gt.npe)nmax=npe
      if(nmin.le.nmax)then
        c4=-vol/((del(1)-del(4))*(del(2)-del(4))*(del(3)-del(4)))
        s4=((a(1)-a(4))/(del(1)-del(4))
     $     +(a(2)-a(4))/(del(2)-del(4))
     $     +(a(3)-a(4))/(del(3)-del(4)))/3.d0
        do k=nmin,nmax
          e=emin+(k-1)*de
          xi(k)=xi(k)+eint(e,del(4),a(4),c4,s4)
        enddo
      endif
      end
