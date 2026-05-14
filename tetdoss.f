      subroutine tetdos(nemax,lmax,lsp,limzn,limallzn
     $     ,emin,de,dos,totdos
     $     ,eigen,almsl,nbmin,nbmax)
c$$$      subroutine tetdos(nemax,lmax,lsp,limzn,limallzn
c$$$     $     ,emin,emax,de,tos,dos,totdos
c$$$     $     ,eigen,almsl,nband,itos,nbmin,nbmax)
************************************************
*                                              *
*  Calculations of l- and sort- projected DOS  *
*                                              *
************************************************
c
c calculation of the number of states with itos=1 does not work.
c
      use m_bzmesh
      INCLUDE 'PREC.FI'

      dimension eigen(limallzn,*),almsl(lmax+1,limzn,*) ! input
      dimension dos(nemax,lsp),totdos(nemax) ! output
      dimension a(4,lsp),e(4)           ! local

c$$$      if (ndimspin.eq.2.or.iall_bands.eq.1) then
c$$$        omg48s=omg48/2.d0
c$$$      else
c$$$        omg48s=omg48
c$$$      endif
      volomg=vol1t*omg48s
      epsde=de*1.d-2
      eps=epsilon(epsde)
c
      totdos(:)=0.d0
      dos(:,:)=0.d0

!$omp parallel default(shared) private( ot , ib , i , ipnt , e , l , a )
!$omp do schedule(dynamic,1) reduction(+:totdos,dos)
      do lk=1,ntibz                     !ntet(1)
        ot=itetr(0,lk)*volomg
        do ib=nbmin,nbmax
          do i=1,4
            ipnt=itetr(i,lk)
            e(i)=eigen(ib,ipnt)
            do l=1,lsp
              a(i,l)=almsl(l,ib,ipnt)
            enddo
          enddo
c$$$          write(*,*)'ib',ib,e(1:4)
c$$$          if(itos.eq.0)then
          call totos(ot,e,a,lsp,nemax,emin,de,dos,totdos,eps,epsde)
c$$$          else
c$$$            call tdostos(ot,e,a,lsp,nemax,emin,de,tos,dos,totdos
c$$$     $           ,0,eps,epsde)
c$$$          endif
        enddo                           ! ib
      enddo                             ! lk
!$omp end do
!$omp end parallel
c too small dos produces an error in orx
      do i=1,nemax
        if(abs(totdos(i)).lt.1.d-15)totdos(i)=0.d0
      enddo
      do l=1,lsp
        do i=1,nemax
          if(abs(dos(i,l)).lt.1.d-15)dos(i,l)=0.d0
        enddo
      enddo
      end

! old tetdos
      subroutine tetdostos(nemax,lmax,lsp,limzn,limallzn
     $     ,emin,de,tos,dos,totdos
     $     ,eigen,almsl,itos,nbmin,nbmax)
c$$$      subroutine tetdos(nemax,lmax,lsp,limzn,limallzn
c$$$     $     ,emin,emax,de,tos,dos,totdos
c$$$     $     ,eigen,almsl,nband,itos,nbmin,nbmax)

!  Calculations of l- and sort- projected DOS and NOS(uncalled)

      use m_bzmesh
      INCLUDE 'PREC.FI'

      dimension eigen(limallzn,*),almsl(lmax+1,limzn,*) ! input
      dimension tos(nemax,lsp)          ! is not used (itos=0)
      dimension dos(nemax,lsp),totdos(nemax) ! output
      dimension a(4,lsp),e(4)           ! local

c$$$      if (ndimspin.eq.2.or.iall_bands.eq.1) then
c$$$        omg48s=omg48/2.d0
c$$$      else
c$$$        omg48s=omg48
c$$$      endif
      volomg=vol1t*omg48s
      epsde=de*1.d-2
      eps=epsilon(epsde)
c
      totdos(:)=0.d0
      dos(:,:)=0.d0
      do lk=1,ntibz                     !ntet(1)
        ot=itetr(0,lk)*volomg
        do ib=nbmin,nbmax
          do i=1,4
            ipnt=itetr(i,lk)
            e(i)=eigen(ib,ipnt)
            do l=1,lsp
              a(i,l)=almsl(l,ib,ipnt)
            enddo
          enddo
c$$$          write(*,*)'ib',ib,e(1:4)
          if(itos.eq.0)then
            call totos(ot,e,a,lsp,nemax,emin,de,dos,totdos
     $           ,eps,epsde)
          else
            call tdostos(ot,e,a,lsp,nemax,emin,de,tos,dos,totdos
     $           ,0,eps,epsde)
          endif
        enddo                           ! ib
      enddo                             ! lk
c too small dos produces an error in orx
      do i=1,nemax
        if(abs(totdos(i)).lt.1.d-15)totdos(i)=0.d0
      enddo
      do l=1,lsp
        do i=1,nemax
          if(abs(dos(i,l)).lt.1.d-15)dos(i,l)=0.d0
        enddo
      enddo
      end

      subroutine totos(ot,e0,a,lsp,nemax,emin0,de,dos,totdos,eps,epsde)

! Integration of DOS over a microtetrahedron
! Only dos is calculated. Calculation of tos is removed

      INCLUDE 'PREC.FI'
c
      parameter (d13=1.d0/3.d0,d14=0.25d0)

      real(8),intent(in) :: e0(4),a(4,*) ! input
      integer,intent(in) :: lsp,nemax
      real(8),intent(in) :: emin0,de,eps,epsde
      real(8),intent(out) :: dos(nemax,*),totdos(*) ! output
      dimension et(4)                   ! local
      dimension ind(4),j(4)
c
      do i=1,4
        ind(i)=1
      enddo
      do i=1,3
        do k=i+1,4
          if(e0(i).le.e0(k))then
            ind(i)=ind(i)+1
          else
            ind(k)=ind(k)+1
          endif
        enddo
      enddo
      do i=1,4
        j(ind(i))=i
        et(ind(i))=e0(i)
      enddo
c descending order: et(4) < et(3) < et(2) < et(1)
c there may be points between e(4) <emin0< e(1)
      if(et(1).lt.emin0)return

      iemin=(et(4)-emin0-epsde)/de+2    ! e(iemin)>=et(4) 
      iemin=max(iemin,1)
      ot3=ot*3.d0
      iemax=(et(1)-emin0)/de+1          !  e(iemax)<et(1) 
      iemax=min(iemax,nemax)
      do ie=iemin,iemax
        ei=emin0+(ie-1)*de
        if(ei.lt.et(4)+eps)then
          ds=0.d0
        elseif(ei.le.et(3)+eps)then
          ei4=ei-et(4)
          e14=ei4/(et(1)-et(4))
          e24=ei4/(et(2)-et(4))
          e34=ei4/(et(3)-et(4))
          ds=ot3*e24*e14*e34/ei4
          do l=1,lsp
            a14=e14*(a(j(1),l)-a(j(4),l))
            a24=e24*(a(j(2),l)-a(j(4),l))
            a34=e34*(a(j(3),l)-a(j(4),l))
            dos(ie,l)=dos(ie,l)+ds*(a(j(4),l)+(a24+a14+a34)*d13)
          enddo
        elseif(ei.lt.et(2)-eps)then     ! section is a quadrangle
          e23=et(2)-et(3)
c dividing into two tetrahedra
          ei4=ei-et(4)                  ! the first tetrahedron
          e14=et(1)-et(4)
          e24=et(2)-et(4)
          f1=ot3*ei4*(et(2)-ei)/e23/e14/e24
          ei1=et(1)-ei                  ! the second tetrahedron
          e13=et(1)-et(3)
          f2=ot3*ei1*(ei-et(3))/e23/e13/e14
          ei3=(ei-et(3))/e23
          ds=f1+f2
          do l=1,lsp
            a13=(a(j(1),l)-a(j(3),l))/e13
            a14=(a(j(1),l)-a(j(4),l))/e14
            a24=(a(j(2),l)-a(j(4),l))/e24
            a3=a(j(3),l)+(a(j(2),l)-a(j(3),l))*ei3
            ds1=a(j(4),l)+(ei4*(a14+a24)+a3-a(j(4),l))*d13
            ds2=a(j(1),l)-(a(j(1),l)-a3+ei1*(a13+a14))*d13
            dos(ie,l)=dos(ie,l)+f1*ds1+f2*ds2
          enddo
        elseif(ei.lt.et(1)-eps)then
          ei1=ei-et(1)                  ! ei1<0
          e14=ei1/(et(1)-et(4))
          e12=ei1/(et(1)-et(2))
          e13=ei1/(et(1)-et(3))
          ds=ot3*e12*e13*e14/ei1
          do l=1,lsp
            a14=e14*(a(j(1),l)-a(j(4),l))
            a13=e13*(a(j(1),l)-a(j(3),l))
            a12=e12*(a(j(1),l)-a(j(2),l))
            dos(ie,l)=dos(ie,l)+
     $           ds*(a(j(1),l)+(a12+a14+a13)*d13)
          enddo
        else
          ds=0.d0
        endif
        totdos(ie)=totdos(ie)+ds
      enddo                             ! ie
      end

      subroutine tdostos(ot,e0,a,lsp,nemax,emin0,de,tos,dos,totdos
     $     ,itos,eps,epsde)

! Integration of DOS and number of states (tos) over a microtetrahedron

      INCLUDE 'PREC.FI'
c
      parameter (d13=1.d0/3.d0,d14=0.25d0)
c
      dimension e0(4),et(4),a(4,*)
      dimension tos(nemax,*),dos(nemax,*),totdos(nemax)
      dimension ind(4),j(4)
c
      do i=1,4
        ind(i)=1
      enddo
      do i=1,3
        do k=i+1,4
          if(e0(i).le.e0(k))then
            ind(i)=ind(i)+1
          else
            ind(k)=ind(k)+1
          endif
        enddo
      enddo
      do i=1,4
        j(ind(i))=i
        et(ind(i))=e0(i)
      enddo
c descending order: et(4) < et(3) < et(2) < et(1)
c$$$      if(et(4).ge.emin0)then
c there may be points between e(4) <emin0< e(1)
      if(et(1).ge.emin0)then
        iemin=(et(4)-emin0-epsde)/de+2  ! e(iemin)>=et(4) 
        iemin=max(iemin,1)
        if(itos.eq.0)then
          ot3=ot*3.d0
          iemax=(et(1)-emin0)/de+1      !  e(iemax)<et(1) 
          iemax=min(iemax,nemax)
          do ie=iemin,iemax
            ei=emin0+(ie-1)*de
            if(ei.lt.et(4)+eps)then
              ds=0.d0
            elseif(ei.le.et(3)+eps)then
              ei4=ei-et(4)
              e14=ei4/(et(1)-et(4))
              e24=ei4/(et(2)-et(4))
              e34=ei4/(et(3)-et(4))
              ds=ot3*e24*e14*e34/ei4
              do l=1,lsp
                a14=e14*(a(j(1),l)-a(j(4),l))
                a24=e24*(a(j(2),l)-a(j(4),l))
                a34=e34*(a(j(3),l)-a(j(4),l))
                dos(ie,l)=dos(ie,l)+ds*(a(j(4),l)+(a24+a14+a34)*d13)
              enddo
            elseif(ei.lt.et(2)-eps)then ! section is a quadrangle
              e23=et(2)-et(3)
c dividing into two tetrahedra
              ei4=ei-et(4)              ! the first tetrahedron
              e14=et(1)-et(4)
              e24=et(2)-et(4)
              f1=ot3*ei4*(et(2)-ei)/e23/e14/e24
              ei1=et(1)-ei              ! the second tetrahedron
              e13=et(1)-et(3)
              f2=ot3*ei1*(ei-et(3))/e23/e13/e14
              ei3=(ei-et(3))/e23
              ds=f1+f2
              do l=1,lsp
                a13=(a(j(1),l)-a(j(3),l))/e13
                a14=(a(j(1),l)-a(j(4),l))/e14
                a24=(a(j(2),l)-a(j(4),l))/e24
                a3=a(j(3),l)+(a(j(2),l)-a(j(3),l))*ei3
                ds1=a(j(4),l)+(ei4*(a14+a24)+a3-a(j(4),l))*d13
                ds2=a(j(1),l)-(a(j(1),l)-a3+ei1*(a13+a14))*d13
                dos(ie,l)=dos(ie,l)+f1*ds1+f2*ds2
              enddo
            elseif(ei.lt.et(1)-eps)then
              ei1=ei-et(1)              ! ei1<0
              e14=ei1/(et(1)-et(4))
              e12=ei1/(et(1)-et(2))
              e13=ei1/(et(1)-et(3))
              ds=ot3*e12*e13*e14/ei1
              do l=1,lsp
                a14=e14*(a(j(1),l)-a(j(4),l))
                a13=e13*(a(j(1),l)-a(j(3),l))
                a12=e12*(a(j(1),l)-a(j(2),l))
                dos(ie,l)=dos(ie,l)+
     $               ds*(a(j(1),l)+(a12+a14+a13)*d13)
              enddo
            else
              ds=0.d0
            endif
            totdos(ie)=totdos(ie)+ds
          enddo
        else
c$$$          do ie=1,nemax
          do ie=iemin,nemax
            ei=emin0+(ie-1)*de
            if(ei.le.et(4))then
              ts=0.d0
              ds=0.d0
            elseif(ei.le.et(3))then
              ei4=ei-et(4)
              e24=ei4/(et(2)-et(4))
              e14=ei4/(et(1)-et(4))
              e34=ei4/(et(3)-et(4))
              ts=e24*e14*e34*ot
              ds=ts*3.d0/ei4
              do l=1,lsp
                a14=e14*(a(j(1),l)-a(j(4),l))
                a34=e34*(a(j(3),l)-a(j(4),l))
                a24=e24*(a(j(2),l)-a(j(4),l))
                dos(ie,l)=dos(ie,l)+
     $               ds*(a(j(4),l)+(a24+a14+a34)*d13)
                tos(ie,l)=tos(ie,l)+
     $               ts*(a(j(4),l)+(a24+a14+a34)*d14)
              enddo
            elseif(ei.le.et(2))then     ! section is a quadrangle
              om23=et(2)-et(3)
c dividing into two tetrahedra
              dom1=ei-et(4)             ! the first tetrahedron
              dom11=et(1)-et(4)
              dom12=et(2)-et(4)
              dom13=dom1
              f1=1.d0/dom11/dom12*(et(2)-ei)/om23
              dom2=-(ei-et(1))          ! the second tetrahedron
              dom21=dom2
              dom22=et(1)-et(3)
              dom23=et(1)-et(4)
              f2=1.d0/dom22/dom23
              g2=(ei-et(3))/om23
c
              ds=(dom1*f1+dom2*f2*g2)*3.d0*ot
              ts=(dom1**2*f1+(1.d0-dom2**2*f2)*g2)*ot
              do l=1,lsp
                a4=a(j(4),l)
                a1=a(j(1),l)
                pmdl=a(j(3),l)+(a(j(2),l)-
     $               a(j(3),l))*(ei-et(3))/(et(2)-et(3))
                pi11=(a1-a4)/dom11
                pi12=(a(j(2),l)-a4)/dom12
                pi13=(pmdl-a4)/dom13
                pids1=(pi11+pi12+pi13)*d13
                pits1=(pi11+pi12+pi13)*d14
                pi21=(a1-pmdl)/dom21
                pi22=(a1-a(j(3),l))/dom22
                pi23=(a1-a4)/dom23
                pids2=-(pi21+pi22+pi23)*d13
                pits2=-(pi21+pi22+pi23)*d14
                dos(ie,l)=dos(ie,l)+(dom1*f1*(a4+dom1*pids1)
     $               +dom2*f2*g2*(a1+dom2*pids2))*3.d0*ot
                tos(ie,l)=tos(ie,l)+(dom1**2*f1*(a4+dom1*pits1)
     $               +((pmdl+a4+a1+a(j(3),l))*d14
     $               -dom2**2*f2*(a1+dom2*pits2))*g2)*ot
              enddo
            elseif(ei.lt.et(1))then
              ei1=ei-et(1)              ! ei1<0
              e14=ei1/(et(1)-et(4))
              e12=ei1/(et(1)-et(2))
              e13=ei1/(et(1)-et(3))
              ts=e12*e13*e14*ot
              ds=ts*3.d0/ei1
              do l=1,lsp
                a14=e14*(a(j(1),l)-a(j(4),l))
                a13=e13*(a(j(1),l)-a(j(3),l))
                a12=e12*(a(j(1),l)-a(j(2),l))
                dos(ie,l)=dos(ie,l)+ds*(a(j(1),l)+(a12+a14+a13)*d13)
                tos(ie,l)=tos(ie,l)+
     $               (a(1,l)+a(2,l)+a(3,l)+a(4,l))*d14*ot
     $               +ts*(a(j(1),l)+(a12+a14+a13)*d14)
              enddo
              ts=ot+ts
            else
              ds=0.d0
              ts=ot
              do l=1,lsp
                tos(ie,l)=tos(ie,l)+(a(1,l)+a(2,l)+a(3,l)+a(4,l))*d14*ot
              enddo
            endif
            totdos(ie)=totdos(ie)+ds
          enddo
        endif
      endif
      end
