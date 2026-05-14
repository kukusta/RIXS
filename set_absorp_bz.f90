SUBROUTINE SET_ABSORP_BZ(ia,kpk)
use m_aux, only: lmtindex
use m_bnd
use m_bz
use m_bzmesh, only: ndxyz,nkbz,ig4q,nkibz
use m_functions
use m_rixs
implicit none
! input
integer, intent(in) :: ia
! output
real(REAL64) kpk(nb,nkbz)
! local vars
complex(REAL64), parameter :: z0=(0._REAL64,0._REAL64)
complex(REAL64) mmek(3),I(3,3)
integer ib0,kbz,j3,j2,j1,kibz,igk,ib,iat,iat1,iatn

if(ia < 1 .or. ia > n_mme) stop ' In SET_ABSORP_BZ: wrong ia value.'

! write(*,*)'ia=',ia,' nkibz=',nkibz
! do j3=1,ndxyz(3)
!   do j2=1,ndxyz(2)
!     write(*,*)'newline'
!     write(*,*)ipqibz(:,j2,j3)
!   enddo
! enddo

if(.not.allocated(lmtdata))then
  write(*,*)'SET_ABSORP_BZ: lmtdata is not allocated, ia=',ia
  return
endif
if(.not.allocated(lmtindex))then
  write(*,*)'SET_ABSORP_BZ: lmtindex is not allocated, ia=',ia
  return
endif
if(.not.allocated(ipqibz))then
  write(*,*)'SET_ABSORP_BZ: ipqibz is not allocated, ia=',ia
  return
endif
if(.not.allocated(ik2ibz))then
  write(*,*)'SET_ABSORP_BZ: ik2ibz is not allocated, ia=',ia
  return
endif
if(.not.allocated(ig4q))then
  write(*,*)'SET_ABSORP_BZ: ig4q is not allocated, ia=',ia
  return
endif
if(.not.allocated(lmtdata(lmtindex(ia))%mme))then
  write(*,*)'SET_ABSORP_BZ: lmtdata(lmtindex(ia))%mme is not allocated, ia=',ia
endif

iat1=lmtdata(lmtindex(ia))%iat1
iatn=lmtdata(lmtindex(ia))%iatn
! write(*,*)'iat1,iatn',iat1,iatn
ib0=0
kbz=0
do j3=1,ndxyz(3)
  do j2=1,ndxyz(2)
    do j1=1,ndxyz(1)
      kbz=kbz+1
! here kibz can be referenced as kibz=ik2ibz(kbz), but I use ipq to:
! 1. double check data in ipqibz and ik2ibz
! 2. preserve the same code style for kibz and kqibz
      kibz=ipqibz(j1,j2,j3)
      if(kibz <= 0 .or. kibz > nkibz)stop' In SET_ABSORP_BZ: wrong kibz.'
      if(kibz /= ik2ibz(kbz))stop ' In SET_ABSORP_BZ: kibz /= ik2ibz(kbz).'
      igk=ig4q(j1,j2,j3)
      if(igk <= 0 .or. igk > nopused)stop' In SET_ABSORP_BZ: wrong igk.'
! kibz is transformed into kbz using g=iopnum(ig): kbz=matmul(g,kibz)
      ! if(igk > 1) opk(1:3,1:3)=transpose(g(1:3,1:3,igk))
      do ib=1,nb
        do iat=iat1,iatn
          mmek(:)=lmtdata(lmtindex(ia))%mme(:,iat,ia,ib,kibz)
          if(igk /= 1) mmek(:)=mme_trans(igk,mmek)
          call tensor_from_mme(mmek,I)
          kpk(ib0+ib,kbz)=real(I(1,1)+I(2,2),kind=REAL64)
        enddo         ! iat
      enddo           ! ib
    enddo             ! j1
  enddo               ! j2
enddo                 ! j3

kpk(:,:)=kpk(:,:)/cmplx(8*(iatn-iat1+1),kind=REAL64)

END SUBROUTINE SET_ABSORP_BZ

! OLD CODE WAS SUPPOSED TO INTEGRATE OVER BZ USING k+q
! INSTEAD OF k. SO NOW I REMOVED IT.

! SUBROUTINE SET_ABSORP_BZ(ia,kpk,kpkq)
! use m_aux, only: ispec,i1,i2,i3
! use m_bnd
! use m_bz
! use m_bzmesh, only: ndxyz,nkbz,ig4q
! use m_rixs
! implicit none
! ! input
! integer, intent(in) :: ia
! ! output
! real(REAL64) kpk(nb,nkbz)
! real(REAL64), optional :: kpkq(nb,nkbz)
! ! local vars
! complex(REAL64), parameter :: z0=(0._REAL64,0._REAL64)
! complex(REAL64) mmek(3),mmekq(3),I(3,3)
! integer ib0,kbz,j3,j2,j1,kibz,igk,kqbz,kqibz,ib,igkq,iat,iat1,iatn
! real(REAL64) opk(3,3),opkq(3,3)

! if(ia < 1 .or. ia > na) stop' In SET_ABSORP_BZ: wrong ia value.'
! iat1=rixs(ispec(ia))%iat1
! iatn=rixs(ispec(ia))%iatn
! ib0=0
! kbz=0
! do j3=1,ndxyz(3)
!   do j2=1,ndxyz(2)
!     do j1=1,ndxyz(1)
!       kbz=kbz+1
! ! here kibz can be referenced as kibz=ik2ibz(kbz), but I use ipq to:
! ! 1. double check data in ipqibz and ik2ibz
! ! 2. preserve the same code style for kibz and kqibz
!       kibz=ipqibz(j1,j2,j3)
!       if(kibz /= ik2ibz(kbz))stop ' In SET_ABSORP_BZ: kibz /= ik2ibz(kbz).'
!       igk=ig4q(j1,j2,j3)
! ! kibz is transformed into kbz using g=iopnum(ig): kbz=matmul(g,kibz)
!       if(igk > 1) opk(1:3,1:3)=transpose(g(1:3,1:3,igk))
!       if(present(kpkq))then
!         kqbz=ik2kq(kbz,ispec(ia))
!         kqibz=ipqibz(i1(kqbz),i2(kqbz),i3(kqbz))

!         ! if(kqibz /= ik2kqibz(kqbz,ispec(ia)) .and. ia == 1) &
!         !   write(*,*)'kqibz,ik2kqibz',kqibz,ik2kqibz(kqbz,ispec(ia))
!         if(kqibz /= ik2kqibz(kqbz,ispec(ia))) stop ' kqibz /= ik2kqibz'

!         kqibz=ik2kqibz(kqbz,ispec(ia))
!         igkq=ig4q(i1(kqbz),i2(kqbz),i3(kqbz))
!         if(igkq > 1) opkq(1:3,1:3)=transpose(g(1:3,1:3,igkq))
!       endif
!       do ib=1,nb
!         do iat=iat1,iatn
! ! Absorption using BZ and ik2ibz
!           mmek(:)=rixs(ispec(ia))%mme(:,iat,ia,ib,kibz)
!           if(igk > 1) mmek(:)=matmul(opk,mmek)
!           call tensor_from_mme(mmek,I)
!           kpk(ib0+ib,kbz)=real(I(1,1)+I(2,2),kind=REAL64)
! ! Absorption using BZ and ik2kqibz
!           if(present(kpkq))then
!             mmekq(:)=rixs(ispec(ia))%mme(:,iat,ia,ib,kqibz)
! !                if(igkq > 1) call mme_op(g(:,:,igkq),mme,mme)
!             if(igkq > 1)mmekq(:)=matmul(opkq,mmekq)
!             call tensor_from_mme(mmekq,I)
!             kpkq(ib0+ib,kqbz)=real(I(1,1)+I(2,2),kind=REAL64)
!           endif
!         enddo         ! iat
!       enddo           ! ib
!     enddo             ! j1
!   enddo               ! j2
! enddo                 ! j3
! kpk(:,:)=kpk(:,:)/cmplx(8*(iatn-iat1+1),kind=REAL64)
! if(present(kpkq)) kpkq(:,:)=kpkq(:,:)/cmplx(8*(iatn-iat1+1),kind=REAL64)

! END SUBROUTINE SET_ABSORP_BZ

