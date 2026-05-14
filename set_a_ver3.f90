SUBROUTINE SET_A(itr,a)
use m_aux
use m_bnd
use m_bz
use m_bzmesh, only : ndxyz
use m_functions
use m_params
use m_rixs
implicit none
integer :: itr
real a(rixs(ispec(itr))%nkrbz)
! local vars
character(9), parameter :: srcname=' in SET_A'
complex(REAL64), parameter :: z0 = ( 0._REAL64 , 0._REAL64 )
complex(REAL64), parameter :: zi = ( 0._REAL64 , 1._REAL64 )
complex(REAL64) dpk, dpkq, denom, za, mme_kibz(3), mme_kqibz(3), &
                mme_krbz(3), mme_kqrbz(3), dza
integer krbz, kqibz, kibz, jat, jat1, jatn, i, igk2ibz, igkq2ibz, &
        igk2rbz, igkq2rbz, kbz, kqbz , k1 , k2 , k3
real akrbz
real(REAL64), parameter :: eps = 1.e-6_REAL64
real(REAL64) phase, q(3), rat(3), dpi, op(3,3)

! set all elements to -1. to check if all of them were calculated
a(:) = -1.
q(:) = rixs( ispec(itr) )%q(:)
jat1 = lmtdata( lmtindex( ia(itr) ) )%iat1
jatn = lmtdata( lmtindex( ia(itr) ) )%iatn

kbz = 0
do k3 = 1 , ndxyz(3)
  do k2 = 1 , ndxyz(2)
    do k1 = 1 , ndxyz(1)
      kbz = kbz + 1
      krbz = 0 ! to check
      if ( nopused == 1 ) then
        kibz = kbz
        krbz = kbz
        kqbz = rixs( ispec(itr) )%ik2kq(kbz)
        kqibz = kbz
        igk2ibz = 1
        igkq2ibz = 1
      else
        kibz = ipqibz( k1 , k2 , k3 ) !ik2ibz(kbz)
        kqbz = rixs( ispec(itr) )%ik2kq(kbz)
        kqibz = ipqibz( i1(kqbz) , i2(kqbz) , i3(kqbz) ) !ik2ibz(kqbz)
        igk2ibz = ig4qibz( k1 , k2 , k3 )
        igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
        if ( rixs( ispec(itr) )%ng == 1 ) then
          krbz = kbz
        else if ( rixs( ispec(itr) )%ng == nopused ) then
          krbz = kibz
          ! if ( rixs( ispec(itr) )%ik2rbz(kbz) /= kibz ) &
          !    stop '   why ik2rbz(kbz) /= kibz?'
        else if ( rixs( ispec(itr) )%ig4q( k1 , k2 , k3 ) == 1 ) then
            krbz = rixs( ispec(itr) )%ik2rbz(kbz)
        else
          cycle
        endif
      endif

      za = z0
      do jat = jat1, jatn
        mme_kibz(:) = lmtdata( lmtindex( ia(itr) ) )%mme( :, jat, ia(itr), &
                      ibi(itr), kibz )
        mme_kqibz(:) = lmtdata( lmtindex( ia(itr) ) )%mme( :, jat, ia(itr), &
                       ibf(itr), kqibz )

      if ( igk2ibz == 1 ) then
        mme_krbz(:) = mme_kibz(:)
      else
        op( : , : ) = g( : , : , igk2ibz )
        if ( iopnum(igk2ibz) > 0 ) then
          mme_krbz(:) = matmul( op , mme_kibz )
        else
          mme_krbz(:) = matmul( op , conjg( mme_kibz ) )
        endif
      endif
      if ( igkq2ibz == 1 ) then
        mme_kqrbz(:) = mme_kqibz(:)
      else
        op( : , : ) = g( : , : , igkq2ibz )
        if ( iopnum(igkq2ibz) > 0 ) then
          mme_kqrbz(:) = matmul( op , mme_kqibz )
        else
          mme_kqrbz(:) = matmul( op , conjg( mme_kqibz ) )
        endif
      endif

      denom = rixs( ispec(itr) )%en + &
              lmtdata( lmtindex( ia(itr) ) )%esp( ia(itr) ) - &
              e( ibf(itr) , kibz ) + zi * rixs( ispec(itr) )%Gamma
      rat(:) = lmtdata( lmtindex( ia(itr) ) )%rat( : , jat )
      phase = 2._REAL64 * dpi() * dot_product( q, rat )
      select case ( abs(irixs) )
        case( 1 )
          dpk = dot_product( mme_krbz, rixs( ispec(itr) )%e_in )
          dpkq = dot_product( mme_kqrbz, rixs( ispec(itr) )%e_out )
          za = za + exp( -zi * phase ) * conjg(dpkq) * dpk / denom
        case ( 2 )
          dza = z0
          do i = 1, 3
            dza = dza + conjg( mme_kqrbz(i) ) * mme_krbz(i) / denom
          enddo
          za = za + exp( -zi * phase ) * dza
        case ( 3 )
          za = za + 1._REAL64 / denom
        case ( 4 )
          dpk = dot_product( mme_krbz, rixs( ispec(itr) )%e_in )
          dpkq = dot_product( mme_kqrbz, rixs( ispec(itr) )%e_out )
          za = za + exp( -zi * phase ) * conjg( dpkq ) * dpk
        case default
      end select
    enddo          ! jat
    akrbz = real( abs(za)**2 )
    if ( a(krbz) > -0.5 ) then
      if ( abs( a(krbz) - akrbz ) > eps ) then
!$omp critical
        write( unit = output_unit , fmt = "( 4x , a , i0 , a ,2f15.10)" ) &
          'For krbz=' , krbz , ' a(krbz) .ne. real(abs(za)**2): ' , &
          a(krbz) , akrbz
        call exit_on_error( srcname )
!$omp end critical
      endif
    else 
      a(krbz) = akrbz
    endif

    enddo             ! k3
  enddo               ! k2
enddo                 ! k1

if(abs(irixs) == 2) a(:)=a(:)/3.

!$$$ do krbz=1,rixs(ispec(itr))%nkrbz
!$$$   if(a(krbz) < 0.)write(unit=output_unit,fmt=55)krbz
!$$$ enddo
!$$$ if ( minval( a ) < 0. ) then
!$$$   write( unit = output_unit , fmt = 60 )
!$$$   call exit_on_error( srcname )
!$$$ endif

 15 format(4x,"dza=z0 for transition with ibi= ",i3,", ibf= ",i3&
          ,", ia= ",i2,", iat= ",i3,", k= ",i0)
 20 format(4x,"Contributions from different atoms:")
 25 format(4x,"itr=",i8," iat=",i3," phase/2/pi=",f12.6&
          ," exp(phase)= (",f12.6,",",f12.6,")"&
          ," conjg(dpkq)=(",f12.6,",",f12.6,")"&
          ," dpk/denom=(",f12.6,",",f12.6,")")
 55 format(4x,"ERROR: mme < 0 for krbz=",i0)
 60 format(4x,"ERROR: a(krbz) < 0 for some point in RBZ")

rixs(ispec(itr))%max_mme=max(rixs(ispec(itr))%max_mme,maxval(a))
rixs(ispec(itr))%min_mme=min(rixs(ispec(itr))%min_mme,minval(a))
rixs(ispec(itr))%av_mme =    rixs(ispec(itr))%av_mme +   sum(a)

END SUBROUTINE SET_A

! do krbz = 1 , rixs( ispec(itr) )%nkrbz
! !$$$do kbz = 1 , nkbz

!   ! kbz = rixs( ispec(itr) )%ikbz2rbz(krbz)
!   ! kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!   ! kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
!   ! kqibz = ipqibz( i1(kqbz), i2(kqbz), i3(kqbz) )
!   ! igkibz = ig4qibz( i1(kbz), i2(kbz), i3(kbz) )
!   ! igkqibz = ig4qibz( i1(kqbz), i2(kqbz), i3(kqbz) )

! !   if ( rixs( ispec(itr) )%ng == 1 ) then
! ! ! no symmetry
! !     krbz = kbz
! !     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
! !     ! kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
! !     ! kqibz = ipqibz( i1(kqbz), i2(kqbz), i3(kqbz) )
! !     kibz = ik2ibz( kbz )
! !     kqibz = ik2ibz( kqbz )
! !     igk2ibz = 1
! !     igkq2ibz = 1
! !   elseif ( rixs( ispec(itr) )%ng == nopused ) then
! ! ! full symmetry
! !     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
! !     kibz = ik2ibz( kbz )
! !     kqibz = ik2ibz( kqbz )
! !     krbz = rixs( ispec(itr) )%ik2rbz(kbz)
! !     if ( krbz /= kibz ) stop '   krbz /= kibz '
! !     igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
! !     igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !   else
! ! ! partial symmetry
! !     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
! !     kibz = ik2ibz( kbz )
! !     kqibz = ik2ibz( kqbz )
! !     krbz = rixs( ispec(itr) )%ik2rbz(kbz)
! !     igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
! !     igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !   endif

! !   kqbz = rixs( ispec(itr) )%ik2kq(kbz)
! !   if ( nopused == 1 ) then
! !     krbz = kbz
! !     kibz = kbz
! !     kqibz = kqbz
! !     igk2ibz = 1
! !     igkq2ibz = 1
! !     igk2rbz = 1
! !     igkq2rbz = 1
! !   else
! !     if ( rixs( ispec(itr) )%ng == 1 ) then
! !       krbz = kbz
! !       kibz = ik2ibz(kbz)
! !       kqibz = ik2ibz(kqbz)
! !       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
! !       igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !       igk2rbz = 1
! !       igkq2rbz = 1
! !     else if ( rixs( ispec(itr) )%ng == nopused ) then
! !       kibz = ik2ibz(kbz)
! !       kqibz = ik2ibz(kqbz)
! !       krbz = kibz
! !       igk2ibz = 1 ! ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
! !       igkq2ibz = 1 ! ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !       igk2rbz = 1
! !       igkq2rbz = 1   
! !     else
! !       krbz = rixs( ispec(itr) )%ik2rbz(kbz)
! ! !$$$      kqbz = rixs( ispec(itr) )%ik2kq(kbz)
! !       kibz = ik2ibz(kbz)
! !       kqibz = ik2ibz(kqbz)
! !       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
! !       igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !       igk2rbz = rixs( ispec(itr) )%ig4q( i1(kbz) , i2(kbz) , i3(kbz) )
! !       igkq2rbz = rixs( ispec(itr) )%ig4q( i1(kqbz) , i2(kqbz) , i3(kqbz) )
! !     endif
! !   endif

!   if ( nopused == 1 ) then
!     kbz = krbz
!     kibz = kbz
!     igk2ibz = 1
!     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!     kqibz = kqbz
!     igkq2ibz = 1
!   else
!     if ( rixs( ispec(itr) )%ng == 1 ) then
!       kbz = krbz
!       kibz = ik2ibz(kbz)
!       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!       kqibz = ik2ibz(kqbz)
!       igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
!     else
!       kbz = rixs( ispec(itr) )%ikbz2rbz(krbz)
!       kibz = ipqibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!       kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!       kqibz = ipqibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
!       igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
!     endif
!   endif

! enddo
