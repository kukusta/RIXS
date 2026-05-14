SUBROUTINE SET_A(itr,a)
use m_aux
use m_bnd
use m_bz
use m_bzmesh, only : nkbz
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
        igk2rbz, igkq2rbz, kbz, kqbz
real akrbz
real(REAL64), parameter :: eps = 1.e-6_REAL64
real(REAL64) phase, q(3), rat(3), dpi, op(3,3)

! set all elements to -1. to check if all of them were calculated
a(:) = -1.
q(:) = rixs( ispec(itr) )%q(:)
jat1 = lmtdata( lmtindex( ia(itr) ) )%iat1
jatn = lmtdata( lmtindex( ia(itr) ) )%iatn

!$$$do krbz = 1 , rixs( ispec(itr) )%nkrbz
do kbz = 1 , nkbz

  ! kbz = rixs( ispec(itr) )%ikbz2rbz(krbz)
  ! kqbz = rixs( ispec(itr) )%ik2kq(kbz)
  ! kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
  ! kqibz = ipqibz( i1(kqbz), i2(kqbz), i3(kqbz) )
  ! igkibz = ig4qibz( i1(kbz), i2(kbz), i3(kbz) )
  ! igkqibz = ig4qibz( i1(kqbz), i2(kqbz), i3(kqbz) )

!   if ( rixs( ispec(itr) )%ng == 1 ) then
! ! no symmetry
!     krbz = kbz
!     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!     ! kibz = ipqibz( i1(kbz), i2(kbz), i3(kbz) )
!     ! kqibz = ipqibz( i1(kqbz), i2(kqbz), i3(kqbz) )
!     kibz = ik2ibz( kbz )
!     kqibz = ik2ibz( kqbz )
!     igk2ibz = 1
!     igkq2ibz = 1
!   elseif ( rixs( ispec(itr) )%ng == nopused ) then
! ! full symmetry
!     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!     kibz = ik2ibz( kbz )
!     kqibz = ik2ibz( kqbz )
!     krbz = rixs( ispec(itr) )%ik2rbz(kbz)
!     if ( krbz /= kibz ) stop '   krbz /= kibz '
!     igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!     igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
!   else
! ! partial symmetry
!     kqbz = rixs( ispec(itr) )%ik2kq(kbz)
!     kibz = ik2ibz( kbz )
!     kqibz = ik2ibz( kqbz )
!     krbz = rixs( ispec(itr) )%ik2rbz(kbz)
!     igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
!     igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
!   endif

  kqbz = rixs( ispec(itr) )%ik2kq(kbz)
  if ( nopused == 1 ) then
    krbz = kbz
    kibz = kbz
    kqibz = kqbz
    igk2ibz = 1
    igkq2ibz = 1
    igk2rbz = 1
    igkq2rbz = 1
  else
    if ( rixs( ispec(itr) )%ng == 1 ) then
      krbz = kbz
      kibz = ik2ibz(kbz)
      kqibz = ik2ibz(kqbz)
      igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
      igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
      igk2rbz = 1
      igkq2rbz = 1
    else
      krbz = rixs( ispec(itr) )%ik2rbz(kbz)
!$$$      kqbz = rixs( ispec(itr) )%ik2kq(kbz)
      kibz = ik2ibz(kbz)
      kqibz = ik2ibz(kqbz)
      igk2ibz = ig4qibz( i1(kbz) , i2(kbz) , i3(kbz) )
      igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
      igk2rbz = rixs( ispec(itr) )%ig4q( i1(kbz) , i2(kbz) , i3(kbz) )
      igkq2rbz = rixs( ispec(itr) )%ig4q( i1(kqbz) , i2(kqbz) , i3(kqbz) )
    endif
  endif

!$$$!$omp ordered
!$$$  write(*,*)'itr=',itr,' ispec=',ispec(itr),' igk2ibz=',igk2ibz&
!$$$       ,' igkq2ibz=',igkq2ibz,' igk2rbz=',igk2rbz,' igkq2rbz=',igkq2rbz
!$$$!$omp end ordered

  za = z0
  do jat = jat1, jatn
    mme_kibz(:) = lmtdata( lmtindex( ia(itr) ) )%mme( :, jat, ia(itr), &
                  ibi(itr), kibz )
    mme_kqibz(:) = lmtdata( lmtindex( ia(itr) ) )%mme( :, jat, ia(itr), &
                   ibf(itr), kqibz )
    ! if ( igk2ibz == 1 ) then
    !   mme_krbz(:) = mme_kibz(:)
    ! else
    !   mme_krbz(:) = mme_trans( igk2ibz, mme_kibz )
    ! endif
    ! if ( igkq2ibz == 1 ) then
    !   mme_kqrbz(:) = mme_kqibz(:)
    ! else
    !   mme_kqrbz(:) = mme_trans( igkq2ibz, mme_kqibz )
    ! endif
    ! if ( igk2rbz == 1 ) then
    !   mme_krbz(:) = mme_krbz(:)
    ! else
    !   mme_krbz(:) = mme_trans_inv( igk2rbz, mme_krbz )
    ! endif
    ! if ( igkq2rbz == 1 ) then
    !   mme_kqrbz(:) = mme_kqrbz(:)
    ! else
    !   mme_kqrbz(:) = mme_trans_inv( igkq2rbz, mme_kqrbz )
    ! endif

    if ( igk2ibz == 1 ) then
      mme_krbz(:) = mme_kibz(:)
    else
      op( : , : ) = g( : , : , igk2ibz )
      if ( iopnum(igk2ibz) > 0 ) then
!$$$        mme_krbz(:) = mme_trans_new( op , mme_kibz )
        mme_krbz(:) = matmul( op , mme_kibz )
      else
!$$$        mme_krbz(:) = mme_trans_new( op , conjg(mme_kibz) )
        mme_krbz(:) = matmul( op , conjg( mme_kibz ) )
      endif
    endif
    if ( igkq2ibz == 1 ) then
      mme_kqrbz(:) = mme_kqibz(:)
    else
      op( : , : ) = g( : , : , igkq2ibz )
      if ( iopnum(igkq2ibz) > 0 ) then
!$$$        mme_kqrbz(:) = mme_trans_new( op , mme_kqibz )
        mme_kqrbz(:) = matmul( op , mme_kqibz )
      else
!$$$        mme_kqrbz(:) = mme_trans_new( op , conjg(mme_kqibz) )
        mme_kqrbz(:) = matmul( op , conjg( mme_kqibz ) )
      endif
    endif

    if ( igk2rbz /= 1 ) then
      op(:,:) = transpose( rixs( ispec(itr) )%g( : , : , igk2rbz ) )
      if ( rixs( ispec(itr) )%iopnum( igk2rbz ) > 0 ) then
!$$$        mme_krbz(:) = mme_trans_new( op , mme_krbz )
        mme_krbz(:) = matmul( op , mme_krbz )
      else
!$$$        mme_krbz(:) = mme_trans_new( op , conjg(mme_krbz) )
        mme_krbz(:) = matmul( op , conjg(mme_krbz) )
      endif
    endif
    if ( igkq2rbz /= 1 ) then
      op(:,:) = transpose( rixs( ispec(itr) )%g(:,:,igkq2rbz) )
      if ( rixs(ispec(itr))%iopnum(igkq2rbz) > 0 ) then
!$$$        mme_kqrbz(:) = mme_trans_new( op , mme_kqrbz )
        mme_kqrbz(:) = matmul( op , mme_kqrbz )
      else
!$$$        mme_kqrbz(:) = mme_trans_new( op , conjg(mme_kqrbz) )
        mme_kqrbz(:) = matmul( op , conjg(mme_kqrbz) )
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
enddo

if(abs(irixs) == 2) a(:)=a(:)/3.

!$$$ do krbz=1,rixs(ispec(itr))%nkrbz
!$$$   if(a(krbz) < 0.)write(unit=output_unit,fmt=55)krbz
!$$$ enddo
if ( minval( a ) < 0. ) then
  write( unit = output_unit , fmt = 60 )
  call exit_on_error( srcname )
endif

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
