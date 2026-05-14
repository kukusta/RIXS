SUBROUTINE SET_A( itr, a )
use m_aux
use m_bnd
use m_bz
use m_bzmesh, only: ndxyz
use m_functions, only: exit_on_error
use m_params
use m_rixs
implicit none
integer :: itr
real a( rixs( ispec(itr) )%nkrbz )
! local vars
character(9)   , parameter :: srcname = ' in SET_A'
complex(REAL64), parameter :: z0 = ( 0._REAL64, 0._REAL64 )
complex(REAL64), parameter :: zi = ( 0._REAL64, 1._REAL64 )
complex(REAL64) dpk, dpkq, denom, za, mme_kbz(3), mme_kqbz(3), dza
integer krbz, kibz, jat, jat1, jatn, i, kbz, kqbz, k1, k2, k3, ja
real(REAL64), parameter :: eps = 1.e-6_REAL64
real(REAL64) phase, q(3), rat(3), dpi, akrbz

! set all elements to -1. to check if all of them were calculated
a( : ) = -1.
q( : ) = rixs( ispec(itr) )%q( : )
kbz = 0
krbz = 0
do k3 = 1, ndxyz(3)
  do k2 = 1, ndxyz(2)
    do k1 = 1, ndxyz(1)
      kbz = kbz + 1

!$$$! if ( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 ) the point ( k1, k2, k3 )
!$$$! can be reduced
!$$$      if(allocated( rixs( ispec(itr) )%ig4q ) )then
!$$$        ; !write(*,*)shape( rixs( ispec(itr) )%ig4q )
!$$$      else
!$$$        write(*,*)'not allocated',rixs( ispec(itr) )%ng
!$$$        call exit_on_error(' in set_a')
!$$$      endif

      if( rixs( ispec(itr) )%ng > 1 ) then
        if( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 ) then
          cycle
        endif
      endif

      kqbz = rixs( ispec(itr) )%ik2kq( kbz )
!$$$      if ( nopused == 1 ) then
!$$$        kibz = kbz
!$$$        krbz = krbz + 1
!$$$      else

!$$$      kibz = ik2ibz( kbz )
!$$$      if( rixs( ispec(itr) )%ng == 1 ) then
!$$$        krbz = krbz + 1
!$$$      else
!$$$! if ( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 ) the point ( k1, k2, k3 )
!$$$! can be reduced
!$$$        if( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 ) cycle
!$$$        krbz = krbz + 1
!$$$        if( krbz /= rixs( ispec(itr) )%ik2rbz( kbz ) ) then
!$$$          write( unit = output_unit, fmt = "( 4x, a )") &
!$$$            'krbz /= rixs( ispec(itr) ) %ik2rbz( kbz )'
!$$$           call exit_on_error( srcname )
!$$$        endif
!$$$      endif
!$$$!$$$      endif

      if ( rixs( ispec(itr) )%ng == 1 ) then
        krbz = kbz
        kibz = ik2ibz(kbz)
        ! kqibz = ik2ibz(kqbz)
      else
        kibz = ik2ibz(kbz)
        krbz = rixs( ispec(itr) )%ik2rbz(kbz)
        ! kqibz = ik2ibz(kqbz)
      endif

      za = z0
!$$$      akrbz = 0._REAL64
      do ja = 1, n_mme
        if( map( ja, ispec( itr ) ) == 0 ) then
!$$$          if ( print ) then
!$$$!$omp critical
!$$$            write( unit = output_unit, fmt = 10 ) ja, ispec(itr)
!$$$!$omp end critical
!$$$          endif
          cycle
        endif
        jat1 = lmtdata( lmtindex( ja ) )%iat1
        jatn = lmtdata( lmtindex( ja ) )%iatn

        do jat = jat1, jatn
! if ( nopused == 1 ) then
!   lmtdata(isp)%mmebz is NOT neither allocated or referenced
!   use lmtdata(isp)%mme instead
! else
!   lmtdata(isp)%mmebz IS allocated and referenced
! endif
!$$$          if ( nopused == 1 ) then
!$$$            mme_kbz(:) = lmtdata( lmtindex( ja ) )%mme( : , jat , ja , &
!$$$                         ibi(itr) , kbz )
!$$$            mme_kqbz(:) = lmtdata( lmtindex( ja ) )%mme( : , jat , ja , &
!$$$                          ibf(itr) , kqbz )
!$$$          else

          mme_kbz(:) = lmtdata( lmtindex( ja ) )%mmebz( : , jat , ja , &
                       ibi(itr) , kbz ) ! ibi
          mme_kqbz(:) = lmtdata( lmtindex( ja ) )%mmebz( : , jat , ja , &
                        ibf(itr) , kqbz ) ! ibf

!$$$          endif

          denom = rixs( ispec(itr) )%en + &
                  lmtdata( lmtindex( ja ) )%esp( ja ) - &
                  e( ibf(itr), kibz ) + zi * rixs( ispec(itr) )%Gamma

          rat(:) = lmtdata( lmtindex( ja ) )%rat( : , jat )
          phase = 2._REAL64 * dpi() * dot_product( q, rat )

          select case ( abs(irixs) )
            case( 1 )
              dpk = dot_product( mme_kbz, rixs( ispec(itr) )%e_in )
              dpkq = dot_product( mme_kqbz, rixs( ispec(itr) )%e_out )
              za = za + exp( -zi * phase ) * conjg(dpkq) * dpk / denom
!$$$              akrbz = akrbz + abs( conjg(dpkq) * dpk / denom )
            case ( 2 )
              dza = z0
              do i = 1, 3
                dza = dza + conjg( mme_kqbz(i) ) * mme_kbz(i) / denom
              enddo
              za = za + exp( -zi * phase ) * dza
!$$$              akrbz = akrbz + abs( dza )
            case ( 3 )
              za = za + 1._REAL64 / denom
!$$$              akrbz = akrbz + 1._REAL64 / abs( denom )
            case ( 4 )
              dpk = dot_product( mme_kbz, rixs( ispec(itr) )%e_in )
              dpkq = dot_product( mme_kqbz, rixs( ispec(itr) )%e_out )
              za = za + exp( -zi * phase ) * conjg( dpkq ) * dpk
!$$$              akrbz = akrbz + abs( conjg( dpkq ) * dpk )
            case default
              call exit_on_error( ' in SET_A: unknown irixs value' )
          end select
        enddo          ! jat
      enddo          ! ja
      akrbz = real( abs( za ) ** 2 )
      if ( a(krbz) > -0.5 ) then
        if ( abs( a(krbz) - akrbz ) > eps ) then
!$omp critical
          if( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 )then
            write( unit = output_unit, fmt = 70 ) kbz, k1, k2, k3, krbz, &
                 rixs(ispec(itr))%ig4q( k1, k2, k3 ), a(krbz), akrbz
          endif
!$$$          write( unit = output_unit, fmt = "( 4x, a, i0, a, 2f15.10 )" ) &
!$$$            'For krbz=', krbz, ' a(krbz) .ne. real(abs(za)**2): ', &
!$$$            a(krbz), akrbz
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

 10 format( 4x, 'Skipping ia=', i0, ' for ispec=', i0 )
 15 format( 4x, "dza=z0 for transition with ibi= ", i3, ", ibf= ", i3&
          , ", ia= ", i2, ", iat= ", i3, ", k= ", i0 )
 20 format( 4x, "Contributions from different atoms:" )
 25 format( 4x, "itr=", i8, " iat=", i3, " phase/2/pi=", f12.6&
          , " exp(phase)= (", f12.6, ",", f12.6, ")"&
          , " conjg(dpkq)=(", f12.6, ",", f12.6, ")"&
          , " dpk/denom=(", f12.6, ",", f12.6, ")" )
 55 format( 4x, "ERROR: mme < 0 for krbz=", i0 )
 60 format( 4x, "ERROR: a(krbz) < 0 for some point in RBZ" )
 70 format( /, 4x, 'For kbz=', i0, ': i1=', i0, ', i2=', i0, ', i3=', &
            i0, ', krbz=', i0, ', ig4q(kbz -> krbz)=', i0, &
            ' AND a(krbz) != real(abs(za)**2):', f10.7, ' !=', f10.7)

rixs(ispec(itr))%max_mme=max(rixs(ispec(itr))%max_mme,maxval(a))
rixs(ispec(itr))%min_mme=min(rixs(ispec(itr))%min_mme,minval(a))
rixs(ispec(itr))%av_mme =    rixs(ispec(itr))%av_mme +   sum(a)

END SUBROUTINE SET_A

      ! if ( nopused == 1 ) then
      !   kibz = kbz
      !   krbz = kbz
      !   kqbz = rixs( ispec(itr) )%ik2kq(kbz)
      !   kqibz = kbz
      !   igk2ibz = 1
      !   igkq2ibz = 1
      ! else
      !   kibz = ipqibz( k1 , k2 , k3 ) !ik2ibz(kbz)
      !   kqbz = rixs( ispec(itr) )%ik2kq(kbz)
      !   kqibz = ipqibz( i1(kqbz) , i2(kqbz) , i3(kqbz) ) !ik2ibz(kqbz)
      !   igk2ibz = ig4qibz( k1 , k2 , k3 )
      !   igkq2ibz = ig4qibz( i1(kqbz) , i2(kqbz) , i3(kqbz) )
      !   if ( rixs( ispec(itr) )%ng == 1 ) then
      !     krbz = kbz
      !   else if ( rixs( ispec(itr) )%ng == nopused ) then
      !     krbz = kibz
      !     ! if ( rixs( ispec(itr) )%ik2rbz(kbz) /= kibz ) &
      !     !    stop '   why ik2rbz(kbz) /= kibz?'
      !   else if ( rixs( ispec(itr) )%ig4q( k1 , k2 , k3 ) == 1 ) then
      !       krbz = rixs( ispec(itr) )%ik2rbz(kbz)
      !   else
      !     cycle
      !   endif
      ! endif

      ! if ( igk2ibz == 1 ) then
      !   mme_krbz(:) = mme_kibz(:)
      ! else
      !   op( : , : ) = g( : , : , igk2ibz )
      !   if ( iopnum(igk2ibz) > 0 ) then
      !     mme_krbz(:) = matmul( op , mme_kibz )
      !   else
      !     mme_krbz(:) = matmul( op , conjg( mme_kibz ) )
      !   endif
      ! endif
      ! if ( igkq2ibz == 1 ) then
      !   mme_kqrbz(:) = mme_kqibz(:)
      ! else
      !   op( : , : ) = g( : , : , igkq2ibz )
      !   if ( iopnum(igkq2ibz) > 0 ) then
      !     mme_kqrbz(:) = matmul( op , mme_kqibz )
      !   else
      !     mme_kqrbz(:) = matmul( op , conjg( mme_kqibz ) )
      !   endif
      ! endif
