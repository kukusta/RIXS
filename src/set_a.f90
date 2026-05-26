SUBROUTINE set_a(itr, a)
use m_aux
use m_bnd
use m_bz
use m_bzmesh, only: nkbz, ndxyz
use m_constants, only: z0, zi, eps6
use m_functions, only: exit_on_error, int2string
use m_params
use m_rixs
implicit none
integer :: itr
!vvv real a(rixs(ispec(itr))%nkrbz)
real a(nkbz)
! local vars
character(9), parameter :: srcname = ' in SET_A'
complex(REAL64) dpk, dpkq, denom, za, mme_kbz(3), mme_kqbz(3), dza
integer kibz, jat, jat1, jatn, i, kbz, kqbz, k1, k2, k3, ja
real(REAL64) phase, q(3), rat(3), dpi, akbz

! set all elements to -1. to check if all of them were calculated
a(:) = -1.
q(:) = rixs(ispec(itr))%q(:)
kbz  = 0
!vvv krbz = 0
do k3 = 1, ndxyz(3)
  do k2 = 1, ndxyz(2)
    do k1 = 1, ndxyz(1)
      kbz = kbz + 1
!vvv      if (rixs(ispec(itr))%ng > 1) then
!vvv        if( rixs( ispec(itr) )%ig4q( k1, k2, k3 ) /= 1 ) then
!vvv          cycle
!vvv        endif
!vvv      endif

      kibz = ik2ibz(kbz)
      kqbz = rixs(ispec(itr))%ik2kq(kbz)

      if (rixs(ispec(itr))%Reg_BZ_states(:) == 'above ' .and. &
        rixs(ispec(itr))%map_region(kbz) == 0) then
        a(kbz) = 0.
        cycle
      endif

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

!vvv      if (rixs(ispec(itr))%ng == 1) then
!vvv        krbz = kbz
!vvv        kibz = ik2ibz(kbz)
!vvv        ! kqibz = ik2ibz(kqbz)
!vvv      else
!vvv        kibz = ik2ibz(kbz)
!vvv        krbz = rixs(ispec(itr))%ik2rbz(kbz)
!vvv        ! kqibz = ik2ibz(kqbz)
!vvv      endif

      za = z0
!$$$      akrbz = 0._REAL64
      do ja = 1, n_mme
        if (map(ja, ispec(itr)) == 0) cycle
        jat1 = lmtdata(lmtindex(ja))%iat1
        jatn = lmtdata(lmtindex(ja))%iatn

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

          if (switch_ini_final) then
! old version
            mme_kbz(:)  = lmtdata(lmtindex(ja))%mmebz(:, jat, ja, ibi(itr), kbz)  ! ibi
            mme_kqbz(:) = lmtdata(lmtindex(ja))%mmebz(:, jat, ja, ibf(itr), kqbz) ! ibf
          else
! default version. I think it is correct. ini states with k lie below E_F
! final states with k+q lie above E_F
            mme_kbz(:)  = lmtdata(lmtindex(ja))%mmebz(:, jat, ja, ibf(itr), kbz)  ! ibf
            mme_kqbz(:) = lmtdata(lmtindex(ja))%mmebz(:, jat, ja, ibi(itr), kqbz) ! ibi
          endif

!$$$          endif

          denom = rixs(ispec(itr))%en + lmtdata(lmtindex(ja))%esp(ja) - &
               e(ibf(itr), kibz) + zi * rixs(ispec(itr))%Gamma

          rat(:) = lmtdata(lmtindex(ja))%rat(:, jat)
          phase  = 2._REAL64 * dpi() * dot_product(q, rat)

          ! dpk  = dot_product(mme_kbz,  rixs(ispec(itr))%e_in )
          ! dpkq = dot_product(mme_kqbz, rixs(ispec(itr))%e_out)
          dpk  = z0
          dpkq = z0
          do i = 1, 3
            dpk  = dpk  + mme_kbz(i)  * rixs(ispec(itr))%e_in(i)
            dpkq = dpkq + mme_kqbz(i) * rixs(ispec(itr))%e_out(i)
          enddo

          if (irixs == 1) then
             za   = za + exp(-zi * phase) * conjg(dpkq) * dpk / denom
!$$$              akrbz = akrbz + abs( conjg(dpkq) * dpk / denom )
          else if (irixs == 2) then
            dza = z0
            do i = 1, 3
              dza = dza + conjg(mme_kqbz(i)) * mme_kbz(i) / denom
            enddo
            za = za + exp(-zi * phase) * dza / 3._REAL64
!$$$              akrbz = akrbz + abs( dza )
          else if (irixs == 3) then
            za = za + 1._REAL64 / denom
!$$$              akrbz = akrbz + 1._REAL64 / abs( denom )
          else if (irixs == 4) then
            za   = za + exp(-zi * phase) * conjg(dpkq) * dpk
!$$$              akrbz = akrbz + abs( conjg( dpkq ) * dpk )
          else
            call exit_on_error(' in SET_A: unknown irixs=' // &
                               trim(int2string(irixs)) // ' value')
          endif
        enddo
      enddo
!vvv      akrbz = real(abs(za) ** 2)
      akbz = real(abs(za) ** 2)
      if (.false. .and. a(kbz) > -0.5) then
        if (abs( a(kbz) - akbz ) > eps6 ) then
!$omp critical
          if (rixs(ispec(itr))%ig4q(k1, k2, k3) /= 1) then
            write(unit = output_unit, fmt = 70) kbz, k1, k2, k3, kbz, &
                 rixs(ispec(itr))%ig4q(k1, k2, k3), a(kbz), akbz
          endif
          call exit_on_error(srcname)
!$omp end critical
        endif
      else 
!vvv      a(krbz) = akrbz
        a(kbz) = akbz
      endif
    enddo
  enddo
enddo

if (rixs(ispec(itr))%Reg_BZ_states(:) == 'below ') then
  do kqbz = 1, nkbz
    if (rixs(ispec(itr))%map_region(kqbz) == 0) a(kqbz) = 0.
  enddo
endif

!$omp critical
if (minval(a(:)) < -eps6) then
  write(*,*) 'min(a) =', minval(a)
  call exit_on_error(srcname)
endif
!$omp end critical

!$$$ do krbz=1,rixs(ispec(itr))%nkrbz
!$$$   if(a(krbz) < 0.)write(unit=output_unit,fmt=55)krbz
!$$$ enddo
!$$$ if ( minval( a ) < 0. ) then
!$$$   write( unit = output_unit , fmt = 60 )
!$$$   call exit_on_error( srcname )
!$$$ endif

 55 format(   4x, "ERROR: mme < 0 for krbz=", i0)
 60 format(   4x, "ERROR: a(krbz) < 0 for some point in RBZ")
 70 format(/, 4x, 'For kbz=', i0, ': i1=', i0, ', i2=', i0, ', i3=', &
              i0, ', krbz=', i0, ', ig4q(kbz -> krbz)=', i0,         &
              ' AND a(krbz) != real(abs(za)**2):', f10.7, ' !=', f10.7)

rixs(ispec(itr))%max_mme=max(rixs(ispec(itr))%max_mme,maxval(a))
rixs(ispec(itr))%min_mme=min(rixs(ispec(itr))%min_mme,minval(a))
rixs(ispec(itr))%av_mme =    rixs(ispec(itr))%av_mme +   sum(a)

END SUBROUTINE set_a

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
