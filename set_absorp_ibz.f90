SUBROUTINE SET_ABSORP_IBZ(ia,kpk)
use m_aux , only : lmtindex
use m_bnd
use m_bz
use m_bzmesh , only : nkibz
use m_functions
use m_sdt , only : igtta
use m_rixs
implicit none
integer , intent(in) :: ia
real(REAL64) kpk( nb, nkibz )
! local vars
character(18), parameter :: srcname = ' in SET_ABSORP_IBZ'
complex(REAL64) , allocatable :: Isym( : , : , : )
complex(REAL64) , parameter :: z0 = ( 0._REAL64 , 0._REAL64 )
complex(REAL64) mmek(3), mmek2(3), I( 3, 3 )
integer ib0, ib, k, ind, iop, iat, iat1, iatn
real(REAL64) op( 3, 3 )

! if(ia < 1 .or. ia > 2*sum(abs(rixs(:)%nk))) &
!   stop ' In SET_ABSORP_IBZ: wrong ia value.'

if ( ia < 1 .or. ia > n_mme ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) &
    'In SET_ABSORP_BZ: wrong ia value.'
  call exit_on_error( srcname )
endif

! if(.not.allocated(lmtdata))then
!   write(unit=output_unit,fmt=5)'lmtdata',ia
!   return
! endif
! if(.not.allocated(lmtindex))then
!   write(unit=output_unit,fmt=5)'lmtindex',ia
!   return
! endif
! if(.not.allocated(lmtdata(lmtindex(ia))%mme))then
!   write(unit=output_unit,fmt=5)'lmtdata(lmtindex(ia))%mme',ia
! endif

  5 format( 4x, 'SET_ABSORP_IBZ: ', a, ' is not allocated, ia=', i0 )

iat1 = lmtdata( lmtindex(ia) )%iat1
iatn = lmtdata( lmtindex(ia) )%iatn
! allocate(Isym(3,3,iat1:iatn),stat=iok,errmsg=msg)
! if(iok /= 0)call print_allocation_error('Isym',iok,msg)

call allocate( Isym, 'Isym' // srcname, udim1 = 3, udim2 = 3, &
     ldim3 = iat1, udim3 = iatn )

! now always ib0 = nbmin - 1 = 0. It's just kept for future as possible way
! of generalization
ib0 = 0
do k = 1, nkibz                    ! over all k from IBZ
  do ib = 1, nb                    ! over all bands
    Isym( :, :, : ) = z0
    do iat = iat1, iatn
      mmek(:) = lmtdata( lmtindex(ia) ) % mme( :, iat, ia, ib, k )
! symmetrize atomic tensor in global system

      ! call tensor_from_mme( mmek , I )
      ! Isym( : , : , iat ) = Isym( : , : , iat ) + I( : , : )

      do iop = 1, nopused
        ind = igtta( iop, iat )
        op( : , : ) = g( : , : , iop )
!$$$        mmek2(:) = mme_trans( iop , mmek )
        if ( iop == 1 ) then
          mmek2(:) = mmek(:)
        else
          if ( iopnum(iop) > 0 ) then
            mmek2(:) = matmul( op, mmek )
          else
            mmek2(:) = matmul( op , conjg( mmek ) )
          endif
        endif
        call tensor_from_mme( mmek2 , I )
        Isym( : , : , ind ) = Isym( : , : , ind ) + I( : , : )
      enddo                       ! iop
    enddo                         ! iat
!    Isym(:,:,:)=Isym(:,:,:)/cmplx(nopused*8*(iatn-iat1+1),kind=REAL64)
    Isym( : , : , : ) = Isym( : , : , : ) / &
                        cmplx( nopused * 8 * nfu , kind = REAL64)
    kpk( ib0 + ib , k ) = real( sum( Isym( 1 , 1 , iat1 : iatn ) ) &
                        + sum( Isym( 2 , 2 , iat1 : iatn ) ) , kind = REAL64 )
  enddo                           ! ib
enddo                             ! k    

END SUBROUTINE SET_ABSORP_IBZ

