SUBROUTINE GEOMETRY
use, intrinsic :: iso_fortran_env
use m_files
use m_functions, only: exit_on_error
use m_params, only: iprint, irixs, check_kstar, debug_mode, read_rixfile
use m_rixs
use m_sdt, only: ut1
implicit none
! local vars
character(12), parameter :: srcname=' in GEOMETRY'
complex(REAL64), parameter :: zi=(0._REAL64,1._REAL64)
integer i,isp
real(REAL64) sigma(3),pi_in(3),pi_out(3),norm,ang_in,lenq,q(3)&
     ,ang_out,ang_inout,Qglobal(3),n(3),nlen,R(3,3),nDual(3,3),alpha&
     ,AbasR(3,3),AbasRinv(3,3),det,in(3),out(3),len_e_in,len_e_out
real(REAL64), parameter :: eps=1.e-5_REAL64,s13=sqrt(1._REAL64/3._REAL64)
real(REAL64), parameter :: r2d=45._REAL64/atan(1._REAL64)

if ( check_kstar .and. debug_mode ) return
if ( read_rixfile ) return

  5 format(4x,'Can not build scattering plane for SPEC=',i0&
          ,/,4x,'Vectors k_(in,out) and q are collinear.')
  7 format('')
  8 format(4x,'WARNING: for ISPEC=',i0,' q=0 so polarization vectors '&
          ,'must be set manually.')
 10 format(4x,'GEOMETRICAL DETAILS (angles and rotational axis).',/&
          ,4x,'ISPEC',3x,'(q,k_in)',5x,'(q,k_out)',3x,'(k_in,k_out)'&
          ,13x,'rotational axis',14x,'angle')
 11 format(4x,'Error: for isp=',i2,' rotational matrix is not orthogonal'&
          ,', it"s det=',f0.10)
 ! 13 format(4x,'For isp=',i2,' rotational axis is ( ',3f12.8,' ), rotational'&
 !          ,' angle is ',f10.4,' deg.')
 15 format(4x,i3,3f13.4)
 16 format(4x,i3,3f13.4,4x,'(',3f12.8,')',f9.3)

if(iprint < 0)write(unit=output_unit,fmt=7)
! Transform q vectors into local Cartesian system. q vector is set as a
! fractions of G1, G2, G3 in units of 2pi/a, 2pi/b, 2pi/c. After matmul
! q is expressed in units of 2pi/a. So I need to apply 2pi multiplier
! while calculating RIXS mme.
do isp=1,nrixs

! shift q from G point (where it's not possible to build scattering plane)
! in z direction so (0,0,0) -> (0,0,1/ndiv(z))
!  forall(i=1:3) dndxyz(i)=1._REAL64/real(ndxyz(i),kind=REAL64)
!  if( dot_product(rixs(isp)%q(:),rixs(isp)%q(:)) < &
!      min(dndxyz(1),dndxyz(2),dndxyz(3)) &
!    )then
!    rixs(isp)%q(:)=(/0._REAL64,0._REAL64,dndxyz(3)/)
!    rixs(isp)%q(:)=matmul(ut1,rixs(isp)%q)
!    cycle
!  endif

  rixs(isp)%q(:)=matmul(ut1,rixs(isp)%q)

! length of q vector is calculated here to print a message if necessary
  q(:)=rixs(isp)%q(:)
  lenq=sqrt(dot_product(q,q))
  if(lenq < eps .and. abs(irixs) == 1) write(unit=output_unit,fmt=8)isp

enddo
write( unit = output_unit, fmt = 7 )

if(irixs == 0 .or. irixs == 2)then
! Joint DOS or average over polarizations
! e_in and e_out are these variables really needed?
  do isp=1,nrixs
    rixs(isp)%e_in(:)=0._REAL64
    rixs(isp)%e_out(:)=0._REAL64
    rixs(isp)%k_out(:)=k_in(:)  !  rixs(isp)%k_in(:)
  enddo
  return
endif
norm = sqrt( dot_product( Qaxis, Qaxis ) )
if( norm < eps ) then
  write( unit = output_unit, fmt = "( 4x, a )" ) 'ERROR: Length of Qaxis is zero'
  call exit_on_error( srcname )
endif
Qaxis(:) = Qaxis(:) / norm

if( iprint > 0 ) write( unit = output_unit, fmt = 10 )
do isp = 1, nrixs
  q(:) = rixs(isp)%q(:)
  lenq = sqrt( dot_product( q, q ) )
  if(lenq < eps) then
    ! rixs(isp)%e_in(:)=(/s13,s13,s13/)
    ! rixs(isp)%e_out(:)=(/s13,s13,s13/)
! for q=0 it's not possible to build scattering plane. If one wants to
! calculate spectra for q=0 it's necessary to specify polarization vectors
! manually
    rixs(isp)%channel(:) = 'man '
    call polarization_vectors_for_man(isp)
    rixs(isp)%k_in(:) = k_in(:)
    rixs(isp)%k_out(:) = k_in(:)
  else
! Find components of input vector in local Cartesian frame
    Qglobal(:) = matmul( Abas, q )
    norm = sqrt( dot_product( Qglobal, Qglobal ) )
    Qglobal(:) = Qglobal(:) / norm
! Rotational axis
    call cross( n, Qglobal, Qaxis )
    nlen = sqrt( dot_product( n, n ) )
! defaults: rotational tensor, tensor dual to n(:), and rotational angle
    R( : , : ) = 0._REAL64
    forall( i = 1 : 3 ) R( i, i ) = 1._REAL64
    nDual( :, : ) = 0._REAL64
    alpha = 0._REAL64
    rixs(isp)%k_in(:) = k_in(:)
    call dinv33( Abas, 0, AbasRinv, det )
    if( nlen > eps ) then
      n(:) = n(:) / nlen
      nDual( 1, 2 )  = -n(3)
      nDual( 1, 3 )  =  n(2)
      nDual( 2, 3 )  = -n(1)
      nDual( : , : ) = nDual( : , : ) - transpose( nDual )
      alpha = acos( dot_product( Qglobal, Qaxis ) )
      R( : , :) = R( : , : ) + nDual( : , : ) * sin(alpha) + &
                  matmul( nDual, nDual ) * ( 1 - cos(alpha) )
      AbasR( : , : ) = matmul( R, Abas )
      call dinv33( AbasR, 0, AbasRinv, det )
      if( det < 1._REAL64 - eps ) then
        write( unit = output_unit, fmt = 11 ) isp, det
        call exit_on_error( srcname )
      endif
!      if(iprint > 0)write(unit=output_unit,fmt=13)isp,n(:),r2d*alpha
    endif
    rixs(isp)%k_in(:) = matmul( AbasRinv, k_in )
! Calculate k_out in such a way that abs(klocal)=abs(k_out)
! in and out are local normalized copies of rixs(isp)%k_in and rixs(isp)%k_out
    norm = sqrt( dot_product( rixs(isp)%k_in, rixs(isp)%k_in ) )
    if( norm > eps ) then
      in(:) = rixs(isp)%k_in(:) / norm
    else
      stop '    Input vector k_in(:) == 0.'
    endif
    rixs(isp)%k_out(:) = rixs(isp)%k_in(:) - &
         2._REAL64 * q(:) * dot_product( rixs(isp)%k_in, q ) / lenq ** 2
    norm = sqrt( dot_product( rixs(isp)%k_out, rixs(isp)%k_out ) )
    out(:) = rixs(isp)%k_out(:) / norm
    ang_in = r2d * acos( dot_product( in, q ) / lenq )
    ang_out = r2d * acos( dot_product( out, q ) / lenq )
    ang_inout = r2d * acos( dot_product( in, out ) )
    if( iprint > 0 ) then
      if( nlen <= eps ) then
        write( unit = output_unit, fmt = 15 ) isp, ang_in, ang_out, &
             ang_inout
      else
        write( unit = output_unit, fmt = 16 ) isp, ang_in, ang_out, &
             ang_inout, n(:), r2d*alpha
      endif
    endif
    call cross( sigma, in, out )
    norm = sqrt( dot_product( sigma, sigma ) )
    if( norm < eps ) then
      write( unit = output_unit, fmt = 5 ) isp
      call exit_on_error( srcname )
    endif
    sigma(:) = sigma(:) / norm
    call cross( pi_in, in, sigma )
    norm = sqrt( dot_product( pi_in, pi_in ) )
    if( norm < eps ) then
      write( unit = output_unit, fmt = 5 ) isp
      call exit_on_error( srcname )
    endif
    pi_in(:) = pi_in(:) / norm
    call cross( pi_out, out, sigma )
    norm = sqrt( dot_product( pi_out, pi_out ) )
    if( norm < eps ) then
      write( unit = output_unit, fmt = 5 ) isp
      call exit_on_error( srcname )
    endif
    pi_out(:) = pi_out(:) / norm
    len_e_in = sqrt( real( dot_product( rixs(isp)%e_in, rixs(isp)%e_in ), &
                           kind = REAL64 ) )
    len_e_out = sqrt( real( dot_product( rixs(isp)%e_out, rixs(isp)%e_out ), &
                          kind = REAL64 ) )
    select case(rixs(isp)%channel)
      case ( 's-s ' )
        if( len_e_in < eps ) then
          rixs(isp)%e_in(:) = cmplx( sigma(:), kind = REAL64 )
        else
          rixs(isp)%e_in(:) = rixs(isp)%e_in(:) / len_e_in
        endif
        if( len_e_out < eps ) then
          rixs(isp)%e_out(:) = cmplx( sigma(:), kind = REAL64 )
        else
          rixs(isp)%e_out(:) = rixs(isp)%e_out(:) / len_e_out
        endif
      case ( 's-p ' )
        if( len_e_in < eps ) then
          rixs(isp)%e_in(:) = cmplx( sigma(:), kind = REAL64 )
        else
          rixs(isp)%e_in(:) = rixs(isp)%e_in(:) / len_e_in
        endif
        if( len_e_out < eps ) then
          rixs(isp)%e_out(:) = cmplx( pi_out(:), kind = REAL64 )
        else
          rixs(isp)%e_out(:) = rixs(isp)%e_out(:) / len_e_out
        endif
      case ( 'p-s ' )
        if( len_e_in < eps ) then
          rixs(isp)%e_in(:) = cmplx( pi_in(:), kind = REAL64 )
        else
          rixs(isp)%e_in(:) = rixs(isp)%e_in(:) / len_e_in
        endif
        if( len_e_out < eps ) then
          rixs(isp)%e_out(:) = cmplx( sigma(:), kind = REAL64 )
        else
          rixs(isp)%e_out(:) = rixs(isp)%e_out(:) / len_e_out
        endif
      case ( 'p-p ' )
        if( len_e_in < eps ) then
          rixs(isp)%e_in(:) = cmplx( pi_in(:), kind = REAL64 )
        else
          rixs(isp)%e_in(:) = rixs(isp)%e_in(:) / len_e_in
        endif
        if( len_e_out < eps ) then
          rixs(isp)%e_out(:) = cmplx( pi_out(:), kind = REAL64 )
        else
          rixs(isp)%e_out(:) = rixs(isp)%e_out(:) / len_e_out
        endif
! possible choice for circular polarization
      case ( 'r-l ' )
        rixs(isp)%e_in(:)= ( cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_in(:),kind=REAL64))/sqrt(2._REAL64)
        rixs(isp)%e_out(:)=(-cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_out(:),kind=REAL64))/sqrt(2._REAL64)
      case ( 'r-r ' )
        rixs(isp)%e_in(:)= (cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_in(:),kind=REAL64))/sqrt(2._REAL64)
        rixs(isp)%e_out(:)=(cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_out(:),kind=REAL64))/sqrt(2._REAL64)
      case ( 'l-l ' )
        rixs(isp)%e_in(:)= (-cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_in(:),kind=REAL64))/sqrt(2._REAL64)
        rixs(isp)%e_out(:)=(-cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_out(:),kind=REAL64))/sqrt(2._REAL64)
      case ( 'l-r ' )
        rixs(isp)%e_in(:)= (-cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_in(:),kind=REAL64))/sqrt(2._REAL64)
        rixs(isp)%e_out(:)=( cmplx(sigma(:),kind=REAL64)&
                           +zi*cmplx(pi_out(:),kind=REAL64))/sqrt(2._REAL64)
! if it's required to set the polarization vector manually
      case ( 'man ' )
        call polarization_vectors_for_man(isp)
    end select
  endif
enddo

END SUBROUTINE GEOMETRY

SUBROUTINE POLARIZATION_VECTORS_FOR_MAN(isp)
use, intrinsic :: iso_fortran_env
use m_functions, only : exit_on_error
use m_rixs, only: rixs
implicit none
integer isp
! local vars
character(32), parameter :: srcname=' in POLARIZATION_VECTORS_FOR_MAN'
real(REAL64) norm
real(REAL64), parameter :: eps=1.e-5_REAL64

norm=sqrt(dot_product(rixs(isp)%e_in,rixs(isp)%e_in))
if(norm < eps)then
  write(unit=output_unit,fmt=25)isp,'e_in'
  call exit_on_error( srcname )
else
  rixs(isp)%e_in(:)=rixs(isp)%e_in(:)/norm
endif
norm=sqrt(dot_product(rixs(isp)%e_out,rixs(isp)%e_out))
if(norm < eps)then
  write(unit=output_unit,fmt=25)isp,'e_out'
  call exit_on_error( srcname )
else
  rixs(isp)%e_out(:)=rixs(isp)%e_out(:)/norm
endif

 25 format(/,4x,'Error for SPEC=',i2,': polarization vector ',a &
          ,' equals 0.')

END SUBROUTINE POLARIZATION_VECTORS_FOR_MAN
