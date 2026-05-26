SUBROUTINE GEOMETRY
use, intrinsic :: iso_fortran_env
use m_files
use m_functions, only: exit_on_error
use m_params
use m_rixs
implicit none
! local vars
character(*), parameter :: srcname=' in GEOMETRY'
integer isp
real(REAL64) sigma(3), pi_in(3), pi_out(3), q(3), alpha, inv_R(3, 3)

if (read_rixfile) return
if (check_kstar) return

 10 format(4x, 'GEOMETRICAL DETAILS (angles and rotation axis)', /, &
           4x, 'ISPEC', 11x, 'rotational axis', 12x, 'angle', &
           4x, '(q,k_in)', 4x, '(q,k_out)', 4x, '(k_in,k_out)')

call get_q_in_cartesian_frame

if (irixs == 0 .or. irixs == 2) then
! Joint DOS or average over polarizations
! e_in and e_out are these variables really needed?
  do isp = 1, nrixs
    rixs(isp)%e_in(:)  = sqrt(1._REAL64/3._REAL64)
    rixs(isp)%e_out(:) = sqrt(1._REAL64/3._REAL64)
    rixs(isp)%k_out(:) = k_in(:)   !  rixs(isp)%k_in(:)
  enddo
  return
endif

write(unit = output_unit, fmt = '(4x, a)') 'LMTO cartesian coordinate system is used'
if (iprint > 0) write(unit = output_unit, fmt = 10)
do isp = 1, nrixs
  q(:) = rixs(isp)%q(:)
  if (norm2(q) < eps5) then
! for q=0 it's not possible to build scattering plane. If one wants to
! calculate spectra for q=0 it's necessary to specify polarization vectors
! manually
    rixs(isp)%channel(:) = 'user'
    call manual_polarization_vectors(isp)
    rixs(isp)%k_in(:)  = k_in(:)
    rixs(isp)%k_out(:) = k_in(:)
  else
! Find components of input vector in local Cartesian frame
! Calculate k_out in such a way that abs(klocal)=abs(k_out)
! in and out are local normalized copies of rixs(isp)%k_in and rixs(isp)%k_out
    rixs(isp)%k_in(:) = k_in(:)
    call get_rotation_angle_and_matrix(isp, alpha, inv_R)
    rixs(isp)%k_in(:) = matmul(inv_R, k_in)
    rixs(isp)%k_out(:) = rixs(isp)%k_in(:) - 2._REAL64 * q(:) * &
                         dot_product(rixs(isp)%k_in, q) / norm2(q) ** 2
    call get_basic_polarization_vectors(isp, sigma, pi_in, pi_out)
    call set_polarization_vectors(isp, sigma, pi_in, pi_out)
  endif
enddo

END SUBROUTINE GEOMETRY

SUBROUTINE get_q_in_cartesian_frame

use, intrinsic :: iso_fortran_env
use m_constants, only: eps5
use m_params, only: iprint, irixs
use m_rixs, only: nrixs, rixs
use m_sdt, only: ut1
implicit none
integer isp

if (iprint < 0) write(unit = output_unit, fmt = 7)
! Transform q vectors into local Cartesian system. In inputfile
! q vector is set as a fractions of G1, G2, G3 in units of 2pi/a, 2pi/b,
! 2pi/c. After matmul q is expressed in units of 2pi/a. So I need to apply
! 2pi multiplier while calculating RIXS mme.
do isp = 1, nrixs
  rixs(isp)%q(:) = matmul(ut1, rixs(isp)%q)
  if (norm2(rixs(isp)%q) < eps5 .and. irixs == 1) then
    write(unit = output_unit, fmt = 8) isp
  endif
enddo
write(unit = output_unit, fmt = 7)

 7 format('')
 8 format(4x, 'WARNING: for ISPEC=', i0, ' q=0 so polarization vectors ', &
        'must be set manually.')

END SUBROUTINE get_q_in_cartesian_frame

SUBROUTINE get_rotation_angle_and_matrix(isp, alpha, inv_R)
use, intrinsic :: iso_fortran_env
use m_constants, only: eps5, r2d
use m_functions, only: exit_on_error
use m_params
use m_rixs
implicit none
character(*), parameter :: srcname=' in GEOMETRY'
integer i, isp
real(REAL64) q(3), n(3), R(3, 3), nDual(3, 3), alpha, inv_R(3, 3), det

q(:) = rixs(isp)%q
if (norm2(q) > eps5) q(:) = q(:) / norm2(q)
! q_global(:) = matmul(A_bas, rixs(isp)%q)
! if (norm2(q_global) > eps5) q_global(:) = q_global(:) / norm2(q_global)
! call cross(n, q_global, q_axis)
call cross(n, q, q_axis)
if (norm2(n) > eps5) then
!  write(unit = output_unit, fmt = 10) 'norm2(n) = 0.'
!  call exit_on_error(srcname)
!else
  n(:) = n(:) / norm2(n)
endif
! defaults: rotation tensor, tensor dual to n(:), and rotation angle
R(:,:) = 0._REAL64
forall (i = 1:3) R(i, i) = 1._REAL64
nDual(:,:)   = 0._REAL64
nDual(1, 2)  = -n(3)
nDual(1, 3)  =  n(2)
nDual(2, 3)  = -n(1)
nDual(:, :)  =  nDual(:,:) - transpose(nDual)
! alpha        =  acos(dot_product(q_global, q_axis))
alpha        =  acos(dot_product(q, q_axis))
R(:,:)       =  R(:,:) + nDual(:,:) * sin(alpha) + &
                matmul(nDual, nDual) * (1 - cos(alpha))
! R_A_bas(:,:) =  matmul(R, A_bas)
! call dinv33(R_A_bas, 0, inv_R_A_bas, det)
call dinv33(R, 0, inv_R, det)
if (det < 1._REAL64 - eps5) then
  write(unit = output_unit, fmt = 11) isp, det
  call exit_on_error(srcname)
endif
if (iprint > 0) then
  write(unit = output_unit, fmt = 16, advance = 'no') isp, n(:), r2d * alpha
endif

 10 format(4x, a)
 11 format(4x, 'Error: for ISPEC=', i2, &
          ' rotational matrix is not orthogonal', ", it's det=", f8.5)
 16 format(4x, i3, 4x, '(', 3f10.6, ')', f9.3)

END SUBROUTINE get_rotation_angle_and_matrix

SUBROUTINE get_basic_polarization_vectors(isp, sigma, pi_in, pi_out)
use, intrinsic :: iso_fortran_env
use m_constants, only: eps5, r2d
use m_functions, only: exit_on_error
use m_params, only: iprint, k_in
use m_rixs
implicit none
character(*), parameter :: srcname=' in GEOMETRY'
integer isp
real(REAL64) sigma(3), pi_in(3), pi_out(3), in(3), out(3), q(3), ang_in
real(REAL64) ang_out, ang_inout

q(:) = rixs(isp)%q(:)
in(:) = rixs(isp)%k_in(:) / norm2(rixs(isp)%k_in)
! Why do I need to calculate rixs(isp)%theta_k_in and rixs(isp)%phi_k_in?
rixs(isp)%theta_k_in = acos(in(3))
if (abs(k_in(1)) < eps5 .and. abs(k_in(2)) < eps5) then
  rixs(isp)%phi_k_in = 0._REAL64
else
  rixs(isp)%phi_k_in = atan2(in(2), in(1))
endif
out(:)    = rixs(isp)%k_out(:) / norm2(rixs(isp)%k_out)
ang_in    = r2d * acos(dot_product(in, q)  / norm2(q))
ang_out   = r2d * acos(dot_product(out, q) / norm2(q))
ang_inout = r2d * acos(dot_product(in, out))
if (iprint > 0) then
  write(unit = output_unit, fmt = 16) ang_in, ang_out, ang_inout
endif
call cross(sigma, in, out)
if (norm2(sigma) < eps5) then
  write(unit = output_unit, fmt = 5) isp
  call exit_on_error(srcname)
else
  sigma(:) = sigma(:) / norm2(sigma)
endif
call cross(pi_in, in, sigma)
if (norm2(pi_in) < eps5) then
  write(unit = output_unit, fmt = 5) isp
  call exit_on_error(srcname)
else
  pi_in(:) = pi_in(:) / norm2(pi_in)
endif
call cross(pi_out, out, sigma)
if (norm2(pi_out) < eps5) then
  write(unit = output_unit, fmt = 5) isp
  call exit_on_error(srcname)
else
  pi_out(:) = pi_out(:) / norm2(pi_out)
endif
! Here I try to get all the components positive
pi_in(:) = -pi_in(:)
pi_out(:) = -pi_out(:)
sigma(:) = - sigma(:)

  5 format(4x, 'Can not build scattering plane for SPEC=', i0, &
           /, 4x, 'Vectors k_(in,out) and q are collinear.')
 16 format(4x, f8.4, 4x, f8.4, 7x, f8.4)

END SUBROUTINE get_basic_polarization_vectors

SUBROUTINE set_polarization_vectors(isp, sigma, pi_in, pi_out)
use, intrinsic :: iso_fortran_env
use m_constants, only: zi
use m_rixs, only: rixs
implicit none
complex(REAL64) csigma(3), cpi_in(3), cpi_out(3)
integer isp
real(REAL64) sigma(3), pi_in(3), pi_out(3)

csigma(:)  = cmplx(sigma(:),  kind = REAL64) ! complex sigma
cpi_in(:)  = cmplx(pi_in(:),  kind = REAL64) ! complex pi_in
cpi_out(:) = cmplx(pi_out(:), kind = REAL64) ! complex pi_out
if (rixs(isp)%channel(:) == 's-s ') then
  rixs(isp)%e_in(:)  = csigma(:)
  rixs(isp)%e_out(:) = csigma(:)
  return
endif
if (rixs(isp)%channel(:) == 's-p ') then
  rixs(isp)%e_in(:)  = csigma(:)
  rixs(isp)%e_out(:) = cpi_out(:)
  return
endif
if (rixs(isp)%channel(:) == 'p-s ') then
  rixs(isp)%e_in(:)  = cpi_in(:)
  rixs(isp)%e_out(:) = csigma(:)
  return
endif
if (rixs(isp)%channel(:) == 'p-p ') then
  rixs(isp)%e_in(:)  = cpi_in(:)
  rixs(isp)%e_out(:) = cpi_out(:)
  return
endif
csigma(:)  = csigma(:)  / sqrt(2._REAL64)
cpi_in(:)  = cpi_in(:)  / sqrt(2._REAL64)
cpi_out(:) = cpi_out(:) / sqrt(2._REAL64)
if (rixs(isp)%channel(:) == 'r-l ') then
  rixs(isp)%e_in(:)  =  csigma(:) + zi * cpi_in(:)
  rixs(isp)%e_out(:) = -csigma(:) + zi * cpi_out(:)
  return
endif
if (rixs(isp)%channel(:) == 'r-r ') then
  rixs(isp)%e_in(:)  = csigma(:) + zi * cpi_in(:)
  rixs(isp)%e_out(:) = csigma(:) + zi * cpi_out(:)
  return
endif
if (rixs(isp)%channel(:) == 'l-l ') then
  rixs(isp)%e_in(:)  = -csigma(:) + zi * cpi_in(:)
  rixs(isp)%e_out(:) = -csigma(:) + zi * cpi_out(:)
  return
endif
if (rixs(isp)%channel(:) == 'l-r ') then
  rixs(isp)%e_in(:)  = -csigma(:) + zi * cpi_in(:)
  rixs(isp)%e_out(:) =  csigma(:) + zi * cpi_out(:)
  return
endif
if (rixs(isp)%channel(:) == 'user') then
  call manual_polarization_vectors(isp)
endif

! write(unit = output_unit, fmt = 20) 'For ISPEC=', isp, &
!      ' unknown scattering channel: ' // trim(rixs(isp)%channel)
! 20 format(4x, i2, a)

END SUBROUTINE set_polarization_vectors

SUBROUTINE manual_polarization_vectors(isp)
use, intrinsic :: iso_fortran_env
use m_functions, only : exit_on_error
use m_constants, only: eps5
use m_rixs, only: rixs
implicit none
integer isp
! local vars
character(31), parameter :: srcname=' in MANUAL_POLARIZATION_VECTORS'
real(REAL64) norm

norm = sqrt(dot_product(rixs(isp)%e_in, rixs(isp)%e_in))
if (norm < eps5) then
  write(unit = output_unit, fmt = 25) isp, 'e_in'
  call exit_on_error(srcname)
endif
rixs(isp)%e_in(:) = rixs(isp)%e_in(:) / norm
norm = sqrt(dot_product(rixs(isp)%e_out, rixs(isp)%e_out))
if (norm < eps5) then
  write(unit = output_unit, fmt = 25) isp, 'e_out'
  call exit_on_error(srcname)
endif
rixs(isp)%e_out(:) = rixs(isp)%e_out(:) / norm

25 format(/, 4x, 'Error for SPEC=', i2, ': polarization vector ', a, &
         ' equals 0.')

END SUBROUTINE manual_polarization_vectors

