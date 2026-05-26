SUBROUTINE KPMESH
use m_bnd
use m_bz
use m_bzmesh
use m_functions
use m_params
use m_rixs, only: nrixs, rixs
use omp_lib
implicit none
character(10), parameter :: srcname=' in KPMESH'
integer ib, iop, is_g4q, ng, k1, k2, k3, iprint_local, isp, ik
real(REAL64) op(3, 3), dpi, start_ref, finish_ref

if (read_rixfile) return
if (check_kstar) return

 10 format(4x, a)

if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Symmetry for RIXS.'
else
  write(unit = output_unit, fmt = "()")
endif
iprint_local = -1
! get matrices for symmetry operations
ng = nopused
call allocate(g, 'g' // srcname, udim1 = 3, udim2 = 3, udim3 = ng)
do iop = 1, ng
  call g4ig(iopnum(iop), op)
  g(:,:, iop) = op(:,:)
enddo
if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'Symmetry for RIXS has been defined.'
  write(unit = output_unit, fmt = 10) 'K-mesh for IBZ.'
endif
! ipq should be allocated outside of kmesh
if (.not. allocated(ipq)) then
  call allocate(ipq, 'ipq' // srcname, udim1 = ndxyz(1), udim2 = ndxyz(2), &
       udim3 = ndxyz(3))
else
  write(unit = output_unit, fmt = 10) 'In KPMESH: ipq array is allocated. Why?'
  call exit_on_error(srcname)
endif
! get ipq and ig4q (set is_g4q = 1) for IBZ, maximal symmetry from LMTO
is_g4q = 1 ! is_g4q = 0 for not referencing ig4q
if (.not. allocated(ig4q)) then
  call allocate(ig4q, 'ig4q' // srcname, udim1 = ndxyz(1), &
          udim2 = ndxyz(2), udim3 = ndxyz(3))
else
  write(unit = output_unit, fmt = 10) 'ig4q array is allocated. Why?'
  call exit_on_error(srcname)
endif
if (.not. allocated(ipqibz)) then
  call allocate(ipqibz, 'ipqibz ' // srcname, udim1 = ndxyz(1), &
       udim2 = ndxyz(2), udim3 = ndxyz(3))
else
  write(unit = output_unit, fmt = 10) 'ipqibz array is allocated. Why?'
  call exit_on_error(srcname)
endif
if (.not. allocated(ig4qibz)) then
  call allocate(ig4qibz, 'ig4qibz' // srcname, udim1 = ndxyz(1), &
       udim2 = ndxyz(2), udim3 = ndxyz(3))
else
  write(unit = output_unit, fmt = 10) 'ig4qibz array is allocated. Why?'
  call exit_on_error(srcname)
endif

call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
       ibzext, is_idold, is_avtet, iprint_local)
ipqibz(:,:,:)  = ipq(:,:,:)
ig4qibz(:,:,:) = ig4q(:,:,:)

! several extra checks
if (nkbz /= ndxyz(1) * ndxyz(2) * ndxyz(3)) then
  write(unit = output_unit, fmt = 10) 'In KPMESH: nkbz != n1 * n2 * n3'
  call exit_on_error(srcname)
endif
if (nkibz /= npnt) then
  write(unit = output_unit, fmt = 10) 'In KPMESH: nkibz != npnt.'
  call exit_on_error(srcname)
endif

call allocate(ik2ibz, 'ik2ibz' // srcname, udim1 = nkbz)
do k3 = 1, ndxyz(3)
  do k2 = 1, ndxyz(2)
    do k1 = 1, ndxyz(1)
      ik = (k3 - 1) * ndxyz(2) * ndxyz(1) + (k2 - 1) * ndxyz(1) + k1
      ik2ibz(ik) = ipq(k1, k2, k3)
    enddo
  enddo
enddo

if (iprint > 0) then
  write(unit = output_unit, fmt = 10) 'K-mesh for IBZ has been set.'
endif

  5 format(4x, 'ERROR in KPMESH: for iop=', i0, ' op. matrix is degenerate.')

call cpu_time(start_ref)
if ( iprint > 0 ) &
  write( unit = output_unit, fmt = "( 4x, a )" ) &
       'Start referencing arrays for all RIXS spectra.'

do isp = 1, nrixs
   rixs(isp)%fnorm = 2._REAL64 * pcellvol * dpi()
   call allocate(rixs(isp)%ik2kq, 'rixs(' // trim(int2string(isp)) // &
        ')%ik2kq' // srcname, udim1 = nkbz)
enddo

call cpu_time(finish_ref)

! check if all q vectors lie on q mesh
call check_q
! get arrays for referencing band energies. ik2kq(:,:) and ik2kqibz(:,:)
! contain k-point numbers for all q-vectors and are used in OpenMP loop.
call get_ik2kq

! Restore data for full symmetry from LMTO. It is used in set_ncf_ncl.f90
is_g4q = 0
    ng = nopused
call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
     ibzext, is_idold, is_avtet, iprint_local)
call set_nbi_nbf  ! set_ncf_ncl is called here
is_g4q = 0
    ng = 1
call kmesh(rbas, ndxyz(1), ndxyz(2), ndxyz(3), g, ng, is_g4q, is_wgt, &
     ibzext, is_idold, is_avtet, iprint_local)
if (iprint > 0) write(unit = output_unit, fmt = 10) 'K-mesh has been defined.'

do isp = 1, nrixs
  call allocate(rixs(isp)%map_region, 'rixs(' // trim(int2string(isp)) // &
       ')%map_region' // srcname, udim1 = nkbz)
enddo
call check_point_in_BZ
call data_region_in_BZ
call allocate(ebz, 'ebz' // srcname, udim1 = nb, udim2 = nkbz)
forall(ib = 1:nb, k1 = 1:nkbz) ebz(ib, k1) = e(ib, ik2ibz(k1))

! call check_sym

END SUBROUTINE KPMESH

SUBROUTINE CHECK_SYM
use, intrinsic :: ISO_FORTRAN_ENV
use m_bnd
implicit none

complex(REAL64), parameter :: z0 = ( 0._REAL64, 0._REAL64 )
complex(REAL64), parameter :: zi = ( 0._REAL64, 1._REAL64 )
complex(REAL64) vec(3), vec2(3), rm(3,3), rm2(3,3), rm_iop(3,3)
integer iop, ix, iy
real m( 3, 3 ), vec_re(3), vec_im(3)
real(REAL64) op(3,3)

! symmetry check: initilaze random matrix and apply symmetry operations
! from iopnum. r_matr is random complex matrix with real numbers on its
! main diagonal
call random_seed()
call random_number( vec_re )
call random_number( vec_im )
vec( : ) = cmplx( vec_re, kind = REAL64 ) + zi * cmplx( vec_im, kind = REAL64 )
call tensor_from_mme( vec, rm )

write( unit = output_unit, fmt = "( /, 4x, a )" ) &
  'The first algorithm for symmetrization'
rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
    rm2( : , : ) = rm2( : , : ) + matmul( transpose(op), matmul( rm, op ) )
  else
    rm2( : , : ) = rm2( : , : ) + &
                   matmul( transpose(op), matmul( conjg(rm), op ) )
  endif
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit, &
       fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit , fmt = "( / )")

write( unit = output_unit , fmt = "( 4x, a )") &
  'The second algorithm for symmetrization'
rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
!$$$    vec2( : ) = matmul( transpose(op) , vec )
    vec2( : ) = matmul( op, vec )
  else
!$$$    vec2( : ) = matmul( transpose(op) , conjg( vec ) )
    vec2( : ) = matmul( op, conjg( vec ) )
  endif
  call tensor_from_mme( vec2, rm_iop )
  rm2( : , : ) = rm2( : , : ) + rm_iop( : , : )
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3 
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit &
       , fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit, fmt = "( / )")

write( unit = output_unit, fmt = "( 4x, a )") &
  'The third algorithm for symmetrization'
call random_number( m )
rm( : , : ) = m( : , : )
call random_number( m )
rm( : , : ) = rm( : , : ) + zi * m( : , : )

rm2( : , : ) = z0
do iop = 1, nopused
  call g4ig( iopnum(iop), op )
  if ( iopnum(iop) > 0 ) then
    rm2( : , : ) = rm2( : , : ) + matmul( transpose(op), matmul( rm, op ) )
  else
    rm2( : , : ) = rm2( : , : ) &
                 + matmul( transpose(op), matmul( conjg(rm), op ) )
  endif
enddo

! normalize and make all the elements with abs(...) < z0 equal to z0.
rm2( : , : ) = rm2( : , : ) / nopused
do ix = 1, 3
  do iy = 1, 3 
    if ( abs( rm2( ix, iy ) ) < 1.e-12_REAL64 ) rm2( ix, iy ) = z0
  enddo
enddo

! find nonzero elements
write( unit = output_unit, fmt = "( 4x, a )", advance = 'no' ) &
  'Absorption tensor: nonzero elements are'
do ix = 1, 3
  do iy = 1, 3
    if ( abs( rm2( ix, iy ) ) > 1.e-12_REAL64 ) write( unit = output_unit &
       , fmt = "( ' ', i1, i1 )", advance = 'no' ) ix, iy
  enddo
enddo
write( unit = output_unit, fmt = "( / )")

END SUBROUTINE CHECK_SYM
