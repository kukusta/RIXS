SUBROUTINE SET_ABSORP_IBZ(ia, kpk)
use m_aux, only: lmtindex
use m_bnd
use m_bz
use m_bzmesh, only: nkibz
use m_functions
use m_params
use m_sdt, only: igtta
use m_rixs
implicit none
integer, intent(in) :: ia
real(REAL64) kpk(nb, nkibz)
! local vars
character(18), parameter :: srcname = ' in SET_ABSORP_IBZ'
complex(REAL64), allocatable :: Isym(:, :, :)
complex(REAL64) mmek(3), mmek2(3), I(3, 3)
integer ib0, ib, k, ind, iop, iat, iat1, iatn
logical rot2qph
real(REAL64) ang(3), op(3, 3), opqph(3, 3), theta, phi

 1 format(4x, a, i0, a)

if (ia < 1 .or. ia > n_mme) then
  write(unit = output_unit, fmt = 1)'In SET_ABSORP_BZ: wrong ia=', ia, ' value.'
  call exit_on_error(srcname)
endif
iat1 = lmtdata(lmtindex(ia))%iat1
iatn = lmtdata(lmtindex(ia))%iatn
theta = rixs(lmtindex(ia))%theta_k_in
phi = rixs(lmtindex(ia))%phi_k_in
rot2qph = sqrt(theta**2 + phi**2) > eps2
if (rot2qph) then
  ang(:) = (/ 0._REAL64, -theta, -phi /)
  call turnm(ang, opqph)
endif
call allocate(Isym, 'Isym' // srcname, udim1 = 3, udim2 = 3, ldim3 = iat1, &
     udim3 = iatn)
! now always ib0 = nbmin - 1 = 0. It's just kept for future as possible way
! of generalization
ib0 = 0
do k = 1, nkibz                    ! over all k from IBZ
  do ib = 1, nb                    ! over all bands
    Isym(:, :, :) = z0
    do iat = iat1, iatn
      mmek(:) = lmtdata(lmtindex(ia))%mme(:, iat, ia, ib, k)
! symmetrize atomic tensor in global system
      do iop = 1, nopused
        ind = igtta(iop, iat)
        op(:, :) = g(:, :, iop)
!$$$        mmek2(:) = mme_trans( iop , mmek )
        if (iop == 1) then
          mmek2(:) = mmek(:)
        else
          if (iopnum(iop) > 0) then
            mmek2(:) = matmul(op, mmek)
          else
            mmek2(:) = matmul(op, conjg(mmek))
          endif
        endif
        call tensor_from_mme(mmek2, I)
! rotate tensor to photon local system
        if(rot2qph) I(:,:) = matmul(opqph, matmul(I, transpose(opqph)))
        Isym(:, :, ind) = Isym(:, :, ind) + I(:, :)
      enddo
    enddo
    Isym(:, :, :) = Isym(:, :, :) / cmplx(nopused * 8 * nfu, kind = REAL64)
    kpk(ib0 + ib, k) = real( sum(Isym(1, 1, iat1 : iatn)) &
                           + sum(Isym(2, 2, iat1 : iatn)), kind = REAL64)
  enddo
enddo

END SUBROUTINE SET_ABSORP_IBZ

