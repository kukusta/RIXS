SUBROUTINE PRINTDATA
use, intrinsic :: iso_fortran_env
use m_aux, only: map
use m_bnd, only: ef
! use m_bzmesh, only: nkbz, ntbz
use m_files
use m_functions
use m_params
use m_rixs
use m_sdt, only: ut
use units, only: dry2ev
implicit none
! local vars
character(5) where_states_lie
integer iat,isp
real(REAL64) min_esp, max_esp
real(REAL64), external :: dpi

if (read_rixfile) return
if (check_kstar ) return

if (iabsorp /= 0) then
  write(unit = output_unit, fmt = 5) wmin, wmax, dw, nw, emin, emax, de, &
       int((emax - emin) / de + 1)
else
  write(unit = output_unit, fmt = 6) wmin, wmax, dw, nw
endif

5 format(/, 4x, 'ENERGY RANGE', &
         /, 4x, 'wmin= ', f8.5, ' Ry, wmax= ', f8.5, ' Ry, dw= ', f8.5, &
                ' Ry, nw= ', i0, &
         /, 4x, 'emin= ', f8.5, ' Ry, emax= ', f8.5, ' Ry, de= ', f8.5, &
                ' Ry, ne= ', i0)
6 format(/, 4x, 'ENERGY RANGE', &
         /, 4x, 'wmin= ', f8.5, ' Ry, wmax= ', f8.5, ' Ry, dw= ', f8.5, &
                ' Ry, nw= ', i0)

if (irixs == 0) then
  write(unit = output_unit, fmt = 10) 'JOINT DOS SPECTRA'
elseif (irixs == 2) then
  write(unit = output_unit, fmt = 10) 'RIXS SPECTRA WITH MME AVERAGED' // &
       ' OVER POLARIZATIONS'
elseif (irixs == 3) then
  write(unit = output_unit, fmt = 10) 'RIXS LORENTZIAN'
else
  write(unit = output_unit, fmt = 10)'RIXS SPECTRA'
endif  
do isp = 1, nrixs
  call get_min_max_esp(isp, min_esp, max_esp)
  if (irixs == 0 .or. irixs == 10) then
    write(unit = output_unit, fmt = 14) isp, rixs(isp)%txtel(:2),    &
         rixs(isp)%sname(:), matmul(ut,rixs(isp)%q), rixs(isp)%q(:), &
         rixs(isp)%nbi(:), rixs(isp)%nbf(:)
  else if (irixs == 1) then
    write(unit = output_unit, fmt = 15) isp, rixs(isp)%txtel(:2),           &
         rixs(isp)%sname(:), rixs(isp)%channel(:), matmul(ut, rixs(isp)%q), &
         rixs(isp)%q(:), rixs(isp)%en, rixs(isp)%en*dry2ev,                 &
         rixs(isp)%Gamma * dry2ev, dry2ev*(rixs(isp)%en + min_esp - ef),    &
         dry2ev*(rixs(isp)%en + max_esp - ef),                              &
         rixs(isp)%nbi(:), rixs(isp)%nbf(:), rixs(isp)%Reg_BZ(:),           &
         rixs(isp)%Reg_BZ_states(:), rixs(isp)%Reg_BZ_size(:),              &
         rixs(isp)%Reg_BZ_center(:)
    if (iprint > 10) call print_phases(isp)
  elseif (irixs == 2) then
    write(unit = output_unit, fmt = 17) isp, rixs(isp)%txtel(:2),        &
         rixs(isp)%sname(:), rixs(isp)%en, rixs(isp)%en*dry2ev,          &
         rixs(isp)%Gamma * dry2ev, dry2ev*(rixs(isp)%en + min_esp - ef), &
         dry2ev*(rixs(isp)%en + max_esp - ef), matmul(ut,rixs(isp)%q),   &
         rixs(isp)%q(:), rixs(isp)%nbi(:), rixs(isp)%nbf(:)
    if (iprint > 10) call print_phases(isp)
  elseif(irixs == 3) then
    write(unit = output_unit, fmt = 16) isp, rixs(isp)%txtel(:2),        &
         rixs(isp)%sname(:), rixs(isp)%en, rixs(isp)%en*dry2ev,          &
         rixs(isp)%Gamma * dry2ev, dry2ev*(rixs(isp)%en + min_esp - ef), &
         dry2ev*(rixs(isp)%en + max_esp - ef), rixs(isp)%nbi(:),         &
         rixs(isp)%nbf(:)
  endif
  if(iprint > 0 .and. abs(irixs) == 1)then
    write(unit=output_unit,fmt=22) rixs(isp)%k_in(:) &
         /sqrt(dot_product(rixs(isp)%k_in, rixs(isp)%k_in)) &
         ,rixs(isp)%theta_k_in * r2d, rixs(isp)%theta_k_in * r2d &
         ,rixs(isp)%k_out(:)&
         /sqrt(dot_product(rixs(isp)%k_out,rixs(isp)%k_out))
    select case(rixs(isp)%channel)
      case('s-s ','s-p ','p-s ','p-p ','user')
        write(unit=output_unit,fmt=25)real(rixs(isp)%e_in(:))&
             ,real(rixs(isp)%e_out(:))
      case('r-l ','r-r ','l-r ','l-l ')
        write(unit=output_unit,fmt=30)rixs(isp)%e_in(:),rixs(isp)%e_out(:)
      case default
        write(unit=output_unit,fmt=33)rixs(isp)%channel,isp
        call deallocate_global_arrays
        stop
    end select
  endif
enddo

 10 format(/, 4x, a)
 14 format(4x, 'ISPEC=', i2, ": atom= '", a, "'", 4x, "edge= '", a, "'", &
           /, 14x, 'q= (', 3f10.5, ' ) in G1, G2, G3', &
           /, 14x, 'q= (', 3f10.5, ' ) in 2pi/a in Cart. CS', &
           /, 14x, 'nbi(1:2)= ', 2(i0, 2x), ' nbf(1:2)= ', 2(i0, 2x))
 15 format(4x, 'ISPEC=', i2, ": atom= '", a, "'", 4x, "edge= '", a, "'", 4x, &
           "channel='", a, "'", &
           /, 14x, 'q= (', 3f10.5, ' ) in G1, G2, G3', &
           /, 14x, 'q= (', 3f10.5, ' ) in 2pi/a in Cart. CS', &
           /, 14x, 'input_energy= ', f9.4, ' Ry = ', f11.4, ' eV', &
               4x, 'Gamma= ', f6.4, ' eV', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (min)', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (max)', &
           /, 14x, 'nbi(1:2)= ', 2(i0, 2x), ' nbf(1:2)= ', 2(i0, 2x), &
           /, 14x, "reg_BZ= '", a, "'", 4x, "reg_BZ_states= '", &
                a, "'", &
           /, 14x, 'reg_BZ_size(:)= (', 3i4, ' )', &
           /, 14x, 'reg_BZ_center(:)= (', 3f11.6, ")")
 16 format(4x, 'ISPEC=', i2, ": atom= '", a, "'", 4x, "edge= '", a, "'", &
           /, 14x, 'input_energy= ', f9.4, ' Ry = ', f11.4, ' eV', &
               4x, 'Gamma= ', f6.4, ' eV', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (min)', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (max)', &
           /, 14x, 'nbi(1:2)= ', 2(i0, 2x), ' nbf(1:2)= ', 2(i0, 2x))
 17 format(4x, 'ISPEC=', i2, ": atom= '", a, "'", 4x, "edge= '", a, "'", 4x, &
           /, 14x, 'input_energy= ', f9.4, ' Ry = ', f11.4, ' eV', &
               4x, 'Gamma= ', f6.4, ' eV', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (min)', &
           /, 14x, 'input_energy relative to Ef= ', f9.4, ' eV (max)', &
           /, 14x, 'q= (', 3f10.5, ' ) in G1, G2, G3', &
           /, 14x, 'q= (', 3f10.5, ' ) in 2pi/a in Cart. CS', &
           /, 14x, 'nbi(1:2)= ', 2(i0, 2x), ' nbf(1:2)= ', 2(i0, 2x))
 18 format(4x, 'Unknown irixs=', i0, ' value.')
 ! 20 format(13x,'ng=',i2,' from ',i2,4x,'ntrbz= ',i0,' from ',i0&
 !          ,4x,'nkrbz= ',i0,' from ',i0,/,13x,'iopnum= ',48(i0,2x))
 22 format(14x,'input  k-vector= (',3f10.6,')',/ &
          ,14x,'for input k-vector theta=',f10.5,' deg, phi=',f10.5,' deg',/ &
          ,14x,'output k-vector= (',3f10.6,')')
 25 format(14x,'input  polarization= ','(',3f10.6,')',/&
          ,14x,'output polarization= ','(',3f10.6,')')
 30 format(14x,'input  polarization= ',2('(',3f10.6,')'),/&
          ,14x,'output polarization= ',2('(',3f10.6,')'))
 33 format(4x, "Error: unknow scattering channel '", a, "' for ispec=", i0)
 40 format(/, 4x, 'FILES')

write(unit = output_unit, fmt = 40) ! ,(achar(9),i=1,1)
write(unit = output_unit, fmt = 45) trim(inptfile), trim(sdtfile)
write(unit = output_unit, fmt = 50) trim(rixsfile), trim(bndfile)
if ((irixs == 0 .or. abs(irixs) == 10) .and. iabsorp == 0) then
  write(unit = output_unit, fmt = 55) trim(datafile)
  return
else
  write(unit = output_unit, fmt = 56) trim(datafile), trim(mmefile)
endif

 45 format(4x,"inptfile= '",a,"'",4x,"sdtfile= '",a,"'")
 50 format(4x,"rixsfile= '",a,"'",4x,"bndfile= '",a,"'")
 55 format(4x,"datafile= '",a,"'")
 56 format(4x, "datafile= '", a, "'", /, 4x, "mmefile= '", a, "'")

write(unit = output_unit, fmt = 10) 'STATES'
where_states_lie = 'below'
if (switch_ini_final) where_states_lie = 'above'
write(unit = output_unit, fmt = '(4x, a)') 'States with k+q lie ' // &
     where_states_lie // ' Fermi level'

write(unit = output_unit, fmt = 10) 'LMTO DATA USED IN RIXS'
do isp = 1, nlmtdata
  if (sum( map(lmtdata(isp)%ia1, : ) ) == 0) cycle
  write(unit = output_unit, fmt = 65) isp, lmtdata(isp)%sname(:), &
       lmtdata(isp)%txtel(:), (ef - lmtdata(isp)%esp(lmtdata(isp)%ian)) * dry2ev
  if (iprint > 0) then
    write(unit = output_unit, fmt = 70) lmtdata(isp)%isort, lmtdata(isp)%ia1, &
         lmtdata(isp)%ian, lmtdata(isp)%iat1, lmtdata(isp)%iatn, &
         lmtdata(isp)%label
    write(unit = output_unit, fmt = 75) lmtdata(isp)%ec, &
         lbound(lmtdata(isp)%esp), ubound(lmtdata(isp)%esp), &
         lmtdata(isp)%esp(:)
    do iat = lmtdata(isp)%iat1, lmtdata(isp)%iatn
      write(unit = output_unit, fmt = 80) iat, lmtdata(isp)%rat(:,iat)
    enddo
  endif
enddo

 65 format(4x, 'ISPEC=', i2, ": edge='", a, "'", 4x, "txtel='", a, "'", 4x, &
          'w_edge= ', f13.6, ' eV')
 70 format(14x, 'isort=', i0, 4x, 'ia1=', i0, 4x, 'ian=', i0, 4x, 'iat1=', i0, &
          4x, 'iatn=', i0, 4x, 'label=', i0)
 75 format( 14x, 'ec= ', f12.8, ' Ry', &
          /,14x, 'esp(', i0, ':', i0, ')=', 20(2x, f9.5,' Ry'))
 80 format(14x, 'rat(1:3,', i0, ')=', 3(4x, f12.8))

END SUBROUTINE PRINTDATA

SUBROUTINE PRINT_PHASES(isp)
use, intrinsic :: iso_fortran_env
use m_params
use m_rixs
implicit none
integer isp, iat, isplmt

write(unit = output_unit, fmt = 18) 'Phases (in 2*pi)'
do isplmt = 1, nlmtdata
  if (rixs(isp)%txtel(:2) == lmtdata(isplmt)%txtel(:2) .and. &
      rixs(isp)%sname(:) == lmtdata(isplmt)%sname(:)) then
    if (lmtdata(isplmt)%iatn - lmtdata(isplmt)%iat1 > 0) then
      do iat = lmtdata(isplmt)%iat1, lmtdata(isplmt)%iatn, 2
        write(unit = output_unit, fmt = 19) iat, &
             dot_product(lmtdata(isplmt)%rat(:, iat), rixs(isp)%q), &
             iat + 1, dot_product(lmtdata(isplmt)%rat(:, iat + 1), rixs(isp)%q)
      enddo
    else
!     one atom of this sort
      write(unit = output_unit, fmt = 20) lmtdata(isplmt)%iat1, &
           dot_product(lmtdata(isplmt)%rat(:,lmtdata(isplmt)%iat1),rixs(isp)%q)
    endif
  endif
enddo

 18 format(14x, a)
 19 format(10x, 2(4x, 'iat=', i4, ' phase= ', f12.8))
 20 format(10x, 4x, 'iat=', i4, ' phase= ', f12.8)

END SUBROUTINE PRINT_PHASES

SUBROUTINE get_min_max_esp(isp, min_esp, max_esp)
  use m_aux
  use m_rixs, only: lmtdata, n_mme
  implicit none
  integer ia, isp
  real(REAL64) esp_ia, min_esp, max_esp

  min_esp = 0.
  max_esp = -huge(0.)
  do ia = 1, n_mme
    if (map(ia, isp) == 0) cycle
    esp_ia = lmtdata(lmtindex(ia))%esp(ia)
    if (esp_ia > max_esp) max_esp = esp_ia
    if (esp_ia < min_esp) min_esp = esp_ia
  enddo
  
END SUBROUTINE get_min_max_esp
