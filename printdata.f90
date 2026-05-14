SUBROUTINE PRINTDATA
use, intrinsic :: iso_fortran_env
use m_bnd, only: ef, nopused
use m_bzmesh, only: nkbz, ntbz
use m_files
use m_functions
use m_params
use m_rixs
use m_sdt, only: ut
use units, only: dry2ev
implicit none
! local vars
integer iat,isp
real(REAL64), external :: dpi

if ( check_kstar .and. debug_mode ) return
if ( read_rixfile ) return

if(iabsorp /= 0)then
  write(unit=output_unit,fmt=5)wmin,wmax,dw,nw,emin,emax,de&
       ,int((emax-emin)/de+1)
else
  write(unit=output_unit,fmt=6)wmin,wmax,dw,nw
endif

  5 format(/,4x,'ENERGY RANGE',/,4x,'wmin= ',f0.8,' Ry, wmax= '&
          ,f0.8,' Ry, dw= ',f0.8,' Ry, nw= ',i0,/,4x,'emin= ',f0.8&
          ,' Ry, emax= ',f0.8,' Ry, de= ',f0.8,' Ry, ne= ',i0)
  6 format(/,4x,'ENERGY RANGE',/,4x,'wmin= ',f0.8,' Ry, wmax= '&
          ,f0.8,' Ry, dw= ',f0.8,' Ry, nw= ',i0)

select case(abs(irixs))
  case(0)
    write(unit=output_unit,fmt=10)'JOINT DOS SPECTRA'
  case(2)
    write(unit=output_unit,fmt=10)&
         'RIXS SPECTRA WITH MME AVERAGED OVER POLARIZATIONS'
  case(3)
    write(unit=output_unit,fmt=10)'RIXS LORENTZIAN'
  case default
    write(unit=output_unit,fmt=10)'RIXS SPECTRA'
end select

do isp=1,nrixs
  select case(abs(irixs))
    case(0,10)
      write(unit=output_unit,fmt=14)isp,rixs(isp)%txtel(:2)&
           ,rixs(isp)%sname(:),matmul(ut,rixs(isp)%q),rixs(isp)%q(:)&
           ,sqrt(dot_product(rixs(isp)%q(:),rixs(isp)%q(:)))&
           ,rixs(isp)%nbi(:),rixs(isp)%nbf(:)
    case(1)
      write(unit=output_unit,fmt=15)isp,rixs(isp)%txtel(:2)&
           ,rixs(isp)%sname(:),rixs(isp)%channel(:),matmul(ut,rixs(isp)%q),&
            rixs(isp)%q(:)&
           ,rixs(isp)%en,rixs(isp)%en*dry2ev&
           ,dry2ev*(rixs(isp)%en+rixs(isp)%ec-ef)&
           ,sqrt(dot_product(rixs(isp)%q,rixs(isp)%q))&
           ,rixs(isp)%nbi(:),rixs(isp)%nbf(:)
      call print_phases(isp)
    case(2)
      write(unit=output_unit,fmt=17)isp,rixs(isp)%txtel(:2)&
           ,rixs(isp)%sname(:),rixs(isp)%en,rixs(isp)%en*dry2ev&
           ,dry2ev*(rixs(isp)%en+rixs(isp)%ec-ef),matmul(ut,rixs(isp)%q),&
            rixs(isp)%q(:)&
           ,sqrt(dot_product(rixs(isp)%q,rixs(isp)%q))&
           ,rixs(isp)%nbi(:),rixs(isp)%nbf(:)
      call print_phases(isp)
    case(3)
      write(unit=output_unit,fmt=16)isp,rixs(isp)%txtel(:2)&
           ,rixs(isp)%sname(:),rixs(isp)%en,rixs(isp)%en*dry2ev&
           ,dry2ev*(rixs(isp)%en+rixs(isp)%ec-ef)&
           ,rixs(isp)%nbi(:),rixs(isp)%nbf(:)
    case default
      write(unit=output_unit,fmt=18)irixs
      call deallocate_global_arrays
      stop
  end select
  write(unit=output_unit,fmt=20)rixs(isp)%ng,nopused,rixs(isp)%ntrbz&
       ,ntbz,rixs(isp)%nkrbz,nkbz,rixs(isp)%iopnum(:rixs(isp)%ng)

  if(iprint > 0 .and. abs(irixs) == 1)then
    write(unit=output_unit,fmt=22)rixs(isp)%k_in(:)&
         /sqrt(dot_product(rixs(isp)%k_in,rixs(isp)%k_in))&
         ,rixs(isp)%k_out(:)&
         /sqrt(dot_product(rixs(isp)%k_out,rixs(isp)%k_out))
    select case(rixs(isp)%channel)
      case('s-s ','s-p ','p-s ','p-p ','man ')
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

 10 format(/,4x,a)
 14 format(4x,'SPEC=',i2,": atom= '",a,"'",4x,"edge= '",a,"'"&
          ,/,13x,'q= (',3f17.12,' ) in r. l. u.'&
          ,/,13x,'q= (',3f17.12,' ) in 2pi/a'&
          ,/,13x,'qlen=',f12.7,4x,'nbi(1:2)= ',2(i0,2x),' nbf(1:2)= ',2(i0,2x))
 15 format(4x,'SPEC=',i2,": atom= '",a,"'",4x,"edge= '",a,"'",4x&
          ,"channel='",a,"'"&
          ,/,13x,'q= (',3f17.12,' ) in r. l. u.'&
          ,/,13x,'q= (',3f17.12,' ) in 2pi/a'&
          ,/,13x,'input_energy= ',f0.7,' Ry = ',f0.7,' eV'&
          ,/,13x,'input_energy relative to Ef= ',f8.4,' eV'&
          ,/,13x,'qlen=',f12.7,4x,'nbi(1:2)= ',2(i0,2x),' nbf(1:2)= ',2(i0,2x))
 16 format(4x,'SPEC=',i2,": atom= '",a,"'",4x,"edge= '",a,"'"&
          ,/,13x,'input_energy= ',f0.7,' Ry = ',f0.7,' eV'&
          ,/,13x,'input_energy relative to Ef= ',f8.4,' eV'&
          ,/,13x,'nbi(1:2)= ',2(i0,2x),' nbf(1:2)= ',2(i0,2x))
 17 format(4x,'SPEC=',i2,": atom= '",a,"'",4x,"edge= '",a,"'",4x&
          ,/,13x,'input_energy= ',f0.7,' Ry = ',f0.7,' eV'&
          ,/,13x,'input_energy relative to Ef= ',f8.4,' eV'&
          ,/,13x,'q= (',3f17.12,' ) in r. l. u.'&
          ,/,13x,'q= (',3f17.12,' ) in 2pi/a'&
          ,/,13x,'qlen=',f12.7,4x,'nbi(1:2)= ',2(i0,2x),' nbf(1:2)= '&
          ,2(i0,2x))
 18 format(4x,'Unknown irixs=',i0,' value.')
 20 format(13x,'ng=',i2,' from ',i2,4x,'ntrbz= ',i0,' from ',i0&
          ,4x,'nkrbz= ',i0,' from ',i0,/,13x,'iopnum= ',48(i0,2x))
 22 format(13x,'input  k-vector (',2(f15.12,','),f15.12,')',/&
          ,13x,'output k-vector (',2(f15.12,','),f15.12,')')
 25 format(13x,'input  polarization ','(',2(f15.12,','),f15.12,')',/&
          ,13x,'output polarization ','(',2(f15.12,','),f15.12,')')
 30 format(13x,'input  polarization ',2('(',2(f15.12,','),f15.12,')'),/&
          ,13x,'output polarization ',2('(',2(f15.12,','),f15.12,')'))
 33 format(4x,"Error: unknow scattering channel '",a,"' for ispec=",i0)
 40 format(/,4x,'FILES')

write(unit=output_unit,fmt=40)
write(unit=output_unit,fmt=45)trim(inptfile),trim(sdtfile) ! ,(achar(9),i=1,1)
write(unit=output_unit,fmt=50)trim(rixsfile),trim(bndfile)
if((irixs == 0 .or. abs(irixs) == 10 ) .and. iabsorp == 0)then
  write(unit=output_unit,fmt=55)trim(datafile)
  return
else
  write(unit=output_unit,fmt=56)trim(datafile),trim(mmefile)
endif

 45 format(4x,"inptfile= '",a,"'",4x,"sdtfile= '",a,"'") ! 1(a)
 50 format(4x,"rixsfile= '",a,"'",4x,"bndfile= '",a,"'")
 55 format(4x,"datafile= '",a,"'")
 56 format(4x,"datafile= '",a,"'",4x,"mmefile= '",a,"'")

write(unit=output_unit,fmt=60)
do isp=1,nlmtdata
  write(unit=output_unit,fmt=65)isp,lmtdata(isp)%sname(:)&
       ,lmtdata(isp)%txtel(:),(ef-lmtdata(isp)%esp(lmtdata(isp)%ian))*dry2ev
  if(iprint > 0)then
    write(unit=output_unit,fmt=70)lmtdata(isp)%isort,lmtdata(isp)%ia1&
         ,lmtdata(isp)%ian,lmtdata(isp)%iat1,lmtdata(isp)%iatn&
         ,lmtdata(isp)%label
    write(unit=output_unit,fmt=75)lmtdata(isp)%ec&
         ,lbound(lmtdata(isp)%esp),ubound(lmtdata(isp)%esp)&
         ,lmtdata(isp)%esp(:)
    do iat=lmtdata(isp)%iat1,lmtdata(isp)%iatn
      write(unit=output_unit,fmt=80)iat,lmtdata(isp)%rat(:,iat)
    enddo
  endif
enddo

 60 format(/,4x,'LMTO DATA USED IN RIXS')
 65 format(4x,'ISPEC=',i2,": edge='",a,"'",4x,"txtel='",a,"'",4x&
          ,'w_edge= ',f0.6,' eV')
 70 format(14x,'isort=',i0,4x,'ia1=',i0,4x,'ian=',i0,4x,'iat1=',i0&
          ,4x,'iatn=',i0,4x,'label=',i0)
 75 format(14x,'ec= ',f0.8,' Ry,',4x,'esp(',i0,':',i0,')=',20(2x,f0.8,' Ry'))

 80 format(14x,'rat(1:3,',i0,')=',3(4x,f12.8))
 90 format(14x,'')

END SUBROUTINE PRINTDATA

SUBROUTINE PRINT_PHASES(isp)

use, intrinsic :: iso_fortran_env
use m_params
use m_rixs

implicit none
! input vars
integer isp
! local vars
integer iat,isplmt

write(unit=output_unit,fmt=18)'Phases (in 2*pi)'
do isplmt=1,nlmtdata
  if(rixs(isp)%txtel(:2) == lmtdata(isplmt)%txtel(:2) .and. &
      rixs(isp)%sname(:) == lmtdata(isplmt)%sname(:))then
    if(lmtdata(isplmt)%iatn-lmtdata(isplmt)%iat1 > 0)then
      do iat=lmtdata(isplmt)%iat1,lmtdata(isplmt)%iatn,2
        write(unit=output_unit,fmt=19)iat&
             ,dot_product(lmtdata(isplmt)%rat(:,iat),rixs(isp)%q)&
             ,iat+1,dot_product(lmtdata(isplmt)%rat(:,iat+1),rixs(isp)%q)
      enddo
    else
!     one atom of this sort
      write(unit=output_unit,fmt=20)lmtdata(isplmt)%iat1&
           ,dot_product(lmtdata(isplmt)%rat(:,lmtdata(isplmt)%iat1),rixs(isp)%q)
    endif
  endif
enddo

 18 format(13x,a)
 19 format(9x,2(4x,'iat=',i4,' phase= ',f12.8))
 20 format(9x,4x,'iat=',i4,' phase= ',f12.8)

END SUBROUTINE PRINT_PHASES
