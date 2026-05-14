SUBROUTINE SET_NCF_NCL(ncf,ncl)
use, intrinsic :: iso_fortran_env
use m_bzmesh
use m_functions
use m_params, only: iprint
use m_bnd, only: e, nb, ef
implicit none
! output
integer ncf,ncl
! local vars
integer ib,k,istop
real emin(nb),emax(nb)

! nspin==1 and ndimspin==1 is assumed
! ncf is the lowest NOT fully occupied band
! ncl the last NOT completely empty band

istop=0
emin(:)= huge(1.)
emax(:)=-huge(1.)
do k=1,nkibz
  do ib=1,nb
    if(real(e(ib,k)) < emin(ib))emin(ib)=real(e(ib,k))
    if(real(e(ib,k)) > emax(ib))emax(ib)=real(e(ib,k))
  enddo          ! ib
enddo            ! k
ncf=0
ncl=nb+1
do ib=1,nb
  if(emax(ib) < real(ef)+epsilon(1.))ncf=ib
  if(emin(ib) >= real(ef))then
    ncl=ib
    exit
  endif
enddo            ! ib
ncf=ncf+1        ! the first not fully occupied
ncl=ncl-1        ! the last not completely empty
if(ncf > nb .or. ncl > nb)then
  write(unit=output_unit,fmt="(4x,a)") 'In READBND: too few bands in bndfile.'
  call deallocate_global_arrays
  stop
endif
if(iprint > 10) write(unit=output_unit,fmt=5)ncl,ncf,nb

 5 format(4x,'In SET_NBI_NBF: ncl=',i0,' ncf=',i0,' nbands=',i0)
 
END SUBROUTINE SET_NCF_NCL

