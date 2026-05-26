INTEGER FUNCTION MEMORY_USAGE()
use, intrinsic :: iso_fortran_env
implicit none
! local vars
character(256) line,units,file
integer id,iok
logical exist

memory_usage=0_INT64
file(:)="/proc/self/status"

inquire(file=file,exist=exist)
if(.not. exist) return

open(newunit=id,file=file,action='read',status='unknown',iostat=iok)
if(iok /= 0)return
do
  read(unit=id,fmt='(a)',iostat=iok)line
  if(iok /= 0)exit
  if(line(1:7) == 'VmPeak:')then
     read(line(8:),*)memory_usage,units(:)
     if(units(:2) == 'GB')memory_usage=memory_usage*1024
     if(units(:2) == 'kB')memory_usage=memory_usage/1024
     if(units(:1) ==  'B')memory_usage=memory_usage/(1024**2)
     exit
  endif
enddo
close(unit=id,iostat=iok)

END FUNCTION MEMORY_USAGE
