subroutine run
use common_V
USE OMP_LIB
use rng, only: rng_seed, rng_uniform
implicit none
integer :: i, i2, iph, nph=100000
real :: rv(3), kv(3), dva, dva_kmps, dva_off_kmps
real, external :: vmag
integer,allocatable :: nscat(:)
character(len=5) :: sdva
integer :: icount
real :: dummy
character(len=500) :: outputfile, outputPOfile
real :: kvsq

outputfile=trim(outputbase)//'_1e6_Photons_irun'//trim(crun)//'.bin'; nph=1e6
!outputfile=trim(outputbase)//'_Nph1e5_irun'//trim(crun)//'.bin'; nph=1e5
allocate(rngstate(nph))
allocate(rvo(3,nph))
allocate(kvo(3,nph))
allocate(dvao(nph))
allocate(nscat(nph))

allocate(output_data(nodata,nph))
!odata(1)    ! Number of scattering
!odata(2:4)  ! Last scattering location 
!odata(5:7)  ! Final scattering direction
!odata(8)    ! Final scattering wavelength
!odata(9:11) ! First scattering location 
!odata(12:14)! Initial location
!odata(15:17)! Initial direction
!odata(18)   ! Initial wavelength


write(sdva,'(I4.4)')nint(abs(dva_kmps))
if(dva_kmps.lt.0)sdva='m'//trim(adjustl(sdva))
if(dva_kmps.ge.0)sdva='p'//trim(adjustl(sdva))
dva_off_kmps=Rvir_cMpcoh*dble(Hhub)*dble(Mpc_in_cm)/h/(1+zs)/1e5


print*,'# Initial wavelength in km/s:',dva_kmps+dva_off_kmps
print*, '# Initial wavelength of photon:',dva_off_kmps,'(km/s)'

print*,'# output file:'//trim(outputfile)
open(1,file=trim(outputfile),form='unformatted',access='stream',action='write')
write(1)nodata,nph
icount=0
!$OMP PARALLEL PRIVATE(i,dummy,iph,rv,kv,dva)
!$OMP DO
do iph=1, nph ! Loop for photons
 call rng_seed(rngstate(iph), iph+irun*nph)
 do i=1,100; dummy=rng_uniform(rngstate(iph)); enddo

 kvsq = 4
 do while (kvsq.ge.1)
  do i=1,3
   kv(i)=2*rng_uniform(rngstate(iph))-1.
  enddo
  kvsq = sum(kv**2)
 enddo
 if(vmag(kv).eq.0.) then
 print*,'# vmag(kv)=0! at ',iph,'th.'
 endif
 kv=kv/vmag(kv); rv=kv*Rvir_cMpcoh/Lbox
 dva_kmps=rng_uniform(rngstate(iph))*2.5e3 - 2.2e3
 dva=(dva_kmps+dva_off_kmps)*1e5/c_sl
 !print*, '# Initial location of ',iph,'th photon:', rv*Lbox,'(cMpc/h)'

 output_data(12:14,iph)=rv                ! Initial location
 output_data(15:17,iph)=kv                ! Initial direction
 output_data(18,iph)   =dva_kmps*1e5/c_sl ! Initial wavelength
 call LyART(rv,kv,dva,iph,output_data(:,iph))
 !read(*,*)
 write(1)output_data(:,iph)
 icount=icount+1
 if(mod(icount,100).eq.0) print*,'#',icount,'photons done.'
enddo! iph
!$OMP END DO
!$OMP END PARALLEL

close(1)


if(PO_plane_flag) then
outputPOfile=trim(outputbase)//trim(sdva)//'_PO_plane.bin'
open(iu_po,file=trim(outputPOfile),form='unformatted',access='stream',action='write')
write(iu_po)6,ndim,ndim;write(iu_po)SB_PO; close(iu_po)
print*,'# Peeling-off plane data saved: '//trim(outputPOfile)
endif


end subroutine run
