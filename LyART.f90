subroutine LyART(rv_in,kv_in,dva_in,iph,odata)
use rng, only: rng_uniform
use common_v
implicit none
real, intent(in) :: rv_in(3), kv_in(3), dva_in
integer, intent(IN) :: iph
!integer, intent(out) :: nscat
!real, intent(out) :: rv_out(3), kv_out(3), dva_out
real, intent(inout) :: odata(nodata)

double precision :: rv(3), rv_pre(3)
real :: kv(3), kphi, kmu, kv_pr(3)
integer :: i, iloop
real :: rind(3), dummy
double precision :: tau_f, tau_t, tau, Deltau
double precision, external :: tau_btw_2p
real, external :: vmag
double precision, external :: vmagd
double precision :: dri=0.1, dr
real :: temp_in, dva, dva_p
double precision :: ar, vth, xf, vpar, vperp
real :: harv
real :: vpf(3), vphi
real :: vpe_los, vpe_in(3), ct_costh,ct_sinth,ct_cosphi,ct_sinphi
integer :: iscatter
logical :: print_flag=.false.,exit_flag, db_flag=.false.
logical :: print_flag2=.false., print_flag3=.false.
!integer, save :: iphoton=0
integer :: iphoton=0
real :: rbfac=1.2, ran
integer :: nseed
integer, allocatable :: seed_arr(:)
real :: g
real :: dlambda_A, dmfp, nHI_in
double precision, external :: sig_a_app
logical :: HD_flag=.false.

!print_flag3=.true.
!print_flag2=.true.
!print_flag=.true.
if (print_flag)  write(42,*) iph

!iphoton=iphoton+1
if(print_flag2) then
print*, '##############################'
print*, '######## ',iph,'th photon.'
print*, '##############################'
endif

odata(1:11)=0.
rv=rv_in; kv=kv_in/vmag(kv_in); dva=dva_in
rv_pre=rv
iscatter=0

odata(2:4) = rv ! Last scattering location 
odata(5:7) = kv ! Final scattering direction
odata(8) = dva  ! Final scattering wavelength

do !iscatter


exit_flag=.false.

if(print_flag2) then
print*, '#### ',iscatter+1,'th path'
endif
!!!! Propagate photon until scattered
13 continue
tau_f=rng_uniform(rngstate(iph))
if (tau_f.eq.0d0) go to 13 ! To prevent tau_f=NaN
tau_f=-log(tau_f)
if (print_flag)  write(42,*) real(rv*Lbox),dva*c_sl/1e5, tau_f

if(print_flag2)print'(X,A,f6.3,A)','# tau_f: ',tau_f
tau=0.; dr=dri; iloop=0

!!! in case of very high optical depth !!!
rind=(rv+0.5)*ndim
vpe_in  = vpe(int(rind(1)),int(rind(2)),int(rind(3)),:)
temp_in = temp(int(rind(1)),int(rind(2)),int(rind(3)))
nHI_in = nHI(int(rind(1)),int(rind(2)),int(rind(3)))
vpe_los = IP(vpe_in,kv)
dlambda_A = (dva + vpe_los)*lambda_lya
dmfp = (sig_a_app(dble(dlambda_A),dble(Temp_in))*nHI_in)**(-1)
HD_flag=.false.
if(dmfp/cell_cgs.lt.1e-1) then
 HD_flag=.true.
 dr = tau_f*dmfp/cell_cgs/dble(ndim)
 rv = rv+dr*kv
 dva=dva+dr*Hhub*cell_cgs*ndim/c_sl
 iloop = -1
 go to 11
endif

do while(abs(tau-tau_f).gt.min(0.01,tau_f*0.01))

 call exit_condition(real(rv+dr*kv),exit_flag)
 do while(exit_flag.and.dr.gt.1./ndim)
 dr=dr/2.
 call exit_condition(real(rv+dr*kv),exit_flag)
 enddo
 if(exit_flag.and.dr.lt.1./ndim) goto 10
 Deltau=tau_btw_2p(rv,rv+dr*kv,dva,iph)
 tau_t=tau+Deltau


 if(tau_t.lt.tau_f) then
  rv=rv+dr*kv
  dva=dva+dr*Hhub*cell_cgs*ndim/c_sl
  tau=tau_t
  dr=dr*rbfac
 else
  dr=dr/2.
 endif
 iloop=iloop+1
enddo ! do while

! Final second order correction
dr=dr/rbfac
dr=dr*(tau_f-tau)/Deltau
if(vmagd(dr*kv).gt.1e-1/real(ndim)) then
 call exit_condition(real(rv+dr*kv),exit_flag); if(exit_flag) goto 10
 Deltau=tau_btw_2p(rv,rv+dr*kv,dva,iph)
 if(print_flag2)print*,'# r & tau correction:',dr,Deltau
 tau_t=tau+Deltau
 rv=rv+dr*kv
 dva=dva+dr*Hhub*cell_cgs*ndim/c_sl
 tau=tau_t
endif
!print_flag=.true.

11 continue
if(print_flag) then
 print*,'# Scattered after ',iloop,' iterations.'
 print*,'# tau error ',real(tau_f-tau),' out of ',real(tau_f),'.'
 print*,'# Wavelength at scatter:', dva*c_sl/1e5,'(km/s)'
 print*,'# Location at scatter: ',real(rv*Lbox),'(cMpc/h)'
 !write(42,*) real(rv*Lbox),dva*c_sl/1e5
endif
if(print_flag2)print'(X,A,f6.2,A)','# ds at scatter: ',real(vmagd(rv)*Lbox),' (cMpc/h)'
if(print_flag2)print'(X,A,f10.4,A)','# Moved ',real(vmagd(rv-rv_pre)*Lbox*1e3),' (ckpc/h)'
rv_pre=rv
if(print_flag2)print'(X,A,ES11.3,A)','# nHI = ',real(nHI(int(rind(1)),int(rind(2)),int(rind(3)))),' /cm^3'
if(print_flag2)print*,'# HD flag: ',HD_flag

rind=(rv+0.5)*ndim
if(int(rind(1)).lt.0.or.int(rind(1)).ge.ndim) then
print*,'# Array index limit exceeded. Discard photon.'
print*,'!1',iph,iscatter,rv,rind
iscatter=-1; odata(2:4) = 0.; odata(5:7) = 0.; odata(8) = 0.; go to 10
endif
if(int(rind(2)).lt.0.or.int(rind(2)).ge.ndim) then
print*,'# Array index limit exceeded. Discard photon.'
print*,'!2',iph,iscatter,rv,rind
iscatter=-1; odata(2:4) = 0.; odata(5:7) = 0.; odata(8) = 0.; go to 10
endif
if(int(rind(3)).lt.0.or.int(rind(3)).ge.ndim) then
print*,'# Array index limit exceeded. Discard photon.'
print*,'!3',iph,iscatter,rv,rind
iscatter=-1; odata(2:4) = 0.; odata(5:7) = 0.; odata(8) = 0.; go to 10
endif
!!!! Assign direction and frequency to the scattered photon (Hyojeong's contribution)
temp_in=temp(int(rind(1)),int(rind(2)),int(rind(3)))
temp_in=max(temp_in,10.)
ar=4.702d-4*(temp_in/1d4)**(-0.5) ! (Natural line width)/(Doppler line width) 
vth=(2d0*k_b*temp_in/m_p)**0.5/c_sl ! Thermal velocity
if(print_flag2)print'(X,A,f6.0,A,f6.1,A)','# T = ',temp_in,' K, vth = ',vth*c_sl/1e5,' km/s'


!!!! Calculate scatterer velocity vpf(3)
vpe_in=vpe(int(rind(1)),int(rind(2)),int(rind(3)),:) ! Peculiar velocity at the scattering location
vpe_los=sum(kv*vpe_in) ! LOS component of vpe_in
xf=-(dva+vpe_los)/vth ! D-less freq in the gas frame. Note that positive LOS v add to dva.
call scat_vpar(ar, xf, vpar, iph); vpar=vpar*vth ! Photon direction component
call gasdev_s(harv, iph);          vperp=harv/SQRT(2.)*vth ! Perpendicular component
vphi=2*pi*rng_uniform(rngstate(iph)) ! Azimuthal direction of the perp component
vpf=(/vperp*cos(vphi),vperp*sin(vphi),vpar/) ! Atom velocity in gas frame w.r.t. kv

!!! Use two rotations get vpf into the global frame
ct_costh=kv(3);ct_sinth=(kv(1)**2+kv(2)**2)**0.5 
if(ct_sinth.eq.0.) then
ct_cosphi=1.; ct_sinphi=0.
else
ct_cosphi=kv(1)/ct_sinth; ct_sinphi=kv(2)/ct_sinth
endif
vpf=(/vpf(1)*ct_costh+vpf(3)*ct_sinth,vpf(2),-vpf(1)*ct_sinth+vpf(3)*ct_costh/) ! Rotate w.r.t. the y axis
vpf=(/ct_cosphi*vpf(1)-ct_sinphi*vpf(2),ct_sinphi*vpf(1)+ct_cosphi*vpf(2),vpf(3)/) ! Rotate w.r.t. the z axis
vpf=vpf+vpe_in ! Atom velocity in the global frame


!!!! New scattered photon direction from isotropic random distribution
kv_pr = kv
kphi=2.*pi*rng_uniform(rngstate(iph))
kmu=2.*rng_uniform(rngstate(iph))-1.
kv=(/(1.-kmu**2)**0.5*cos(kphi),(1.-kmu**2)**0.5*sin(kphi),kmu/) ! New photon direction

if(print_flag2)print'(X,A,f11.4)','# xf = ',real(xf)
if(print_flag2)print*,'# Wavelength before scatter:',(dva+IP(kv_pr,vpe_in))*c_sl/1e5,'(km/s)'
!!!!!!!!!!!!!!!!!! New photon frequency !!!!!!!!!!!!!!!!!!!!
g = 2.6e-4*(temp_in/1e4)**(-0.5)
dva=dva - IP((kv-kv_pr),vpf) -g*vth*(IP(kv,kv_pr)-1) ! Scattered photon frequency
!!!dva=-(xf*vth+sum(kv*vpf)) ! Scattered photon frequency old version
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


dva_p=dva-sum(kv*rv)*Hhub*cell_cgs*ndim/c_sl ! in the source plane
if(PO_plane_flag) call PO_plane(iph,rv,vpf,real(-xf*vth))


if(db_flag)stop

if(print_flag) then
 print*,iph,'# Coordinate:',int(rind(1)),int(rind(2)),int(rind(3))
 print*,iph,'# HI density:',real(nHI(int(rind(1)),int(rind(2)),int(rind(3))))
 print*,iph,'# New direction:',kv
 print*,iph,'# Wavelength after scatter:',dva*c_sl/1e5,'(km/s)'
endif
if(print_flag2)print*,'# New wavelength after scatter:',(dva+IP(kv,vpe_in))*c_sl/1e5,'(km/s)'
if(print_flag2)print'(X,A,f7.3)','# Direction:',IP(kv,real(rv))/vmagd(rv)


iscatter=iscatter+1
if(iscatter.eq.1)odata(9:11) = rv ! First scattering location 
if(iscatter.ge.100000000)then
print*,'# Photon No.',iph,' scattered 100,000,000 times. Discard photon.'
odata(2:4) = 0.; odata(5:7) = 0.; odata(8) = 0.; go to 10
endif


odata(2:4) = rv ! Last scattering location 
odata(5:7) = kv ! Final scattering direction
odata(8) = dva  ! Final scattering wavelength

enddo !iscatter

10 continue
if(print_flag2.or.print_flag3) then
print*,'# Photon No.',iph,' escaped after ',iscatter,'scatters.'
!print*,'# Final location: ', real(rv*Lbox),' (cMpc/h)'
kmu=sum(kv*rv)/vmagd(rv)
!print*,'# Final direction cosine: ', kmu
!print*,'# Final wavelength: ',dva*c_sl/1e5,' (km/s)'
dva_p=dva-sum(kv*rv)*Hhub*cell_cgs*ndim/c_sl
print*,'# Final wavelength at the source plane: ',dva_p*c_sl/1e5,' (km/s)'
print*,'# Apprent distance from the source: ',real((1.-kmu**2.)**0.5*vmagd(rv)*Lbox),' (cMpc/h)'
endif

!write(42,*) real(rv*Lbox),dva*c_sl/1e5, 0

if(iscatter.gt.10000000) print*,'# Photon ',iph,' scattered ',iscatter,' times.'
!if(db_flag)stop
!nscat=iscatter
!rv_out=rv; kv_out=kv; dva_out=dva
odata(1) = real(iscatter) ! Number of scattering

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!
real function IP(vec1,vec2)
real,intent(IN) :: vec1(3),vec2(3)
IP=sum(vec1*vec2)
end function IP
end subroutine LyART
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function tau_btw_2p(rvi,rvf,dva,iph)
use common_V, only : rvi_g, kv_g, dva_g
implicit none
double precision, intent(IN) :: rvi(3), rvf(3)
real, intent(IN) :: dva
integer, intent(IN) :: iph
double precision, external :: dtau
double precision :: dr, vmagd
real, external :: vmag

dva_g=dva
rvi_g=rvi
kv_g=rvf-rvi; kv_g=kv_g/vmag(kv_g)
dr=vmagd(rvf-rvi)
if(dr.eq.0d0.or.dva.ne.dva) then
print*,'# dr=0 or dva=NaN found!'
print*,'# iph,rvi,rvf,dva=',iph,rvi,rvf,dva
stop
endif
call dqsimp(dtau,0d0,dr,tau_btw_2p)

end function tau_btw_2p
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function dtau(dr)
use common_V
implicit none
double precision, intent(IN) :: dr
double precision :: rv(3)
real :: rind(3)
real :: nHI_in, vpe_in(3), vplos_in, Temp_in, v_Hub, dlambda_A
integer :: i
double precision, external :: sig_a_app

rv=rvi_g+dr*kv_g

rind= (rv+0.5)*ndim

if(int(rind(1)).lt.0.or.int(rind(1)).ge.ndim) then
print*,'# Array index limit exceeded. (2)'
print*,rind,dr,kv_g
stop
endif
if(int(rind(2)).lt.0.or.int(rind(2)).ge.ndim) then
print*,'# Array index limit exceeded. (2)'
print*,rind,dr,kv_g
stop
endif
if(int(rind(3)).lt.0.or.int(rind(3)).ge.ndim) then
print*,'# Array index limit exceeded. (2)'
print*,rind,dr,kv_g
stop
endif

nHI_in  = nHI(int(rind(1)),int(rind(2)),int(rind(3)))
do i=1,3
 vpe_in(i)=vpe(int(rind(1)),int(rind(2)),int(rind(3)),i)
enddo
vplos_in= sum(vpe_in*kv_g)
Temp_in = Temp(int(rind(1)),int(rind(2)),int(rind(3)))
temp_in=max(temp_in,10.)

v_Hub   = Hhub*cell_cgs*dr*ndim/c_sl
dlambda_A = (dva_g+vplos_in+v_Hub)*lambda_lya
dtau = sig_a_app(dble(dlambda_A),dble(Temp_in))*nHI_in*cell_cgs*ndim 


end function dtau
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
real function vmag(dv)
real, intent(IN) :: dv(3)
vmag=sum(dv**2)**0.5
end function vmag
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function vmagd(dv)
double precision, intent(IN) :: dv(3)
vmagd=sum(dv**2)**0.5
end function vmagd
