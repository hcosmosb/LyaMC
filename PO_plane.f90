subroutine PO_plane(iph,rv,vp,dva)
use common_v, only : iu_PO, Lbox, SB_PO, ndim
implicit none
double precision, intent(IN) :: rv(3)
real, intent(IN) :: dva, vp(3)
integer, intent(IN) :: iph
real :: dvaf, kv(3)
double precision :: rvf(3), tau, tau_lim=10.
double precision, external :: tau_btw_2p
double precision :: blim=0.495d0
integer :: ind(3)

ind=int((rv+0.5)*ndim)

kv=(/1.,0.,0./); dvaf=dva-sum(kv*vp); rvf=(/blim,rv(2),rv(3)/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(1,ind(2),ind(3))=SB_PO(1,ind(2),ind(3))+exp(-tau)

kv=(/-1.,0.,0./); dvaf=dva-sum(kv*vp); rvf=(/-blim,rv(2),rv(3)/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(2,ndim-1-ind(2),ind(3))=SB_PO(2,ndim-1-ind(2),ind(3))+exp(-tau)

kv=(/0.,1.,0./); dvaf=dva-sum(kv*vp); rvf=(/rv(1),blim,rv(3)/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(3,ind(1),ind(3))=SB_PO(3,ind(1),ind(3))+exp(-tau)

kv=(/0.,-1.,0./); dvaf=dva-sum(kv*vp); rvf=(/rv(1),-blim,rv(3)/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(4,ndim-1-ind(1),ind(3))=SB_PO(4,ndim-1-ind(1),ind(3))+exp(-tau)

kv=(/0.,0.,1./); dvaf=dva-sum(kv*vp); rvf=(/rv(1),rv(2),blim/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(5,ind(1),ind(2))=SB_PO(5,ind(1),ind(2))+exp(-tau)

kv=(/0.,0.,-1./); dvaf=dva-sum(kv*vp); rvf=(/rv(1),rv(2),-blim/)
tau=tau_btw_2p(rv,rvf,dvaf,iph)
if(tau.lt.tau_lim)SB_PO(6,ndim-1-ind(1),ind(2))=SB_PO(6,ndim-1-ind(1),ind(2))+exp(-tau)

end subroutine PO_plane
