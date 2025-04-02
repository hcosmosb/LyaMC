!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
module common_variables
  real :: T, x
end module common_variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine calcTr(ncell,nH,Temp,vp,nnu,lambda_arr,z_red,tau_arr,nRcell,cell_size)
implicit none
integer, intent(IN) :: ncell, nnu, nRcell
real, dimension(ncell), intent(IN) :: nH, Temp, vp
real, dimension(nnu), intent(IN) :: lambda_arr
real, dimension(nnu), intent(OUT) :: tau_arr
real, intent(IN) :: z_red, cell_size
integer :: i
real :: &
k_boltz = 1.380648e-16, &
hhubpar = 0.6776, &
Mpcincm = 3.086e24, &
omega_m = 0.307115e+00, &
omega_b = 0.45e-01, &
m_p = 1.6726219e-24, &
lambda_lya=1215.67, &
c_sl = 2.99792e10
real :: v_hub, vp_rel, v_th, a_rlw
real :: Hhub
double precision :: dl_a_Angst, dlambda
real :: dtau,  cell_Mpcoh
double precision, external :: sig_a
integer :: inu

Hhub = hhubpar * 100d0 * 1d5/Mpcincm &
     * dsqrt(omega_m*(1d0 + z_red)**3 + (1d0 - omega_m))
cell_Mpcoh=cell_size/(Mpcincm/hhubpar)

tau_arr=0d0
write(99,'(A)')' d [Mpc/h],  dLamda [A],  tau,     nH [cm^-3],    Temp [K],    vpe-vhub [km/s]'
write(88,'(A)')' d [Mpc/h],  dLamda [A],  tau,     nH [cm^-3],    Temp [K],    vpe-vhub [km/s]'
write(77,'(A)')' d [Mpc/h],  dLamda [A],  tau,     nH [cm^-3],    Temp [K],    vpe-vhub [km/s]'
do i=nRcell+1,ncell
!do i=1,256
vp_rel  = vp(i)!-vp(1)

v_hub = Hhub*cell_size*(i-1)/c_sl
v_th  = (2d0*k_boltz*Temp(i)/m_p)**0.5
!print*, i, vp_rel*c_sl/1e5, v_hub*c_sl/1e5, (vp_rel+v_hub)*c_sl/1e5
a_rlw = 4.7d-4*(Temp(i)/1.d4)**(-0.5)
   do inu=1,nnu
      dl_a_Angst=(lambda_arr(inu)-lambda_lya)+(vp_rel+v_hub)*lambda_lya
      dtau = sig_a(dl_a_Angst,dble(Temp(i)))*nH(i)*cell_size !!!
      tau_arr(inu)  = tau_arr(inu) + dtau
      dlambda = lambda_arr(inu)-lambda_lya
      !if(abs(dlambda-0.1).lt.1e-5)write(1,*)i,tau_arr(inu),nH(i),Temp(i)
      if(abs(dlambda).lt.1e-3)print'(2f9.3,f15.5,ES13.5,2f13.2,I4)',&
       (i-0.5)*cell_Mpcoh,dl_a_Angst,tau_arr(inu),nH(i),Temp(i),vp_rel*c_sl/1e5+v_hub*c_sl/1e5,i
      if(abs(dlambda).lt.1e-3)write(99,'(2f9.3,f15.5,ES13.5,2f12.2)')&
      (i-0.5)*cell_Mpcoh,dl_a_Angst,tau_arr(inu),nH(i),Temp(i),vp_rel*c_sl/1e5+v_hub*c_sl/1e5
      if(abs(dlambda-0.8166).lt.1e-3)write(88,'(2f9.3,f15.5,ES13.5,2f12.2)')&
      (i-0.5)*cell_Mpcoh,dl_a_Angst,tau_arr(inu),nH(i),Temp(i),vp_rel*c_sl/1e5+v_hub*c_sl/1e5
      if(abs(dlambda+0.8166).lt.1e-3)write(77,'(2f9.3,f15.5,ES13.5,2f12.2)')&
      (i-0.5)*cell_Mpcoh,dl_a_Angst,tau_arr(inu),nH(i),Temp(i),vp_rel*c_sl/1e5+v_hub*c_sl/1e5

   enddo
enddo
!stop

end subroutine calcTr
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function sig_a_dw(dl_a_Angst,T_in)
  use common_variables
  implicit none
  double precision, intent(IN) :: dl_a_Angst, T_in
  double precision :: Delnu_D, dnu_a, l_a_Angst=1216d0
  double precision ::  phi,  a_V
  real :: c_sl = 2.99792e10, k_b = 1.380648e-16, m_p = 1.6726219e-24
  T=T_in
  dnu_a = c_sl/((dl_a_Angst+l_a_Angst)*1d-8) - c_sl/(l_a_Angst*1d-8)
  Delnu_D=2.46d15*(2d0*k_b*T/(m_p*c_sl**2))**0.5
  x=abs(dnu_a/Delnu_D)
  !print*,'# !!!',(2d0*k_b*T/(m_p*c_sl**2))**0.5
  a_V = 4.7d-4
  phi=4.7d-4/dsqrt(3.141592d0)/x**2
  sig_a_dw=phi*5.889d-14

end function sig_a_dw
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function sig_a(dl_a_Angst,T_in)
  use common_variables
  implicit none
  double precision, intent(IN) :: dl_a_Angst, T_in
  double precision, external :: dphi
  double precision :: Delnu_D, dnu_a, l_a_Angst=1216d0
  double precision :: phi1, phi2, phi3, phi
  real :: & 
  c_sl = 2.99792e10, &
  k_boltz = 1.380648e-16, &
  m_p = 1.6726219e-24

  T=T_in
  dnu_a = c_sl/((dl_a_Angst+l_a_Angst)*1d-8) - c_sl/(l_a_Angst*1d-8)
  Delnu_D=2.46d15*(2d0*k_boltz*T/(m_p*c_sl**2))**0.5
  x=abs(dnu_a/Delnu_D)

  if(x.le.10d0) then
    call dqsimp(dphi, -4d0,     x-0.01d0, phi1)
    call dqsimp(dphi, x-0.01d0, x+0.01d0, phi2)
    call dqsimp(dphi, x+0.01d0, 11d0,     phi3)
    phi=phi1+phi2+phi3
    !print*,phi1,phi2,phi3
  else if(x.gt.10d0) then
    call dqsimp(dphi, -6d0,6d0,     phi)
  endif
  sig_a=phi*5.889d-14*(T/1d4)**(-0.5)

end function sig_a
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function sig_a_app(dl_a_Angst,T_in)
 implicit none
 double precision, intent(IN) :: dl_a_Angst, T_in
 double precision :: x,x2,z,q, a_V, phi_app, dnu_a, Delnu_D, l_a_Angst=1216d0
 real :: c_sl = 2.99792e10, k_boltz = 1.380648e-16, m_p = 1.6726219e-24,&
  pi=3.141592

 a_V = 4.7d-4*(T_in/1d4)**(-0.5d0)
 dnu_a = c_sl/((dl_a_Angst+l_a_Angst)*1d-8) - c_sl/(l_a_Angst*1d-8)
 Delnu_D=2.46d15*(2d0*k_boltz*T_in/(m_p*c_sl**2))**0.5
 x=abs(dnu_a/Delnu_D);x2=x**2
 z=(x2-0.855)/(x2+3.42)
 if(z.le.0)  then
  q=0d0
 else
  q=z*(1.+21./x2)*a_V/pi/(x2+1)*(0.1117+z*(4.421+z*(-9.207+5.674*z)))
 endif
 phi_app=(q+exp(-x2)/1.77245385)*sqrt(pi)
 !phi_app=(q+exp(-x2))
 sig_a_app = phi_app*5.889d-14*(T_in/1d4)**(-0.5)
end function sig_a_app
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
double precision function dphi(y)
  use common_variables
  implicit none
  double precision, intent(IN) :: y
  double precision :: a_V

  a_V = 4.7d-4*(T/1d4)**(-0.5d0)
  dphi = dexp(-y**2)/((y-x)**2+a_V**2)
  dphi = dphi*a_V/3.141592d0


end function dphi
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE dqsimp(func,a,b,s)
      INTEGER JMAX
      REAL*8 a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=2d-3, JMAX=63)
!U    USES dtrapzd
      INTEGER j
      REAL*8 os,ost,st
      ost=-1.d30
      os= -1.d30
      do 11 j=1,JMAX
        call dtrapzd(func,a,b,st,j)
        s=(4.d0*st-ost)/3.d0
!       print*,'#!1',j,a,b,real(st)
        if (abs(s-os).le.EPS*abs(os).and.j.gt.6.or.j.gt.20) then
        return
        endif
        os=s
        ost=st
11    continue
!      pause 'too many steps in qsimp'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE dtrapzd(func,a,b,s,n)
      INTEGER n
      REAL*8 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*8 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5d0*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE qsimp(func,a,b,s)
      INTEGER JMAX
      REAL*4 a,b,func,s,EPS
      EXTERNAL func
      PARAMETER (EPS=2e-3, JMAX=63)
!U    USES trapzd
      INTEGER j
      REAL*4 os,ost,st
      ost=-1.e30
      os= -1.e30
      do 11 j=1,JMAX
        call trapzd(func,a,b,st,j)
        s=(4.*st-ost)/3.
        !print*,'#!1',j,a,b,real(st)
        if (abs(s-os).le.EPS*abs(os).and.j.gt.9.or.j.gt.20) then
        return
        endif
        os=s
        ost=st
11    continue
!      pause 'too many steps in qsimp'
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE trapzd(func,a,b,s,n)
      INTEGER n
      REAL*4 a,b,s,func
      EXTERNAL func
      INTEGER it,j
      REAL*4 del,sum,tnm,x
      if (n.eq.1) then
        s=0.5*(b-a)*(func(a)+func(b))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5*del
        sum=0.
        do 11 j=1,it
          sum=sum+func(x)
          x=x+del
11      continue
        s=0.5*(s+(b-a)*sum/tnm)
      endif
      return
      END
!  (C) Copr. 1986-92 Numerical Recipes Software :)z%+.

