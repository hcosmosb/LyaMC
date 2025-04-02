FUNCTION voi_func(x, x_freqD, a_rlw)
  IMPLICIT NONE
  double precision, INTENT(IN):: x, x_freqD, a_rlw
  double precision, PARAMETER :: c_sl = 2.9979245800d10, & !Speed of light in cm 
       pi=3.141592653d0
  double precision :: voi_func, ans, nc
  INTEGER :: n, it, j
  
  call voigt(a_rlw, x_freqD, ans)
  nc =  a_rlw/(pi*ans)
  voi_func = nc * exp(-x * x)/((x_freqD-x)**2+a_rlw**2)
               
END FUNCTION voi_func

SUBROUTINE scat_vpar(a, v, vparallel, iph)
  use rng, only: rng_uniform
  use common_v, only: rngstate
  IMPLICIT NONE

  double precision, INTENT(IN):: a, v
  double precision, INTENT(OUT) :: vparallel
  integer, intent(IN) :: iph
  double precision, PARAMETER :: &
  c_sl = 29979245800.d0, & !Speed of light in cm 
  pi=3.141592653d0, &
  temperat=1d4, &
  TAUES = 1d5, &
  nu_lya=2.47d15, &
  nu_D=1.06d11, &
  mass_hydro=1.6737236d-24, &      !hydrogen mass in g
  k_boltz = 1.380648d-16, &        !Boltzmann constant in erg/K
  cell_size = 1.96298685d24, &     !cell size in cm
  Mpc = 3.08568025d24, &
  Ol=0.73, Om=1.-Ol, hhubpar = 0.7, & !Hubble constant (h) at today in 100 km/s/Mpc
  Hhub = hhubpar * 100.0 * 1.e5/Mpc, & !Hubble constant in s^-1
  gbbox = 114d0, &
  n_box = 256*24d0
  INTEGER, parameter :: tbl_size=201
        
  INTEGER :: i1, j, k, k1, k2, k3(1), ii, jj, kk
  real :: harvest, harvest2 
  double precision :: a_rlw, x_freqD, &
              vki_pe1_1, vki_pe1_2, vki_pe1_3, vki_pe2_1, vki_pe2_2, vki_pe2_3, &
              v_thermal, vperpendicular, &
              x_temp, y_temp, y_random, x_init, &
              ans, ans_i, ans_f, dansda, dansdv, x_freqD_sign, &
              x_freqD_min, x_freqD_max, x_freqD_int, integ_voigt, &
              x_temp_p, x_1, y_1, det_1, z_1, & 
              p_max, p_max1, p_max2, prob1, prob2, prob3, voifunc, &
              num1, num2, num3, num4, num5, x_interp, p_interp1, p_interp2, &
              P_resonant, pdf_max, &
              tau, delta_Crd, tau_crx, delta_tau, tau_proceed, tau_part, density, &
              voi_func, s
              
     
a_rlw = a
x_freqD = v

        
IF (x_freqD < 0d0) then
x_freqD = abs(x_freqD)
x_freqD_sign = -1d0
else
x_freqD = abs(x_freqD)
x_freqD_sign = 1d0
end if

IF (x_freqD < 1.5d0) THEN

     y_random = 1.d0; y_temp = 0.d0

     Do while (y_random .gt. y_temp)

        det_1=2
            Do while (det_1 > 1)
                x_1=2*rng_uniform(rngstate(iph))-1
                y_1=2*rng_uniform(rngstate(iph))-1
                det_1=x_1*x_1+y_1*y_1
            END DO
        z_1=x_1/y_1
        x_temp = x_freqD+a_rlw*z_1
        y_temp = exp(-1d0*x_temp*x_temp)
        y_random = rng_uniform(rngstate(iph))

     END DO
        
ELSE 

!print*, '2-1'
!       call FIND1D(x_freqD, p_interp1, x_freqd_tbl, p1_tbl, tbl_size)
!       call FIND1D(x_freqD, p_interp2, x_freqd_tbl, p2_tbl, tbl_size)
!       prob1=10**(p_interp1)
!       prob2=10**(p_interp2)
!       prob3=1-prob1-prob2

! print*, prob1, prob2, prob3

! print*, a_rlw, x_freqD, 'here is the code'
 ! call qsimp(voi_func, -8.d0, 1.d0,s, a_rlw, x_freqD)
! print*, voi_func()
! print*,'!',a_rlw, x_freqD
 call voi_int(a_rlw, x_freqD, -10d0, 1d0, s)
 prob1=s
 
 call voi_int(a_rlw, x_freqD, 1d0, x_freqD-0.25d0, s)
 prob2=s
 prob3=1-prob1-prob2

        !call ran1_s(harvest)
        x_temp_p=rng_uniform(rngstate(iph))
   IF (x_temp_p < prob1) Then
        y_random = 1d0
        y_temp = 0d0

                Do while (y_random .gt. y_temp)
                x_temp=2.d0
                    Do while (x_temp > 1.d0)
                    x_temp=dsqrt(-1d0*LOG(rng_uniform(rngstate(iph))))*&
                              cos(2d0*pi*rng_uniform(rngstate(iph)))
                    END DO
                y_random=rng_uniform(rngstate(iph))
                y_temp=(x_freqD-1)**2/(x_freqD-x_temp)**2 
        END DO

   ELSE IF (x_temp_p < (prob1+prob2)) THEN

!      call voigt(a_rlw, x_freqD, ans, dansda, dansdv)
        call voigt(a_rlw, x_freqD, ans)
!print*, '2-2'
        y_random = 1.d0
        y_temp = 0.d0

        p_max1=exp(-1.D0)*a_rlw/(((x_freqD-1.d0)**2+a_rlw**2)*pi*ans)
        p_max2=exp(-1.D0*(x_freqD-0.25)**2)*a_rlw/((0.25**2+a_rlw**2)*pi*ans)

        if (p_max1 > p_max2) then
                p_max=p_max1
        else
                p_max=p_max2
        end if

        Do while (y_random .gt. y_temp)
                !call ran1_s(harvest)
                x_temp=1.d0+rng_uniform(rngstate(iph))*(x_freqD-1.25)
                y_temp=exp(-1.D0*x_temp*x_temp)*a_rlw/( &
                ans*p_max*pi*(a_rlw**2+(x_freqD-x_temp)**2))
                !call ran1_s(harvest)
                y_random=rng_uniform(rngstate(iph))

        END DO

   ELSE

        y_random = 1.d0
        y_temp = 0.d0

        Do while (y_random .gt. y_temp)

         x_temp=0
         Do while (x_temp .le. x_freqD-0.25)
            call gasdev_s(harvest,iph)
            x_1=harvest
            !if(iph.eq.10) print*, '$$1$',x_temp,harvest
            call gasdev_s(harvest,iph)
            !if(iph.eq.10) print*, '$$2$',x_temp,harvest
            y_1=harvest
            z_1=x_1/y_1
            x_temp = x_freqD-a_rlw*z_1
         END DO

            y_temp=exp((x_freqD-0.25)**2-x_temp**2)
            !call ran1_s(harvest)
            y_random=rng_uniform(rngstate(iph))
        ENd DO

   END IF

END IF

vparallel = x_temp*x_freqD_sign
!x_freqD = x_freqD*x_freqD_sign

!if(iph.eq.10) print*, '$$$',x_temp,x_freqD_sign

  END SUBROUTINE scat_vpar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
!call gasdev_s(harv, iph)
        SUBROUTINE gasdev_s(harvest,iph)
        use rng, only: rng_uniform
        use common_v, only: rngstate
        IMPLICIT NONE
        integer, intent(IN) :: iph
        REAL, INTENT(OUT) :: harvest
        REAL :: rsq,v1,v2
        !REAL, SAVE :: g
        !LOGICAL, SAVE :: gaus_stored=.false.
        !if (gaus_stored) then
        !        harvest=g
        !        gaus_stored=.false.
        !else
                do
                     v1=rng_uniform(rngstate(iph))
                     v2=rng_uniform(rngstate(iph))
                     v1=2.0*v1-1.0
                     v2=2.0*v2-1.0
                     rsq=v1**2+v2**2
                     if (rsq > 0.0 .and. rsq < 1.0) exit
                end do
                rsq=sqrt(-2.0*log(rsq)/rsq)
                harvest=v1*rsq
        !        g=v2*rsq
        !        gaus_stored=.true.
        !end if
        END SUBROUTINE gasdev_s
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
SUBROUTINE voi_int (a_rlw, x_freqD, a,b,s)
   IMPLICIT NONE
  double precision, INTENT(IN):: a, b, x_freqD, a_rlw
  double precision, INTENT(OUT) :: s
  double precision :: voi_func,  os, ost, st
  INTEGER ::  j
  double precision, PARAMETER :: EPS=1d-3
  INTEGER, PARAMETER :: JMAX=63
      ost=-1.d30
      os= -1.d30
      do j=1, JMAX
        call vtrapzd(a_rlw, x_freqD, a, b, st, j)
        s=(4.d0*st-ost)/3.d0
        if (abs(s-os).le.EPS*abs(os).and.j.gt.7.or.j.gt.20) then
        return
        endif
        os=s
        ost=st
      end do

END SUBROUTINE voi_int
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
SUBROUTINE voigt(a, v, ans)
  IMPLICIT NONE

      double precision, INTENT(IN):: a, v
      double precision, INTENT(OUT) :: ans
      double precision, PARAMETER :: zero   = 0d0
      double precision, PARAMETER ::    pi=3.141592653d0
      double precision :: a_rlw, x_freqD, xx, z, q, H_a_x

     a_rlw = a
     x_freqD = v
     
     xx = x_freqD*x_freqD
     z = (xx-0.855)/(xx+3.42)
     
     if (z .le. zero) then
        q = zero
     else
        q = z*(1+(21/xx))*(a_rlw/(pi*(xx+1)))*(0.1117+z*(4.421+z*(5.674*z-9.207)))
     end if
     
     H_a_x = sqrt(pi)*(q+(exp(-xx)/1.77245385))

     ans = H_a_x
     
     END SUBROUTINE voigt
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
SUBROUTINE vtrapzd(a_rlw, x_freqD, a,b,s,n)

  IMPLICIT NONE
  double precision, INTENT(IN):: a, b, x_freqD, a_rlw
  double precision, INTENT(INOUT) :: s
  INTEGER, INTENT(In) :: n
  double precision :: voi_func, del,sum,tnm,x
  INTEGER :: it, j
  
  if (n.eq.1) then
     
        s=0.5d0*(b-a)*(voi_func(a, x_freqD, a_rlw)+voi_func(b, x_freqD, a_rlw))
      else
        it=2**(n-2)
        tnm=it
        del=(b-a)/tnm
        x=a+0.5d0*del
        sum=0.d0
        do 11 j=1,it

          sum=sum+voi_func(x, x_freqD, a_rlw)
          x=x+del
11      continue
        s=0.5d0*(s+(b-a)*sum/tnm)
      endif
      return
      
END SUBROUTINE vtrapzd


