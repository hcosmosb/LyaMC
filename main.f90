module common_V
use rng, only:rng_t

real, parameter :: &
Msun=2e33, pi=3.141592, Grav=6.67e-8,Mpc_in_cm = 3.086e24, c_sl=2.99792e10,&
lambda_lya=1215.67, k_b = 1.380648e-16, m_p = 1.6726219e-24
real, allocatable :: nHI(:,:,:), vpe(:,:,:,:), Temp(:,:,:) 
real :: zs, Lbox, Omega_M, Omega_B, h, Rvir_cMpcoh, cell_cgs, Hhub
integer :: ndim
double precision :: rvi_g(3)
real :: kv_g(3), dva_g
type(rng_t), dimension(:), allocatable :: rngstate
character(len=500) :: outputbase
logical :: dbprint_flag=.false., PO_plane_flag=.false.
integer, parameter :: iu_po=10
real, allocatable :: SB_PO(:,:,:)
real,allocatable :: rvo(:,:),kvo(:,:),dvao(:)
real, allocatable :: output_data(:,:)
integer, parameter :: nodata=18
integer :: irun
character(len=200) :: crun
!$OMP THREADPRIVATE(rvi_g,kv_g,dva_g,dbprint_flag)
end module common_V
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program main
use common_v, only: PO_plane_flag
implicit none

!PO_plane_flag=.true.
PO_plane_flag=.false.
call check_input_arguments
call initial_setup
call run

end program main
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine check_input_arguments
use common_V, only: irun, crun
implicit none
integer :: cnt, len, status

cnt = command_argument_count ()
if (cnt.eq.0) then
print*,'# No argument found. Setting irun=1.'
irun=0
else if (cnt.eq.1) then
call get_command_argument (1, crun, len, status)
write (*,*) '# irun = ', crun(1:len)
read(crun(1:len),*)irun
else
print*,'# Too many arguments.(',cnt,') Pls, input only one.'
stop
endif

end subroutine check_input_arguments
