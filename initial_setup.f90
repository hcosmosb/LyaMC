subroutine initial_setup
use common_V
implicit none
character(len=1000) :: SCRATCH, infilebase
double precision :: MhDM, rhocM
integer :: i,j,k
real:: rcell
character(len=5) :: sh
real :: nHI_in

print*,'#########################'
print*,'#### Initial setup ######'
print*,'#########################'

zs=3.; print*, '# Redshift at the source:', real(zs)

call get_environment_variable('SCRATCH',SCRATCH)!; print*,SCRATCH; stop
infilebase='(empty)'
outputbase='results'; Lbox=32.; nHI_in=1e-8
Lbox=32.; nHI_in=1e-8; ndim = 256

print*, '# Output file base: '//trim(outputbase)
print*, '# Box size:', real(Lbox), ' (cMpc/h)'

allocate(nHI(0:ndim-1,0:ndim-1,0:ndim-1)); nHI=nHI_in
allocate(Temp(0:ndim-1,0:ndim-1,0:ndim-1)); Temp=2e4
allocate(vpe(0:ndim-1,0:ndim-1,0:ndim-1,3)); vpe=0.
print*,'# Data reading done.'

print*, '# Mesh size:', ndim
Omega_M=0.31; Omega_B=0.048; h=0.68
Hhub = h * 100. * 1d5/Mpc_in_cm &
     * sqrt(Omega_m*(1. + zs)**3 + (1. - Omega_M))
cell_cgs=Lbox*Mpc_in_cm/h/(1.+zs)/ndim
print*, '# Cell size:', real(cell_cgs), '(cm)'

print*,'#########################'
end subroutine initial_setup
