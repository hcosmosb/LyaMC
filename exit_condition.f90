subroutine exit_condition(rv,exit_flag)
use common_v, only: Lbox
implicit none
real, intent(IN) :: rv(3)
logical, intent(OUT) :: exit_flag
real :: ebfac=0.49

exit_flag=.false.
if (sum(rv**2).gt.ebfac**2)exit_flag=.true.

!print*,'!',sum(rv**2)**0.5
end subroutine exit_condition
