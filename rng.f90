module rng
  implicit none

  private
  public :: rng_t, rng_seed, rng_uniform

  ! Dimension of the state
  integer, parameter :: ns = 4

  ! Default seed vector
  integer, parameter, dimension(ns) :: default_seed &
       = (/ 521288629, 362436069, 16163801, 1131199299 /)

  ! A data type for storing the state of the RNG
  type :: rng_t
     integer, dimension(ns) :: state = default_seed
  end type rng_t

contains

  ! Seeds the RNG using a single integer and a default seed vector.
  subroutine rng_seed(self, seed)
    type(rng_t), intent(inout) :: self
    integer, intent(in) :: seed
    self%state(ns) = seed
    self%state(1:ns-1) = default_seed(2:ns)
  end subroutine rng_seed

  ! Draws a uniform real number on [0,1].
  function rng_uniform_old(self) result(u)
    type(rng_t), intent(inout) :: self
    real :: u
    integer :: imz

    imz = self%state(1) - self%state(3)

    if (imz < 0) imz = imz + 2147483579

    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz
    self%state(4) = 69069 * self%state(4) + 1013904243
    !print*,self%state,imz
    imz = imz + self%state(4)
    u = 0.5d0 + 0.23283064e-9 * imz
  end function rng_uniform_old

  function rng_uniform(self) result(u)
    type(rng_t), intent(inout) :: self
    real :: u
    integer :: imz

    ! Xorshift algorithm to improve randomness
    imz = self%state(4)
    imz = ieor(imz, ishft(imz, 13))
    imz = ieor(imz, ishft(imz, -17))
    imz = ieor(imz, ishft(imz, 5))
    self%state(4) = imz

    ! Combine with subtractive generator for better mixing
    imz = self%state(1) - self%state(3)
    if (imz < 0) imz = imz + 2147483579

    self%state(1) = self%state(2)
    self%state(2) = self%state(3)
    self%state(3) = imz

    ! Combine Xorshift and subtractive generator results
    imz = ieor(imz, self%state(4))

    ! Map to [0, 1)
    u = 0.5d0 + 0.23283064e-9 * imz
  end function rng_uniform

end module rng
