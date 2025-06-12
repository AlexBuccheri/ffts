module fft_utils_m
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: fftfreq

contains

    subroutine fftfreq(n_points, spacing, x_frequencies)
    implicit none
    integer, intent(in)       :: n_points
    real(real64), intent(in)  :: spacing
    real(real64), intent(out) :: x_frequencies(:)
    integer      :: i, midpoint
    real(real64) :: one_over_n_dx, n_points_real

    if (n_points == 1) then
        x_frequencies(1) = 0
        return
    endif

    n_points_real = real(n_points, real64)
    one_over_n_dx = 1.0_real64 / ( n_points_real * spacing )
    midpoint = (n_points + 1) / 2

    do i = 1, n_points
        x_frequencies(i) = real(i-1, real64) * one_over_n_dx
    end do

    ! Wrap the upper half into negatives
    do i = midpoint + 1, n_points
        x_frequencies(i) = x_frequencies(i) - (n_points_real * one_over_n_dx)
    end do

    end subroutine fftfreq

end module fft_utils_m