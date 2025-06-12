module fft_utils_m
    use, intrinsic :: iso_fortran_env
    implicit none
    private
    public :: fftfreq

contains

    subroutine fftfreq(n_points, spacing, x_frequencies, order)
        integer,          intent(in ) :: n_points
        real(real64),     intent(in ) :: spacing
        real(real64),     intent(out) :: x_frequencies(:)
        character(len=*), intent(in ), optional :: order !< Default is consistent with scipy ordering

        real(real64)      ::  one_over_n_dx
        integer           :: i, j, midpoint, offset
        character(len=10) :: ordering
 
        if (present(order)) then
            ordering = order
        else
            ordering = ''
        endif

        one_over_n_dx = 1.0_real64 / (real(n_points, real64) * spacing)

        ! Note, this is probably over-engineered, as I did not bother changing
        ! real(midpoint +1 ...) for odd N
        if (mod(n_points, 2) == 0) then
            ! Even
            midpoint = n_points / 2
            offset = 0
            ! Could add check on size of x_frequencies
        else
            ! Odd
            midpoint = (n_points -1) / 2
            offset = 1
            ! Could add check on size of x_frequencies
        endif

        if (trim(ordering) == 'size') then
            ! Negative frequencies
            do i = 1, midpoint + offset
                x_frequencies(i) = -real(midpoint + 1 - i , real64) * one_over_N_dx
            enddo 
            ! Followed by positive frequencies
            do i = 1, midpoint
                j = i + (midpoint + offset)
                x_frequencies(j) = real(i - 1, real64) * one_over_N_dx
            enddo 
        else
            ! Positive frequencies
            do i = 1, midpoint + offset
                x_frequencies(i) = real(i - 1, real64) * one_over_N_dx
            enddo 
            ! Followed by negative frequencies
            do i = 1, midpoint
                j = i + (midpoint + offset)
                x_frequencies(j) = -real(midpoint + 1 - i , real64) * one_over_N_dx
            enddo 
        endif

    end subroutine fftfreq

end module fft_utils_m