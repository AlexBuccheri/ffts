!> Perform a 1D forward FFT on a real function, then backward transform the result to recover it.
!!
!! From https://numpy.org/doc/stable/reference/generated/numpy.fft.rfft.html#numpy.fft.rfft
!! When the DFT is computed for purely real input, the output is Hermitian-symmetric, 
!! i.e. the negative frequency terms are just the complex conjugates of the corresponding 
!! positive-frequency terms, and the negative-frequency terms are therefore redundant. 
!! This function does not compute the negative frequency terms, and the length of the 
!! transformed axis of the output is therefore n//2 + 1.
!!
!! If n is even, yc[n] contains the term representing both positive and negative Nyquist frequency 
!! (+fs/2 and -fs/2), and must also be purely real.
!
program fftw_real_example
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use fftw3_interface
  use constants_m
  implicit none

  integer(int32),  parameter        :: n = 600            !< Grid size
  real(real64),    dimension(n)     :: y                  !< Input transform    
  real(real64),    dimension(n)     :: y_recov            !< Result of backward transform     
  complex(real64), dimension(n/2+1) :: yc                 !< Result of forward transform
  type(c_ptr)                       :: plan_fwd, plan_bwd !< FFTW3 plans are c pointers
  integer                           :: i
  real(real64)                      :: dx, x_i, one_over_N_dx, phase, rad_to_degrees

  ! Initial real 1D function
  ! Using the function: https://docs.scipy.org/doc/scipy/tutorial/fft.html
  dx = 1.0_real64 / 800.0_real64
  rad_to_degrees = 180.0_real64 / acos(-1.0_real64)

  do i = 1, n
    x_i = 0.0_real64 + (i - 1) * dx
    y(i) = sin(50.0_real64 * 2.* pi * x_i) + (0.5 * sin(80.0_real64 * 2.* pi * x_i))
    write(100, '(F12.5, X, F12.5)') x_i, y(i)
  enddo
  
  ! Create FFTW plans
  plan_fwd = fftw_plan_dft_r2c_1d(n, y, yc, FFTW_ESTIMATE)
  plan_bwd = fftw_plan_dft_c2r_1d(n, yc, y_recov, FFTW_ESTIMATE)

  ! Execute forward real -> complex FFT
  call fftw_execute_dft_r2c(plan_fwd, y, yc)

  one_over_N_dx = 1.0_real64 / (real(n, real64) * dx)

  ! yc is a complex amplitude whose real part is the cosine‐coefficient and 
  ! whose imaginary part is the sine‐coefficient at frequencies defined below
  write(*, *) 'Forward FFT (R2C):'
  do i = 1, size(yc)
     ! Frequencies
     x_i = real(i - 1, real64) * one_over_N_dx
     ! Phase
     phase = atan2(aimag(yc(i)), real(yc(i))) !* rad_to_degrees
     ! Magnitude abs(yc(i)) == sqrt(re^2 + im^2)
     ! Where the scaling gives the true amplitude of your sinusoid components (rather than the raw FFT “counts”)
     ! Note that one should apply wrapping to remove the saw-tooth noise
     write(200, '(F12.5, X, F12.5, X, F12.5)') x_i, (2.0_real64 / real(n, real64)) * abs(yc(i)), phase
  end do

  ! Execute inverse complex -> real FFT
  call fftw_execute_dft_c2r(plan_bwd, yc, y_recov)

  ! Normalise and print the reconstructed real signal
  write(*, *) 'Inverse FFT (C2R) [normalised]:'
  do i = 1, n
     x_i = 0.0_real64 + (i - 1) * dx
     write(101, '(F12.5, X, F12.5)') x_i, y_recov(i) / real(n, real64)
  end do

  ! Clean up
  call fftw_destroy_plan(plan_fwd)
  call fftw_destroy_plan(plan_bwd)
  call fftw_cleanup()

end program fftw_real_example
