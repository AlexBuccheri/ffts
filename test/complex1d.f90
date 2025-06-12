program complex1d
  use, intrinsic :: iso_c_binding
  use, intrinsic :: iso_fortran_env

  use fftw3_interface
  use constants_m, only: pi
  use fft_utils_m
  implicit none

  complex(real64), parameter :: i_img = cmplx(0.0_real64, 1.0_real64)
  integer,         parameter :: N = 400

  ! TODO On should strictly be passing C types from iso_c_binding
  ! however, I appear to be getting away with it (most types map one to one). 
  ! One could also write an interface layer that maps fortran to C types using c_loc and c_f_pointer
  real(real64)    :: dx
  real(real64)    :: x(N), xf(N)
  complex(real64) :: y(N), yf(N), y_recovered(N)
  integer         :: i
  type(c_ptr)     :: plan_fwd, plan_bwd !< FFTW3 plans are c pointers

  ! Initialise complex function
  dx = 1.0_real64 / 800.0_real64
  do i = 1, N
    x(i) = 0.0 + (i-1) * dx
    y(i) = exp(50.0_real64 * i_img * 2.0_real64 * pi * x(i)) + &
      0.5 * exp(-80.0_real64 * i_img * 2.0_real64 * pi * x(i))
  enddo

! Create FFTW plans
  plan_fwd = fftw_plan_dft_1d(N, y, yf, FFTW_FORWARD, FFTW_PATIENT)
  plan_bwd = fftw_plan_dft_1d(N, yf, y_recovered, FFTW_BACKWARD, FFTW_PATIENT)

  ! Execute forward complex -> complex FFT
  ! Note, the default index-ordering for the FFT is consistent 
  ! between scipy and FFTW3, so fftfreq returns frequencies in the correct order
  call fftfreq(N, dx, xf)
  call fftw_execute_dft(plan_fwd, y, yf)

  ! Execute backward complex -> complex FFT
  call fftw_execute_dft(plan_bwd, yf, y_recovered)

  ! Output initial vs reconstructed spectrum (real part)
  ! y is asymmetric, returning over all frequencies, so the normalisation is 1/N
  do i = 1, n
     write(101, '(3(F12.5, X))') x(i), real(y(i)), real(y_recovered(i)) / real(N, real64)
  end do

  ! Output FFT amplitude spectrum
  do i = 1, n
     write(201, '(F12.5, X, F12.5)') xf(i), abs(yf(i)) / real(N, real64)
  end do
  
  call fftw_destroy_plan(plan_fwd)
  call fftw_destroy_plan(plan_bwd)
  call fftw_cleanup()

end program complex1d
