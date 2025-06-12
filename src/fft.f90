! module fft_m
!     implicit none

! contains    

! !> @brief Forward 3D FFT using FFTW3 
! subroutine fft_forward_3d(fft, in, out, norm)
!     type(fft_t),                 intent(in)    :: fft
!     R_TYPE,          contiguous, intent(inout) :: in(:,:,:)
!     complex(real64), contiguous, intent(out)   :: out(:,:,:)
!     real(real64),    optional,   intent(out)   :: norm

!     integer :: ii, jj, kk, slot, n1, n2, n3

!     if (all(fft_array(slot)%rs_n(1:3) >= 1)) then
!         ! Real    
!         call fftw_execute_dft_r2c(fft_array(slot)%planf, in(:,:,:), out(:,:,:))
!         ! Complex
!         !call fftw_execute_dft(fft_array(slot)%planf, in(:,:,:), out(:,:,:))
!     else
!         ii = min(1, fft_array(slot)%rs_n(1))
!         jj = min(1, fft_array(slot)%rs_n(2))
!         kk = min(1, fft_array(slot)%rs_n(3))
!         ! Real
!         call fftw_execute_dft_r2c(fft_array(slot)%planf, in(ii:,jj:,kk:), out(ii:,jj:,kk:))
!         ! Complex
!         !call fftw_execute_dft(fft_array(slot)%planf, in(ii:,jj:,kk:), out(ii:,jj:,kk:))
!     end if
   
!     call fft_operation_count(fft)

!     call profiling_out(TOSTRING(X(FFT_FORWARD)))

! end subroutine X(fft_forward_3d)

! end module fft_m