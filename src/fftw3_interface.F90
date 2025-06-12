!> @brief Full module layer to FFTW3 include
!! TODO Add the comprehensive list of MPI routines
module fftw3_interface
  use, intrinsic :: iso_c_binding
  implicit none
  private

! Constants/enums
public :: &
  C_FFTW_R2R_KIND, &
  FFTW_R2HC, &
  FFTW_HC2R, &
  FFTW_DHT, &
  FFTW_REDFT00, &
  FFTW_REDFT01, &
  FFTW_REDFT10, &
  FFTW_REDFT11, &
  FFTW_RODFT00, &
  FFTW_RODFT01, &
  FFTW_RODFT10, &
  FFTW_RODFT11, &
  FFTW_FORWARD, &
  FFTW_BACKWARD, &
  FFTW_MEASURE, &
  FFTW_DESTROY_INPUT, &
  FFTW_UNALIGNED, &
  FFTW_CONSERVE_MEMORY, &
  FFTW_EXHAUSTIVE, &
  FFTW_PRESERVE_INPUT, &
  FFTW_PATIENT, &
  FFTW_ESTIMATE, &
  FFTW_WISDOM_ONLY, &
  FFTW_ESTIMATE_PATIENT, &
  FFTW_BELIEVE_PCOST, &
  FFTW_NO_DFT_R2HC, &
  FFTW_NO_NONTHREADED, &
  FFTW_NO_BUFFERING, &
  FFTW_NO_INDIRECT_OP, &
  FFTW_ALLOW_LARGE_GENERIC, &
  FFTW_NO_RANK_SPLITS, &
  FFTW_NO_VRANK_SPLITS, &
  FFTW_NO_VRECURSE, &
  FFTW_NO_SIMD, &
  FFTW_NO_SLOW, &
  FFTW_NO_FIXED_RADIX_LARGE_N, &
  FFTW_ALLOW_PRUNING

#ifdef HAVE_FFTW3_MPI
  public :: FFTW_MPI_DEFAULT_BLOCK
#endif

! Types  
public :: &
  fftw_iodim, &
  fftw_iodim64

! Subroutines
public :: &
  fftw_execute_dft, &  
  fftw_execute_split_dft, &  
  fftw_execute_dft_r2c, &  
  fftw_execute_dft_c2r, &  
  fftw_execute_split_dft_r2c, &  
  fftw_execute_split_dft_c2r, &  
  fftw_execute_r2r, &  
  fftw_destroy_plan, &  
  fftw_forget_wisdom, &  
  fftw_cleanup, &  
  fftw_set_timelimit, &  
  fftw_plan_with_nthreads, &  
  fftw_cleanup_threads, &  
  fftw_make_planner_thread_safe, &  
  fftw_export_wisdom_to_file, &  
  fftw_export_wisdom, &  
  fftw_fprint_plan, &  
  fftw_print_plan, &  
  fftw_free, &  
  fftw_flops, &  
  fftwf_execute_dft, &  
  fftwf_execute_split_dft, &  
  fftwf_execute_dft_r2c, &  
  fftwf_execute_dft_c2r, &  
  fftwf_execute_split_dft_r2c, &  
  fftwf_execute_split_dft_c2r, &  
  fftwf_execute_r2r, &  
  fftwf_destroy_plan, &  
  fftwf_forget_wisdom, &  
  fftwf_cleanup, &  
  fftwf_set_timelimit, &  
  fftwf_export_wisdom_to_file, &  
  fftwf_export_wisdom, &  
  fftwf_fprint_plan, &  
  fftwf_print_plan, &  
  fftwf_free, &  
  fftwf_flops

#if defined(HAVE_FFTW3_THREADS)
  public ::                  &
    fftw_init_threads,       &
    fftw_plan_with_nthreads, &
    fftw_cleanup_threads,    &
    fftwf_plan_with_nthreads, &
    fftwf_cleanup_threads,    &
    fftwf_make_planner_thread_safe
#endif

! Functions
public :: &
  fftw_plan_dft, &
  fftw_plan_dft_1d, &
  fftw_plan_dft_2d, &
  fftw_plan_dft_3d, &
  fftw_plan_many_dft, &
  fftw_plan_guru_dft, &
  fftw_plan_guru_split_dft, &
  fftw_plan_guru64_dft, &
  fftw_plan_guru64_split_dft, &
  fftw_plan_many_dft_r2c, &
  fftw_plan_dft_r2c, &
  fftw_plan_dft_r2c_1d, &
  fftw_plan_dft_r2c_2d, &
  fftw_plan_dft_r2c_3d, &
  fftw_plan_many_dft_c2r, &
  fftw_plan_dft_c2r, &
  fftw_plan_dft_c2r_1d, &
  fftw_plan_dft_c2r_2d, &
  fftw_plan_dft_c2r_3d, &
  fftw_plan_guru_dft_r2c, &
  fftw_plan_guru_dft_c2r, &
  fftw_plan_guru_split_dft_r2c, &
  fftw_plan_guru_split_dft_c2r, &
  fftw_plan_guru64_dft_r2c, &
  fftw_plan_guru64_dft_c2r, &
  fftw_plan_guru64_split_dft_r2c, &
  fftw_plan_guru64_split_dft_c2r, &
  fftw_plan_many_r2r, &
  fftw_plan_r2r, &
  fftw_plan_r2r_1d, &
  fftw_plan_r2r_2d, &
  fftw_plan_r2r_3d, &
  fftw_plan_guru_r2r, &
  fftw_plan_guru64_r2r, &
  fftw_export_wisdom_to_string, &
  fftw_sprint_plan, &
  fftw_malloc, &
  fftw_alloc_real, &
  fftw_alloc_complex, &
  fftwf_plan_dft, &
  fftwf_plan_dft_1d, &
  fftwf_plan_dft_2d, &
  fftwf_plan_dft_3d, &
  fftwf_plan_many_dft, &
  fftwf_plan_guru_dft, &
  fftwf_plan_guru_split_dft, &
  fftwf_plan_guru64_dft, &
  fftwf_plan_guru64_split_dft, &
  fftwf_plan_many_dft_r2c, &
  fftwf_plan_dft_r2c, &
  fftwf_plan_dft_r2c_1d, &
  fftwf_plan_dft_r2c_2d, &
  fftwf_plan_dft_r2c_3d, &
  fftwf_plan_many_dft_c2r, &
  fftwf_plan_dft_c2r, &
  fftwf_plan_dft_c2r_1d, &
  fftwf_plan_dft_c2r_2d, &
  fftwf_plan_dft_c2r_3d, &
  fftwf_plan_guru_dft_r2c, &
  fftwf_plan_guru_dft_c2r, &
  fftwf_plan_guru_split_dft_r2c, &
  fftwf_plan_guru_split_dft_c2r, &
  fftwf_plan_guru64_dft_r2c, &
  fftwf_plan_guru64_dft_c2r, &
  fftwf_plan_guru64_split_dft_r2c, &
  fftwf_plan_guru64_split_dft_c2r, &
  fftwf_plan_many_r2r, &
  fftwf_plan_r2r, &
  fftwf_plan_r2r_1d, &
  fftwf_plan_r2r_2d, &
  fftwf_plan_r2r_3d, &
  fftwf_plan_guru_r2r, &
  fftwf_plan_guru64_r2r, &
  fftwf_export_wisdom_to_string, &
  fftwf_sprint_plan, &
  fftwf_malloc, &
  fftwf_alloc_real, &
  fftwf_alloc_complex

#ifdef HAVE_FFTW3_MPI
  ! The MPI header already includes the serial header below
  include "fftw3-mpi.f03"
#else
  include "fftw3.f03"
#endif

end module fftw3_interface
