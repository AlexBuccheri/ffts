cmake_minimum_required(VERSION 3.25)

project(fft_prototyping VERSION 1.0
		DESCRIPTION "Examples of FFTs, and prototyped batching"
		LANGUAGES C Fortran)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_COLOR_DIAGNOSTICS ON)
set(CMAKE_MESSAGE_CONTEXT_SHOW ON)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif ()

# Build all fortran modules in a common module directory
if (NOT CMAKE_Fortran_MODULE_DIRECTORY)
	set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
endif ()
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

# Load and execute CMake code from provided modules (probs not needed)
include(GNUInstallDirs)

# Add custom modules in cmake/
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

#[==============================================================================================[
#                                 Compiler Flags and Options                                             #
#          Note, flags "should" be specified in cmake presets rather in CMakeLists.txt                   #
#]==============================================================================================]

# Not portable. Added for convenience
# Avoid picking up Apple Clang
set(CMAKE_C_COMPILER "gcc-14")
# Brew path to FFTW3
# set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/opt/homebrew/Cellar/fftw/3.3.10_2" CACHE PATH "Prefix path" FORCE)

# gfortran flags
set(GCC_DEBUG
     -O0                                          # 0 optimisation. Probably implicit
     -fbacktrace                                  # symbolic stack traceback
     -fcheck=all                                  # Enable all run-time test of -fcheck: array-temps, bits, bounds, do, mem, pointer, recursion
     -Wall                                        # Enable all warnings
 #   -ffpe-trap=invalid,zero,overflow,underflow   # control over floating-point exception
 #   -finit-real=nan                              # All real scalars are initialised to NaN
    )
                    
set(GCC_RELEASE
    -O3
   )

# Interface to control scope of compiler flags and options
add_library(common_flags INTERFACE)

target_compile_options(common_flags INTERFACE 
   $<$<AND:$<CONFIG:Debug>,$<Fortran_COMPILER_ID:GNU>>:${GCC_DEBUG}>
   $<$<AND:$<CONFIG:Release>,$<Fortran_COMPILER_ID:GNU>>:${GCC_RELEASE}>
   )

# --------------------------------------------
# External Dependencies
# --------------------------------------------

# BLAS/lapack
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# OpenMP
find_package(OpenMP COMPONENTS Fortran)

# FFTW3
find_package(FFTW)

# MPI-dependent settings
option(MPI "Use MPI" OFF)
if(MPI)
   find_package(MPI REQUIRED Fortran C)
   # Set preprocessor variable 
   target_compile_definitions(common_flags INTERFACE USE_MPI)
else()
   message("-- Not using MPI")
endif()

# --------------------------------------------
# Project Targets
# --------------------------------------------
option(BUILD_SHARED_LIBS "Build using shared libraries" ON)

# Library
add_library(fft_proto_lib)

set_target_properties(fft_proto_lib PROPERTIES
		VERSION ${PROJECT_VERSION}
		SOVERSION ${PROJECT_VERSION_MAJOR}
		EXPORT_NAME fft_proto
		OUTPUT_NAME fft_proto
		)

# Location of source - added after defining library variables,
# as they're referenced in the source-level CMakeLists.txt
add_subdirectory(src)

# Link external dependencies to target
# Fine to have all these dependencies (and flags) propagate to the tests
target_link_libraries(fft_proto_lib PUBLIC 
    common_flags
    BLAS::BLAS
    LAPACK::LAPACK)

# Neither are direct dependencies of the caller
# OpenMP::OpenMP_Fortran 
# $<$<BOOL:${MPI}>:MPI::MPI_Fortran> 

# FFTW::OpenMP/Threads is a wrapper lib. Linking to the main library is still necessary
if (TARGET FFTW::DoubleOpenMP)
	target_link_libraries(fft_proto_lib PUBLIC FFTW::DoubleOpenMP)
    option(HAVE_FFTW3_THREADS "Use FFTW3 threading" ON)
elseif ()
    # TODO CHECK THIS
	target_link_libraries(fft_proto_lib PUBLIC FFTW::DoubleThreads)
    option(HAVE_FFTW3_THREADS "Use FFTW3 threading" OFF)
endif ()
target_link_libraries(fft_proto_lib PUBLIC FFTW::Double)

option(HAVE_FFTW3_MPI "Use FFTW3 MPI" OFF)

target_compile_definitions(common_flags INTERFACE HAVE_FFTW3_THREADS)
target_compile_definitions(common_flags INTERFACE HAVE_FFTW3_MPI)

# Tests
add_executable(real1d)
add_executable(complex1d)

add_subdirectory(test)

target_link_libraries(real1d PRIVATE fft_proto_lib)
target_link_libraries(complex1d PRIVATE fft_proto_lib)
