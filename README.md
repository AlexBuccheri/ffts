# FFTs and Poisson Prototyping

Building the fortran code:

```shell
export CMAKE_PREFIX_PATH="/opt/homebrew/Cellar/fftw/3.3.10_2/:${CMAKE_PREFIX_PATH}"
$ cmake -B ./build -G Ninja --fresh
$ cmake --build ./build
```

<!-- TODO: Configure such that the test binaries are run by ctest
$ ctest --test-dir ./build -->

## Part 1. Fast Fourier Transforms Theory and Code Examples

Examples of 1D, 2D and 3D FFTs using scipy in python, and FFTW3 called from fortran, respectively.

### 1D Real Example

Python:

```shell
python python/real1d.py
```

Fortran:

```shell
./build/real1d
```

### 2D Example

### 3D Real Example

### 1D Complex Example

### 3D Complex Example

### Batched 3D Complex Example
