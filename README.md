# FFTs and Poisson Prototyping

**Disclaimer:** This repository currently makes no effort to support any environment other than
the my local machine.

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

### Normalisation and Amplitudes with Scipy

Scipy defines their forward transform as:
$$
Y[k] = \sum_{n=0}^{N-1} x[n]\;e^{-2\pi i\,k n/N}.
$$  
Given a sin wave at $k_0$:
$$
x[n] = A\,\sin\!\Bigl(\frac{2\pi\,k_0\,n}{N}\Bigr),
$$
which can equivalently be expressed in terms of exponentials: 
$$
x[n] = \frac{A}{2i}\Bigl(e^{2\pi i\,k_0 n/N} - e^{-2\pi i\,k_0 n/N}\Bigr).
$$

Substituting this into the definition of the forward transform gives:

\[
\begin{aligned}
  Y[k]
  &= \sum_{n=0}^{N-1}
     \frac{A}{2i}\bigl(e^{2\pi i\,k_0 n/N}-e^{-2\pi i\,k_0 n/N}\bigr)\,
     e^{-2\pi i\,k n/N},\\
  &= \frac{A}{2i}\Bigl[
      \sum_{n=0}^{N-1} e^{-2\pi i\,n\,(k - k_0)/N}
     -\sum_{n=0}^{N-1} e^{-2\pi i\,n\,(k + k_0)/N}
    \Bigr].
\end{aligned}
\]

Using the geometric‐series identity:
   $$
     \sum_{n=0}^{N-1} e^{-2\pi i\,n\,m/N}
     =
     \begin{cases}
       N, & m \equiv 0 \pmod N,\\
       0, & \text{otherwise}
     \end{cases}
   $$
one can rewrite the forward transform in terms of delta functions (note $k_0 \equiv N-k_0$):

   $$
     Y[k]
     = \frac{A}{2i}\bigl(N\,\delta_{k,k_0} - N\,\delta_{k,N-k_0}\bigr)
     = \frac{N A}{2i}\,\bigl(\delta_{k,k_0} - \delta_{k,N-k_0}\bigr).
   $$

The magnitudes at $k_0$ and $-k_0$ are therefore: 

- At $k = k_0$:
  $$
    Y[k_0] = \frac{NA}{2i}
    \quad\Longrightarrow\quad
    \bigl|Y[k_0]\bigr| = \frac{NA}{2},
  $$

- At $k = N - k_0$:
  $$
    Y[N-k_0] = -\frac{NA}{2i}
    \quad\Longrightarrow\quad
    \bigl|Y[N-k_0]\bigr| = \frac{NA}{2},
  $$

hence the need to multiply by $2/N$ to correctly normalise with $A$. 

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
