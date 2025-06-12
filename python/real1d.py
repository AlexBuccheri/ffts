""" 1D FFT using Scipy
Reference: https://docs.scipy.org/doc/scipy/tutorial/fft.html
"""
from scipy.fft import fft, ifft, fftfreq
import numpy as np
import matplotlib.pyplot as plt

# Number of sample points
N = 600
# sample spacing
T = 1.0 / 800.0

x = np.linspace(0.0, N*T, N, endpoint=False)
y = np.sin(50.0 * 2.0*np.pi*x) + 0.5*np.sin(80.0 * 2.0*np.pi*x)

# Forward transform
# Note, could be using rfft, which would only output complex
# coefficients for the positive frequencies
yf = fft(y)

# Frequencies for the forward transform
xf = fftfreq(N, T)

# which are explicitly defined for the positive frequencies as:
xf_explicit = np.empty(shape=(N//2))
for i in range(0, N//2):
    xf_explicit[i] = float(i) / (N * T)

np.testing.assert_allclose(xf[0:N//2], xf_explicit)

# Negative frequencies are returned from most negative to least negative i.e. -400.0 to -1.33
# However xf_explicit contains the first frequency (0) and misses the last (400)
negative_xf = np.empty_like(xf_explicit)
negative_xf[0] = xf[N//2]
negative_xf[1:] = -xf_explicit[::-1][:-1]
np.testing.assert_allclose(xf[N//2:N], negative_xf)

# Single-sided spectrum
# For a real input x[n], the FFT is conjugate‐symmetric:
#  Y[N-k] = Y[k]^* i.e. Y[0] = Y[N]^*, Y[1] = Y[N-1]^*, etc 
#
# Normalisation
#  2/N is required because of the definition of the raw FFT.
#  This can be seen when passing x[n] = A sin(2 pi k_0 n/N). The resulting FFT has
#  a magnitude NA/2, so one must multiply by 2/N, where
#  the 2 then puts both the positive- and negative-frequency halves back together.
#
# Plotting the norm (with abs) i.e. |yf|
#  Complex coefficients, yf, encode the amplitude and the phase
plt.plot(xf[0:N//2], 2.0/N * np.abs(yf[0:N//2]), label='Real frequencies')
plt.plot(xf[N//2:N], 2.0/N * np.abs(yf[N//2:N]), label='Imaginary frequencies')
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude |Y|")
plt.title("Single-sided Amplitude Spectrum")
plt.legend()
plt.grid()
plt.show()

# Phase
# Return the angle of the complex argument: Pass the imaginary and real parts of the argument to arctan2
phase = np.angle(yf)
# removes the ±2π discontinuities
phase_unwrapped = np.unwrap(phase)
plt.figure()
plt.plot(xf[0:N//2], phase[0:N//2], label="Raw Phase")
plt.plot(xf[0:N//2], phase_unwrapped[0:N//2], label="Unwrapped Phase")
plt.xlabel('Frequency (Hz)')
plt.ylabel('Phase (radians)')
plt.title('Single-Sided Phase Spectrum')
plt.legend()
plt.grid(True)
plt.show()

# FFT back and plot against original
y_recovered = ifft(yf)
plt.plot(x, y, label='Original')
plt.plot(x, np.real(y_recovered), 'r--', label='Recovered')
plt.xlabel("Time")
plt.ylabel("y = sin(50 2pi x) + 0.5 sin(80 2pi x)")
plt.legend()
plt.grid()
plt.show()
