""" 1D FFT using Scipy
Reference: https://docs.scipy.org/doc/scipy/tutorial/fft.html
"""
from scipy.fft import fft, ifft, fftfreq, fftshift
import numpy as np
import matplotlib.pyplot as plt

# An asymmetric spectrum using a complex function
# Number of signal points
N = 400
# Sample spacing
T = 1.0 / 800.0
x = np.linspace(0.0, N*T, N, endpoint=False)
y = np.exp(50.0 * 1.j * 2.0*np.pi*x) + 0.5*np.exp(-80.0 * 1.j * 2.0*np.pi*x)

# Forward FFT
yf = fft(y)

# Positive frequency range, followed by negative frequency range
#[0, 400) then |-400, 0)
xf = fftfreq(N, T)

# Order from smallest number to largest number (for plotting)
# [-400, 400)]
xf = fftshift(xf)
yplot = fftshift(yf)

# Notice that there is no need for factor 2 in normalisation
plt.plot(xf, 1.0/N * np.abs(yplot))
plt.xlabel("Frequency (Hz)")
plt.ylabel("Amplitude |Y|")
plt.title("Amplitude Spectrum")
plt.grid()
plt.show()

# FFT back and plot against original
y_recovered = ifft(yf)
plt.plot(x, np.real(y), label='Original')
plt.plot(x, np.real(y_recovered), 'r--', label='Recovered')
plt.legend()
plt.show()
