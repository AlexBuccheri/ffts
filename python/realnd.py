""" Test 2D and 3D FFTs
"""
from scipy.fft import fftfreq, fftshift, fftn, ifftn
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import windows
import plotly.graph_objects as go


def gaussian(XYZ:list | tuple, sigma:list | tuple):
    """ND Gaussian

    FFT expects an array (function) with the same shape as the grid i.e. (Nx, Ny, Nz)

    Args:
        XYZ (list): List or array of arrays produced by np.meshgrid
        sigma (_type_): Broadening per dimension
    Return:
       Gaussian: shape XYZ[i].shape == (Nx, Ny, Nz) (depending on the dimensionality of 
       the input)
    """
    assert len(XYZ) == len(sigma)
    exponent = (XYZ[0] / sigma[0])**2
    for i in range(1, len(XYZ)):
        exponent += (XYZ[i] / sigma[i])**2
    return np.exp(-0.5 * exponent)


def grid_wrapper(n_points, limits):

    assert len(n_points) == len(limits)

    xyz_values = []
    spacings = []
    for i in range(len(n_points)):
        vals, spacing = np.linspace(limits[i][0], limits[i][1], n_points[i], 
                                    endpoint=False, retstep=True)
        xyz_values.append(vals)
        spacings.append(spacing)

    # Produce a tuple with each entry of shape (n_points[0], n_points[1], ..., n_points[n_dim])
    XYZ = np.meshgrid(*xyz_values, indexing='xy')
    return XYZ, spacings


def gaussian_fft_2d():
    n_points = [400, 400]
    limits = [[-1.0, 1.0], [-1.0, 1.0]]
    XYZ, spacings = grid_wrapper(n_points, limits)
    dr = np.prod(spacings)

    # No need for a window
    sigma = np.array([0.1, 0.1])

    # Using a window
    # sigma = np.array([1.0, 1.0])
    # wy = windows.hann(n_points[1])
    # wx = windows.hann(n_points[0])
    # W2d = np.outer(wy, wx)

    g_real = gaussian(XYZ, sigma) #* W2d
    assert g_real.shape == (n_points[0], n_points[1])

    # Plot function in real space
    extent = []
    for limit in limits:
        extent += limit
    plt.imshow(g_real, extent=extent, origin='lower')
    plt.colorbar()
    plt.show()

    # Forward transform
    # xyz_f = []
    # for i in range(0, len(n_points)):
    #     xf = fftfreq(n_points[i], spacings[i])
    #     xyz_f.append(fftshift(xf))

    g_forward = fftn(g_real)
    assert g_forward.shape == (n_points[0], n_points[1])

    # Plot function in Fourier space
    # Note, one can multiply fftfreq by 2*np.pi to convert to angular frequency
    extent = []
    for i in range(0, len(n_points)):
        xf = fftfreq(n_points[i], spacings[i])
        # Reorder from negative frequencies to positive frequencies
        xf = fftshift(xf)
        extent += [xf[0], xf[1]]
    
    g_plot = np.abs(fftshift(g_forward)) / (np.prod(n_points))

    plt.figure(figsize=(6,5))
    plt.title('Log-scale magnitude of FFT')
    plt.imshow(np.log10(g_plot),
        extent=extent,
        origin='lower',
        aspect='equal'
    )
    plt.show()

    # FFT back and plot against original
    g_recovered = ifftn(g_forward)

    g_diff = g_real - np.real(g_recovered)
    assert np.all(g_diff < 1.e-10), "Back-transformed Gaussian disagrees with original"

    extent = []
    for limit in limits:
        extent += limit
    plt.imshow(np.real(g_recovered), extent=extent, origin='lower')
    plt.colorbar()
    plt.show()


def parseval_theorem(y, y_forward) -> bool:
    lhs = np.sum(np.abs(y)**2)
    rhs = np.sum(np.abs(y_forward)**2) / y_forward.size
    return np.isclose(lhs, rhs)


def gaussian_fft_3d():
    plotting = True
    n_points = [100, 100, 100]
    limits = [[-1.0, 1.0], [-1.0, 1.0], [-1.0, 1.0]]
    XYZ, spacings = grid_wrapper(n_points, limits)
    sigma = np.array([0.1, 0.1, 0.1])
    
    g_real = gaussian(XYZ, sigma) 
    assert g_real.shape == (n_points[0], n_points[1], n_points[2])

    # Forward transform
    g_forward = fftn(g_real)
    assert g_forward.shape == (n_points[0], n_points[1], n_points[2])
    assert parseval_theorem(g_real, g_forward), "Energy not conserved"
    assert np.isclose(np.sum(g_real), g_forward[0, 0, 0]), "Sum of G should equal the 0-frequency component of the forward transform"

    # Backward transform
    g_recovered = ifftn(g_forward)

    g_diff = g_real - np.real(g_recovered)
    assert np.all(g_diff < 1.e-10), "Back-transformed Gaussian disagrees with original"

    if plotting:
        # Plot initial function (quite slow if a larger grid is used)
        # fig = go.Figure(data=go.Isosurface(
        #     x=XYZ[0].ravel(), y=XYZ[1].ravel(), z=XYZ[2].ravel(),
        #     value=g_real.ravel(),
        #     isomin=0.5 * g_real.max(),
        #     isomax=g_real.max(),
        #     surface_count=2,   # number of nested shells
        #     caps=dict(x_show=False, y_show=False)
        # ))
        # fig.update_layout(title="3D Gaussian Isosurface")
        # fig.show()

        # Build k-axes for plotting, using angular frequency
        kxyz = []
        for i in range(0, len(n_points)):
            kxyz.append(2*np.pi * fftfreq(n_points[i], d=spacings[i]))
        KX, KY, KZ = np.meshgrid(*kxyz, indexing='ij')

        # # Fourier transform plot
        # # Note: Does not show anything  - probably the scale
        # gf_mag = np.abs(g_forward) / np.prod(n_points)
        # fig = go.Figure(data=go.Isosurface(
        #     x=KX.ravel(), y=KY.ravel(), z=KZ.ravel(),
        #     value=gf_mag.ravel(),
        #     isomin= 0.5 * gf_mag.max(),
        #     isomax=gf_mag.max(),
        #     surface_count=2
        # ))
        # fig.update_layout(title="3D FFT Magnitude Isosurface")
        # fig.show()

        # Radial Plot
        r  = np.sqrt(KX**2 + KY**2 + KZ**2).ravel()
        gf_mag = np.abs(g_forward) / np.prod(n_points)

        rbins = np.linspace(0, r.max(), 100)
        which = np.digitize(r, rbins)
        radial_mean = [gf_mag.ravel()[which==i].mean() for i in range(1,len(rbins))]

        plt.plot(rbins[1:], radial_mean)
        plt.xlabel('||k||')
        plt.ylabel('⟨|Ĝ|⟩')
        plt.title('Radial Profile of 3D FFT')
        plt.show()


def sin_cos_fft():
    n_points = [64, 64, 64]
    limits = [[0.0, 1.0], [0.0, 1.0], [0.0, 1.0]]
    XYZ, spacings = grid_wrapper(n_points, limits)

    # Frequencies (in cycles per unit length)
    fx, fy, fz = 5, 10, 15

    n_peaks = 6
    y = np.sin(2 *np.pi * fx *XYZ[0]) + np.cos(2 *np.pi * fy *XYZ[1]) + np.sin(2 *np.pi * fz *XYZ[2])

    # Three cos/sin waves, each with unit amplitude
    assert np.isclose(np.max(y), 3) and np.isclose(np.min(y), -3), "Max amplitude of y"
    # Integrals of cos(2pi n) and sin(2pi n) over [0, 1] = 0   
    assert np.isclose(np.mean(y), 0), "Mean of y"

    # Forward transform, with normalisation
    y_f = fftn(y) / np.prod(n_points)

    # This component represents the average value of the signal over the entire domain:
    # i.e. k = 0, so the transform is just the sum of the real-space signal
    # If the real-space mean is zero, so is the reciprocal space mean.
    assert np.isclose(y_f[0,0,0], 0.0), "y has a mean of zero, implying the zero-frequency component should also be zero"

    # --------------------------------------
    # Frequency points in reciprocal space
    # y is a sum of 3 pure harmonics => each sinusoid corresponds to a sharp peak in the frequency domain
    # at ± fx, fy, fz
    # --------------------------------------
    kxyz = [fftfreq(n, d=s) for n, s in zip(n_points, spacings)]
    KX, KY, KZ = np.meshgrid(*kxyz, indexing='ij')

    y_f_abs = np.abs(y_f)
    # Flatten from (n_points[0], n_points[1], n_points[2]) to 1D product(n_points)
    y_f_flat = y_f_abs.ravel()
    assert y_f_flat.shape == np.prod(n_points)

    # Partially sort the array
    inds = np.argpartition(y_f_flat, -n_peaks)[-n_peaks:]

    # Given indices from a ravelled (flattened array), get their corresponding (ix, iy, iz) indices
    ix, iy, iz = np.unravel_index(inds, y_f_abs.shape)
    for i, j, k in zip(ix, iy, iz):
        print(f"peak @ grid ({i},{j},{k}) → "
            f"k = ({kxyz[0][i]:.3f}, {kxyz[1][j]:.3f}, {kxyz[2][k]:.3f})")

    # ----------------------------------------
    # 3D scatter plot of FFT peaks
    # ----------------------------------------
    scatter_plot = False
    if scatter_plot:
        # Flatten and mask top peaks
        kx_flat = KX.ravel()
        ky_flat = KY.ravel()
        kz_flat = KZ.ravel()

        # Plotting 0.1% of the 64*64*64 points = 262
        # Should contain n_peaks
        threshold = np.percentile(y_f_flat, 99.9)
        mask = y_f_flat > threshold
        kx_strong = kx_flat[mask]
        ky_strong = ky_flat[mask]
        kz_strong = kz_flat[mask]
        y_strong = y_f_flat[mask]

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        p = ax.scatter(kx_strong, ky_strong, kz_strong, c=y_strong, s=20, cmap='viridis')
        fig.colorbar(p, label='FFT Magnitude')

        ax.set_title("3D Scatter of Strongest FFT Peaks")
        ax.set_xlabel("KX (cycles/unit)")
        ax.set_ylabel("KY (cycles/unit)")
        ax.set_zlabel("KZ (cycles/unit)")
        plt.tight_layout()
        plt.show()

    # 1D radial plot
    radial_plot = False
    if radial_plot:
        k_norm = np.sqrt(KX**2 + KY**2 + KZ**2).ravel()
        energy = (np.abs(y_f_abs)**2).ravel()

        # Radial bins
        num_bins = 100
        bins = np.linspace(0, np.max(k_norm), num_bins + 1)
        bin_centers = 0.5 * (bins[:-1] + bins[1:])
        bin_indices = np.digitize(k_norm, bins)

        # Sum FFT energy in each radial bin
        radial_energy = np.zeros(num_bins)
        # counts = np.zeros(num_bins)

        for i in range(1, num_bins + 1):
            # Index mask for bin i
            in_bin = bin_indices == i
            radial_energy[i - 1] = np.sum(energy[in_bin])
            # counts[i - 1] = np.sum(in_bin)

        plt.figure(figsize=(7, 4))
        plt.plot(bin_centers, radial_energy)
        plt.xlabel("Frequency Norm ||k|| (cycles/unit)")
        plt.ylabel("Total FFT Energy (|FFT|²)")
        plt.title("Radial Energy Spectrum")
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    # Slice not particularly informative
    # mag_max_z = np.max(y_f_abs, axis=2)  # Collapse along Z

    # plt.figure(figsize=(6, 5))
    # plt.imshow(np.fft.fftshift(mag_max_z), origin='lower', extent=[
    #     kxyz[1][0], kxyz[1][-1], kxyz[0][0], kxyz[0][-1]
    # ])
    # plt.colorbar(label='Max FFT Magnitude (over KZ)')
    # plt.title("Max-Intensity Projection (KZ)")
    # plt.xlabel("KY (cycles/unit)")
    # plt.ylabel("KX (cycles/unit)")
    # plt.tight_layout()
    # plt.show()

    # Transform back and show backward transform agrees with the original signal
    y_recovered = ifftn(y_f * np.prod(n_points))
    diff = y - np.real(y_recovered)
    assert np.all(diff < 1.e-10), "Back-transformed function disagrees with original signal"
    # Note, forgetting the normalisation can lead to errors ~e-6. Look small, but really are errors

if __name__ == '__main__':
    sin_cos_fft()
