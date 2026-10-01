
import numpy as np


import numba
from numba import jit

# def jit(**k):
#    return lambda x: x  # deactivate numba


# This computes the reflection for a stack of multiple incoherent layers.
# A convenience function is provided to compute the reflection of a stack of pair of layers

# References:
# https://hal.science/hal-01155614/document

# https://arxiv.org/pdf/1603.02720.pdf

VPOL = 0
HPOL = 1


def compute_multi_alternating_layer_reflection2(eps_0, eps_1, eps_2, npair, kd, mu):
    """the snowpack is composed of a top layer (usually air, eps_0=1) and a number (npair) of pairs of layers with
       permittivities eps_1 and eps_2. The thickness of each layer is 'thickness'. Hence, the total thickness of 
       the snowpack is 2 * npair * thickness. Absorption is taken into account through the imaginary part of the permittivities

    :param eps_0: permittivity of the top layer (usually air=1)
    :param eps_1: first permittivity of the alternative layers
    :param eps_2: second permittivity of the alternative layers
    :param npair: number of pairs of layers
    :param kd: thickness of the layers relative * wavenumber (k = 2 pi / wavelength)
    :param mu: cosine of the zenith angle (incidence angle)
"""
    mu = np.atleast_1d(mu)

    mu1 = snell_angle(eps_0, eps_1, mu)
    # kd0 = 1 * 2 * np.pi * 1.4e9 / 3e8       # For SMOS (1.4 GHz) with a 1 m thick top layer

    L = combine([forward_matrix(eps_0, eps_1, mu, 0), matrix_power(forward_matrix(eps_1, eps_2, mu1, kd), 2 * npair)])
    # L = combine([forward_matrix(eps_0, eps_1, mu, 0), forward_matrix(eps_0, eps_1, mu1, kd0),
    #               matrix_power(forward_matrix(eps_2, eps_3, mu1, kd), 2 * npair)])
    
    Rv, Rh = compute_reflection(L)
    return Rv, Rh


@jit(nopython=True, cache=True)
def forward_matrix(eps1, eps2, mu, kd, remove_transmission=True):
    """compute the operator to go from layer 1 to 2.

    :param remove_transmission: the scaling by 1 / transmission is removed because it is infinite for total refractions,
        causing +inf. This coefficient is only necessary to compute the total transmission through the slab, the reason
        why it is necessary to exclude it from the calculation when only the reflection coefficient is needed.

    """

    mu = np.atleast_1d(mu)

    rv, rh, mu2 = fresnel_coefficients_maezawa09_rigorous(eps1, eps2, mu)

    r = np.stack((rv.real**2 + rv.imag**2,
                  rh.real**2 + rh.imag**2))

    trans_v = np.exp(-2 * np.sqrt(eps2).imag * kd / mu)  # power attenuation
    trans_v = np.stack((trans_v, trans_v))
    # trans_v = trans_v[np.newaxis, :]  # P = np.stack((P, P))

    # F = np.array([[np.ones_like(r) / P, -r / P],
    #               [r * P, (1 - 2 * r) * P]])

    F = np.empty((2, 2, 2, len(mu)))
    F[0, 0] = 1 / trans_v
    F[0, 1] = -r * trans_v
    F[1, 0] = r / trans_v
    F[1, 1] = (1 - 2 * r) * trans_v

    # if not remove_transmission:
    #     F /= (1 - r)[None, None, :, :]

    # save the layer transmission

    # return F, trans_v
    return F


@jit(nopython=True, cache=True)
def combine(Fs, axis=-1):

    Fout = Fs[0]

    for F in Fs[1:]:
        # axes = [(0, 1), (0, 1), (0, 1)]  # assume the first two dimensions are the matrix dimension.
        # Fout = np.matmul(Fout, F, axes=axes)
        Fout = matmul2(Fout, F)
    return Fout


@jit(nopython=True, cache=True)
def matrix_power(a, n, axis=-1):
    # taken from matrix_power in numpy but use the broadcast

    # axes = [(0, 1), (0, 1), (0, 1)]  # assume the first two dimensions are the matrix dimension.
    assert n >= 3

    z = result = None

    while n > 0:
        # z = a if z is None else np.matmul(z, z, axes=axes)
        z = a if z is None else matmul2(z, z)
        n, bit = divmod(n, 2)
        if bit:
            # result = z if result is None else np.matmul(result, z, axes=axes)
            result = z if result is None else matmul2(result, z)

    return result


def compute_reflection(L):

    # L, trans_v = L
    return L[1, 0, VPOL, :] / L[0, 0, VPOL, :], L[1, 0, HPOL, :] / L[0, 0, HPOL, :]


# extracted from SMRT
@jit(nopython=True, cache=True)
def fresnel_coefficients_maezawa09_rigorous(eps_1, eps_2, mu1, full_output=False):
    """compute the reflection in two polarizations (H and V) for lossly media with the "rigorous Fresnel" based
    on Maezawa, H., & Miyauchi, H. (2009). Rigorous expressions for the Fresnel equations at interfaces between absorbing media. 
    Journal of the Optical Society of America A, 26(2), 330. https://doi.org/10.1364/josaa.26.000330

    The 'rigorous' derivation respect the energy conservation even for strongly loosly media.

    :param eps_1: permittivity of medium 1.
    :param eps_2: permittivity of medium 2.
    :param mu1: cosine zenith angle in medium 1.

    :returns: rv, rh, mu2 the cosine of the angle in medium 2
"""
    # y is the axis normal to the interface (usually z, but here it is y!)

    # incident wavenumber
    n1 = np.sqrt(eps_1)
    kiz2 = n1.real**2 * (1 - mu1**2)   # this the square of kiz = n1 * sin(theta)
    kyi = - np.sqrt(eps_1 - kiz2)                  # Eq 8 for i

    ktz2 = kiz2   # unnumbered equation before 22  -> tangential k is conserved throught the interface (=Snell law)
    kyt = - np.sqrt(complex(eps_2) - ktz2)                  # Eq 8 for t

    rh = (kyi - kyt) / (kyi.conjugate() + kyt)  # Eq 59

    rv = n1.conjugate() * (eps_2 * kyi - eps_1 * kyt) / (n1 * (eps_2 * kyi.conjugate() + eps_1.conjugate() * kyt))  # Eq 61

    mu2 = - kyt.real / np.sqrt(eps_2).real  # by definition of kyt

    return rv, rh, mu2


@jit(nopython=True, cache=True)
def snell_angle(eps_1, eps_2, mu1):
    """compute mu2 the cos(angle) in the second medium according to Snell's law."""

    # incident wavenumber
    n1 = np.sqrt(eps_1)
    kiz2 = n1.real**2 * (1 - mu1**2)   # this the square of kiz = n1 * sin(theta)

    ktz2 = kiz2   # unnumbered equation before 22  -> tangential k is conserved throught the interface (=Snell law)
    kyt = - np.sqrt(complex(eps_2) - ktz2)                  # Eq 8 for t

    mu2 = - kyt.real / np.sqrt(eps_2).real  # by definition of kyt

    return mu2


# @numba.vectorize([numba.float64(numba.complex128), numba.float32(numba.complex64)], cache=True)
# #@jit(nopython=True, cache=True)
# def abs2(x):
#     return x.real**2 + x.imag**2


@jit(nopython=True, cache=True)
def matmul2(a, b):

    # a, trans_v_a = a
    # b, trans_v_b = b

    assert a.shape == b.shape

    c = np.empty((2, 2, *a.shape[2:]), dtype=a.dtype)

    # c[0, 0] = a[0, 0] * b[0, 0] + a[0, 1] * b[1, 0]
    # c[0, 1] = a[0, 0] * b[0, 1] + a[0, 1] * b[1, 1]
    # c[1, 0] = a[1, 0] * b[0, 0] + a[1, 1] * b[1, 0]
    # c[1, 1] = a[1, 0] * b[0, 1] + a[1, 1] * b[1, 1]

    # row-
    # c[:, 0] = a[:, 0] * b[0, 0] + a[:, 1] * b[1, 0]
    # c[:, 1] = a[:, 0] * b[0, 1] + a[:, 1] * b[1, 1]

    c[0, :] = a[0, 0] * b[0, :] + a[0, 1] * b[1, :]
    c[1, :] = a[1, 0] * b[0, :] + a[1, 1] * b[1, :]
    # return c, trans_v_a * trans_v_b
    return c
