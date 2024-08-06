import numpy as np


# radial unit conversion functions

def twotheta_to_q(twotheta, wavelength):
    """
    Converts two-theta value to units of q, based on the wavelength

    Args:
        twotheta (array_like): float or float array of two-theta values
        wavelength (float): X-ray beam wavelength value

    Returns:
        (array_like): twotheta converted to units of q
    """
    twotheta = np.asarray(twotheta)
    wavelength = float(wavelength)
    pre_factor = ((4 * np.pi) / wavelength)
    return pre_factor * np.sin(twotheta / 2)

def q_to_twotheta(q, wavelength):
    """
    Converts q value to units of two-theta, based on the wavelength

    Args:
        q (array_like):  float or float array of q values
        wavelength (float): X-ray beam wavelength value

    Returns:
        (array_like): q converted to units of two-theta
    """
    q = np.asarray(q)
    wavelength = float(wavelength)
    pre_factor = wavelength / (4 * np.pi)
    return 2 * np.arcsin(q * pre_factor)

def q_to_d(q):
    """
    Converts q to units of d

    Args:
        q (array_like): float or float array of q values

    Returns:
        array_like: q converted to units of d
    """
    return (2 * np.pi) / np.asarray(q)

def d_to_q(d):
    """
    Converts d to units of q

    Args:
        d (array_like): float or float array of d values

    Returns:
        array_like: d converted to units of q
    """
    return (2 * np.pi) / np.asarray(d)

def twotheta_to_d(twotheta, wavelength):
    """
    Converts two-theta value to units of d, based on the wavelength

    Args:
        twotheta (array_like): float or float array of two-theta values
        wavelength (float): X-ray beam wavelength value

    Returns:
        (array_like): twotheta converted to units of d
    """
    th = np.asarray(twotheta)/2
    rad = np.radians(th)
    t = 2*np.sin(rad)
    d = (wavelength)/t
    return d

def tth_wl1_to_wl2(tth1,wl1=0.187,wl2=0.4592):
    """
    Converts two-theta value with one wavelength, to two-theta value 
    with another wavelength.

    Args:
        tth1 (array_like): float or float array of two-theta value
        wl1 (float, optional): initial wavelength. Defaults to 0.187.
        wl2 (float, optional): new wavelength. Defaults to 0.4592.

    Returns:
        (array_like): two-theta converted to correspond to new wavelength
    """
    q = twotheta_to_q(np.deg2rad(tth1), wl1)
    return np.rad2deg(q_to_twotheta(q,wl2))





























