import numpy as np
from ctypes import c_float, c_double
from numpy.polynomial.hermite import Hermite
from typing import Union
import math

def numpy_eval_hermite(n : int, x : Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    """This function duplicates what scipy.special.eval_hermite does:
    
    Evaluate physicist's Hermite polynomial at a point.
    
    Have found that installing `scipy` on some clusters can be problematic
    (looking at you Pawsey) so implement a numpy version to sidestep needed
    scipy.
    
    See https://docs.scipy.org/doc/scipy/reference/generated/scipy.special.eval_hermite.html#scipy-special-eval-hermite for more detail

    Parameters
    ----------
    n : int
        Order of hermite polynomial
    x : Union[float, np.ndarray]
        Coordinate(s) to evaluate the polynomial at

    Returns
    -------
    Union[float, np.ndarray]
        The evaulated polynomial values
    """
    
    ##`Hermite` below will sum all orders of polynomials, with given
    ##coefficients. We just want 
    herm_coeffs = np.zeros(n + 1)
    herm_coeffs[n] = 1
    
    herm_func = Hermite(herm_coeffs)
    
    return herm_func(x)


def calc_basis_func_1D_numpy(n : int, x : Union[float, np.ndarray], beta=1) -> np.ndarray:
    """Calculate shapelet basis function values to use in WODEN. These
    values work with models fitted using SHAMFI https://github.com/JLBLine/SHAMFI
    
    See Line et al. 2020 https://doi.org/10.1017/pasa.2020.18 for more
    information on shapelet fitting

    Parameters
    ----------
    n : int
        Order of shapelet basis function to calculate
    x : Union[float, np.ndarray]
        Coordinate(s) to evaluate the basis function at
    beta : int, optional
        The x-coord scaling factor, should be 1 when used with WODEN, by default 1

    Returns
    -------
    np.ndarray
        The basis function values
    """

    ##this is the fourier transform of image-based basis functions written into
    ##Line et al. 2020, so there is a factor of pi difference. Something
    ##something fourier transform convention something something
    norm = np.sqrt(beta)*np.sqrt((2**n*float(math.factorial(n))))
    hermite = numpy_eval_hermite(n, x)
    gauss = np.exp(-0.5*((x*beta)**2))
    
    return (hermite*gauss) / norm

def create_sbf(precision = "double", sbf_N = 101, sbf_c = 5000, sbf_dx = 0.01):
    """Creates a ctypes array of float or doubles, containing shapelet basis
    functions to feed into `run_woden`.
    
    All settings default to settings hard-coded into `libwoden_double.so`,
    so only mess with them if you know what you're doing. A future iteration
    of WODEN could have `sbf_N`, `sbf_c`, `sbf_dx` put into `woden_settings`
    to allow for on-the-fly basis function adjustment.

    With the defaults listed below, each basis function will be calculate over
    a range of x = -50 to x = 50, with a resolution of x = 0.01 (for a total
    of 10001 x values per basis function).
    
    The final array is 1D, populated by basis function n = 0 for all x values,
    up to basis function n = (sbf_N - 1).

    Parameters
    ----------
    precision : str, optional
        Precision to return. If set to 'float', uses `c_float`, if set to
        `double` uses `c_double`. By default set to "double".
    sbf_N : int, optional
        Sets what order of basis function to generate up to, by default 101
    sbf_c : int, optional
        Sets the central array value, where x = 0 will be, by default 5000
    sbf_dx : float, optional
        The `x` resolution of each array element, by default 0.01

    Returns
    -------
    c_float_Array or c_double_Array
        All of the basis functions stored in a 1D array of len = 2*sbf_c + 1
    """
    
    x_range = np.arange(-sbf_c*sbf_dx, sbf_c*sbf_dx + sbf_dx, sbf_dx)
    sbf_L = len(x_range)
    
    if precision == 'float':
        sbf = (c_float*(sbf_L*sbf_N))()
    else:
        sbf = (c_double*(sbf_L*sbf_N))()
    
    for n in range(sbf_N):
        low_ind = n*sbf_L
        basis = calc_basis_func_1D_numpy(n, x_range)
        
        for basis_ind, sbf_ind in enumerate(range(low_ind, low_ind + sbf_L)):
            sbf[sbf_ind] = basis[basis_ind]
    
    return sbf