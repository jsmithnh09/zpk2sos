import pytest
from ctypes import CDLL, c_size_t, POINTER, c_double
from pathlib import Path
import numpy as np
import sys
import os
import atexit
from scipy.signal import butter, cheby1, cheby2, ellip, zpk2sos, sosfreqz
from typing import Tuple

# CMake ought to configure this.
_libenv = os.environ.get("ZSOS_LIB_PATH")
if not _libenv:
    raise RuntimeError("ZSOS_LIB_PATH environment variable not set by Cmake.")
else:
    _libpath = Path(_libenv)
    if not _libpath.is_file():
        raise RuntimeError("Specified library path does not exist.")
_libstr = str(_libpath.resolve())


# putting in the unload function for the DLL.
def unload_library() -> int:
    if sys.platform == "win32":
        import kernel32

        close = kernel32.FreeLibrary
    else:
        import _ctypes

        close = _ctypes.dlclose
    close(lib._handle)
    return 0


# Load the library using the path from Cmake/Ctest.
lib = CDLL(_libstr)
atexit.register(unload_library)

# exposing the zpk2sos signature for interleaved complex values.
lib.soscount.argtypes = [c_size_t, c_size_t]
lib.soscount.restype = c_size_t
lib.zpk2sos.argtypes = [
    POINTER(c_double),
    c_size_t,
    POINTER(c_double),
    c_size_t,
    c_double,
    POINTER(c_double),
]
lib.zpk2sos.restype = c_size_t


def soscmp(sos1, sos2, num=1024, atol=1e-08) -> bool:
    """
    Compare SOS responses.

    Parameters
    ----------
    sos1 : ndarray
           the reference SOS matrix.
    sos2 : ndarray
           the test SOS matrix.
    num  : Union[float,int]
           an indicator of the number of bins to sample.
    fs   : float
           the sampling rate (since this is discrete.)

    Returns
    -------
    flag : bool
           indicates True if the respones are within
           default precision.
    """
    (wref, href) = sosfreqz(sos1, worN=num)
    (wtest, htest) = sosfreqz(sos2, worN=num)
    return bool(np.all(np.isclose(href, htest, atol=atol)))


def zpk2sos_wrapper(z, p, k) -> Tuple[np.ndarray, int]:
    """
    ZPK2SOS C-wrapper.

    Parameters
    ----------
    z : ndarray
        complex array of zeroes on z-plane.
    p : ndarray
        complex array of poles on the z-plane.
    k : double
        real scalar gain of the system.

    Returns
    -------
    sos, err : ndarray, int
        Second-order sections of the IIR filter, along with
        an indicator if an internal error code (must be zero.)
    """
    # need to enforce complex values.
    z = np.complex128(z)
    p = np.complex128(p)

    z_ptr = z.ctypes.data_as(POINTER(c_double))
    p_ptr = p.ctypes.data_as(POINTER(c_double))
    kc = c_double(k)

    numz = c_size_t(len(z))
    nump = c_size_t(len(p))
    nstages = lib.soscount(numz, nump)
    sosmat = np.zeros((nstages, 6), dtype=np.float64)

    sos_ptr = sosmat.ctypes.data_as(POINTER(c_double))
    err = lib.zpk2sos(z_ptr, numz, p_ptr, nump, kc, sos_ptr)
    return (sosmat, err)


# defining the orders and corner frequencies.
_fc = np.logspace(np.log10(20), np.log10(18e3), num=10)
_ord = range(2, 13, 1)


# packing multiple of the same parametrizations into a single decorator.
def common_setup(func):
    func = pytest.mark.parametrize("fc", _fc)(func)
    func = pytest.mark.parametrize("btype", ["lowpass", "highpass"])(func)
    func = pytest.mark.parametrize("order", _ord)(func)
    return func


@common_setup
def test_butterworth(fc, btype, order):
    (z, p, k) = butter(order, fc, fs=48e3, btype=btype, analog=False, output="zpk")
    sosref = zpk2sos(z, p, k, pairing="nearest")
    (sostest, err) = zpk2sos_wrapper(z, p, k)
    assert soscmp(sosref, sostest, num=1024)


@common_setup
def test_ellip(fc, btype, order):
    rp, rs = 40, 100
    (z, p, k) = ellip(
        order, rp, rs, fc, fs=48e3, btype=btype, analog=False, output="zpk"
    )
    sosref = zpk2sos(z, p, k, pairing="nearest")
    (sostest, err) = zpk2sos_wrapper(z, p, k)
    assert soscmp(sosref, sostest, num=1024, atol=1e-9)


@common_setup
def test_cheby1(fc, btype, order):
    rp = 40
    (z, p, k) = cheby1(order, rp, fc, fs=48e3, btype=btype, analog=False, output="zpk")
    sosref = zpk2sos(z, p, k, pairing="nearest")
    (sostest, err) = zpk2sos_wrapper(z, p, k)
    assert soscmp(sosref, sostest, num=1024, atol=1e-9)


@common_setup
def test_cheby2(fc, btype, order):
    rs = 100
    (z, p, k) = cheby2(order, rs, fc, fs=48e3, btype=btype, analog=False, output="zpk")
    sosref = zpk2sos(z, p, k, pairing="nearest")
    (sostest, err) = zpk2sos_wrapper(z, p, k)
    assert soscmp(sosref, sostest, num=1024, atol=1e-9)
