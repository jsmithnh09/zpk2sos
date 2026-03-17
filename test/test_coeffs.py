import pytest
from ctypes import CDLL, c_size_t, POINTER, c_double
from pathlib import Path
import numpy as np
import sys
import atexit
from scipy.signal import butter, cheby1, cheby2, ellip, zpk2sos

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
def unload_library():
    if sys.platform == "win32":
        import kernel32

        close = kernel32.FreeLibrary
    else:
        import _ctypes

        close = _ctypes.dlclose
    close(lib._handle)


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


def zpk2sos_wrapper(z, p, k):
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

