# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import ctypes
import pathlib

_fastjet_core = pathlib.Path(__file__).parent.resolve() / "_fastjet_core" / "lib"
_libfastjet = ctypes.cdll.LoadLibrary(_fastjet_core / "libfastjet.so")
_libfastjettools = ctypes.cdll.LoadLibrary(_fastjet_core / "libfastjettools.so")
_libsiscone = ctypes.cdll.LoadLibrary(_fastjet_core / "libsiscone.so")
_libsiscone_spherical = ctypes.cdll.LoadLibrary(
    _fastjet_core / "libsiscone_spherical.so"
)
_libfastjetplugins = ctypes.cdll.LoadLibrary(_fastjet_core / "libfastjetplugins.so")


import fastjet._ext  # noqa: F401, E402
import fastjet._swig  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)
