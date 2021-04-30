# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE
import numpy as np

import fastjet
import fastjet._ext
import fastjet._swig


def test_swig_bindings():
    fastjet._swig.PseudoJet(1, 2, 3, 4)


def test_pybind11():
    px = np.ndarray([1.12, 1.15, 11.1])
    py = np.ndarray([2.11, 2.41, 24.1])
    pz = np.ndarray([3.71, 3.18, 35.11])
    E = np.ndarray([5.1, 52.1, 51.11])
    fastjet._ext.interface(px, py, pz, E, 0.8, "antikt_algorithm")
