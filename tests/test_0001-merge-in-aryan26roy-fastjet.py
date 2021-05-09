# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE
import numpy as np

import fastjet
import fastjet._ext
import fastjet._swig


def test_swig_bindings():
    fastjet._swig.PseudoJet(1, 2, 3, 4)


def test_pybind11():
    px = np.array([1.12, 1.15, 11.1], dtype="float32")
    py = np.array([2.11, 2.41, 24.1], dtype="float32")
    pz = np.array([3.71, 3.18, 35.11], dtype="float32")
    E = np.array([5.1, 52.1, 51.11], dtype="float32")
    jetdef = fastjet._swig.JetDefinition(fastjet._swig.antikt_algorithm, 0.6)
    inps = {"algor": "anti-kt"}
    inpf = {"R": 0.6}
    fastjet._ext.interface(px, py, pz, E, inps, inpf, jetdef)
