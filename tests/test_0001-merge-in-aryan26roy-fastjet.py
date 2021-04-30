# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import fastjet
import fastjet._ext
import fastjet._swig
import numpy as np
import awkward as ak


def test_swig_bindings():
    fastjet._swig.PseudoJet(1, 2, 3, 4)


def test_pybind11():

    array = ak.Array(
        [
            {"px": 1.12, "py": 2.11, "pz": 3.71, "E": 5.1},
            {"px": 1.15, "py": 2.41, "pz": 3.18, "E": 52.1},
            {"px": 11.1, "py": 24.1, "pz": 35.11, "E": 51.11},
        ]
    )
    aa = fastjet.pyjet.AwkwardClusterSequence(array, 0.8, "antikt_algorithm")
