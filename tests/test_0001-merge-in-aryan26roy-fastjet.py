# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import fastjet
import fastjet._ext
import fastjet._swig


def test_swig_bindings():
    fastjet._swig.PseudoJet(1, 2, 3, 4)


def test_pybind11():
    fastjet._ext.interface()
