# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import fastjet
import fastjet._ext


def test_pybind11():
    assert fastjet._ext.add(1, 2) == 3
    assert fastjet._ext.subtract(1, 2) == -1


def test_version():
    assert fastjet.__version__
