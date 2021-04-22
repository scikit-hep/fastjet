# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import fastjet
import fastjet._ext
import fastjet._swig
import numpy as np

def test_swig_bindings():
    fastjet._swig.PseudoJet(1, 2, 3, 4)


def test_pybind11():
    data = np.array([(6, 7, 8, 9), (6.1, 7.2, 8.3, 9.4), (4.1, 2.2, 8.1, 9.4), (6.1, 71.2, 18.3, 9.41), (6.1, 7.211, 128.3, 91.4), (6.41, 27.2, 8.3,19.4)],dtype=[("px", np.float64), ("py", np.float32), ("pz", np.int32), ("E", np.int8)])
    data = (
            data.astype(
                [
                    ("px", np.float32),
                    ("py", np.float32),
                    ("pz", np.float32),
                    ("E", np.float32),
                ]
            )
            .view(np.float32)
            .reshape(-1, 4)
        )
    data = data.tolist()
    R = 0.8
    algor = "antikt_algorithm"
    fastjet._ext.interface(data, R, algor)
