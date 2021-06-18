import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet.pyjet  # noqa: F401

vector = pytest.importorskip("vector")  # noqa: F401


def test_unique_history_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    out = [[0, 5, 1, 2, 3, 4], [0, 3, 1, 2], [0, 5, 1, 2, 3, 4]]
    assert out == cluster.unique_history_order().to_list()


def test_unique_history_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    out = [0, 5, 1, 2, 3, 4]
    assert out == cluster.unique_history_order().to_list()
