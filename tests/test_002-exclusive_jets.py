import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet.pyjet  # noqa: F401

vector = pytest.importorskip("vector")  # noqa: F401


def test_exclusive_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ]
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    single_exclusive_dcut = [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}]
    assert single_exclusive_dcut == cluster.exclusive_jets(dcut=0.0001).to_list()
    single_exclusive_njets = [
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
    ]
    assert single_exclusive_njets == cluster.exclusive_jets(n_jets=3).to_list()


def test_exclusive_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ]
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    multi_exclusive_njets = [
        [
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        ],
        [
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        ],
    ]
    assert multi_exclusive_njets == cluster.exclusive_jets(n_jets=3).to_list()
    multi_exclusive_dcut = [
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
    ]
    assert multi_exclusive_dcut == cluster.exclusive_jets(dcut=0.0001).to_list()