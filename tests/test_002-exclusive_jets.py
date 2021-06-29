import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401

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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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


def test_exclusive_ycut():
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    multi_exclusive_ycut = [
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        ],
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        ],
    ]
    assert (
        multi_exclusive_ycut == cluster.exclusive_jets_ycut(ycut=0.0000000001).to_list()
    )


def test_exclusive_ycut_multi():
    array2 = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef2 = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array2, jetdef2)
    single_exclusive_ycut = []
    assert single_exclusive_ycut == cluster.exclusive_jets_ycut(ycut=0.01).to_list()
