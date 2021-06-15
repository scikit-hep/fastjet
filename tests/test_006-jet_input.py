import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet
import fastjet.pyjet  # noqa: F401


def test_exclusive_subjets_multi():

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
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    dcut_out = [
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}],
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
    ]
    assert dcut_out == cluster.exclusive_subjets(test_input, dcut=0.0000003).to_list()
    nsub_out = [
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}],
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
    ]
    assert nsub_out == cluster.exclusive_subjets(test_input, nsub=1).to_list()


def test_exclusive_subjets_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    dcut_out = [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}]
    assert dcut_out == cluster.exclusive_subjets(test_input, dcut=0.0000003).to_list()
    nsub_out = [
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
    ]
    assert nsub_out == cluster.exclusive_subjets(test_input, nsub=2).to_list()


def test_exclusive_subjets_up_to_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = [
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
    ]
    assert nsub_out == cluster.exclusive_subjets_up_to(test_input, nsub=2).to_list()


def test_exclusive_subjets_up_to_multi():

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
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}],
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
    ]
    assert nsub_out == cluster.exclusive_subjets_up_to(test_input, nsub=1).to_list()


def test_exclusive_subdmerge_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}]
    )
    nsub_out = 1.1713099970894107e-07
    assert nsub_out == cluster.exclusive_subdmerge(test_input, nsub=1)


def test_exclusive_subdmerge_multi():

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
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [1.1713099970894107e-07, 0.0]
    assert nsub_out == cluster.exclusive_subdmerge(test_input, nsub=1).to_list()


def test_exclusive_subdmerge_max_multi():

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
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [1.1713099970894107e-07, 0.0]
    assert nsub_out == cluster.exclusive_subdmerge_max(test_input, nsub=1).to_list()


def test_exclusive_subdmerge_max_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}]
    )
    nsub_out = 1.1713099970894107e-07
    assert nsub_out == cluster.exclusive_subdmerge_max(test_input, nsub=1)


def test_n_exclusive_subjets_multi():

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
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = [1, 1]
    assert nsub_out == cluster.n_exclusive_subjets(test_input, dcut=0.0001).to_list()


def test_n_exclusive_subjets_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = 1
    assert nsub_out == cluster.n_exclusive_subjets(test_input, dcut=0.0001)
