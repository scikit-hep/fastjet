import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet
import fastjet._pyjet  # noqa: F401


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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}]
    )
    nsub_out = 1.1713099970894107e-07
    assert nsub_out == pytest.approx(cluster.exclusive_subdmerge(test_input, nsub=1))


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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [1.1713099970894107e-07, 0.0]
    assert nsub_out == pytest.approx(
        cluster.exclusive_subdmerge(test_input, nsub=1).to_list()
    )


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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [1.1713099970894107e-07, 0.0]
    assert nsub_out == pytest.approx(
        cluster.exclusive_subdmerge_max(test_input, nsub=1).to_list()
    )


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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [{"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68}]
    )
    nsub_out = 1.1713099970894107e-07
    assert nsub_out == pytest.approx(
        cluster.exclusive_subdmerge_max(test_input, nsub=1)
    )


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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = 1
    assert nsub_out == cluster.n_exclusive_subjets(test_input, dcut=0.0001)


def test_has_parents_multi():

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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = [False, True]
    assert nsub_out == cluster.has_parents(test_input).to_list()


def test_has_parents_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = True
    assert nsub_out == cluster.has_parents(test_input)


def test_has_child_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = False
    assert nsub_out == cluster.has_child(test_input)


def test_has_child_multi():

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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    nsub_out = [False, False]
    assert nsub_out == cluster.has_child(test_input).to_list()


def test_jet_scale_for_algorithm_multi():

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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = [4.898244145524307e-05, 0.08561643835616438]
    assert nsub_out == pytest.approx(
        cluster.jet_scale_for_algorithm(test_input).to_list()
    )


def test_jet_scale_for_algorithm_single():

    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ]
    )
    nsub_out = 0.08561643835616438
    assert nsub_out == pytest.approx(cluster.jet_scale_for_algorithm(test_input))


def test_get_parents_multi():
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    out = [
        [],
        [
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        ],
    ]
    assert out == cluster.get_parents(test_input).to_list()


def test_get_parents_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    out = [
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
    ]
    assert out == cluster.get_parents(test_input).to_list()


def test_get_child_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    out = []
    assert out == cluster.get_child(test_input).to_list()


def test_get_child_multi():
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
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    test_input = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ]
    )
    out = [
        [],
        [],
    ]
    assert out == cluster.get_child(test_input).to_list()
