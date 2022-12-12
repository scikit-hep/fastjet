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


def test_exclusive_constituents_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ]
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)

    constituent_output = [
        [
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        ],
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
    ]
    assert constituent_output == cluster.exclusive_jets_constituents(2).to_list()
    constituent_index_output = [[1, 2], [0]]
    assert (
        constituent_index_output
        == cluster.exclusive_jets_constituent_index(2).to_list()
    )


def test_exclusive_lund_declustering_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
            {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )

    jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, 0.8)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)

    lds = cluster.exclusive_jets_lund_declusterings(2)

    lund_declustering_output = [
        [
            {"Delta": 0.08755181299980186, "kt": 0.30179987478357223},
            {"Delta": 0.019481226884377707, "kt": 0.06602095529127928},
        ],
        [{"Delta": 0.014750342295225208, "kt": 1.0480537658466145}],
    ]

    assert lund_declustering_output == lds.to_list()


def test_exclusive_lund_declustering_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
                {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
                {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
        with_name="Momentum4D",
    )

    jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, 0.8)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)

    lds = cluster.exclusive_jets_lund_declusterings(2)

    lund_declustering_output = [
        [
            [
                {"Delta": 0.08755181299980186, "kt": 0.30179987478357223},
                {"Delta": 0.019481226884377707, "kt": 0.06602095529127928},
            ],
            [{"Delta": 0.014750342295225208, "kt": 1.0480537658466145}],
        ],
        [
            [
                {"Delta": 0.08755181299980186, "kt": 0.30179987478357223},
                {"Delta": 0.019481226884377707, "kt": 0.06602095529127928},
            ],
            [{"Delta": 0.014750342295225208, "kt": 1.0480537658466145}],
        ],
    ]

    assert lund_declustering_output == lds.to_list()


def test_exclusive_constituents_multi():
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

    constituents_output = [
        [
            [
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78}],
        ],
        [
            [
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78}],
        ],
    ]
    assert cluster.exclusive_jets_constituents(2).to_list() == constituents_output
    constituent_index_out = [[[1, 2], [0]], [[1, 2], [0]]]
    assert (
        constituent_index_out == cluster.exclusive_jets_constituent_index(2).to_list()
    )


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


def test_listoffset_indexed_input():
    inputs = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
            ],
        ],
        with_name="Momentum4D",
    )
    out = ak.Array(
        ak.contents.RecordArray(
            [
                ak.contents.NumpyArray(np.asarray(ak.Array(inputs.layout.content).px)),
                ak.contents.NumpyArray(np.asarray(ak.Array(inputs.layout.content).py)),
                ak.contents.NumpyArray(np.asarray(ak.Array(inputs.layout.content).pz)),
                ak.contents.NumpyArray(np.asarray(ak.Array(inputs.layout.content).E)),
            ],
            ["px", "py", "pz", "E"],
        )
    )
    out = ak.Array(
        ak.contents.IndexedArray(
            ak.index.Index64([7, 2, 3, 1, 0, 5, 4, 6, 8]), out.layout
        )
    )
    out = ak.Array(ak.contents.ListOffsetArray(inputs.layout.offsets, out.layout))
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(out, jetdef)
    inclusive_jets = [
        [
            {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
        [
            {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
