import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401


def test_listoffset_inclusive():
    array = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    inclusive_jets_out = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]

    assert inclusive_jets_out == cluster.inclusive_jets().to_list()
    constituent_index_out = [[[[[0], [1, 2]]]], [[[[0], [1, 2]]]]]
    assert constituent_index_out == cluster.constituent_index().to_list()


def test_bytemasked_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = inputs.mask[[True, False, True]]
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        None,
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], None, [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_regulararray_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(ak.contents.RegularArray(inputs.layout, 1))
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[[0], [1, 2]]]]], [[[[[0], [1, 2]]]]], [[[[[0], [1, 2]]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_listarray_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.ListArray(
            ak.index.Index64([0, 1]), ak.index.Index64([1, 2]), inputs.layout
        )
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[[0], [1, 2]]]]], [[[[[0], [1, 2]]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_listarray32_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.ListArray(
            ak.index.Index32([0, 1]), ak.index.Index32([1, 2]), inputs.layout
        )
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[[0], [1, 2]]]]], [[[[[0], [1, 2]]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_listarrayU32_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.ListArray(
            ak.index.IndexU32([0, 1]), ak.index.IndexU32([1, 2]), inputs.layout
        )
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[[0], [1, 2]]]]], [[[[[0], [1, 2]]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexedarray_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.IndexedArray(ak.index.Index64([0, 2, 1]), inputs.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], [[[[0], [1, 2]]]], [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexedarray32_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.IndexedArray(ak.index.Index32([0, 2, 1]), inputs.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], [[[[0], [1, 2]]]], [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexedarrayU32_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.IndexedArray(ak.index.IndexU32([0, 2, 1]), inputs.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], [[[[0], [1, 2]]]], [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexedoptionarray32_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.IndexedOptionArray(ak.index.Index32([0, -2, 1]), inputs.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        None,
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], None, [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexedoptionarray64_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(
        ak.contents.IndexedOptionArray(ak.index.Index64([0, -2, 1]), inputs.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        None,
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], None, [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_unmaskedarray_input():
    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    inputs = ak.Array(ak.contents.UnmaskedArray(inputs.layout))
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], [[[[0], [1, 2]]]], [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_virtualarray_input_NOT_REALLY():
    def generate():
        inclusive_jets = ak.Array(
            [
                [
                    [
                        [
                            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                            {
                                "px": 32.2,
                                "py": 64.21,
                                "pz": 543.34,
                                "E": 24.12,
                                "ex": 0.35,
                            },
                            {
                                "px": 32.45,
                                "py": 63.21,
                                "pz": 543.14,
                                "E": 24.56,
                                "ex": 0.0,
                            },
                        ]
                    ]
                ],
                [
                    [
                        [
                            {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                            {
                                "px": 32.2,
                                "py": 64.21,
                                "pz": 543.34,
                                "E": 24.12,
                                "ex": 0.35,
                            },
                            {
                                "px": 32.45,
                                "py": 63.21,
                                "pz": 543.14,
                                "E": 24.56,
                                "ex": 0.0,
                            },
                        ]
                    ]
                ],
                [
                    [
                        [
                            {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                            {
                                "px": 32.2,
                                "py": 64.21,
                                "pz": 543.34,
                                "E": 24.12,
                                "ex": 0.35,
                            },
                            {
                                "px": 32.45,
                                "py": 63.21,
                                "pz": 543.14,
                                "E": 24.56,
                                "ex": 0.1,
                            },
                        ]
                    ]
                ],
            ],
            with_name="Momentum4D",
        )
        inclusive_jets = inclusive_jets.mask[[True, False, True]]
        return inclusive_jets

    inputs = generate()
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(inputs, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        None,
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[0], [1, 2]]]], None, [[[[0], [1, 2]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_indexed_subtree_input():
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
    out = ak.Array(
        ak.contents.ListOffsetArray(ak.index.Index64([0, 1, 2, 3]), out.layout)
    )
    out = ak.Array(
        ak.contents.ListOffsetArray(ak.index.Index64([0, 1, 2, 3]), out.layout)
    )
    out = ak.Array(
        ak.contents.ListOffsetArray(ak.index.Index64([0, 1, 2, 3]), out.layout)
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(out, jetdef)
    inclusive_jets = [
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
        [
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [[[[[[2], [0, 1]]]]], [[[[[1], [0, 2]]]]], [[[[[1], [0, 2]]]]]]
    assert const_idx == cluster.constituent_index().to_list()


def test_union_8_64_input():

    inputs = ak.Array(
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    input2 = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    out = ak.concatenate((input2, inputs), axis=0)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(out, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
        [
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
        [
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [
        [[[[0], [1, 2]]]],
        [[[[0], [1, 2]]]],
        [[[[0], [1, 2]]]],
        [[[0], [1, 2]]],
        [[[0], [1, 2]]],
        [[[0], [1, 2]]],
    ]
    assert const_idx == cluster.constituent_index().to_list()


def test_partitioned_input_NOT_REALLY():

    inputs = ak.Array(
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    input2 = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    out = ak.concatenate((input2, inputs), axis=0, highlevel=True)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(out, jetdef)
    inclusive_jets = [
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ],
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
        [
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
        [
            [
                {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
            ]
        ],
    ]
    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [
        [[[[0], [1, 2]]]],
        [[[[0], [1, 2]]]],
        [[[[0], [1, 2]]]],
        [[[0], [1, 2]]],
        [[[0], [1, 2]]],
        [[[0], [1, 2]]],
    ]
    assert const_idx == cluster.constituent_index().to_list()


def test_record_input():

    inputs = ak.Array(
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                ]
            ],
            [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    input2 = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    out = ak.zip({"a": inputs, "b": input2}, depth_limit=1)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(out, jetdef)
    inclusive_jets = [
        {
            "a": [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ],
            "b": [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ],
        },
        {
            "a": [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ],
            "b": [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ],
        },
        {
            "a": [
                [
                    {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ],
            "b": [
                [
                    [
                        {"px": 11.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                        {
                            "px": 64.65,
                            "py": 127.41999999999999,
                            "pz": 1086.48,
                            "E": 48.68,
                        },
                    ]
                ]
            ],
        },
    ]

    assert inclusive_jets == cluster.inclusive_jets().to_list()
    const_idx = [
        {"a": [[[0], [1, 2]]], "b": [[[[0], [1, 2]]]]},
        {"a": [[[0], [1, 2]]], "b": [[[[0], [1, 2]]]]},
        {"a": [[[0], [1, 2]]], "b": [[[[0], [1, 2]]]]},
    ]

    assert const_idx == cluster.constituent_index().to_list()


def test_func_call__input():

    inputs = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    input2 = ak.Array(
        [
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
                    ]
                ]
            ],
            [
                [
                    [
                        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.1},
                    ]
                ]
            ],
        ],
        with_name="Momentum4D",
    )
    array1 = ak.Array(
        [
            [
                [
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                ]
            ]
        ]
    )
    array2 = ak.Array(
        [
            [
                [
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                    {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                    {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
                ]
            ]
        ]
    )
    inn = ak.zip({"a": array1, "b": array2}, depth_limit=1)
    out = ak.zip({"a": inputs, "b": input2}, depth_limit=1)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet.ClusterSequence(out, jetdef)
    get_parents = [
        {
            "a": [
                [
                    [],
                    [
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
                    ],
                    [],
                ]
            ],
            "b": [
                [
                    [
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
                    ],
                    [],
                    [
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
                    ],
                ]
            ],
        }
    ]
    assert get_parents == cluster.get_parents(inn).to_list()
