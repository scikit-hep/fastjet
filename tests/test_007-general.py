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
    inputs = ak.Array(ak.layout.RegularArray(inputs.layout, 1))
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
        ak.layout.ListArray64(
            ak.layout.Index64([0, 1]), ak.layout.Index64([1, 2]), inputs.layout
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
