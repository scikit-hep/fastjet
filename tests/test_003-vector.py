import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401

vector = pytest.importorskip("vector")  # noqa: F401


def test_vector_single():
    array = ak.Array(
        [
            {"phi": 1.2, "eta": 3.2, "rho": 5.4, "E": 2.5, "ex": 0.78},
            {"phi": 32.2, "eta": 64.21, "rho": 543.34, "E": 24.12, "ex": 0.35},
            {"phi": 32.45, "eta": 63.21, "rho": 543.14, "E": 24.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
        behavior=vector.backends.awkward.behavior,
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    exclusive_dcut = [
        {
            "px": 1.9567318741740376,
            "py": 5.033011064223023,
            "pz": 66.12777358145367,
            "E": 2.5,
        }
    ]
    assert exclusive_dcut == [
        pytest.approx(x) for x in cluster.exclusive_jets(dcut=0.0001).to_list()
    ]
    exclusive_njets = [
        {
            "px": 1.9567318741740376,
            "py": 5.033011064223023,
            "pz": 66.12777358145367,
            "E": 2.5,
        }
    ]
    assert exclusive_njets == [
        pytest.approx(x) for x in cluster.exclusive_jets(n_jets=1).to_list()
    ]
    inclusive_jets = [
        {
            "px": 1.9567318741740376,
            "py": 5.033011064223023,
            "pz": 66.12777358145367,
            "E": 2.5,
        },
        {
            "px": 277.71965200410807,
            "py": 466.76852345752394,
            "pz": 7.684860589309717e29,
            "E": 24.56,
        },
        {
            "px": 384.7080099115854,
            "py": 383.69011286436285,
            "pz": 2.0897309060783084e30,
            "E": 24.12,
        },
    ]
    assert inclusive_jets == [
        pytest.approx(x) for x in cluster.inclusive_jets().to_list()
    ]
    constituents = [
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            }
        ],
        [
            {
                "px": 277.71965200410807,
                "py": 466.76852345752394,
                "pz": 7.684860589309717e29,
                "E": 24.56,
            }
        ],
        [
            {
                "px": 384.7080099115854,
                "py": 383.69011286436285,
                "pz": 2.0897309060783084e30,
                "E": 24.12,
            }
        ],
    ]
    assert constituents == [
        [pytest.approx(y) for y in x] for x in cluster.constituents().to_list()
    ]
    constituent_index = [[0], [2], [1]]
    assert constituent_index == [
        pytest.approx(x) for x in cluster.constituent_index().to_list()
    ]


def test_vector_multi():
    array = ak.Array(
        [
            [
                {"phi": 1.2, "eta": 3.2, "rho": 5.4, "E": 2.5, "ex": 0.78},
                {"phi": 32.2, "eta": 64.21, "rho": 543.34, "E": 24.12, "ex": 0.35},
                {"phi": 32.45, "eta": 63.21, "rho": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"phi": 1.2, "eta": 3.2, "rho": 5.4, "E": 2.5, "ex": 0.78},
                {"phi": 32.2, "eta": 64.21, "rho": 543.34, "E": 24.12, "ex": 0.35},
                {"phi": 32.45, "eta": 63.21, "rho": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
        with_name="Momentum4D",
        behavior=vector.backends.awkward.behavior,
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    exclusive_dcut = [
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            }
        ],
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            }
        ],
    ]

    assert exclusive_dcut == [
        [pytest.approx(y) for y in x]
        for x in cluster.exclusive_jets(dcut=0.0001).to_list()
    ]
    exclusive_njets = [
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            }
        ],
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            }
        ],
    ]

    assert exclusive_njets == [
        [pytest.approx(y) for y in x]
        for x in cluster.exclusive_jets(n_jets=1).to_list()
    ]
    inclusive_jets = [
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            },
            {
                "px": 277.71965200410807,
                "py": 466.76852345752394,
                "pz": 7.684860589309717e29,
                "E": 24.56,
            },
            {
                "px": 384.7080099115854,
                "py": 383.69011286436285,
                "pz": 2.0897309060783084e30,
                "E": 24.12,
            },
        ],
        [
            {
                "px": 1.9567318741740376,
                "py": 5.033011064223023,
                "pz": 66.12777358145367,
                "E": 2.5,
            },
            {
                "px": 277.71965200410807,
                "py": 466.76852345752394,
                "pz": 7.684860589309717e29,
                "E": 24.56,
            },
            {
                "px": 384.7080099115854,
                "py": 383.69011286436285,
                "pz": 2.0897309060783084e30,
                "E": 24.12,
            },
        ],
    ]
    assert inclusive_jets == [
        [pytest.approx(y) for y in x] for x in cluster.inclusive_jets().to_list()
    ]
    constituents = [
        [
            [{"phi": 1.2, "eta": 3.2, "rho": 5.4, "E": 2.5, "ex": 0.78}],
            [{"phi": 32.45, "eta": 63.21, "rho": 543.14, "E": 24.56, "ex": 0.0}],
            [{"phi": 32.2, "eta": 64.21, "rho": 543.34, "E": 24.12, "ex": 0.35}],
        ],
        [
            [{"phi": 1.2, "eta": 3.2, "rho": 5.4, "E": 2.5, "ex": 0.78}],
            [{"phi": 32.45, "eta": 63.21, "rho": 543.14, "E": 24.56, "ex": 0.0}],
            [{"phi": 32.2, "eta": 64.21, "rho": 543.34, "E": 24.12, "ex": 0.35}],
        ],
    ]

    assert constituents == [
        [[pytest.approx(z) for z in y] for y in x]
        for x in cluster.constituents().to_list()
    ]
    constituent_index = [[[0], [2], [1]], [[0], [2], [1]]]
    assert constituent_index == [
        [pytest.approx(y) for y in x] for x in cluster.constituent_index().to_list()
    ]
