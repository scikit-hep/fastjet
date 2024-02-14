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


def test_exclusive_up_to_multi():
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
    # cluster.exclusive_jets(n_jets=4) raises a RuntimeError
    # because there are only 3 constitutents
    with pytest.raises(RuntimeError):
        cluster.exclusive_jets(n_jets=4)
    # cluster.exclusive_jets_up_to(njets=4) returns the 3 constituents
    assert multi_exclusive_njets == cluster.exclusive_jets_up_to(n_jets=4).to_list()


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

    # expected output from x86_64, aarch64 is slightly different
    lund_declustering_output = ak.Array(
        [
            [
                {"Delta": 0.08755181299980186, "kt": 0.30179987478357223},
                {"Delta": 0.019481226884377707, "kt": 0.06602095529127928},
            ],
            [{"Delta": 0.014750342295225208, "kt": 1.0480537658466145}],
        ]
    )

    is_close = ak.ravel(ak.isclose(lund_declustering_output, lds, rtol=1e-12, atol=0))

    assert ak.all(is_close)


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

    # expected output from x86_64, aarch64 is slightly different
    lund_declustering_output = ak.Array(
        [
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
    )

    is_close = ak.ravel(ak.isclose(lund_declustering_output, lds, rtol=1e-12, atol=0))

    assert ak.all(is_close)


def test_exclusive_jets_softdrop_grooming():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
            {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56, "ex": 0.0},
        ],
        with_name="Momentum4D",
    )

    jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, 0.8)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    softdrop = cluster.exclusive_jets_softdrop_grooming()

    softdrop_output = ak.from_iter(
        {
            "constituents": [
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56},
            ],
            "msoftdrop": 488.2395243115817,
            "ptsoftdrop": 142.88274528437645,
            "etasoftdrop": 2.726117171791057,
            "phisoftdrop": 1.1012644074821902,
            "Esoftdrop": 1199.6799999999998,
            "pzsoftdrop": 1086.48,
        },
    )

    is_close = ak.ravel(
        ak.Array(
            [
                ak.isclose(
                    softdrop_output.constituents.px,
                    softdrop.constituents.px,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.py,
                    softdrop.constituents.py,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.pz,
                    softdrop.constituents.pz,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.E,
                    softdrop.constituents.E,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.msoftdrop]),
                    ak.Array([softdrop.msoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.ptsoftdrop]),
                    ak.Array([softdrop.ptsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.etasoftdrop]),
                    ak.Array([softdrop.etasoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.phisoftdrop]),
                    ak.Array([softdrop.phisoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.Esoftdrop]),
                    ak.Array([softdrop.Esoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.pzsoftdrop]),
                    ak.Array([softdrop.pzsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
            ]
        )
    )

    assert ak.all(is_close)


def test_exclusive_jets_softdrop_grooming_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
                {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 1.25, "py": 3.15, "pz": 5.4, "E": 2.4, "ex": 0.78},
                {"px": 1.4, "py": 3.15, "pz": 5.4, "E": 2.0, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56, "ex": 0.0},
            ],
        ],
        with_name="Momentum4D",
    )

    jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, 0.8)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    softdrop = cluster.exclusive_jets_softdrop_grooming()

    softdrop_output = ak.from_iter(
        {
            "constituents": ak.Array(
                [
                    [
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56},
                    ],
                    [
                        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12},
                        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56},
                    ],
                ],
            ),
            "msoftdrop": ak.Array([488.2395243115817, 488.2395243115817]),
            "ptsoftdrop": ak.Array([142.88274528437645, 142.88274528437645]),
            "etasoftdrop": ak.Array([2.726117171791057, 2.726117171791057]),
            "phisoftdrop": ak.Array([1.1012644074821902, 1.1012644074821902]),
            "Esoftdrop": ak.Array([1199.6799999999998, 1199.6799999999998]),
            "pzsoftdrop": ak.Array([1086.48, 1086.48]),
        }
    )

    is_close = ak.ravel(
        ak.Array(
            [
                ak.isclose(
                    softdrop_output.constituents.px,
                    softdrop.constituents.px,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.py,
                    softdrop.constituents.py,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.pz,
                    softdrop.constituents.pz,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    softdrop_output.constituents.E,
                    softdrop.constituents.E,
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.msoftdrop]),
                    ak.Array([softdrop.msoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.ptsoftdrop]),
                    ak.Array([softdrop.ptsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.etasoftdrop]),
                    ak.Array([softdrop.etasoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.phisoftdrop]),
                    ak.Array([softdrop.phisoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.Esoftdrop]),
                    ak.Array([softdrop.Esoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.pzsoftdrop]),
                    ak.Array([softdrop.pzsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
            ]
        )
    )

    assert ak.all(is_close)


def test_exclusive_energy_correlator():
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

    ec1 = cluster.exclusive_jets_energy_correlator(func="generic", npoint=1)
    ec2 = cluster.exclusive_jets_energy_correlator(func="generic", npoint=2)
    ecg2 = cluster.exclusive_jets_energy_correlator(
        func="generalized", npoint=2, angles=1
    )

    is_close = ak.ravel(
        ak.isclose(ak.Array([ec2 / ec1 / ec1]), ak.Array([ecg2]), rtol=1e-12, atol=0)
    )

    assert ak.all(is_close)


def test_exclusive_energy_correlator_multi():
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

    ec1 = cluster.exclusive_jets_energy_correlator(func="generic", npoint=1)
    ec2 = cluster.exclusive_jets_energy_correlator(func="generic", npoint=2)
    ecg2 = cluster.exclusive_jets_energy_correlator(
        func="generalized", npoint=2, angles=1
    )

    is_close = ak.ravel(ak.isclose((ec2 / ec1 / ec1), ecg2, rtol=1e-12, atol=0))

    assert ak.all(is_close)


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
