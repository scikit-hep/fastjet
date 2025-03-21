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
            "deltaRsoftdrop": 0.009873899817126915,
            "symmetrysoftdrop": 0.49727522889673303,
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
                    ak.Array([softdrop_output.deltaRsoftdrop]),
                    ak.Array([softdrop.deltaRsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.symmetrysoftdrop]),
                    ak.Array([softdrop.symmetrysoftdrop]),
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
            "deltaRsoftdrop": ak.Array([0.009873899817126915, 0.009873899817126915]),
            "symmetrysoftdrop": ak.Array([0.49727522889673303, 0.49727522889673303]),
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
                    ak.Array([softdrop_output.deltaRsoftdrop]),
                    ak.Array([softdrop.deltaRsoftdrop]),
                    rtol=1e-12,
                    atol=0,
                ),
                ak.isclose(
                    ak.Array([softdrop_output.symmetrysoftdrop]),
                    ak.Array([softdrop.symmetrysoftdrop]),
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


def test_exclusive_constituents_masking_multi():
    array = ak.Array(
        [
            [
                {"px": -0.181, "py": -0.798, "pz": -0.652, "E": 1.06},
                {"px": 0.377, "py": 0.116, "pz": -0.0749, "E": 0.425},
            ],
            [
                {"px": 0.54, "py": -0.65, "pz": -0.00527, "E": 0.857},
                {"px": 0.253, "py": -0.0342, "pz": -0.0731, "E": 0.3},
            ],
            [
                {"px": 0.294, "py": 0.254, "pz": -0.259, "E": 0.467},
                {"px": 4.65, "py": -1.88, "pz": -3.29, "E": 6},
            ],
            [
                {"px": 1.45, "py": -0.179, "pz": -0.876, "E": 1.71},
                {"px": 12.8, "py": -1.8, "pz": -7.18, "E": 14.8},
            ],
            [
                {"px": -3.55, "py": -1.64, "pz": -0.0941, "E": 3.91},
                {"px": -1.33, "py": -1.03, "pz": 0.147, "E": 1.7},
            ],
        ],
        with_name="Momentum4D",
    )
    mask = ak.Array([True, True, False, True, True])
    jetdef = fastjet.JetDefinition(fastjet.kt_algorithm, 1)
    result_all = fastjet.ClusterSequence(
        array, jetdef
    ).exclusive_jets_constituent_index(njets=2)
    result_mask = fastjet.ClusterSequence(
        array[mask], jetdef
    ).exclusive_jets_constituent_index(njets=2)
    assert ak.all(result_mask == result_all[mask])


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


def test_njettiness():
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

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8)
    jet_cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)

    constituents = jet_cluster.constituents()

    njettiness_cluster = fastjet._pyjet.AwkwardClusterSequence(constituents, jetdef)

    expected = ak.Array(
        [
            [
                [0.041999587156887445, 0.008082162710689409, 0.0, 0.0],
                [0.00921855860347245, 0.0, 0.0, 0.0],
            ],
            [
                [0.041999587156887445, 0.008082162710689409, 0.0, 0.0],
                [0.00921855860347245, 0.0, 0.0, 0.0],
            ],
        ]
    )

    result = njettiness_cluster.njettiness()

    assert ak.all(ak.isclose(result, expected))


def test_njettiness_real_jet():
    import vector

    vector.register_awkward()
    # this is from event 48092, jet 0 of
    # /RunIII2024Summer24NanoAODv15/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM
    # /150X_mcRun3_2024_realistic_v1-v2/120000/26ed0adf-21ff-4992-8f75-8726fd037d7c.root
    array = ak.Array(
        [
            {
                "pt": 0.8681640625,
                "eta": 1.7900390625,
                "phi": -0.713134765625,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.324462890625,
                "eta": 2.21337890625,
                "phi": -0.53466796875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.37060546875,
                "eta": 2.16748046875,
                "phi": 0.158294677734375,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.416259765625,
                "eta": 1.89306640625,
                "phi": -0.580810546875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.47216796875,
                "eta": 1.88671875,
                "phi": 0.5994873046875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.69287109375,
                "eta": 1.4228515625,
                "phi": -0.4420166015625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.681640625,
                "eta": 1.381103515625,
                "phi": 0.5267333984375,
                "mass": 0.138427734375,
            },
            {
                "pt": 0.489990234375,
                "eta": 1.230224609375,
                "phi": 0.206939697265625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.69970703125,
                "eta": 1.79248046875,
                "phi": 0.40966796875,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.6171875,
                "eta": 1.584228515625,
                "phi": 0.239959716796875,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.5146484375,
                "eta": 1.87353515625,
                "phi": 0.024799346923828125,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.724609375,
                "eta": 1.466552734375,
                "phi": 0.54150390625,
                "mass": 0.138427734375,
            },
            {
                "pt": 0.861328125,
                "eta": 1.75341796875,
                "phi": 0.60888671875,
                "mass": 0.1346435546875,
            },
            {
                "pt": 1.005859375,
                "eta": 2.11376953125,
                "phi": -0.552734375,
                "mass": 0.1395263671875,
            },
            {
                "pt": 16.734375,
                "eta": 1.852294921875,
                "phi": -0.353515625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 2.779296875,
                "eta": 1.8505859375,
                "phi": 0.3648681640625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 4.4765625,
                "eta": 1.789306640625,
                "phi": -0.2935791015625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 4.359375,
                "eta": 1.765625,
                "phi": -0.29522705078125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 3.85546875,
                "eta": 1.75390625,
                "phi": 0.177459716796875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 4.125,
                "eta": 1.667724609375,
                "phi": 0.30877685546875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.1630859375,
                "eta": 1.60986328125,
                "phi": 0.45782470703125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 2.12109375,
                "eta": 1.5986328125,
                "phi": 0.239349365234375,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.892578125,
                "eta": 1.560546875,
                "phi": 0.228424072265625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 4.4609375,
                "eta": 1.5458984375,
                "phi": 0.33154296875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.724609375,
                "eta": 1.47998046875,
                "phi": 0.26092529296875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 3.21875,
                "eta": 1.412841796875,
                "phi": 0.43035888671875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 19.125,
                "eta": 1.37841796875,
                "phi": -0.29901123046875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 2.357421875,
                "eta": 1.279541015625,
                "phi": 0.35614013671875,
                "mass": -5.960464477539063e-08,
            },
            {
                "pt": 1.162109375,
                "eta": 1.268798828125,
                "phi": 0.4775390625,
                "mass": 0.1373291015625,
            },
            {
                "pt": 3.35546875,
                "eta": 2.00439453125,
                "phi": 0.232818603515625,
                "mass": 0.134033203125,
            },
            {
                "pt": 7.69921875,
                "eta": 2.0556640625,
                "phi": 0.151092529296875,
                "mass": -4.2445026338100433e-07,
            },
            {
                "pt": 5.75390625,
                "eta": 1.959716796875,
                "phi": 0.5150146484375,
                "mass": 0.0,
            },
            {
                "pt": 4.16015625,
                "eta": 1.876220703125,
                "phi": -0.2843017578125,
                "mass": -1.1874362826347351e-07,
            },
            {
                "pt": 17.375,
                "eta": 1.798583984375,
                "phi": -0.29998779296875,
                "mass": 0.0,
            },
            {
                "pt": 1.634765625,
                "eta": 1.72314453125,
                "phi": -0.45245361328125,
                "mass": 5.9371814131736755e-08,
            },
            {
                "pt": 1.25390625,
                "eta": 1.6884765625,
                "phi": -0.4105224609375,
                "mass": -5.9371814131736755e-08,
            },
            {
                "pt": 0.95751953125,
                "eta": 1.671875,
                "phi": 0.03612518310546875,
                "mass": -0.0,
            },
            {
                "pt": 3.716796875,
                "eta": 1.64990234375,
                "phi": -0.4093017578125,
                "mass": 0.0,
            },
            {
                "pt": 1.912109375,
                "eta": 1.56494140625,
                "phi": -0.29974365234375,
                "mass": -5.9371814131736755e-08,
            },
            {
                "pt": 1.47265625,
                "eta": 1.547119140625,
                "phi": 0.208709716796875,
                "mass": 0.0,
            },
            {
                "pt": 2.3828125,
                "eta": 1.476806640625,
                "phi": 0.48907470703125,
                "mass": -5.9371814131736755e-08,
            },
            {
                "pt": 1.5654296875,
                "eta": 1.435791015625,
                "phi": 0.238861083984375,
                "mass": -0.0,
            },
            {
                "pt": 0.63623046875,
                "eta": 1.40087890625,
                "phi": 0.114166259765625,
                "mass": 0.0,
            },
            {
                "pt": 1.2880859375,
                "eta": 1.39501953125,
                "phi": 0.37103271484375,
                "mass": 0.0,
            },
            {
                "pt": 1.7900390625,
                "eta": 1.394287109375,
                "phi": 0.2047119140625,
                "mass": -0.0,
            },
            {
                "pt": 6.14453125,
                "eta": 1.387939453125,
                "phi": 0.33538818359375,
                "mass": 0.0,
            },
            {
                "pt": 3.115234375,
                "eta": 1.383544921875,
                "phi": -0.3463134765625,
                "mass": 0.0,
            },
            {
                "pt": 2.564453125,
                "eta": 1.3603515625,
                "phi": 0.2799072265625,
                "mass": -5.9138983488082886e-08,
            },
            {
                "pt": 0.5986328125,
                "eta": 1.332763671875,
                "phi": -0.46795654296875,
                "mass": -0.0,
            },
            {
                "pt": 11.171875,
                "eta": 1.330322265625,
                "phi": -0.49627685546875,
                "mass": -2.968590706586838e-07,
            },
            {
                "pt": 7.11328125,
                "eta": 1.26611328125,
                "phi": 0.222564697265625,
                "mass": -1.1874362826347351e-07,
            },
            {
                "pt": 0.6337890625,
                "eta": 1.058349609375,
                "phi": -0.24041748046875,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.342529296875,
                "eta": 1.0537109375,
                "phi": -0.119140625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.2100830078125,
                "eta": 0.819091796875,
                "phi": -0.5205078125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.52978515625,
                "eta": 1.127197265625,
                "phi": 0.19189453125,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.6533203125,
                "eta": 0.895263671875,
                "phi": -0.48577880859375,
                "mass": 0.1390380859375,
            },
            {
                "pt": 2.06640625,
                "eta": 0.8082275390625,
                "phi": -0.4981689453125,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.646484375,
                "eta": 1.07568359375,
                "phi": 0.6129150390625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.2626953125,
                "eta": 1.053955078125,
                "phi": -0.2158203125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 0.95458984375,
                "eta": 0.9580078125,
                "phi": 0.136993408203125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 25.515625,
                "eta": 0.8839111328125,
                "phi": 0.205963134765625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.32421875,
                "eta": 0.8670654296875,
                "phi": -0.30401611328125,
                "mass": 0.1395263671875,
            },
            {
                "pt": 5.73828125,
                "eta": 0.864990234375,
                "phi": -0.20947265625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 2.837890625,
                "eta": 0.8509521484375,
                "phi": 0.04062652587890625,
                "mass": 0.1390380859375,
            },
            {
                "pt": 14.0625,
                "eta": 0.84326171875,
                "phi": -0.246002197265625,
                "mass": 0.1395263671875,
            },
            {
                "pt": 7.0546875,
                "eta": 0.8408203125,
                "phi": 0.25,
                "mass": 0.1395263671875,
            },
            {
                "pt": 3.181640625,
                "eta": 0.7752685546875,
                "phi": -0.40350341796875,
                "mass": 0.1395263671875,
            },
            {
                "pt": 1.666015625,
                "eta": 0.994873046875,
                "phi": -0.226165771484375,
                "mass": 0.1390380859375,
            },
            {
                "pt": 0.60498046875,
                "eta": 1.120361328125,
                "phi": -0.023822784423828125,
                "mass": 0.0,
            },
            {
                "pt": 1.1396484375,
                "eta": 0.9481201171875,
                "phi": -0.5491943359375,
                "mass": 0.0,
            },
            {
                "pt": 1.06640625,
                "eta": 0.9312744140625,
                "phi": 0.018955230712890625,
                "mass": -0.0,
            },
            {
                "pt": 1.544921875,
                "eta": 0.914306640625,
                "phi": 0.3428955078125,
                "mass": 0.0,
            },
            {
                "pt": 4.33984375,
                "eta": 0.906005859375,
                "phi": 0.26416015625,
                "mass": -5.960464477539063e-08,
            },
            {
                "pt": 4.07421875,
                "eta": 0.8682861328125,
                "phi": -0.2730712890625,
                "mass": -5.9371814131736755e-08,
            },
            {
                "pt": 3.478515625,
                "eta": 0.8414306640625,
                "phi": -0.29803466796875,
                "mass": -5.9371814131736755e-08,
            },
            {
                "pt": 0.8251953125,
                "eta": 0.8304443359375,
                "phi": -0.10791015625,
                "mass": -0.0,
            },
            {
                "pt": 0.51171875,
                "eta": 0.8106689453125,
                "phi": -0.2691650390625,
                "mass": -0.0,
            },
            {
                "pt": 2.908203125,
                "eta": 0.810302734375,
                "phi": -0.35333251953125,
                "mass": -0.0,
            },
            {
                "pt": 1.61328125,
                "eta": 0.6851806640625,
                "phi": 0.13214111328125,
                "mass": 0.0,
            },
            {
                "pt": 0.66064453125,
                "eta": 0.646240234375,
                "phi": 0.131072998046875,
                "mass": -0.0,
            },
            {
                "pt": 1.3525390625,
                "eta": 0.6334228515625,
                "phi": -0.381591796875,
                "mass": -0.0,
            },
            {
                "pt": 0.66162109375,
                "eta": 0.581787109375,
                "phi": -0.127044677734375,
                "mass": 0.0,
            },
        ],
        with_name="Momentum4D",
    )

    # expected taus from NanoAOD
    expected_taus = ak.Array([[0.6298828125, 0.420166015625, 0.302734375, 0.220703125]])

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array[None], jetdef)

    result = cluster.njettiness()

    assert ak.all(ak.isclose(expected_taus, result, atol=1e-3))
