import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401

dak = pytest.importorskip("dask_awkward")  # noqa: F401
vector = pytest.importorskip("vector")  # noqa: F401


def test_multi():
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
    darray = dak.from_awkward(array, 1)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.DaskAwkwardClusterSequence(darray, jetdef)
    inclusive_jets_out = [
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
    ]
    assert inclusive_jets_out == cluster.inclusive_jets().compute().to_list()
    constituents_output = [
        [
            [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78}],
            [
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
        [
            [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78}],
            [
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    ]
    assert cluster.constituents().compute().to_list() == constituents_output
    constituent_index_out = [[[0], [1, 2]], [[0], [1, 2]]]
    assert constituent_index_out == cluster.constituent_index().compute().to_list()


def test_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ]
    )
    darray = dak.from_awkward(array, 1)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.DaskAwkwardClusterSequence(darray, jetdef)
    inclusive_jets_out = [
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
    ]

    assert inclusive_jets_out == cluster.inclusive_jets().compute().to_list()
    constituent_output = [
        [{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5}],
        [
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        ],
    ]
    assert constituent_output == cluster.constituents().compute().to_list()
    constituent_index_output = [[0], [1, 2]]
    assert constituent_index_output == cluster.constituent_index().compute().to_list()


def test_inclusive_from_file():
    import uproot
    from dask_awkward.lib.testutils import assert_eq

    vector = pytest.importorskip("vector")

    testfile = "https://github.com/scikit-hep/fastjet/blob/main/tests/samples/pfnano_skim.root?raw=true"

    devents = uproot.dask({testfile: "Events"})

    dpfcands = dak.zip(
        {
            "pt": devents.PFCands_pt,
            "eta": devents.PFCands_eta,
            "phi": devents.PFCands_phi,
            "mass": devents.PFCands_mass,
        },
        with_name="Momentum4D",
        behavior=vector.backends.awkward.behavior,
    )

    djetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    dcluseq = fastjet.ClusterSequence(dpfcands, djetdef)

    events = uproot.open({testfile: "Events"}).arrays(
        ["PFCands_pt", "PFCands_eta", "PFCands_phi", "PFCands_mass"]
    )

    pfcands = ak.zip(
        {
            "pt": events.PFCands_pt,
            "eta": events.PFCands_eta,
            "phi": events.PFCands_phi,
            "mass": events.PFCands_mass,
        },
        with_name="Momentum4D",
        behavior=vector.backends.awkward.behavior,
    )

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    cluseq = fastjet.ClusterSequence(pfcands, jetdef)

    assert_eq(dcluseq.inclusive_jets(), cluseq.inclusive_jets())
