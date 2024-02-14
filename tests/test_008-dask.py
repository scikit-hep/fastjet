import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401

dak = pytest.importorskip("dask_awkward")  # noqa: F401
distributed = pytest.importorskip("distributed")
vector = pytest.importorskip("vector")  # noqa: F401


def test_multi():
    from distributed import Client

    with Client() as _:
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
    from distributed import Client

    with Client() as _:
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
        assert (
            constituent_index_output == cluster.constituent_index().compute().to_list()
        )


def test_inclusive_from_file():
    from pathlib import Path

    import uproot
    from dask_awkward.lib.testutils import assert_eq
    from distributed import Client

    with Client() as _:
        data_dir = Path(__file__).parent / "samples"
        testfile = data_dir / "pfnano_skim.root"

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


def dask_multi_test_exclusive_jets_softdrop_grooming():
    from distributed import Client

    with Client() as _:
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
        darray = dak.from_awkward(array, 1)
        jetdef = fastjet.JetDefinition(fastjet.cambridge_algorithm, 0.8)
        cluster = fastjet._pyjet.DaskAwkwardClusterSequence(darray, jetdef)
        softdrop = cluster.exclusive_jets_softdrop_grooming().compute()

        softdrop_output = ak.from_iter(
            {
                "constituents": [
                    {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 600.12},
                    {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 599.56},
                ],
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
