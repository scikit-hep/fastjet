import awkward as ak
import pytest

import fastjet

vector = pytest.importorskip("vector")


def _events():
    # every event differs, so clustering the wrong particles cannot go unnoticed
    return ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
            ],
            [
                {"px": -7.1, "py": 2.9, "pz": 11.0, "E": 14.0},
                {"px": 4.4, "py": -8.3, "pz": -2.0, "E": 10.1},
            ],
            [
                {"px": 0.3, "py": -0.2, "pz": 9.9, "E": 10.0},
                {"px": 15.0, "py": 14.0, "pz": -3.0, "E": 21.0},
                {"px": 15.2, "py": 13.7, "pz": -3.1, "E": 21.1},
                {"px": -20.0, "py": 1.0, "pz": 1.0, "E": 20.5},
            ],
            [
                {"px": 5.0, "py": 5.0, "pz": 5.0, "E": 9.0},
            ],
        ],
        with_name="Momentum4D",
    )


def _jets(array):
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    jets = fastjet.ClusterSequence(array, jetdef).inclusive_jets()
    return ak.zip({"px": jets.px, "py": jets.py, "pz": jets.pz, "E": jets.E}).to_list()


def test_lists_not_starting_at_zero():
    array = _events()
    sliced = array[1:]
    assert sliced.layout.starts[0] != 0
    assert _jets(sliced) == _jets(array)[1:]


def test_lists_out_of_order():
    array = _events()
    order = [2, 0, 3, 1]
    picked = array[order]
    assert list(picked.layout.starts) != sorted(picked.layout.starts)
    expected = _jets(array)
    assert _jets(picked) == [expected[i] for i in order]


def test_dask_partitions():
    dak = pytest.importorskip("dask_awkward")
    array = _events()
    lazy = dak.from_awkward(array, npartitions=2)
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    jets = fastjet.ClusterSequence(lazy, jetdef).inclusive_jets().compute()
    got = ak.zip({"px": jets.px, "py": jets.py, "pz": jets.pz, "E": jets.E}).to_list()
    assert got == _jets(array)
