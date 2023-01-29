import awkward as ak  # noqa: F401
import numpy as np  # noqa: F401
import pytest  # noqa: F401

import fastjet._pyjet  # noqa: F401

vector = pytest.importorskip("vector")  # noqa: F401


def test_unique_history_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [[0, 5, 1, 2, 3, 4], [0, 3, 1, 2], [0, 5, 1, 2, 3, 4]]
    assert out == cluster.unique_history_order().to_list()


def test_unique_history_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [0, 5, 1, 2, 3, 4]
    assert out == cluster.unique_history_order().to_list()


def test_n_particles_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = 3
    assert out == cluster.n_particles()


def test_n_particles_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [3, 2, 3]
    assert out == cluster.n_particles().to_list()


def test_n_exclusive_jets_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [1, 2, 1]
    assert out == cluster.n_exclusive_jets(dcut=0.0001).to_list()


def test_n_exclusive_jets_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = 1
    assert out == cluster.n_exclusive_jets(dcut=0.0001)


def test_unclustered_particles_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [[], [], []]
    assert out == cluster.unclustered_particles().to_list()


def test_unclustered_particles_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = []
    assert out == cluster.unclustered_particles().to_list()


def test_childless_pseudojets_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = []
    assert out == cluster.childless_pseudojets().to_list()


def test_childless_pseudojets_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [[], [], []]
    assert out == cluster.childless_pseudojets().to_list()


def test_jets_single():
    array = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
        {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
    ]
    assert out == cluster.jets().to_list()


def test_jets_multi():
    array = ak.Array(
        [
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
            [
                {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
                {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
                {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
            ],
        ],
    )
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    cluster = fastjet._pyjet.AwkwardClusterSequence(array, jetdef)
    out = [
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        ],
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
            {"px": 64.65, "py": 127.41999999999999, "pz": 1086.48, "E": 48.68},
        ],
    ]
    assert out == cluster.jets().to_list()


def test_constituent_index_min_pt():
    # test for https://github.com/scikit-hep/fastjet/issues/174
    array = ak.Array(
        [
            [
                {
                    "px": 0.32786062359809875,
                    "py": -0.8113271594047546,
                    "pz": -0.9425322413444519,
                    "E": 0.1395701766014099,
                },
                {
                    "px": 1.2147544622421265,
                    "py": -6.664909839630127,
                    "pz": -0.09962931275367737,
                    "E": 0.1395701766014099,
                },
                {
                    "px": 0.031138502061367035,
                    "py": -0.03814707323908806,
                    "pz": -0.19269536435604095,
                    "E": 0.0,
                },
                {
                    "px": 0.25616654753685,
                    "py": -2.1868643760681152,
                    "pz": -1.1193745136260986,
                    "E": 0.0,
                },
                {
                    "px": 0.1986897736787796,
                    "py": -0.08637325465679169,
                    "pz": 0.16585613787174225,
                    "E": 0.0,
                },
                {
                    "px": 7.569584846496582,
                    "py": -8.49626636505127,
                    "pz": 8.851030349731445,
                    "E": 0.9395653009414673,
                },
            ],
            [
                {
                    "px": -1.995545744895935,
                    "py": 1.0673128366470337,
                    "pz": -6.9042792320251465,
                    "E": 0.1395701766014099,
                },
                {
                    "px": -9.267221450805664,
                    "py": 2.7426843643188477,
                    "pz": -0.42873385548591614,
                    "E": 0.1395701766014099,
                },
                {
                    "px": -1.5686458349227905,
                    "py": -0.1490965336561203,
                    "pz": -0.23429453372955322,
                    "E": 0.9395653009414673,
                },
                {
                    "px": -1.043339729309082,
                    "py": 0.6251208186149597,
                    "pz": -1.4011896848678589,
                    "E": 0.0,
                },
                {
                    "px": 2.5927963256835938,
                    "py": -4.718469142913818,
                    "pz": -5.502776145935059,
                    "E": 0.0,
                },
            ],
        ]
    )

    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    cluster = fastjet.ClusterSequence(array, jetdef)
    constituent_idx = cluster.constituent_index(min_pt=0.5)

    assert constituent_idx.to_list() == [[[3], [1], [0, 5]], [[3], [2], [0], [4], [1]]]
    constituents_flat = ak.flatten(constituent_idx, axis=-1)
    array_flat = array[constituents_flat].to_list()

    array_flat_expected = ak.Array(
        [
            [
                {
                    "px": 0.25616654753685,
                    "py": -2.1868643760681152,
                    "pz": -1.1193745136260986,
                    "E": 0.0,
                },
                {
                    "px": 1.2147544622421265,
                    "py": -6.664909839630127,
                    "pz": -0.09962931275367737,
                    "E": 0.1395701766014099,
                },
                {
                    "px": 0.32786062359809875,
                    "py": -0.8113271594047546,
                    "pz": -0.9425322413444519,
                    "E": 0.1395701766014099,
                },
                {
                    "px": 7.569584846496582,
                    "py": -8.49626636505127,
                    "pz": 8.851030349731445,
                    "E": 0.9395653009414673,
                },
            ],
            [
                {
                    "px": -1.043339729309082,
                    "py": 0.6251208186149597,
                    "pz": -1.4011896848678589,
                    "E": 0.0,
                },
                {
                    "px": -1.5686458349227905,
                    "py": -0.1490965336561203,
                    "pz": -0.23429453372955322,
                    "E": 0.9395653009414673,
                },
                {
                    "px": -1.995545744895935,
                    "py": 1.0673128366470337,
                    "pz": -6.9042792320251465,
                    "E": 0.1395701766014099,
                },
                {
                    "px": 2.5927963256835938,
                    "py": -4.718469142913818,
                    "pz": -5.502776145935059,
                    "E": 0.0,
                },
                {
                    "px": -9.267221450805664,
                    "py": 2.7426843643188477,
                    "pz": -0.42873385548591614,
                    "E": 0.1395701766014099,
                },
            ],
        ]
    )
    is_close = ak.ravel(ak.isclose(array_flat_expected, array_flat, rtol=1e-12, atol=0))
    assert ak.all(is_close)
