import awkward as ak
import pytest

import fastjet

graphed = pytest.importorskip("graphed")
ga = pytest.importorskip("graphed.awkward")
vector = pytest.importorskip("vector")


def _events():
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


# one jet per event, for the queries that take a second array
_JETS = ak.Array(
    [
        {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
        {"px": -7.1, "py": 2.9, "pz": 11.0, "E": 14.0},
        {"px": 0.3, "py": -0.2, "pz": 9.9, "E": 10.0},
        {"px": 5.0, "py": 5.0, "pz": 5.0, "E": 9.0},
    ]
)

# the queries tests/test_010-graphed.py does not reach, with arguments this sample clusters to
QUERIES = [
    ("Q", {}, False),
    ("Q2", {}, False),
    ("childless_pseudojets", {}, False),
    ("exclusive_dmerge", {"njets": 1}, False),
    ("exclusive_dmerge_max", {"njets": 1}, False),
    ("exclusive_jets_constituent_index", {"njets": 1}, False),
    ("exclusive_jets_constituents", {"njets": 1}, False),
    ("exclusive_jets_energy_correlator", {"njets": 1}, False),
    ("exclusive_jets_lund_declusterings", {"njets": 1}, False),
    ("exclusive_jets_softdrop_grooming", {"njets": 1}, False),
    ("exclusive_jets_up_to", {"n_jets": 1}, False),
    ("exclusive_jets_ycut", {"ycut": 0.0001}, False),
    ("exclusive_ymerge", {"njets": 1}, False),
    ("exclusive_ymerge_max", {"njets": 1}, False),
    ("jets", {}, False),
    ("n_exclusive_jets", {"dcut": 0.0001}, False),
    ("n_particles", {}, False),
    ("njettiness", {}, False),
    ("unique_history_order", {}, False),
    ("exclusive_subdmerge", {"nsub": 1}, True),
    ("exclusive_subdmerge_max", {"nsub": 1}, True),
    ("exclusive_subjets_up_to", {"nsub": 1}, True),
    ("get_child", {}, True),
    ("has_child", {}, True),
    ("has_parents", {}, True),
    ("jet_scale_for_algorithm", {}, True),
    ("n_exclusive_subjets", {"dcut": 0.0001}, True),
]


@pytest.mark.parametrize(("query", "kwargs", "takes_jets"), QUERIES)
def test_query_form_and_value_match_eager(query, kwargs, takes_jets):
    array = _events()
    jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    session = graphed.Session(
        ga.AwkwardBackend(behavior=vector.backends.awkward.behavior)
    )
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", array), jetdef)
    eager = fastjet.ClusterSequence(array, jetdef)
    args = [ga.from_awkward(session, "jets", _JETS)] if takes_jets else []

    out = getattr(cluseq, query)(*args, **kwargs)
    expected = getattr(eager, query)(*([_JETS] if takes_jets else []), **kwargs)

    assert str(session.form(out).tt.type) == str(
        ak.Array(expected.layout.to_typetracer(forget_length=True)).type
    )
    assert session.materialize(out).to_list() == expected.to_list()


def test_a_second_array_must_be_deferred_too():
    session = graphed.Session(ga.AwkwardBackend())
    cluseq = fastjet.ClusterSequence(
        ga.from_awkward(session, "ev", _events()),
        fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6),
    )
    with pytest.raises(TypeError):
        cluseq.get_parents(_JETS)
