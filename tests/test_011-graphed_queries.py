import subprocess
import sys

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


def _one_hard_particle_and_two_soft_ones():
    # a hard particle at eta = 0 and soft ones at eta = 0.5 and 1.0. anti-kt pulls the near
    # soft one into the hard jet at R = 0.6 and leaves all three alone at R = 0.4; kt merges
    # the two soft ones with each other instead, so radius and algorithm each change the jets
    return ak.Array(
        [
            [
                {"px": 100.0, "py": 0.0, "pz": 0.0, "E": 100.0},
                {
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 0.5210953054937474,
                    "E": 1.1276259652063807,
                },
                {
                    "px": 1.0,
                    "py": 0.0,
                    "pz": 1.1752011936438014,
                    "E": 1.5430806348152437,
                },
            ]
        ],
        with_name="Momentum4D",
    )


def test_the_jet_definition_is_part_of_the_recorded_node():
    array = _one_hard_particle_and_two_soft_ones()
    session = graphed.Session(
        ga.AwkwardBackend(behavior=vector.backends.awkward.behavior)
    )
    source = ga.from_awkward(session, "ev", array)

    def deferred(jetdef):
        return fastjet.ClusterSequence(source, jetdef).inclusive_jets()

    def eager(jetdef):
        return fastjet.ClusterSequence(array, jetdef).inclusive_jets().to_list()

    narrow = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.4)
    wide = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
    kt_wide = fastjet.JetDefinition(fastjet.kt_algorithm, 0.6)

    jets = [deferred(jetdef) for jetdef in (narrow, wide, kt_wide)]
    assert len({jet.node_id for jet in jets}) == 3

    # radius alone and algorithm alone each change the answer, so sharing a node would hand
    # back the other definition's jets rather than merely reusing an equivalent one
    assert eager(narrow) != eager(wide) != eager(kt_wide)
    for jetdef, jet in zip((narrow, wide, kt_wide), jets):
        assert session.materialize(jet).to_list() == eager(jetdef)


_AVAILABILITY_PROBE = """
import os
import sys

import fastjet  # the dispatch imports graphed lazily, so only that import is at stake

jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
if sys.argv[1] == "uninstalled":
    for name in [m for m in sys.modules if m.split(".")[0] == "graphed"]:
        del sys.modules[name]
    sys.path[:] = [
        p
        for p in sys.path
        if not os.path.exists(os.path.join(p, "graphed", "__init__.py"))
    ]
    import graphed  # the empty directory in the cwd, as a namespace package

    assert not hasattr(graphed, "Array"), graphed.__path__
    try:
        fastjet.ClusterSequence("not an array", jetdef)
    except TypeError:
        print("TypeError")
else:
    import awkward as ak
    import graphed
    import graphed.awkward as ga

    session = graphed.Session(ga.AwkwardBackend())
    particles = ga.from_awkward(
        session, "ev", ak.Array([[{"px": 1.0, "py": 2.0, "pz": 3.0, "E": 4.0}]])
    )
    print(type(fastjet.ClusterSequence(particles, jetdef)).__name__)
"""


def test_an_empty_graphed_directory_is_not_a_graphed_install(tmp_path):
    (tmp_path / "graphed").mkdir()

    def run(mode):
        probe = subprocess.run(
            [sys.executable, "-c", _AVAILABILITY_PROBE, mode],
            cwd=tmp_path,
            capture_output=True,
            text=True,
            check=True,
        )
        return probe.stdout.splitlines()[-1]

    # an installed graphed shadows the directory, and the deferred arm is still reached
    assert run("installed") == "GraphedClusterSequence"
    # without one, `import graphed` succeeds over that directory and carries no names
    assert run("uninstalled") == "TypeError"


def test_every_universe_of_a_varied_input_must_hold_one_event_per_entry():
    session = graphed.Session(
        ga.AwkwardBackend(behavior=vector.backends.awkward.behavior)
    )
    jagged = ga.from_awkward(session, "ev", _events())
    varied = graphed.vary(
        jagged, "squash", points={"flat": ga.gak.flatten(jagged, axis=1)}
    )
    # the nominal is one list of particles per event and the varied universe is not
    assert [
        session.form(graphed.universe(varied, label)).tt.ndim
        for label in graphed.labels(varied)
    ] == [2, 1]

    before = session.node_count()
    # checking only the nominal would defer a universe silently clustered as one event,
    # which partitioning may tear apart
    with pytest.raises(TypeError, match="per event"):
        fastjet.ClusterSequence(
            varied, fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
        )
    assert session.node_count() == before


def test_a_second_array_must_be_deferred_too():
    session = graphed.Session(ga.AwkwardBackend())
    cluseq = fastjet.ClusterSequence(
        ga.from_awkward(session, "ev", _events()),
        fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6),
    )
    with pytest.raises(TypeError):
        cluseq.get_parents(_JETS)


def test_a_second_array_from_another_session_is_refused():
    session, other = (
        graphed.Session(ga.AwkwardBackend()),
        graphed.Session(ga.AwkwardBackend()),
    )
    cluseq = fastjet.ClusterSequence(
        ga.from_awkward(session, "ev", _events()),
        fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6),
    )
    recorded = session.node_count()
    with pytest.raises(graphed.GraphedTypeError):
        cluseq.get_parents(ga.from_awkward(other, "jets", _JETS))
    assert session.node_count() == recorded


def test_a_recorded_query_points_at_the_analysts_line():
    session = graphed.Session(ga.AwkwardBackend())
    cluseq = fastjet.ClusterSequence(
        ga.from_awkward(session, "ev", _events()),
        fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6),
    )
    assert session.provenance(cluseq.inclusive_jets()).filename == __file__
