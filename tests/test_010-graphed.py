import pickle
import subprocess
import sys

import awkward as ak
import pytest

import fastjet
import fastjet._pyjet

graphed = pytest.importorskip("graphed")
ga = pytest.importorskip("graphed.awkward")
execution = pytest.importorskip("graphed.core.execution")
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


def _jetdef():
    return fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)


def _session():
    return graphed.Session(ga.AwkwardBackend(behavior=vector.backends.awkward.behavior))


def _recorded_type(array):
    return str(array.session.form_of(array.node_id).tt.type)


def _eager_type(array):
    return str(ak.Array(array.layout.to_typetracer(forget_length=True)).type)


def _count_clusterings(monkeypatch):
    """Every arm reaches the compiled clustering through AwkwardClusterSequence."""
    calls = []
    original = fastjet._pyjet.AwkwardClusterSequence.__init__

    def counted(self, *args, **kwargs):
        jetdef = kwargs["jetdef"] if "jetdef" in kwargs else args[1]
        calls.append(jetdef.description())
        return original(self, *args, **kwargs)

    monkeypatch.setattr(fastjet._pyjet.AwkwardClusterSequence, "__init__", counted)
    return calls


def _parquet_source(tmp_path, chunks):
    paths = []
    for i, chunk in enumerate(chunks):
        path = tmp_path / f"part{i}.parquet"
        ak.to_parquet(ak.Array(chunk), path)
        paths.append(str(path))
    session = _session()
    # from_parquet drops the record name, so the momentum behaviour is re-declared here
    raw = ga.from_parquet(session, "pf", paths)
    particles = ga.gak.zip(
        {"px": raw.px, "py": raw.py, "pz": raw.pz, "E": raw.E},
        with_name="Momentum4D",
    )
    return session, particles


def _counts(values):
    return (int(ak.sum(values[0])), int(ak.sum(ak.flatten(values[1]))))


def _add_counts(left, right):
    return (left[0] + right[0], left[1] + right[1])


def _no_counts():
    return (0, 0)


QUERIES = [
    ("inclusive_jets", {}),
    ("inclusive_jets", {"min_pt": 10.0}),
    ("constituents", {}),
    ("constituent_index", {}),
    ("exclusive_jets", {"n_jets": 1}),
    ("exclusive_jets", {"dcut": 0.0001}),
    ("unclustered_particles", {}),
]


def test_dispatch_accepts_a_graphed_array():
    session = _session()
    jetdef = _jetdef()
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", _events()), jetdef)
    assert type(cluseq) is not fastjet.ClusterSequence
    assert cluseq.jet_def() is jetdef


def test_public_method_set_matches_the_awkward_class():
    def public_callables(cls):
        return {
            name
            for name, value in vars(cls).items()
            if not name.startswith("_") and callable(value)
        }

    session = _session()
    graphed_class = type(
        fastjet.ClusterSequence(ga.from_awkward(session, "ev", _events()), _jetdef())
    )
    assert public_callables(graphed_class) == public_callables(
        fastjet._pyjet.AwkwardClusterSequence
    )


@pytest.mark.parametrize(("query", "kwargs"), QUERIES)
def test_query_form_and_value_match_eager(query, kwargs):
    array = _events()
    jetdef = _jetdef()
    session = _session()
    graphed_out = getattr(
        fastjet.ClusterSequence(ga.from_awkward(session, "ev", array), jetdef), query
    )(**kwargs)
    eager_out = getattr(fastjet.ClusterSequence(array, jetdef), query)(**kwargs)

    assert _recorded_type(graphed_out) == _eager_type(eager_out)
    assert session.materialize(graphed_out).to_list() == eager_out.to_list()


def test_min_pt_is_part_of_the_query():
    array = _events()
    jetdef = _jetdef()
    session = _session()
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", array), jetdef)
    all_jets = session.materialize(cluseq.inclusive_jets())
    hard_jets = session.materialize(cluseq.inclusive_jets(min_pt=10.0))
    assert ak.num(all_jets, axis=1).to_list() == [2, 2, 3, 1]
    assert ak.num(hard_jets, axis=1).to_list() == [1, 0, 2, 0]


def test_a_flat_particle_collection_is_refused():
    flat = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56},
        ],
        with_name="Momentum4D",
    )
    jetdef = _jetdef()
    session = _session()
    source = ga.from_awkward(session, "ev", flat)
    before = session.node_count()

    # the first axis of a deferred array is the partition axis, so a flat collection
    # is one event that partitioning may silently tear apart
    with pytest.raises(TypeError, match="per event"):
        fastjet.ClusterSequence(source, jetdef)
    assert session.node_count() == before

    # the refusal is about deferring, not about the particles: they cluster eagerly
    eager_jets = fastjet.ClusterSequence(flat, jetdef).inclusive_jets()
    assert len(eager_jets) == 2
    assert ak.sum(eager_jets.E) == pytest.approx(ak.sum(flat.E))

    # and one list of particles per event is still accepted
    jagged = fastjet.ClusterSequence(ga.from_awkward(session, "evs", _events()), jetdef)
    assert (
        session.materialize(jagged.inclusive_jets()).to_list()
        == fastjet.ClusterSequence(_events(), jetdef).inclusive_jets().to_list()
    )


def test_two_array_queries_take_a_second_graphed_array():
    array = _events()
    jets = ak.Array(
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5},
            {"px": -7.1, "py": 2.9, "pz": 11.0, "E": 14.0},
            {"px": 0.3, "py": -0.2, "pz": 9.9, "E": 10.0},
            {"px": 5.0, "py": 5.0, "pz": 5.0, "E": 9.0},
        ]
    )
    jetdef = _jetdef()
    session = _session()
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", array), jetdef)
    graphed_jets = ga.from_awkward(session, "jets", jets)
    eager = fastjet.ClusterSequence(array, jetdef)

    parents = cluseq.get_parents(graphed_jets)
    eager_parents = eager.get_parents(jets)
    assert _recorded_type(parents) == _eager_type(eager_parents)
    assert session.materialize(parents).to_list() == eager_parents.to_list()

    subjets = cluseq.exclusive_subjets(graphed_jets, nsub=1)
    eager_subjets = eager.exclusive_subjets(jets, nsub=1)
    assert _recorded_type(subjets) == _eager_type(eager_subjets)
    assert session.materialize(subjets).to_list() == eager_subjets.to_list()


def test_one_jetdef_records_one_node_per_distinct_query():
    session = _session()
    cluseq = fastjet.ClusterSequence(
        ga.from_awkward(session, "ev", _events()), _jetdef()
    )
    jets = cluseq.inclusive_jets()
    hard_jets = cluseq.inclusive_jets(min_pt=10.0)
    index = cluseq.constituent_index()
    assert len({jets.node_id, hard_jets.node_id, index.node_id}) == 3
    # equal arguments hash-cons to the same External node
    assert cluseq.inclusive_jets().node_id == jets.node_id
    assert cluseq.inclusive_jets(min_pt=10.0).node_id == hard_jets.node_id


def test_two_queries_in_one_plan_cluster_once_per_partition(tmp_path, monkeypatch):
    chunks = (_events()[:2].to_list(), _events()[2:].to_list())
    _, particles = _parquet_source(tmp_path, chunks)
    jetdef = _jetdef()
    cluseq = fastjet.ClusterSequence(particles, jetdef)
    n_jets = ga.gak.num(cluseq.inclusive_jets(), axis=1)
    n_constituents = ga.gak.num(cluseq.constituent_index(), axis=2)

    plan = graphed.aggregate_plan(
        n_jets,
        n_constituents,
        reduce=_counts,
        combine=_add_counts,
        empty=_no_counts,
    )
    calls = _count_clusterings(monkeypatch)
    result = execution.SequentialRunner().run(plan)

    assert result.value == (8, 10)
    # dask parity: one clustering per query per partition, all off the one JetDefinition
    assert calls == [jetdef.description()] * 4


def test_partitioned_plan_matches_eager_and_survives_a_pickle(tmp_path):
    array = _events()
    chunks = (array[:2].to_list(), array[2:].to_list())
    _, particles = _parquet_source(tmp_path, chunks)
    cluseq = fastjet.ClusterSequence(particles, _jetdef())
    plan = graphed.aggregate_plan(
        ga.gak.num(cluseq.inclusive_jets(), axis=1),
        ga.gak.num(cluseq.constituent_index(), axis=2),
        reduce=_counts,
        combine=_add_counts,
        empty=_no_counts,
    )

    eager = fastjet.ClusterSequence(array, _jetdef())
    expected = (
        int(ak.sum(ak.num(eager.inclusive_jets(), axis=1))),
        int(ak.sum(ak.flatten(ak.num(eager.constituent_index(), axis=2)))),
    )

    assert execution.SequentialRunner().run(plan).value == expected
    shipped = pickle.loads(pickle.dumps(plan))
    assert execution.SequentialRunner().run(shipped).value == expected


def test_sliced_and_reordered_input_clusters_its_own_events():
    array = _events()
    jetdef = _jetdef()
    whole = fastjet.ClusterSequence(array, jetdef).inclusive_jets().to_list()

    sliced = array[1:]
    assert sliced.layout.starts[0] != 0
    session = _session()
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", sliced), jetdef)
    assert session.materialize(cluseq.inclusive_jets()).to_list() == whole[1:]

    order = [2, 0, 3, 1]
    picked = array[order]
    session = _session()
    cluseq = fastjet.ClusterSequence(ga.from_awkward(session, "ev", picked), jetdef)
    assert session.materialize(cluseq.inclusive_jets()).to_list() == [
        whole[i] for i in order
    ]


def test_projection_and_record_time_typing_use_the_declared_form():
    session = _session()
    source = ga.from_awkward(session, "ev", _events())
    cluseq = fastjet.ClusterSequence(source, _jetdef())
    index = cluseq.constituent_index()

    # the declared form is one level deeper than the input, so axis=2 type-checks here
    n_constituents = ga.gak.num(index, axis=2)
    assert _recorded_type(n_constituents) == "## * var * int64"
    with pytest.raises(graphed.GraphedTypeError):
        ga.gak.num(cluseq.inclusive_jets(), axis=2)

    projection = ga.project_buffers(n_constituents, on_fail="pass")
    assert projection.read_buffers == {
        "ev": {
            "px": graphed.BufferNeed.DATA,
            "py": graphed.BufferNeed.DATA,
            "pz": graphed.BufferNeed.DATA,
            "E": graphed.BufferNeed.DATA,
        }
    }


def test_varied_input_fans_the_clustering_over_its_members():
    array = _events()
    jetdef = _jetdef()
    session = _session()
    source = ga.from_awkward(session, "ev", array)

    def scaled(factor):
        return ga.gak.zip(
            {
                "px": source.px * factor,
                "py": source.py * factor,
                "pz": source.pz * factor,
                "E": source.E * factor,
            },
            with_name="Momentum4D",
        )

    varied = graphed.vary(
        source, "scale", points={"up": scaled(1.05), "down": scaled(0.95)}
    )
    jets = fastjet.ClusterSequence(varied, jetdef).inclusive_jets()

    assert graphed.labels(jets) == ("nominal", "scale_up", "scale_down")
    for label, factor in (("nominal", 1.0), ("scale_up", 1.05), ("scale_down", 0.95)):
        member = ak.zip(
            {
                "px": array.px * factor,
                "py": array.py * factor,
                "pz": array.pz * factor,
                "E": array.E * factor,
            },
            with_name="Momentum4D",
        )
        expected = fastjet.ClusterSequence(member, jetdef).inclusive_jets()
        got = session.materialize(graphed.universe(jets, label))
        assert got.to_list() == expected.to_list()


def test_importing_fastjet_does_not_import_graphed():
    code = "import fastjet, sys; print('graphed' in sys.modules)"
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    assert out.stdout.splitlines()[-1] == "False"


def test_other_dispatch_branches_are_unchanged():
    jetdef = _jetdef()
    array = _events()
    assert isinstance(
        fastjet.ClusterSequence(array, jetdef), fastjet._pyjet.AwkwardClusterSequence
    )
    assert isinstance(
        fastjet.ClusterSequence(
            [fastjet.PseudoJet(1.2, 3.2, 5.4, 2.5)],
            jetdef,
        ),
        fastjet._swig.ClusterSequence,
    )
    dak = pytest.importorskip("dask_awkward")
    assert isinstance(
        fastjet.ClusterSequence(dak.from_awkward(array, 2), jetdef),
        fastjet._pyjet.DaskAwkwardClusterSequence,
    )
    with pytest.raises(TypeError):
        fastjet.ClusterSequence("not an array", jetdef)
