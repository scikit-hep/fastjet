import hashlib
import json

import awkward as ak
from graphed import Array, Varied, expand, labels, universe
from graphed.awkward import AwkwardForm
from graphed.core import PayloadDescriptor
from graphed.provenance import register_internal

import fastjet._pyjet
from fastjet.__init__ import ClusterSequence
from fastjet._pyjet import _default_taus_njettiness
from fastjet.version import __version__

__all__ = ("__version__",)

# a recorded query's provenance is the analyst's line, not the fastjet frame that recorded it
register_internal("fastjet")


class _FnGraphedInternalRepCaller:
    """Cluster one partition and answer one query. Module level and closure free, so a plan
    carrying it can be pickled to a worker."""

    def __init__(self, method_name, jetdef, **kwargs):
        self.name = method_name
        self.jetdef = jetdef
        self.kwargs = kwargs

    def __call__(self, array, *arrays):
        seq = fastjet._pyjet.AwkwardClusterSequence(array, self.jetdef)
        return getattr(seq, self.name)(*arrays, **self.kwargs)


def _holds_one_event_per_entry(data):
    """The first axis of a deferred array is the partition axis, so every event's particles have
    to sit in a list of their own; the recorded form answers that without reading any data."""
    members = (
        [universe(data, label) for label in labels(data)]
        if isinstance(data, Varied)
        else [data]
    )
    return all(member.session.form(member).tt.ndim > 1 for member in members)


def _length_zero(array):
    return ak.Array(
        array.session.form(array).tt.layout.form.length_zero_array(highlevel=False)
    )


def _query_form(jetdef, method_name, arrays, kwargs):
    """The recorded output form: the real eager class run on length-zero arrays of the operands'
    forms, the same recipe the dask arm uses for its meta."""
    length_zero = [_length_zero(array) for array in arrays]
    seq = fastjet._pyjet.AwkwardClusterSequence(length_zero[0], jetdef)
    out = getattr(seq, method_name)(*length_zero[1:], **kwargs)
    return AwkwardForm(ak.Array(out.layout.to_typetracer(forget_length=True)))


def _payload_descriptor(jetdef):
    """The clustering payload is the jet definition; its description is what the algorithm,
    radius and recombination scheme are, so hashing it gives every query off one definition the
    same payload identity."""
    digest = hashlib.sha256(jetdef.description().encode("utf-8")).hexdigest()
    return PayloadDescriptor(
        kind="fastjet.cluster_sequence",
        content_hash=f"sha256:{digest}",
        framework="fastjet",
        version=__version__,
        io_schema="Momentum4D->Momentum4D",
        preprocessing_ref=None,
    )


def _record(jetdef, method_name, arrays, kwargs):
    # the query and its arguments are the node's params, so two queries off one JetDefinition are
    # distinct nodes while equal calls intern; the evaluator is registered on the session, and
    # graphed.aggregate_plan wires it into the plan it builds
    return arrays[0].session.record_external(
        f"fastjet.{method_name}",
        _FnGraphedInternalRepCaller(method_name, jetdef, **kwargs),
        arrays,
        {"query": method_name, "arguments": json.dumps(kwargs, sort_keys=True)},
        descriptor=_payload_descriptor(jetdef),
        form=_query_form(jetdef, method_name, arrays, kwargs),
    )


def _graphed_dispatch(cluseq, method_name, *arrays, **kwargs):
    for array in arrays:
        if not isinstance(array, (Array, Varied)):
            raise TypeError("The input data is not a graphed Array")

    def record(*operands):
        return _record(cluseq._jetdef, method_name, list(operands), kwargs)

    # a varied operand records the query once per universe, as graphed.apply does
    return expand(record, (cluseq._data, *arrays), {})


class GraphedClusterSequence(ClusterSequence):
    def __init__(self, data, jetdef):
        if not isinstance(data, (Array, Varied)):
            raise TypeError("The input data is not a graphed Array")
        if not isinstance(jetdef, fastjet._swig.JetDefinition):
            raise TypeError("JetDefinition is not of valid type")
        if not _holds_one_event_per_entry(data):
            raise TypeError(
                "The input must hold one list of particles per event; cluster a single "
                "event eagerly with an awkward array"
            )
        self._jetdef = jetdef
        self._data = data

    def jet_def(self):
        return self._jetdef

    def inclusive_jets(self, min_pt=0):
        return _graphed_dispatch(self, "inclusive_jets", min_pt=min_pt)

    def unclustered_particles(self):
        return _graphed_dispatch(self, "unclustered_particles")

    def exclusive_jets(self, n_jets=-1, dcut=-1):
        return _graphed_dispatch(self, "exclusive_jets", n_jets=n_jets, dcut=dcut)

    def exclusive_jets_up_to(self, n_jets=-1):
        return _graphed_dispatch(self, "exclusive_jets_up_to", n_jets=n_jets)

    def exclusive_jets_ycut(self, ycut=-1):
        return _graphed_dispatch(self, "exclusive_jets_ycut", ycut=ycut)

    def constituent_index(self, min_pt=0):
        return _graphed_dispatch(self, "constituent_index", min_pt=min_pt)

    def constituents(self, min_pt=0):
        return _graphed_dispatch(self, "constituents", min_pt=min_pt)

    def exclusive_jets_constituent_index(self, njets=10):
        return _graphed_dispatch(self, "exclusive_jets_constituent_index", njets=njets)

    def exclusive_jets_constituents(self, njets=10):
        return _graphed_dispatch(self, "exclusive_jets_constituents", njets=njets)

    def exclusive_jets_softdrop_grooming(
        self,
        njets=1,
        beta=0.0,
        symmetry_cut=0.1,
        symmetry_measure="scalar_z",
        R0=0.8,
        recursion_choice="larger_pt",
        # subtractor = 0,
        mu_cut=float("inf"),
    ):
        return _graphed_dispatch(
            self,
            "exclusive_jets_softdrop_grooming",
            njets=njets,
            beta=beta,
            symmetry_cut=symmetry_cut,
            symmetry_measure=symmetry_measure,
            R0=R0,
            recursion_choice=recursion_choice,
            # subtractor=subtractor,
            mu_cut=mu_cut,
        )

    def njettiness(
        self,
        measure_definition="NormalizedMeasure",
        axes_definition="OnePass_KT_Axes",
        njets=_default_taus_njettiness,
        beta=1.0,
        R0=0.8,
        Rcutoff=None,
        nPass=None,
        akAxesR0=None,
    ):
        return _graphed_dispatch(
            self,
            "njettiness",
            measure_definition=measure_definition,
            axes_definition=axes_definition,
            njets=njets,
            beta=beta,
            R0=R0,
            Rcutoff=Rcutoff,
            nPass=nPass,
            akAxesR0=akAxesR0,
        )

    def exclusive_jets_energy_correlator(
        self,
        njets=1,
        beta=1,
        npoint=0,
        angles=-1,
        alpha=0,
        func="generalized",
        normalized=True,
    ):
        return _graphed_dispatch(
            self,
            "exclusive_jets_energy_correlator",
            njets=njets,
            beta=beta,
            npoint=npoint,
            angles=angles,
            alpha=alpha,
            func=func,
            normalized=normalized,
        )

    def exclusive_jets_lund_declusterings(self, njets=10):
        return _graphed_dispatch(self, "exclusive_jets_lund_declusterings", njets=njets)

    def exclusive_dmerge(self, njets=10):
        return _graphed_dispatch(self, "exclusive_dmerge", njets=njets)

    def exclusive_dmerge_max(self, njets=10):
        return _graphed_dispatch(self, "exclusive_dmerge_max", njets=njets)

    def exclusive_ymerge_max(self, njets=10):
        return _graphed_dispatch(self, "exclusive_ymerge_max", njets=njets)

    def exclusive_ymerge(self, njets=10):
        return _graphed_dispatch(self, "exclusive_ymerge", njets=njets)

    def Q(self):
        return _graphed_dispatch(self, "Q")

    def Q2(self):
        return _graphed_dispatch(self, "Q2")

    def exclusive_subjets(self, data, dcut=-1, nsub=-1):
        return _graphed_dispatch(self, "exclusive_subjets", data, dcut=dcut, nsub=nsub)

    def exclusive_subjets_up_to(self, data, nsub=0):
        return _graphed_dispatch(self, "exclusive_subjets_up_to", data, nsub=nsub)

    def exclusive_subdmerge(self, data, nsub=0):
        return _graphed_dispatch(self, "exclusive_subdmerge", data, nsub=nsub)

    def exclusive_subdmerge_max(self, data, nsub=0):
        return _graphed_dispatch(self, "exclusive_subdmerge_max", data, nsub=nsub)

    def n_exclusive_subjets(self, data, dcut=0):
        return _graphed_dispatch(self, "n_exclusive_subjets", data, dcut=dcut)

    def has_parents(self, data):
        return _graphed_dispatch(self, "has_parents", data)

    def has_child(self, data):
        return _graphed_dispatch(self, "has_child", data)

    def jet_scale_for_algorithm(self, data):
        return _graphed_dispatch(self, "jet_scale_for_algorithm", data)

    def unique_history_order(self):
        return _graphed_dispatch(self, "unique_history_order")

    def n_particles(self):
        return _graphed_dispatch(self, "n_particles")

    def n_exclusive_jets(self, dcut=0):
        return _graphed_dispatch(self, "n_exclusive_jets", dcut=dcut)

    def childless_pseudojets(self):
        return _graphed_dispatch(self, "childless_pseudojets")

    def jets(self):
        return _graphed_dispatch(self, "jets")

    def get_parents(self, data):
        return _graphed_dispatch(self, "get_parents", data)

    def get_child(self, data):
        return _graphed_dispatch(self, "get_child", data)
