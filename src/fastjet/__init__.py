# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

import os as _os

if _os.name == "nt":
    module_dir = _os.path.dirname(_os.path.abspath(__file__))

    dll_dir1 = _os.path.abspath(_os.path.join(module_dir, "bin"))

    if _os.path.exists(dll_dir1):
        _os.add_dll_directory(dll_dir1)


import awkward as ak

import fastjet._ext
import fastjet._pyjet
import fastjet._swig
from fastjet._swig import (
    AreaDefinition,  # noqa: F401
    BackgroundEstimatorBase,  # noqa: F401
    BackgroundJetPtDensity,  # noqa: F401
    BackgroundJetPtMDensity,  # noqa: F401
    BackgroundJetScalarPtDensity,  # noqa: F401
    BackgroundRescalingYPolynomial,  # noqa: F401
    Best,  # noqa: F401
    BestFJ30,  # noqa: F401
    BIpt2_scheme,  # noqa: F401
    BIpt_scheme,  # noqa: F401
    Boost,  # noqa: F401
    CASubJetTagger,  # noqa: F401
    CASubJetTaggerStructure,  # noqa: F401
    ClusterSequence1GhostPassiveArea,  # noqa: F401
    ClusterSequenceActiveArea,  # noqa: F401
    ClusterSequenceActiveAreaExplicitGhosts,  # noqa: F401
    ClusterSequenceArea,  # noqa: F401
    ClusterSequenceAreaBase,  # noqa: F401
    ClusterSequencePassiveArea,  # noqa: F401
    ClusterSequenceStructure,  # noqa: F401
    ClusterSequenceVoronoiArea,  # noqa: F401
    CompositeJetStructure,  # noqa: F401
    E_scheme,  # noqa: F401
    Error,  # noqa: F401
    Et2_scheme,  # noqa: F401
    Et_scheme,  # noqa: F401
    Filter,  # noqa: F401
    FilterStructure,  # noqa: F401
    FunctionOfPseudoJetDouble,  # noqa: F401
    FunctionOfPseudoJetPseudoJet,  # noqa: F401
    GhostedAreaSpec,  # noqa: F401
    GridMedianBackgroundEstimator,  # noqa: F401
    IndexedSortHelper,  # noqa: F401
    InternalError,  # noqa: F401
    JetDefinition0Param,  # noqa: F401
    JetDefinition1Param,  # noqa: F401
    JetDefinition2Param,  # noqa: F401
    JetMedianBackgroundEstimator,  # noqa: F401
    JHTopTagger,  # noqa: F401
    JHTopTaggerStructure,  # noqa: F401
    LimitedWarning,  # noqa: F401
    MassDropTagger,  # noqa: F401
    MassDropTaggerStructure,  # noqa: F401
    MaxRap,  # noqa: F401
    N2MHTLazy9,  # noqa: F401
    N2MHTLazy9Alt,  # noqa: F401
    N2MHTLazy9AntiKtSeparateGhosts,  # noqa: F401
    N2MHTLazy25,  # noqa: F401
    N2MinHeapTiled,  # noqa: F401
    N2Plain,  # noqa: F401
    N2PoorTiled,  # noqa: F401
    N2Tiled,  # noqa: F401
    N3Dumb,  # noqa: F401
    NlnN,  # noqa: F401
    NlnN3pi,  # noqa: F401
    NlnN4pi,  # noqa: F401
    NlnNCam,  # noqa: F401
    NlnNCam2pi2R,  # noqa: F401
    NlnNCam4pi,  # noqa: F401
    Pruner,  # noqa: F401
    PrunerStructure,  # noqa: F401
    PruningPlugin,  # noqa: F401
    PruningRecombiner,  # noqa: F401
    PseudoJet,  # noqa: F401
    PseudoJetStructureBase,  # noqa: F401
    PtYPhiM,  # noqa: F401
    RangeDefinition,  # noqa: F401
    Recluster,  # noqa: F401
    RecombinerPython,  # noqa: F401
    RectangularGrid,  # noqa: F401
    RestFrameNSubjettinessTagger,  # noqa: F401
    RestFrameNSubjettinessTaggerStructure,  # noqa: F401
    Selector,  # noqa: F401
    SelectorAbsEtaMax,  # noqa: F401
    SelectorAbsEtaMin,  # noqa: F401
    SelectorAbsEtaRange,  # noqa: F401
    SelectorAbsRapMax,  # noqa: F401
    SelectorAbsRapMin,  # noqa: F401
    SelectorAbsRapRange,  # noqa: F401
    SelectorCircle,  # noqa: F401
    SelectorDoughnut,  # noqa: F401
    SelectorEMax,  # noqa: F401
    SelectorEMin,  # noqa: F401
    SelectorERange,  # noqa: F401
    SelectorEtaMax,  # noqa: F401
    SelectorEtaMin,  # noqa: F401
    SelectorEtaRange,  # noqa: F401
    SelectorEtMax,  # noqa: F401
    SelectorEtMin,  # noqa: F401
    SelectorEtRange,  # noqa: F401
    SelectorIdentity,  # noqa: F401
    SelectorIsPureGhost,  # noqa: F401
    SelectorIsZero,  # noqa: F401
    SelectorMassMax,  # noqa: F401
    SelectorMassMin,  # noqa: F401
    SelectorMassRange,  # noqa: F401
    SelectorNHardest,  # noqa: F401
    SelectorPhiRange,  # noqa: F401
    SelectorPtFractionMin,  # noqa: F401
    SelectorPtMax,  # noqa: F401
    SelectorPtMin,  # noqa: F401
    SelectorPtRange,  # noqa: F401
    SelectorPython,  # noqa: F401
    SelectorRapMax,  # noqa: F401
    SelectorRapMin,  # noqa: F401
    SelectorRapPhiRange,  # noqa: F401
    SelectorRapRange,  # noqa: F401
    SelectorRectangle,  # noqa: F401
    SelectorStrip,  # noqa: F401
    SelectorWorker,  # noqa: F401
    SelectorWorkerPython,  # noqa: F401
    Subtractor,  # noqa: F401
    SwigPyIterator,  # noqa: F401
    TilingBase,  # noqa: F401
    TopTaggerBase,  # noqa: F401
    TopTaggerBaseStructure,  # noqa: F401
    Transformer,  # noqa: F401
    Unboost,  # noqa: F401
    UserInfoPython,  # noqa: F401
    VoronoiAreaSpec,  # noqa: F401
    WTA_modp_scheme,  # noqa: F401
    WTA_pt_scheme,  # noqa: F401
    active_area,  # noqa: F401
    active_area_explicit_ghosts,  # noqa: F401
    antikt_algorithm,  # noqa: F401
    cambridge_aachen_algorithm,  # noqa: F401
    cambridge_algorithm,  # noqa: F401
    cambridge_for_passive_algorithm,  # noqa: F401
    cpp_string_from_name_py_obj,  # noqa: F401
    cpp_string_from_py_str,  # noqa: F401
    cpp_string_from_str_py_obj,  # noqa: F401
    cvar,  # noqa: F401
    def_ghost_area,  # noqa: F401
    def_ghost_maxrap,  # noqa: F401
    def_grid_scatter,  # noqa: F401
    def_mean_ghost_pt,  # noqa: F401
    def_pt_scatter,  # noqa: F401
    def_repeat,  # noqa: F401
    ee_genkt_algorithm,  # noqa: F401
    ee_kt_algorithm,  # noqa: F401
    eulergamma,  # noqa: F401
    external_scheme,  # noqa: F401
    fastjet_version_string,  # noqa: F401
    genkt_algorithm,  # noqa: F401
    genkt_for_passive_algorithm,  # noqa: F401
    invalid_area,  # noqa: F401
    kt_algorithm,  # noqa: F401
    ln2,  # noqa: F401
    one_ghost_passive_area,  # noqa: F401
    passive_area,  # noqa: F401
    pi,  # noqa: F401
    pisq,  # noqa: F401
    plugin_algorithm,  # noqa: F401
    plugin_strategy,  # noqa: F401
    pseudojet_invalid_phi,  # noqa: F401
    pseudojet_invalid_rap,  # noqa: F401
    pt2_scheme,  # noqa: F401
    pt_scheme,  # noqa: F401
    twopi,  # noqa: F401
    undefined_jet_algorithm,  # noqa: F401
    vectorPJ,  # noqa: F401
    voronoi_area,  # noqa: F401
    zeta2,  # noqa: F401
    zeta3,  # noqa: F401
)
from fastjet._swig import JetDefinition as JetDefinitionNoCast
from fastjet._utils import (
    cos_theta,  # noqa: F401
    dot_product,  # noqa: F401
    have_same_momentum,  # noqa: F401
    join,  # noqa: F401
    sort_indices,  # noqa: F401
    sorted_by_E,  # noqa: F401
    sorted_by_pt,  # noqa: F401
    sorted_by_pz,  # noqa: F401
    sorted_by_rapidity,  # noqa: F401
    theta,  # noqa: F401
)
from fastjet.version import __version__

# TODO: everything should be in this list. Except maybe __version__.
__all__ = ("__version__",)


class JetDefinition(JetDefinitionNoCast):
    def __init__(self, *args, **kwargs):
        r"""

        `JetDefinition(JetAlgorithm jet_algorithm_in, double R_in, RecombinationScheme
            recomb_scheme_in, Strategy strategy_in, int nparameters_in)`

        constructor to fully specify a jet-definition (together with information about
        how algorithically to run it).

        """

        R_in = kwargs.pop("R_in", None)
        as_kwargs = False
        if R_in is None:
            R_in = args[1]
        else:
            as_kwargs = True

        if not isinstance(R_in, (float, int)):
            raise TypeError(
                f"R_in should be a real number, got {R_in} of type {type(R_in)}"
            )

        if isinstance(R_in, int):
            R_in = float(R_in)

        new_args = args
        new_kwargs = kwargs
        if as_kwargs:
            new_kwargs = kwargs.copy()
            new_kwargs["R_in"] = R_in
        else:
            new_args = (args[0], R_in, *args[2:])

        self.args = new_args
        self.kwargs = new_kwargs
        super().__init__(*new_args, **kwargs)

    def __setstate__(self, state):
        self.__init__(*state["args"], **state["kwargs"])

    def __getstate__(self):
        return {"args": self.args, "kwargs": self.kwargs}


class ClusterSequence:  # The super class
    """The base class for all clustering.

    Args:
        data(awkward.highlevel.Array): The data for clustering.
        jetdef(fastjet._swig.JetDefinition): The JetDefinition for clustering specification.
    """

    def __init__(self, data, jetdef):
        if not isinstance(jetdef, fastjet._swig.JetDefinition):
            raise TypeError("JetDefinition is not correct") from None
        if isinstance(data, ak.Array):
            self.__class__ = fastjet._pyjet.AwkwardClusterSequence
            fastjet._pyjet.AwkwardClusterSequence.__init__(
                self, data=data, jetdef=jetdef
            )
        elif isinstance(data, list):
            self.__class__ = fastjet._swig.ClusterSequence
            fastjet._swig.ClusterSequence.__init__(self, data, jetdef)
        else:
            try:
                import dask_awkward as dak
            except ImportError:
                dak = None
            if dak is not None and isinstance(data, dak.Array):
                self.__class__ = fastjet._pyjet.DaskAwkwardClusterSequence
                fastjet._pyjet.DaskAwkwardClusterSequence.__init__(
                    self, data=data, jetdef=jetdef
                )
            else:
                raise TypeError(
                    f"{data} must be an awkward.Array, dask_awkward.Array, or list!"
                )

    def jet_def(self) -> JetDefinition:
        """Returns the Jet Definition Object associated with the instance

        Args:
            None

        Returns:
            JetDefinition: Returns the jetdefinition stored as an attribute.
        """
        raise AssertionError()

    def inclusive_jets(self, min_pt: float = 0) -> ak.Array:
        """Returns the inclusive jets after clustering in the same format as the input awkward array

        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input containting inclusive jets.
        """
        raise AssertionError()

    def unclustered_particles(self) -> ak.Array:
        """Returns the unclustered particles after clustering in the same format as the input awkward array

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input contating the unclustered particles.
        """
        raise AssertionError()

    def exclusive_jets(self, n_jets: int = -1, dcut: float = -1) -> ak.Array:
        """Returns the exclusive jets after clustering in the same format as the input awkward array. Either takes njets or dcut as argument.

        Args:
            n_jets (int): The number of jets it was clustered to.
            dcut (float): The dcut for the result.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_jets_ycut(self, ycut: float = -1) -> ak.Array:
        """Returns the exclusive jets after clustering in the same format as the input awkward array.

        Args:
            ycut (float): The dcut for the result.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def constituent_index(self, min_pt: float = 0) -> ak.Array:
        """Returns the index of the constituent of each Jet.

        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_jets_constituents_index(self, njets: int = 10) -> ak.Array:
        """Returns the index of the constituent of each exclusive jet.

        Args:
            njets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """

        raise AssertionError()

    def constituents(self, min_pt: float = 0) -> ak.Array:
        """Returns the particles that make up each Jet.

        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_jets_constituents(self, njets: int = 10) -> ak.Array:
        """Returns the particles that make up each exclusive jet.

        Args:
            njets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """

        raise AssertionError()

    def exclusive_jets_softdrop_grooming(
        self,
        njets: int = 1,
        beta: int = 0,
        symmetry_cut: float = 0.1,
        symmetry_measure="scalar_z",
        R0: float = 0.8,
        recursion_choice="larger_pt",
        # subtractor: int = 0,
        mu_cut: float = float("inf"),
    ) -> ak.Array:
        """Performs softdrop pruning on jets.
        Args:
          n_jets: number of exclusive subjets.
          beta: softdrop beta parameter.
          symmetry_cut: softdrop symmetry cut value.
          symmetry_measure: Which symmetry measure to use, found in RecursiveSymmetryCutBase.hh
          R0: softdrop R0 parameter.
          recursion_choice: Which recursion choice to use, found in RecursiveSymmetryCutBase.hh
          subtractor: an optional pointer to a pileup subtractor (ignored if zero)
        Returns:
          Returns an array of values from the jet after it has been groomed by softdrop.
        """
        raise AssertionError()

    def exclusive_jets_energy_correlator(
        self,
        njets: int = 1,
        beta: int = 1,
        npoint: int = 0,
        angles: int = -1,
        alpha=0,
        func="generalized",
        normalized=True,
    ) -> ak.Array:
        """Returns the energy correlator of each exclusive jet.

        Args:
            njets (int): The number of jets it was clustered to.
            n_point (int): The number of points in the correlator.
            angle: The number of angles to be used in the correlator (if angle != n_point, ECFG is used).
            beta: The beta value for the correlator.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_jets_lund_declusterings(self, njets: int = 10) -> ak.Array:
        """Returns the Lund declustering Delta and k_T parameters from exclusive n_jets.

        Args:
            njets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """

        raise AssertionError()

    def exclusive_dmerge(self, njets: int = 10) -> ak.Array | float:
        """Returns the dmin corresponding to the recombination that went from n+1 to n jets.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_dmerge_max(self, njets: int = 10) -> ak.Array | float:
        """Returns the maximum of the dmin encountered during all recombinations up to the one that led to an n-jet final state.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_ymerge_max(self, njets: int = 10) -> ak.Array | float:
        """Same as exclusive_dmerge_max, but normalised to squared total energy.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_ymerge(self, njets: int = 10) -> ak.Array | float:
        """Returns the ymin corresponding to the recombination that went from n+1 to n jets.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def Q(self) -> ak.Array | float:
        """Returns the sum of all energies in the event (relevant mainly for e+e-)

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def Q2(self) -> ak.Array | float:
        """Return Q()^2

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_subjets(
        self, data: ak.Array, dcut: float = -1, nsub: int = -1
    ) -> ak.Array:
        """Returns an Awkward Array of all subjets of the current jet (in the sense of the exclusive algorithm) that would be obtained when running the algorithm with the given dcut.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
            dcut (float): The dcut for the result.
            n_sub (int): The number of subjets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_subjets_up_to(self, data: ak.Array, nsub: int = 0) -> ak.Array:
        """Returns the list of subjets obtained by unclustering the supplied jet down to nsub subjets (or all constituents if there are fewer than nsub).

        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_subdmerge(self, data: ak.Array, nsub: int = 0) -> ak.Array | float:
        """Returns the dij that was present in the merging nsub+1 -> nsub subjets inside this jet.

        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_subdmerge_max(
        self, data: ak.Array, nsub: int = 0
    ) -> ak.Array | float:
        """Returns the maximum dij that occurred in the whole event at the stage that the nsub+1 -> nsub merge of subjets occurred inside this jet.

        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def n_exclusive_subjets(self, data: ak.Array, dcut: float = 0) -> ak.Array | int:
        """Returns the size of exclusive_subjets(...); still n ln n with same coefficient, but marginally more efficient than manually taking len(exclusive_subjets)

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
            dcut (float): The dcut for the result.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def has_parents(self, data: ak.Array) -> ak.Array | bool:
        """if the jet has parents in the clustering, it returns true.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def has_child(self, data: ak.Array) -> ak.Array | bool:
        """If the jet has children in the clustering, it returns true.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def jet_scale_for_algorithm(self, data: ak.Array) -> ak.Array | float:
        """Returns the scale associated with a jet as required for this clustering algorithm (kt^2 for the kt-algorithm, 1 for the Cambridge algorithm).

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def unique_history_order(self) -> ak.Array:
        """Routine that returns an order in which to read the history such that clusterings that lead to identical jet compositions but different histories (because of degeneracies in the clustering order) will have matching constituents for each matching entry in the unique_history_order.

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def n_particles(self) -> ak.Array | int:
        """Returns the number of particles that were provided to the clustering algorithm.

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def n_exclusive_jets(self, dcut: float = 0) -> ak.Array | int:
        """Returns the number of jets (in the sense of the exclusive algorithm) that would be obtained when running the algorithm with the given dcut.

        Args:
            dcut (float): The dcut for the result.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def childless_pseudojets(self) -> ak.Array:
        """Return the list of pseudojets in the ClusterSequence that do not have children (and are not among the inclusive jets).

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def jets(self) -> ak.Array:
        """Allows the user to access the internally stored _jets() array, which contains both the initial particles and the various intermediate and final stages of recombination.

        Args:
            none

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def get_parents(self, data: ak.Array) -> ak.Array:
        """If the jet has parents in the clustering, it returns them.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def get_child(self, data: ak.Array) -> ak.Array:
        """If the jet has parents in the clustering, it returns them.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()


class multi_inheritor(
    fastjet._swig.ClusterSequence, ClusterSequence
):  # class that inherits both the custom ClusterSequence and swig ClusterSequence and acts as a trampoline
    def __init__(self):
        pass


def formatwarning(message, category, filename, lineno, line=None):
    """Make warnings resemble the ones from fastjet-core"""
    return f"{category.__name__}: {message}\n"
