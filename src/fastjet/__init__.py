# BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

from typing import Union

import awkward as ak

import fastjet._ext  # noqa: F401, E402
import fastjet._pyjet  # noqa: F401, E402
import fastjet._swig  # noqa: F401, E402
from fastjet._swig import AreaDefinition  # noqa: F401, E402
from fastjet._swig import BackgroundEstimatorBase  # noqa: F401, E402
from fastjet._swig import BackgroundJetPtDensity  # noqa: F401, E402
from fastjet._swig import BackgroundJetPtMDensity  # noqa: F401, E402
from fastjet._swig import BackgroundJetScalarPtDensity  # noqa: F401, E402
from fastjet._swig import BackgroundRescalingYPolynomial  # noqa: F401, E402
from fastjet._swig import Best  # noqa: F401, E402
from fastjet._swig import BestFJ30  # noqa: F401, E402
from fastjet._swig import BIpt2_scheme  # noqa: F401, E402
from fastjet._swig import BIpt_scheme  # noqa: F401, E402
from fastjet._swig import Boost  # noqa: F401, E402
from fastjet._swig import CASubJetTagger  # noqa: F401, E402
from fastjet._swig import CASubJetTaggerStructure  # noqa: F401, E402
from fastjet._swig import ClusterSequence1GhostPassiveArea  # noqa: F401, E402
from fastjet._swig import ClusterSequence_fastjet_banner_stream  # noqa: F401, E402
from fastjet._swig import ClusterSequence_print_banner  # noqa: F401, E402
from fastjet._swig import ClusterSequence_set_fastjet_banner_stream  # noqa: F401, E402
from fastjet._swig import ClusterSequenceActiveArea  # noqa: F401, E402
from fastjet._swig import ClusterSequenceActiveAreaExplicitGhosts  # noqa: F401, E402
from fastjet._swig import ClusterSequenceArea  # noqa: F401, E402
from fastjet._swig import ClusterSequenceAreaBase  # noqa: F401, E402
from fastjet._swig import ClusterSequencePassiveArea  # noqa: F401, E402
from fastjet._swig import ClusterSequenceStructure  # noqa: F401, E402
from fastjet._swig import ClusterSequenceVoronoiArea  # noqa: F401, E402
from fastjet._swig import CompositeJetStructure  # noqa: F401, E402
from fastjet._swig import E_scheme  # noqa: F401, E402
from fastjet._swig import Error  # noqa: F401, E402
from fastjet._swig import Error_set_default_stream  # noqa: F401, E402
from fastjet._swig import Error_set_print_backtrace  # noqa: F401, E402
from fastjet._swig import Error_set_print_errors  # noqa: F401, E402
from fastjet._swig import Et2_scheme  # noqa: F401, E402
from fastjet._swig import Et_scheme  # noqa: F401, E402
from fastjet._swig import Filter  # noqa: F401, E402
from fastjet._swig import FilterStructure  # noqa: F401, E402
from fastjet._swig import FunctionOfPseudoJetDouble  # noqa: F401, E402
from fastjet._swig import FunctionOfPseudoJetPseudoJet  # noqa: F401, E402
from fastjet._swig import GhostedAreaSpec  # noqa: F401, E402
from fastjet._swig import GridMedianBackgroundEstimator  # noqa: F401, E402
from fastjet._swig import IndexedSortHelper  # noqa: F401, E402
from fastjet._swig import InternalError  # noqa: F401, E402
from fastjet._swig import JetDefinition  # noqa: F401, E402
from fastjet._swig import JetDefinition0Param  # noqa: F401, E402
from fastjet._swig import JetDefinition1Param  # noqa: F401, E402
from fastjet._swig import JetDefinition2Param  # noqa: F401, E402
from fastjet._swig import JetDefinition_algorithm_description  # noqa: F401, E402
from fastjet._swig import JetDefinition_n_parameters_for_algorithm  # noqa: F401, E402
from fastjet._swig import JetMedianBackgroundEstimator  # noqa: F401, E402
from fastjet._swig import JHTopTagger  # noqa: F401, E402
from fastjet._swig import JHTopTaggerStructure  # noqa: F401, E402
from fastjet._swig import LimitedWarning  # noqa: F401, E402
from fastjet._swig import LimitedWarning_set_default_max_warn  # noqa: F401, E402
from fastjet._swig import LimitedWarning_set_default_stream  # noqa: F401, E402
from fastjet._swig import LimitedWarning_summary  # noqa: F401, E402
from fastjet._swig import MassDropTagger  # noqa: F401, E402
from fastjet._swig import MassDropTaggerStructure  # noqa: F401, E402
from fastjet._swig import MaxRap  # noqa: F401, E402
from fastjet._swig import N2MHTLazy9  # noqa: F401, E402
from fastjet._swig import N2MHTLazy9Alt  # noqa: F401, E402
from fastjet._swig import N2MHTLazy9AntiKtSeparateGhosts  # noqa: F401, E402
from fastjet._swig import N2MHTLazy25  # noqa: F401, E402
from fastjet._swig import N2MinHeapTiled  # noqa: F401, E402
from fastjet._swig import N2Plain  # noqa: F401, E402
from fastjet._swig import N2PoorTiled  # noqa: F401, E402
from fastjet._swig import N2Tiled  # noqa: F401, E402
from fastjet._swig import N3Dumb  # noqa: F401, E402
from fastjet._swig import NlnN  # noqa: F401, E402
from fastjet._swig import NlnN3pi  # noqa: F401, E402
from fastjet._swig import NlnN4pi  # noqa: F401, E402
from fastjet._swig import NlnNCam  # noqa: F401, E402
from fastjet._swig import NlnNCam2pi2R  # noqa: F401, E402
from fastjet._swig import NlnNCam4pi  # noqa: F401, E402
from fastjet._swig import Pruner  # noqa: F401, E402
from fastjet._swig import PrunerStructure  # noqa: F401, E402
from fastjet._swig import PruningPlugin  # noqa: F401, E402
from fastjet._swig import PruningRecombiner  # noqa: F401, E402
from fastjet._swig import PseudoJet  # noqa: F401, E402
from fastjet._swig import PseudoJetStructureBase  # noqa: F401, E402
from fastjet._swig import PtYPhiM  # noqa: F401, E402
from fastjet._swig import RangeDefinition  # noqa: F401, E402
from fastjet._swig import Recluster  # noqa: F401, E402
from fastjet._swig import RecombinerPython  # noqa: F401, E402
from fastjet._swig import RectangularGrid  # noqa: F401, E402
from fastjet._swig import RestFrameNSubjettinessTagger  # noqa: F401, E402
from fastjet._swig import RestFrameNSubjettinessTaggerStructure  # noqa: F401, E402
from fastjet._swig import Selector  # noqa: F401, E402
from fastjet._swig import SelectorAbsEtaMax  # noqa: F401, E402
from fastjet._swig import SelectorAbsEtaMin  # noqa: F401, E402
from fastjet._swig import SelectorAbsEtaRange  # noqa: F401, E402
from fastjet._swig import SelectorAbsRapMax  # noqa: F401, E402
from fastjet._swig import SelectorAbsRapMin  # noqa: F401, E402
from fastjet._swig import SelectorAbsRapRange  # noqa: F401, E402
from fastjet._swig import SelectorCircle  # noqa: F401, E402
from fastjet._swig import SelectorDoughnut  # noqa: F401, E402
from fastjet._swig import SelectorEMax  # noqa: F401, E402
from fastjet._swig import SelectorEMin  # noqa: F401, E402
from fastjet._swig import SelectorERange  # noqa: F401, E402
from fastjet._swig import SelectorEtaMax  # noqa: F401, E402
from fastjet._swig import SelectorEtaMin  # noqa: F401, E402
from fastjet._swig import SelectorEtaRange  # noqa: F401, E402
from fastjet._swig import SelectorEtMax  # noqa: F401, E402
from fastjet._swig import SelectorEtMin  # noqa: F401, E402
from fastjet._swig import SelectorEtRange  # noqa: F401, E402
from fastjet._swig import SelectorIdentity  # noqa: F401, E402
from fastjet._swig import SelectorIsPureGhost  # noqa: F401, E402
from fastjet._swig import SelectorIsZero  # noqa: F401, E402
from fastjet._swig import SelectorMassMax  # noqa: F401, E402
from fastjet._swig import SelectorMassMin  # noqa: F401, E402
from fastjet._swig import SelectorMassRange  # noqa: F401, E402
from fastjet._swig import SelectorNHardest  # noqa: F401, E402
from fastjet._swig import SelectorPhiRange  # noqa: F401, E402
from fastjet._swig import SelectorPtFractionMin  # noqa: F401, E402
from fastjet._swig import SelectorPtMax  # noqa: F401, E402
from fastjet._swig import SelectorPtMin  # noqa: F401, E402
from fastjet._swig import SelectorPtRange  # noqa: F401, E402
from fastjet._swig import SelectorPython  # noqa: F401, E402
from fastjet._swig import SelectorRapMax  # noqa: F401, E402
from fastjet._swig import SelectorRapMin  # noqa: F401, E402
from fastjet._swig import SelectorRapPhiRange  # noqa: F401, E402
from fastjet._swig import SelectorRapRange  # noqa: F401, E402
from fastjet._swig import SelectorRectangle  # noqa: F401, E402
from fastjet._swig import SelectorStrip  # noqa: F401, E402
from fastjet._swig import SelectorWorker  # noqa: F401, E402
from fastjet._swig import SelectorWorkerPython  # noqa: F401, E402
from fastjet._swig import Subtractor  # noqa: F401, E402
from fastjet._swig import SwigPyIterator  # noqa: F401, E402
from fastjet._swig import TilingBase  # noqa: F401, E402
from fastjet._swig import TopTaggerBase  # noqa: F401, E402
from fastjet._swig import TopTaggerBaseStructure  # noqa: F401, E402
from fastjet._swig import Transformer  # noqa: F401, E402
from fastjet._swig import Unboost  # noqa: F401, E402
from fastjet._swig import UserInfoPython  # noqa: F401, E402
from fastjet._swig import VoronoiAreaSpec  # noqa: F401, E402
from fastjet._swig import WTA_modp_scheme  # noqa: F401, E402
from fastjet._swig import WTA_pt_scheme  # noqa: F401, E402
from fastjet._swig import active_area  # noqa: F401, E402
from fastjet._swig import active_area_explicit_ghosts  # noqa: F401, E402
from fastjet._swig import antikt_algorithm  # noqa: F401, E402
from fastjet._swig import cambridge_aachen_algorithm  # noqa: F401, E402
from fastjet._swig import cambridge_algorithm  # noqa: F401, E402
from fastjet._swig import cambridge_for_passive_algorithm  # noqa: F401, E402
from fastjet._swig import cpp_string_from_name_py_obj  # noqa: F401, E402
from fastjet._swig import cpp_string_from_py_str  # noqa: F401, E402
from fastjet._swig import cpp_string_from_str_py_obj  # noqa: F401, E402
from fastjet._swig import cvar  # noqa: F401, E402
from fastjet._swig import def_ghost_area  # noqa: F401, E402
from fastjet._swig import def_ghost_maxrap  # noqa: F401, E402
from fastjet._swig import def_grid_scatter  # noqa: F401, E402
from fastjet._swig import def_mean_ghost_pt  # noqa: F401, E402
from fastjet._swig import def_pt_scatter  # noqa: F401, E402
from fastjet._swig import def_repeat  # noqa: F401, E402
from fastjet._swig import ee_genkt_algorithm  # noqa: F401, E402
from fastjet._swig import ee_kt_algorithm  # noqa: F401, E402
from fastjet._swig import eulergamma  # noqa: F401, E402
from fastjet._swig import external_scheme  # noqa: F401, E402
from fastjet._swig import fastjet_version_string  # noqa: F401, E402
from fastjet._swig import genkt_algorithm  # noqa: F401, E402
from fastjet._swig import genkt_for_passive_algorithm  # noqa: F401, E402
from fastjet._swig import invalid_area  # noqa: F401, E402
from fastjet._swig import kt_algorithm  # noqa: F401, E402
from fastjet._swig import ln2  # noqa: F401, E402
from fastjet._swig import one_ghost_passive_area  # noqa: F401, E402
from fastjet._swig import passive_area  # noqa: F401, E402
from fastjet._swig import pi  # noqa: F401, E402
from fastjet._swig import pisq  # noqa: F401, E402
from fastjet._swig import plugin_algorithm  # noqa: F401, E402
from fastjet._swig import plugin_strategy  # noqa: F401, E402
from fastjet._swig import pseudojet_invalid_phi  # noqa: F401, E402
from fastjet._swig import pseudojet_invalid_rap  # noqa: F401, E402
from fastjet._swig import pt2_scheme  # noqa: F401, E402
from fastjet._swig import pt_scheme  # noqa: F401, E402
from fastjet._swig import twopi  # noqa: F401, E402
from fastjet._swig import undefined_jet_algorithm  # noqa: F401, E402
from fastjet._swig import vectorPJ  # noqa: F401, E402
from fastjet._swig import voronoi_area  # noqa: F401, E402
from fastjet._swig import zeta2  # noqa: F401, E402
from fastjet._swig import zeta3  # noqa: F401, E402
from fastjet._utils import cos_theta  # noqa: F401, E402
from fastjet._utils import dot_product  # noqa: F401, E402
from fastjet._utils import have_same_momentum  # noqa: F401, E402
from fastjet._utils import join  # noqa: F401, E402
from fastjet._utils import sort_indices  # noqa: F401, E402
from fastjet._utils import sorted_by_E  # noqa: F401, E402
from fastjet._utils import sorted_by_pt  # noqa: F401, E402
from fastjet._utils import sorted_by_pz  # noqa: F401, E402
from fastjet._utils import sorted_by_rapidity  # noqa: F401, E402
from fastjet._utils import theta  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

# TODO: everything should be in this list. Except maybe __version__.
__all__ = ("__version__",)


class ClusterSequence:  # The super class
    """The base class for all clustering.

    Args:
        data(awkward.highlevel.Array): The data for clustering.
        jetdef(fastjet._swig.JetDefinition): The JetDefinition for clustering specification.
    """

    def __init__(self, data, jetdef):
        if not isinstance(jetdef, fastjet._swig.JetDefinition):
            raise AttributeError("JetDefinition is not correct") from None
        if isinstance(data, ak.Array):
            self.__class__ = fastjet._pyjet.AwkwardClusterSequence
            fastjet._pyjet.AwkwardClusterSequence.__init__(
                self, data=data, jetdef=jetdef
            )
        if isinstance(data, list):
            self.__class__ = fastjet._swig.ClusterSequence
            fastjet._swig.ClusterSequence.__init__(self, data, jetdef)

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

    def constituents(self, min_p: float = 0) -> ak.Array:
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

    def exclusive_jets_lund_declusterings(self, njets: int = 10) -> ak.Array:
        """Returns the Lund declustering Delta and k_T parameters from exclusive n_jets.

        Args:
            njets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """

        raise AssertionError()

    def exclusive_dmerge(self, njets: int = 10) -> Union[ak.Array, float]:
        """Returns the dmin corresponding to the recombination that went from n+1 to n jets.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_dmerge_max(self, njets: int = 10) -> Union[ak.Array, float]:
        """Returns the maximum of the dmin encountered during all recombinations up to the one that led to an n-jet final state.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_ymerge_max(self, njets: int = 10) -> Union[ak.Array, float]:
        """Same as exclusive_dmerge_max, but normalised to squared total energy.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def exclusive_ymerge(self, njets: int = 10) -> Union[ak.Array, float]:
        """Returns the ymin corresponding to the recombination that went from n+1 to n jets.

        Args:
            n_jets (int): The number of jets it was clustered to.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def Q(self) -> Union[ak.Array, float]:
        """Returns the sum of all energies in the event (relevant mainly for e+e-)

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def Q2(self) -> Union[ak.Array, float]:
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

    def exclusive_subdmerge(
        self, data: ak.Array, nsub: int = 0
    ) -> Union[ak.Array, float]:
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
    ) -> Union[ak.Array, float]:
        """Returns the maximum dij that occurred in the whole event at the stage that the nsub+1 -> nsub merge of subjets occurred inside this jet.

        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def n_exclusive_subjets(
        self, data: ak.Array, dcut: float = 0
    ) -> Union[ak.Array, int]:
        """Returns the size of exclusive_subjets(...); still n ln n with same coefficient, but marginally more efficient than manually taking len(exclusive_subjets)

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
            dcut (float): The dcut for the result.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def has_parents(self, data: ak.Array) -> Union[ak.Array, bool]:
        """if the jet has parents in the clustering, it returns true.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def has_child(self, data: ak.Array) -> Union[ak.Array, bool]:
        """If the jet has children in the clustering, it returns true.

        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def jet_scale_for_algorithm(self, data: ak.Array) -> Union[ak.Array, float]:
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

    def n_particles(self) -> Union[ak.Array, int]:
        """Returns the number of particles that were provided to the clustering algorithm.

        Args:
            None

        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        raise AssertionError()

    def n_exclusive_jets(self, dcut: float = 0) -> Union[ak.Array, int]:
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
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input."""
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
