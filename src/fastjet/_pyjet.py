import awkward as ak

import fastjet._ext  # noqa: F401, E402
import fastjet._generalevent
import fastjet._multievent
import fastjet._singleevent
from fastjet.__init__ import ClusterSequence
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence(ClusterSequence):
    def __init__(self, data, jetdef):
        """The base class for all clustering
        Attributes:
                _jetdef (fastjet._swig.JetDefinition): The jetdefinition that was passed by the user.
                _jaggedness (int): The maximum depth of the Awkward Array (stored for internal use).
                _internalrep (fastjet._[type]event): The internal class which performs the Jet clustering, changes depending on the input type.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        if not isinstance(jetdef, fastjet._swig.JetDefinition):
            raise TypeError("JetDefinition is not of valid type")
        self._jetdef = jetdef
        self._jagedness = self._check_jaggedness(data)
        if self._check_listoffset(data) and self._jagedness == 2:
            self._internalrep = fastjet._multievent._classmultievent(data, self._jetdef)
        if self._jagedness == 1 and isinstance(data.layout, ak.layout.RecordArray):
            self._internalrep = fastjet._singleevent._classsingleevent(
                data, self._jetdef
            )
        if (
            self._jagedness >= 3
            or self._check_general(data)
            or isinstance(data.layout, ak.partition.IrregularlyPartitionedArray)
        ):
            self._internalrep = fastjet._generalevent._classgeneralevent(data, jetdef)

    # else:
    # raise TypeError(
    #    "This kind of Awkward Array is not supported yet. Please contact the maintainers for further action."
    # )

    def _check_jaggedness(self, data):
        """Internal function for checking the jaggedness of awkward array.
        Args:
            data (awkward.highlevel.Array): The input Awkward Array.
        Returns:
            int: The maximum depth of the awkward array.
        """
        if self._check_general_jaggedness(data) or self._check_listoffset(data):
            return 1 + self._check_jaggedness(ak.Array(data.layout.content))
        if isinstance(
            data.layout,
            (
                ak.layout.UnionArray8_32,
                ak.layout.UnionArray8_U32,
                ak.layout.UnionArray8_64,
            ),
        ):
            return 1 + max(
                [self._check_jaggedness(ak.Array(x)) for x in data.layout.contents]
            )
        if isinstance(
            data.layout,
            (ak.layout.RecordArray,),
        ):
            return 1 + max(
                [self._check_jaggedness(ak.Array(x)) for x in data.layout.contents]
            )
        if isinstance(
            data.layout,
            (ak.partition.IrregularlyPartitionedArray),
        ):
            return 1 + max(
                [self._check_jaggedness(ak.Array(x)) for x in data.layout.partitions]
            )
        if isinstance(data.layout, ak.layout.VirtualArray):
            return 1 + self._check_jaggedness(ak.Array(data.layout.array))
        else:
            return 0

    def _check_general(self, data):
        """Internal function for checking whether the given array is a general Awkward Array or not.
        Args:
            data (awkward.highlevel.Array): The input Awkward Array.
        Returns:
            bool: True if it is a general case, False otherwise.
        """
        out = isinstance(
            data.layout,
            (
                ak.layout.IndexedArray64,
                ak.layout.IndexedArray32,
                ak.layout.IndexedArrayU32,
                ak.layout.ByteMaskedArray,
                ak.layout.BitMaskedArray,
                ak.layout.UnmaskedArray,
                ak.layout.IndexedOptionArray64,
                ak.layout.IndexedOptionArray32,
                ak.layout.VirtualArray,
                ak.layout.UnionArray8_32,
                ak.layout.UnionArray8_U32,
                ak.layout.UnionArray8_64,
                ak.layout.Record,
            ),
        )
        return out

    def _check_general_jaggedness(self, data):
        """Internal function for checking whether the given array is a general case or not for depth calculation.
        Args:
                data (awkward.highlevel.Array): The input Awkward Array.
        Returns:
            bool: Returns True if the awkward array is a general case, False otherwise.
        """
        out = isinstance(
            data.layout,
            (
                ak.layout.IndexedArray64,
                ak.layout.IndexedArray32,
                ak.layout.IndexedArrayU32,
                ak.layout.ByteMaskedArray,
                ak.layout.BitMaskedArray,
                ak.layout.UnmaskedArray,
                ak.layout.IndexedOptionArray64,
                ak.layout.IndexedOptionArray32,
                ak.layout.Record,
            ),
        )
        return out

    def _check_listoffset(self, data):
        """Internal function for checking whether the given array is a listoffset array or not
        Args:
            data (awkward.highlevel.Array): The input Awkward Array.
        Returns:
            Bool: Returns True if the awkward array is a listoffset case, False otherwise.
        """
        out = isinstance(
            data.layout,
            (
                ak.layout.ListOffsetArray64,
                ak.layout.ListOffsetArray32,
                ak.layout.ListOffsetArrayU32,
                ak.layout.ListArray64,
                ak.layout.ListArray32,
                ak.layout.ListArrayU32,
                ak.layout.RegularArray,
            ),
        )
        return out

    def jet_def(self):
        """Returns the Jet Definition Object associated with the instance
        Args:
            None
        Returns:
            JetDefinition: Returns the jetdefinition stored as an attribute.
        """
        return self._jetdef

    def inclusive_jets(self, min_pt=0):
        """Returns the inclusive jets after clustering in the same format as the input awkward array
        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input containting inclusive jets.
        """
        return self._internalrep.inclusive_jets(min_pt)

    def unclustered_particles(self):
        """Returns the unclustered particles after clustering in the same format as the input awkward array
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input contating the unclustered particles.
        """
        return self._internalrep.unclustered_particles()

    def exclusive_jets(self, n_jets=-1, dcut=-1):
        """Returns the exclusive jets after clustering in the same format as the input awkward array. Either takes njets or dcut as argument.
        Args:
            n_jets (int): The number of jets it was clustered to.
            dcut (float): The dcut for the result.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_jets(n_jets, dcut)

    def exclusive_jets_ycut(self, ycut=-1):
        """Returns the exclusive jets after clustering in the same format as the input awkward array.
        Args:
            ycut (float): The dcut for the result.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_jets_ycut(ycut)

    def constituent_index(self, min_pt=0):
        """Returns the index of the constituent of each Jet.
        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.constituent_index(min_pt)

    def constituents(self, min_pt=0):
        """Returns the particles that make up each Jet.
        Args:
            min_pt (float): The minimum value of the pt for the inclusive jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.constituents(min_pt)

    def exclusive_dmerge(self, njets=10):
        """Returns the dmin corresponding to the recombination that went from n+1 to n jets.
        Args:
            n_jets (int): The number of jets it was clustered to.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_dmerge(njets)

    def exclusive_dmerge_max(self, njets=10):
        """Returns the maximum of the dmin encountered during all recombinations up to the one that led to an n-jet final state.
        Args:
            n_jets (int): The number of jets it was clustered to.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_dmerge_max(njets)

    def exclusive_ymerge_max(self, njets=10):
        """Same as exclusive_dmerge_max, but normalised to squared total energy.
        Args:
            n_jets (int): The number of jets it was clustered to.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_ymerge_max(njets)

    def exclusive_ymerge(self, njets=10):
        """Returns the ymin corresponding to the recombination that went from n+1 to n jets.
        Args:
            n_jets (int): The number of jets it was clustered to.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.exclusive_ymerge(njets)

    def Q(self):
        """Returns the sum of all energies in the event (relevant mainly for e+e-)
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.Q()

    def Q2(self):
        """Return Q()^2
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.Q2()

    def exclusive_subjets(self, data, dcut=-1, nsub=-1):
        """Returns an Awkward Array of all subjets of the current jet (in the sense of the exclusive algorithm) that would be obtained when running the algorithm with the given dcut.
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
            dcut (float): The dcut for the result.
            n_sub (int): The number of subjets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subjets(data, dcut, nsub)

    def exclusive_subjets_up_to(self, data, nsub=0):
        """Returns the list of subjets obtained by unclustering the supplied jet down to nsub subjets (or all constituents if there are fewer than nsub).
        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subjets_up_to(data, nsub)

    def exclusive_subdmerge(self, data, nsub=0):
        """Returns the dij that was present in the merging nsub+1 -> nsub subjets inside this jet.
        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subdmerge(data, nsub)

    def exclusive_subdmerge_max(self, data, nsub=0):
        """Returns the maximum dij that occurred in the whole event at the stage that the nsub+1 -> nsub merge of subjets occurred inside this jet.
        Args:
            data (awkward.highlevel.Array): An Awkward Array containing the Jets.
            n_sub (int): The number of subjets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subdmerge_max(data, nsub)

    def n_exclusive_subjets(self, data, dcut=0):
        """Returns the size of exclusive_subjets(...); still n ln n with same coefficient, but marginally more efficient than manually taking len(exclusive_subjets)
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
            dcut (float): The dcut for the result.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.n_exclusive_subjets(data, dcut)

    def has_parents(self, data):
        """if the jet has parents in the clustering, it returns true.
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.has_parents(data)

    def has_child(self, data):
        """If the jet has children in the clustering, it returns true.
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.has_child(data)

    def jet_scale_for_algorithm(self, data):
        """Returns the scale associated with a jet as required for this clustering algorithm (kt^2 for the kt-algorithm, 1 for the Cambridge algorithm).
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.jet_scale_for_algorithm(data)

    def unique_history_order(self):
        """Routine that returns an order in which to read the history such that clusterings that lead to identical jet compositions but different histories (because of degeneracies in the clustering order) will have matching constituents for each matching entry in the unique_history_order.
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.unique_history_order()

    def n_particles(self):
        """Returns the number of particles that were provided to the clustering algorithm.
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.n_particles()

    def n_exclusive_jets(self, dcut=0):
        """Returns the number of jets (in the sense of the exclusive algorithm) that would be obtained when running the algorithm with the given dcut.
        Args:
            dcut (float): The dcut for the result.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.n_exclusive_jets(dcut)

    def childless_pseudojets(self):
        """Return the list of pseudojets in the ClusterSequence that do not have children (and are not among the inclusive jets).
        Args:
            None
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.childless_pseudojets()

    def jets(self):
        """Allows the user to access the internally stored _jets() array, which contains both the initial particles and the various intermediate and final stages of recombination.
        Args:
            none
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        return self._internalrep.jets()

    def get_parents(self, data):
        """If the jet has parents in the clustering, it returns them.
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input."""
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.get_parents(data)

    def get_child(self, data):
        """If the jet has parents in the clustering, it returns them.
        Args:
            data (awkward.highlevel.Array): An Array containing the Jets.
        Returns:
            awkward.highlevel.Array: Returns an Awkward Array of the same type as the input.
        """
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.get_child(data)
