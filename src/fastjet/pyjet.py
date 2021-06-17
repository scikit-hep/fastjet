import awkward as ak

import fastjet._ext  # noqa: F401, E402
import fastjet._multievent
import fastjet._singleevent
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        """The base class for all clustering"""
        self._jetdef = jetdef
        self._jagedness = self._check_jaggedness(data)
        if self._check_listoffset(data):
            if self._jagedness >= 1:
                self._internalrep = fastjet._multievent._classmultievent(
                    data, self._jetdef
                )
        if self._jagedness == 0:
            self._internalrep = fastjet._singleevent._classsingleevent(
                data, self._jetdef
            )

    def _check_jaggedness(self, data):
        """Internal function for checking the jaggedness of awkward array"""
        if isinstance(data.layout, ak.layout.ListOffsetArray64):
            return 1 + self._check_jaggedness(ak.Array(data.layout.content))
        else:
            return 0

    def _check_listoffset(self, data):
        """Internal function for checking whether the given array is a listoffset array or not"""
        out = isinstance(
            data.layout,
            (
                ak.layout.ListOffsetArray64,
                ak.layout.ListOffsetArray32,
                ak.layout.ListOffsetArrayU32,
            ),
        )
        return out

    def jet_def(self):
        """Returns the Jet Definition Object associated with the instance"""
        return self._jetdef

    def inclusive_jets(self, min_pt=0):
        """Returns the inclusive jets after clustering in the same format as the input awkward array"""
        return self._internalrep.inclusive_jets(min_pt)

    def unclustered_particles(self):
        """Returns the unclustered particles after clustering in the same format as the input awkward array"""
        return self._internalrep.unclustered_particles()

    def exclusive_jets(self, n_jets=-1, dcut=-1):
        """Returns the exclusive jets after clustering in the same format as the input awkward array"""
        return self._internalrep.exclusive_jets(n_jets, dcut)

    def exclusive_jets_ycut(self, ycut=-1):
        """Returns the exclusive jets after clustering in the same format as the input awkward array"""
        return self._internalrep.exclusive_jets_ycut(ycut)

    def unclustered_parts(self):
        """Returns the particles that were left unclustered"""
        return self._internalrep.unclustered_parts

    def constituent_index(self, min_pt=0):
        """Returns the index of the constituent of each Jet"""
        return self._internalrep.constituent_index(min_pt)

    def constituents(self, min_pt=0):
        """Returns the particles that make up each Jet"""
        return self._internalrep.constituents(min_pt)

    def exclusive_dmerge(self, njets=10):
        """Returns the dmin corresponding to the recombination that went from n+1 to n jets"""
        return self._internalrep.exclusive_dmerge(njets)

    def exclusive_dmerge_max(self, njets=10):
        """Returns the maximum of the dmin encountered during all recombinations up to the one that led to an n-jet final state"""
        return self._internalrep.exclusive_dmerge_max(njets)

    def exclusive_ymerge_max(self, njets=10):
        """Same as exclusive_dmerge_max, but normalised to squared total energy"""
        return self._internalrep.exclusive_ymerge_max(njets)

    def exclusive_ymerge(self, njets=10):
        """Returns the ymin corresponding to the recombination that went from n+1 to n jets"""
        return self._internalrep.exclusive_ymerge(njets)

    def Q(self):
        """Returns the sum of all energies in the event (relevant mainly for e+e-)"""
        return self._internalrep.Q()

    def Q2(self):
        """Return Q()^2"""
        return self._internalrep.Q2()

    def exclusive_subjets(self, data, dcut=-1, nsub=-1):
        """Returns an Awkward Array of all subjets of the current jet (in the sense of the exclusive algorithm) that would be obtained when running the algorithm with the given dcut."""
        return self._internalrep.exclusive_subjets(data, dcut, nsub)

    def exclusive_subjets_up_to(self, data, nsub=0):
        """Returns the list of subjets obtained by unclustering the supplied jet down to nsub subjets (or all constituents if there are fewer than nsub)."""
        return self._internalrep.exclusive_subjets_up_to(data, nsub)

    def exclusive_subdmerge(self, data, nsub=0):
        """Returns the dij that was present in the merging nsub+1 -> nsub subjets inside this jet."""
        return self._internalrep.exclusive_subdmerge(data, nsub)

    def exclusive_subdmerge_max(self, data, nsub=0):
        """Returns the maximum dij that occurred in the whole event at the stage that the nsub+1 -> nsub merge of subjets occurred inside this jet."""
        return self._internalrep.exclusive_subdmerge_max(data, nsub)

    def n_exclusive_subjets(self, data, dcut=0):
        """Returns the size of exclusive_subjets(...); still n ln n with same coefficient, but marginally more efficient than manually taking len(exclusive_subjets)"""
        return self._internalrep.n_exclusive_subjets(data, dcut)

    def has_parents(self, data):
        """if the jet has parents in the clustering, it returns true"""
        return self._internalrep.has_parents(data)

    def has_child(self, data):
        """if the jet has children in the clustering, it returns true"""
        return self._internalrep.has_child(data)

    def jet_scale_for_algorithm(self, data):
        """Returns the scale associated with a jet as required for this clustering algorithm (kt^2 for the kt-algorithm, 1 for the Cambridge algorithm)."""
        return self._internalrep.jet_scale_for_algorithm(data)
