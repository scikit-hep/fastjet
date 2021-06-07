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
