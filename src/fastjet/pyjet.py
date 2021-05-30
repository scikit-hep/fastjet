import awkward as ak

import fastjet._ext  # noqa: F401, E402
import fastjet._multievent
import fastjet._singleevent
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        """The base class for all clustering"""
        self.jagedness = self.check_jaggedness(data)
        if self.check_listoffset(data):
            if self.jagedness >= 1:
                self._internalrep = fastjet._multievent._classmultievent(data, jetdef)
        if self.jagedness == 0:
            self._internalrep = fastjet._singleevent._classsingleevent(data, jetdef)

    def check_jaggedness(self, data):
        """Internal function for checking the jaggedness of awkward array"""
        if isinstance(data.layout, ak.layout.ListOffsetArray64):
            return 1 + self.check_jaggedness(ak.Array(data.layout.content))
        else:
            return 0

    def check_listoffset(self, data):
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

    @property
    def inclusive_jets(self):
        """Returns the inclusive jets after clustering in the same format as the input awkward array"""
        return self._internalrep.inclusive_jets

    @property
    def unclustered_parts(self):
        """Returns the particles that were left unclustered"""
        return self._internalrep.unclustered_parts

    @property
    def constituent_index(self):
        """Returns the index of the constituent of each Jet"""
        return self._internalrep.constituent_index

    @property
    def constituents(self):
        """Returns the particles that make up each Jet"""
        return self._internalrep.constituents
