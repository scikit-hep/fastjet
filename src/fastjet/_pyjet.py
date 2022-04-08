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
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array or Numpy Array")
        if not isinstance(jetdef, fastjet._swig.JetDefinition):
            raise TypeError("JetDefinition is not of valid type")
        self._jetdef = jetdef
        self._jagedness = self._check_jaggedness(data)
        self._flag = 1
        if (self._check_listoffset(data) and self._jagedness == 2) or (
            self._check_listoffset_index(data)
        ):
            self._flag = 0
            self._internalrep = fastjet._multievent._classmultievent(data, self._jetdef)
        elif self._jagedness == 1 and isinstance(data.layout, ak.layout.RecordArray):
            self._internalrep = fastjet._singleevent._classsingleevent(
                data, self._jetdef
            )
        elif (
            self._jagedness >= 3
            or self._check_general(data)
            or isinstance(data.layout, ak.partition.IrregularlyPartitionedArray)
            and self._flag == 1
        ):
            self._internalrep = fastjet._generalevent._classgeneralevent(data, jetdef)

    # else:
    # raise TypeError(
    #    "This kind of Awkward Array is not supported yet. Please contact the maintainers for further action."
    # )

    def _check_jaggedness(self, data):
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
                self._check_jaggedness(ak.Array(x)) for x in data.layout.contents
            )
        if isinstance(
            data.layout,
            (ak.layout.RecordArray,),
        ):
            return 1 + max(
                self._check_jaggedness(ak.Array(x)) for x in data.layout.contents
            )
        if isinstance(
            data.layout,
            (ak.partition.IrregularlyPartitionedArray),
        ):
            return 1 + max(
                self._check_jaggedness(ak.Array(x)) for x in data.layout.partitions
            )
        if isinstance(data.layout, ak.layout.VirtualArray):
            return 1 + self._check_jaggedness(ak.Array(data.layout.array))
        else:
            return 0

    def _check_listoffset_index(self, data):
        if self._check_listoffset_subtree(ak.Array(data.layout)):
            if self._check_record(
                ak.Array(ak.Array(data.layout.content)),
            ):
                return True
            elif self._check_indexed(
                ak.Array(data.layout.content),
            ):
                if self._check_record(
                    ak.Array(ak.Array(data.layout.content).layout.content)
                ):
                    return True
            else:
                return False
        else:
            return False

    def _check_record(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.RecordArray,
                ak.layout.NumpyArray,
            ),
        )

        return out

    def _check_indexed(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.IndexedArray64,
                ak.layout.IndexedArray32,
                ak.layout.IndexedArray32,
                ak.layout.IndexedOptionArray64,
                ak.layout.IndexedOptionArray32,
            ),
        )

        return out

    def _check_listoffset_subtree(self, data):
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

    def _check_general(self, data):
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
        return self._jetdef

    def inclusive_jets(self, min_pt=0):
        return self._internalrep.inclusive_jets(min_pt)

    def unclustered_particles(self):
        return self._internalrep.unclustered_particles()

    def exclusive_jets(self, n_jets=-1, dcut=-1):
        return self._internalrep.exclusive_jets(n_jets, dcut)

    def exclusive_jets_ycut(self, ycut=-1):
        return self._internalrep.exclusive_jets_ycut(ycut)

    def constituent_index(self, min_pt=0):
        return self._internalrep.constituent_index(min_pt)

    def constituents(self, min_pt=0):
        return self._internalrep.constituents(min_pt)

    def exclusive_dmerge(self, njets=10):
        return self._internalrep.exclusive_dmerge(njets)

    def exclusive_dmerge_max(self, njets=10):
        return self._internalrep.exclusive_dmerge_max(njets)

    def exclusive_ymerge_max(self, njets=10):
        return self._internalrep.exclusive_ymerge_max(njets)

    def exclusive_ymerge(self, njets=10):
        return self._internalrep.exclusive_ymerge(njets)

    def Q(self):
        return self._internalrep.Q()

    def Q2(self):
        return self._internalrep.Q2()

    def exclusive_subjets(self, data, dcut=-1, nsub=-1):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subjets(data, dcut, nsub)

    def exclusive_subjets_up_to(self, data, nsub=0):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subjets_up_to(data, nsub)

    def exclusive_subdmerge(self, data, nsub=0):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subdmerge(data, nsub)

    def exclusive_subdmerge_max(self, data, nsub=0):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.exclusive_subdmerge_max(data, nsub)

    def n_exclusive_subjets(self, data, dcut=0):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.n_exclusive_subjets(data, dcut)

    def has_parents(self, data):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.has_parents(data)

    def has_child(self, data):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.has_child(data)

    def jet_scale_for_algorithm(self, data):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.jet_scale_for_algorithm(data)

    def unique_history_order(self):
        return self._internalrep.unique_history_order()

    def n_particles(self):
        return self._internalrep.n_particles()

    def n_exclusive_jets(self, dcut=0):
        return self._internalrep.n_exclusive_jets(dcut)

    def childless_pseudojets(self):
        return self._internalrep.childless_pseudojets()

    def jets(self):
        return self._internalrep.jets()

    def get_parents(self, data):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.get_parents(data)

    def get_child(self, data):
        if not isinstance(data, ak.Array):
            raise TypeError("The input data is not an Awkward Array")
        return self._internalrep.get_child(data)
