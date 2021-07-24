import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402


class _classgeneralevent:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.data = data
        self.multi_layered_listoffset(self.data)
        self._clusterable_level = ak.Array(
            self._clusterable_level.layout.toListOffsetArray64(True)
        )
        px, py, pz, E, offsets = self.extract_cons(self._clusterable_level)
        px = self.correct_byteorder(px)
        py = self.correct_byteorder(py)
        pz = self.correct_byteorder(pz)
        E = self.correct_byteorder(E)
        offsets = self.correct_byteorder(offsets)
        self._results = fastjet._ext.interfacemulti(px, py, pz, E, offsets, jetdef)

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

    def _check_listoffset(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.ListOffsetArray64,
                ak.layout.ListOffsetArray32,
                ak.layout.ListOffsetArrayU32,
            ),
        )
        return out

    def _check_record(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.RecordArray,
                ak.layout.NumpyArray,
            ),
        )

        return out

    def multi_layered_listoffset(self, data):
        if isinstance(data.layout, ak.layout.VirtualArray):
            self.multi_layered_listoffset(ak.Array(data.layout.array))
        elif self._check_listoffset_subtree(ak.Array(data.layout.content)):
            if self._check_record(
                ak.Array(ak.Array(data.layout.content).layout.content),
            ):
                attributes = dir(data)
                if (
                    "px" in attributes
                    and "py" in attributes
                    and "pz" in attributes
                    and "E" in attributes
                ):
                    self._clusterable_level = ak.Array(data.layout.content)
            else:
                self.multi_layered_listoffset(ak.Array(data.layout.content))
        else:
            self.multi_layered_listoffset(ak.Array(data.layout.content))

    def correct_byteorder(self, data):
        if data.dtype.byteorder == "=":
            pass
        else:
            data = data.dtype.newbyteorder("=")
        return data

    def extract_cons(self, array):
        px = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).px)
        py = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).py)
        pz = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).pz)
        E = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).E)
        off = np.asarray(array.layout.stops)
        off = np.insert(off, 0, 0)
        return px, py, pz, E, off

    def replace(self, layout):
        if isinstance(
            layout,
            (
                ak.layout.ListOffsetArray64,
                ak.layout.ListOffsetArray32,
                ak.layout.ListOffsetArrayU32,
                ak.layout.RegularArray,
            ),
        ) and isinstance(
            layout.content,
            (
                ak.layout.RecordArray,
                ak.layout.NumpyArray,
            ),
        ):
            attributes = dir(ak.Array(layout))
            if (
                "px" in attributes
                and "py" in attributes
                and "pz" in attributes
                and "E" in attributes
            ):
                return self.out.layout
        elif self._check_listoffset(ak.Array(layout)):
            if isinstance(layout, ak.layout.ListOffsetArray64):
                return ak.layout.ListOffsetArray64(
                    layout.offsets,
                    self.replace(layout.content),
                    layout.identities,
                    layout.parameters,
                )
            if isinstance(layout, ak.layout.ListOffsetArray32):
                return ak.layout.ListOffsetArray32(
                    layout.offsets,
                    self.replace(layout.content),
                    layout.identities,
                    layout.parameters,
                )
            if isinstance(layout, ak.layout.ListOffsetArrayU32):
                return ak.layout.ListOffsetArrayU32(
                    layout.offsets,
                    self.replace(layout.content),
                    layout.identities,
                    layout.parameters,
                )
        elif isinstance(layout, ak.layout.ByteMaskedArray):
            return ak.layout.ByteMaskedArray(
                layout.mask,
                self.replace(layout.content),
                layout.valid_when,
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.BitMaskedArray):
            return ak.layout.BitMaskedArray(
                layout.mask,
                self.replace(layout.content),
                layout.valid_when,
                len(layout),
                layout.lsb_order,
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.UnmaskedArray):
            return ak.layout.UnmaskedArray(
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.RegularArray):
            return ak.layout.RegularArray(
                self.replace(layout.content),
                layout.size,
                len(layout),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.ListArray64):
            return ak.layout.ListArray64(
                layout.starts,
                layout.stops,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.ListArray32):
            return ak.layout.ListArray32(
                layout.starts,
                layout.stops,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.ListArrayU32):
            return ak.layout.ListArrayU32(
                layout.starts,
                layout.stops,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.IndexedArray64):
            return ak.layout.IndexedArray64(
                layout.index,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.IndexedArray32):
            return ak.layout.IndexedArray32(
                layout.index,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.IndexedArrayU32):
            return ak.layout.IndexedArrayU32(
                layout.index,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.IndexedOptionArray64):
            return ak.layout.IndexedOptionArray64(
                layout.index,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.IndexedOptionArray32):
            return ak.layout.IndexedOptionArray32(
                layout.index,
                self.replace(layout.content),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.RecordArray):
            return ak.layout.RecordArray(
                [self.replace(x) for x in layout.contents],
                layout.recordlookup,
                len(layout),
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.Record):
            return ak.layout.Record(
                self.replace(layout.array),
                layout.at,
            )
        elif isinstance(layout, ak.layout.UnionArray8_32):
            return ak.layout.UnionArray8_32(
                layout.tags,
                layout.index,
                [self.replace(x) for x in layout.contents],
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.UnionArray8_U32):
            return ak.layout.UnionArray8_U32(
                layout.tags,
                layout.index,
                [self.replace(x) for x in layout.contents],
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.UnionArray8_64):
            return ak.layout.UnionArray8_64(
                layout.tags,
                layout.index,
                [self.replace(x) for x in layout.contents],
                layout.identities,
                layout.parameters,
            )
        elif isinstance(layout, ak.layout.VirtualArray):
            return self.replace(layout.array)

        if isinstance(layout, ak.partition.PartitionedArray):
            return ak.partition.IrregularlyPartitionedArray(
                [self.replace(x) for x in layout.partitions]
            )
        else:
            raise AssertionError(layout)

    def inclusive_jets(self, min_pt):
        np_results = self._results.to_numpy(min_pt)
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        self.out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(of),
                ak.layout.RecordArray(
                    (
                        ak.layout.NumpyArray(np_results[0]),
                        ak.layout.NumpyArray(np_results[1]),
                        ak.layout.NumpyArray(np_results[2]),
                        ak.layout.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                ),
            )
        )
        res = ak.Array(self.replace(self.data.layout))
        return res

    def constituent_index(self, min_pt):
        np_results = self._results.to_numpy_with_constituents(min_pt)
        off = np.insert(np_results[-1], 0, 0)
        out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np_results[0]), ak.layout.NumpyArray(np_results[1])
            )
        )
        self.out = ak.Array(
            ak.layout.ListOffsetArray64(ak.layout.Index64(off), out.layout)
        )
        res = ak.Array(self.replace(self.data.layout))
        return res
