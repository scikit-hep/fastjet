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

    def _check_general(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.IndexedArray64,
                ak.layout.IndexedArray32,
                ak.layout.IndexedArrayU32,
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
                ak.layout.RegularArray,
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
        if (
            self._check_record(
                ak.Array(ak.Array(data.layout.content).layout.content),
            )
            and self._check_listoffset(ak.Array(data.layout.content))
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
        elif self._check_general(ak.Array(layout)) or self._check_listoffset(
            ak.Array(layout)
        ):
            return ak.layout.ListOffsetArray64(
                layout.offsets,
                self.replace(layout.content),
            )
        elif isinstance(layout, ak.layout.ByteMaskedArray):
            return ak.layout.ByteMaskedArray(
                layout.bytemask(), self.replace(layout.content), not layout.valid_when
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
