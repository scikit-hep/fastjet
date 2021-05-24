import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.data = data
        self.jagedness = self.check_jaggedness(self.data)
        # inps, inpf = self.swig_to_params(self.jetdef)
        # offsets = data["part0-node0-offsets"]
        if self.jagedness == 0:
            data = self.single_to_jagged(data)
        px, py, pz, E, offsets = self.extract_cons(data)
        px = self.correct_byteorder(px)
        py = self.correct_byteorder(py)
        pz = self.correct_byteorder(pz)
        E = self.correct_byteorder(E)
        offsets = self.correct_byteorder(offsets)
        self._results = fastjet._ext.interfacemulti(px, py, pz, E, offsets, jetdef)

    def correct_byteorder(self, data):
        if data.dtype.byteorder == "=":
            pass
        else:
            data = data.dtype.newbyteorder("=")
        return data

    def check_jaggedness(self, data):
        if isinstance(data.layout, ak.layout.ListOffsetArray64):
            return 1 + self.check_jaggedness(ak.Array(data.layout.content))
        else:
            return 0

    def extract_cons(self, array):
        px = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).px)
        py = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).py)
        pz = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).pz)
        E = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).E)
        off = np.asarray(array.layout.stops)
        off = np.insert(off, 0, 0)
        return px, py, pz, E, off

    def single_to_jagged(self, array):
        single = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np.array([0, len(array)])),
                ak.layout.RecordArray(
                    [
                        ak.layout.NumpyArray(array.px),
                        ak.layout.NumpyArray(array.py),
                        ak.layout.NumpyArray(array.pz),
                        ak.layout.NumpyArray(array.E),
                    ],
                    ["px", "py", "pz", "E"],
                ),
            )
        )
        return single

    @property
    def inclusive_jets(self):
        np_results = self._results.to_numpy()
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        out = ak.Array(
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
        return out

    @property
    def unclustered_parts(self):
        np_results = self._results.to_numpy_unclustered()
        out = ak.Array(
            ak.layout.RecordArray(
                [
                    ak.layout.NumpyArray(np_results[0]),
                    ak.layout.NumpyArray(np_results[1]),
                    ak.layout.NumpyArray(np_results[2]),
                    ak.layout.NumpyArray(np_results[3]),
                ],
                ["px", "py", "pz", "E"],
            )
        )
        return out

    @property
    def constituent_idx(self):
        np_results = self._results.to_numpy_with_constituents()
        off = np.insert(np_results[-1], 0, 0)
        out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np_results[0]), ak.layout.NumpyArray(np_results[1])
            )
        )
        out = ak.Array(ak.layout.ListOffsetArray64(ak.layout.Index64(off), out.layout))
        return out

    @property
    def constituent_particles(self):
        outputs_to_inputs = self.constituent_idx
        shape = ak.num(outputs_to_inputs)
        total = np.sum(shape)
        duplicate = ak.unflatten(np.zeros(total, np.int64), shape)
        prepared = self.data[:, np.newaxis][duplicate]
        return prepared[outputs_to_inputs]
