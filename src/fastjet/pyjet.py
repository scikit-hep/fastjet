import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.jagedness = self.check_jaggedness(data)
        # inps, inpf = self.swig_to_params(self.jetdef)
        # offsets = data["part0-node0-offsets"]
        if self.jagedness == 0:
            container, length, data = ak.to_buffers(data)
            px = data["part0-node1-data"]
            py = data["part0-node2-data"]
            pz = data["part0-node3-data"]
            E = data["part0-node4-data"]
            self._results = fastjet._ext.interface(
                self.px, self.py, self.pz, self.E, jetdef
            )
        px, py, pz, E, offsets = self.extract_cons(data)
        self.px = self.correct_byteorder(px)
        self.py = self.correct_byteorder(py)
        self.pz = self.correct_byteorder(pz)
        self.E = self.correct_byteorder(E)
        self.offsets = self.correct_byteorder(offsets)
        self._results = fastjet._ext.interfacemulti(
            self.px, self.py, self.pz, self.E, self.offsets, self.jetdef
        )

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

    @property
    def inclusive_jets(self):
        if self.jagedness == 0:
            np_results = self._results.to_numpy()
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
        if self.jagedness == 1:
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
    def constituents(self):
        np_results = self._results.to_numpy_with_constituents()
        off = np.insert(np_results[-1], 0, 0)
        out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(off), ak.layout.NumpyArray(np_results[-2])
            )
        )
        return out
