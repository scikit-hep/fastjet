import awkward as ak

import fastjet._ext  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        inps, inpf = self.swig_to_params(self.jetdef)
        container, length, data = ak.to_buffers(data)
        # offsets = data["part0-node0-offsets"]
        data = self.correct_byteorder(data)
        px = data["part0-node1-data"]
        py = data["part0-node2-data"]
        pz = data["part0-node3-data"]
        E = data["part0-node4-data"]
        self._results = fastjet._ext.interface(px, py, pz, E, jetdef)

    def correct_byteorder(self, data):
        for keys in data:
            if data[keys].dtype.byteorder == "=":
                pass
            else:
                data[keys] = data[keys].dtype.nebyteorder("=")
        return data

    @property
    def inclusive_jets(self):
        np_results = self._results.to_numpy
        out = ak.Array(
            ak.layout.RecordArray(
                [
                    ak.layout.NumpyArray(np_results()[0]),
                    ak.layout.NumpyArray(np_results()[1]),
                    ak.layout.NumpyArray(np_results()[2]),
                    ak.layout.NumpyArray(np_results()[3]),
                ],
                ["px", "py", "pz", "E"],
            )
        )

        return out
