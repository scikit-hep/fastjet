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
        self._results = fastjet._ext.interface(px, py, pz, E, inps, inpf, jetdef)
        # self.inclusive_jets

    def swig_to_params(self, jetdef):
        params = jetdef.description().split()
        if params[0] == "e+e-":
            algor = params[0]
            Rv = float(params[7])  # e+e- generalised
        if params[0] == "e+e-" and params[2] == "(Durham)":
            algor = "Durham"  # e+e- durham
            Rv = 0
        if params[2] == "generalised":
            Rv = params[8]
            algor = "genkt"  # genralised kt
        else:
            algor = params[2]
            Rv = float(params[7])
        inps = {"algor": algor}  # kt or antikt
        inpf = {"R": Rv}
        return inps, inpf

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
        # np_results()[0] = ak.layout.NumpyArray(np_results()[0])
        # np_results()[1] = ak.layout.NumpyArray(np_results()[1])
        # np_results()[2] = ak.layout.NumpyArray(np_results()[2])
        # np_results()[3] = ak.layout.NumpyArray(np_results()[3])
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
