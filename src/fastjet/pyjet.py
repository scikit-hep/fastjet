import awkward as ak

import fastjet._ext  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, jetdef):
        self.params = jetdef.description().split()
        algor = self.params[2]
        Rv = float(self.params[7])
        inps = {"algor": algor}
        inpf = {"R": Rv}
        container, length, data = ak.to_buffers(data)
        # offsets = data["part0-node0-offsets"]
        px = data["part0-node1-data"]
        py = data["part0-node2-data"]
        pz = data["part0-node3-data"]
        E = data["part0-node4-data"]
        self.form = """ {
    "class": "RecordArray",
    "contents": {
        "phi": {
            "class": "NumpyArray",
            "itemsize": 8,
            "format": "d",
            "primitive": "float64",
            "form_key": "node1"
        },
        "rap": {
            "class": "NumpyArray",
            "itemsize": 8,
            "format": "d",
            "primitive": "float64",
            "form_key": "node2"
        },
        "pt": {
            "class": "NumpyArray",
            "itemsize": 8,
            "format": "d",
            "primitive": "float64",
            "form_key": "node3"
        }
    },
    "parameters": {
        "__record__": "Momentum4D"
    },
    "form_key": "node0"
}
"""
        results = fastjet._ext.interface(px, py, pz, E, inps, inpf)
        size = len(results["part0-node1-data"])
        self.inclusive_jets = ak.from_buffers(self.form, size, results)
