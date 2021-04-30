import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402
from fastjet.version import __version__  # noqa: E402

__all__ = ("__version__",)


class AwkwardClusterSequence:
    def __init__(self, data, R, algor):
        container, length, data = ak.to_buffers(data)
        # offsets = data["part0-node0-offsets"]
        px = data["part0-node1-data"]
        py = data["part0-node2-data"]
        pz = data["part0-node3-data"]
        E = data["part0-node4-data"]
        self.form = """ {
    "class": "ListOffsetArray64",
    "offsets": "i64",
    "content": {
        "class": "RecordArray",
        "contents": {
            "phi": {
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
            },
            "rap": {
                "class": "NumpyArray",
                "itemsize": 8,
                "format": "d",
                "primitive": "float64",
                "form_key": "node4"
            },
            "constituents": {
                "class": "ListOffsetArray64",
                "offsets": "i64",
                "content": {
                    "class": "NumpyArray",
                    "itemsize": 8,
                    "format": "l",
                    "primitive": "int64",
                    "form_key": "node6"
                },
                "form_key": "node5"
            }
        },
        "parameters": {
            "__record__": "Momentum4D"
        },
        "form_key": "node1"
    },
    "form_key": "node0"
}

"""
        self.out = fastjet._ext.interface(px, py, pz, E, R, algor)
        length = self.out["nevents"]
        # self.out.pop("nevents")
        # for key in self.out:
        # self.out[key] = np.asarray(self.out[key], dtype = "float64")
        # self.final = ak.from_buffers(self.form,length,self.out)
