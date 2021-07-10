import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402


class _classcomplexevent:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.data = data
        self.multi_layered_listoffset(self.data)
        px, py, pz, E, offsets = self.extract_cons(self.data)
        px = self.correct_byteorder(px)
        py = self.correct_byteorder(py)
        pz = self.correct_byteorder(pz)
        E = self.correct_byteorder(E)
        offsets = self.correct_byteorder(offsets)
        self._results = fastjet._ext.interfacemulti(px, py, pz, E, offsets, jetdef)

    def multi_layered_listoffset(self, data):
        if isinstance(
            ak.Array(ak.Array(data.layout.content).layout.content).layout,
            ak.layout.RecordArray,
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
