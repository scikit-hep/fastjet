import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402


class _classmultievent:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        # data = ak.Array(data.layout.toListOffsetArray64(True)) If I do this vector can't convert the coordinates
        self.data = data
        px, py, pz, E, offsets = self.extract_cons(self.data)
        px = self.correct_byteorder(px)
        py = self.correct_byteorder(py)
        pz = self.correct_byteorder(pz)
        E = self.correct_byteorder(E)
        offsets = self.correct_byteorder(offsets)
        self._results = fastjet._ext.interfacemulti(px, py, pz, E, offsets, jetdef)

    def _check_record(self, data):
        out = isinstance(
            data.layout,
            (
                ak.layout.RecordArray,
                ak.layout.NumpyArray,
            ),
        )
        return out

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
                    parameters={"__record__": "Momentum4D"},
                ),
            )
        )
        return single

    def inclusive_jets(self, min_pt):
        np_results = self._results.to_numpy(min_pt)
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def unclustered_particles(self):
        np_results = self._results.to_numpy_unclustered_particles()
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def exclusive_jets(self, n_jets, dcut):
        of = 0
        np_results = 0
        if n_jets == 0:
            raise ValueError("Njets cannot be 0")
        if dcut == -1 and n_jets != -1:
            np_results = self._results.to_numpy_exclusive_njet(n_jets)
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        if n_jets == -1 and dcut != -1:
            np_results = self._results.to_numpy_exclusive_dcut(dcut)
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        if np_results == 0 and of == 0:
            raise ValueError("Either NJets or Dcut sould be entered")
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def exclusive_jets_ycut(self, ycut):
        np_results = self._results.to_numpy_exclusive_ycut(ycut)
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def constituent_index(self, min_pt):
        np_results = self._results.to_numpy_with_constituents(min_pt)
        off = np.insert(np_results[-1], 0, 0)
        out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np_results[0]), ak.layout.NumpyArray(np_results[1])
            )
        )
        out = ak.Array(ak.layout.ListOffsetArray64(ak.layout.Index64(off), out.layout))
        return out

    def unique_history_order(self):
        np_results = self._results.to_numpy_unique_history_order()
        off = np.insert(np_results[-1], 0, 0)
        out = ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(off), ak.layout.NumpyArray(np_results[0])
            )
        )
        return out

    def exclusive_dmerge(self, njets):
        np_results = self._results.to_numpy_exclusive_dmerge(njets)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def exclusive_dmerge_max(self, njets):
        np_results = self._results.to_numpy_exclusive_dmerge_max(njets)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def exclusive_ymerge_max(self, njets):
        np_results = self._results.to_numpy_exclusive_ymerge_max(njets)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def exclusive_ymerge(self, njets):
        np_results = self._results.to_numpy_exclusive_ymerge(njets)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def Q(self):
        np_results = self._results.to_numpy_q()
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def Q2(self):
        np_results = self._results.to_numpy_q2()
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def constituents(self, min_pt):
        outputs_to_inputs = self.constituent_index(min_pt)
        shape = ak.num(outputs_to_inputs)
        total = np.sum(shape)
        duplicate = ak.unflatten(np.zeros(total, np.int64), shape)
        prepared = self.data[:, np.newaxis][duplicate]
        return prepared[outputs_to_inputs]

    def exclusive_subjets(self, data, dcut, nsub):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        of = 0
        np_results = 0
        if nsub == 0:
            raise ValueError("Njets cannot be 0")
        if dcut == -1 and nsub != -1:
            np_results = self._results.to_numpy_exclusive_subjets_nsub(
                px, py, pz, E, nsub
            )
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        if nsub == -1 and dcut != -1:
            np_results = self._results.to_numpy_exclusive_subjets_dcut(
                px, py, pz, E, dcut
            )
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        if np_results == 0 and of == 0:
            raise ValueError("Either NJets or Dcut sould be entered") from None

        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def exclusive_subjets_up_to(self, data, nsub):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_exclusive_subjets_up_to(px, py, pz, E, nsub)
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def exclusive_subdmerge(self, data, nsub):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_exclusive_subdmerge(px, py, pz, E, nsub)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def exclusive_subdmerge_max(self, data, nsub):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_exclusive_subdmerge_max(px, py, pz, E, nsub)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def n_exclusive_subjets(self, data, dcut):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_n_exclusive_subjets(px, py, pz, E, dcut)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def has_parents(self, data):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_has_parents(px, py, pz, E)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def has_child(self, data):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_has_child(px, py, pz, E)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def jet_scale_for_algorithm(self, data):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_jet_scale_for_algorithm(px, py, pz, E)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def n_particles(self):
        np_results = self._results.to_numpy_n_particles()
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def n_exclusive_jets(self, dcut):
        np_results = self._results.to_numpy_n_exclusive_jets(dcut)
        out = ak.Array(ak.layout.NumpyArray(np_results[0]))
        return out

    def childless_pseudojets(self):
        np_results = self._results.to_numpy_childless_pseudojets()
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def jets(self):
        np_results = self._results.to_numpy_jets()
        of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
        return ak.Array(
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
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def get_parents(self, data):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_get_parents(px, py, pz, E)
        return ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np_results[-1]),
                ak.layout.RecordArray(
                    (
                        ak.layout.NumpyArray(np_results[0]),
                        ak.layout.NumpyArray(np_results[1]),
                        ak.layout.NumpyArray(np_results[2]),
                        ak.layout.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )

    def get_child(self, data):
        try:
            px = data.px
            py = data.py
            pz = data.pz
            E = data.E
        except AttributeError:
            raise AttributeError("Lorentz vector not found") from None
        np_results = self._results.to_numpy_get_child(px, py, pz, E)
        return ak.Array(
            ak.layout.ListOffsetArray64(
                ak.layout.Index64(np_results[-1]),
                ak.layout.RecordArray(
                    (
                        ak.layout.NumpyArray(np_results[0]),
                        ak.layout.NumpyArray(np_results[1]),
                        ak.layout.NumpyArray(np_results[2]),
                        ak.layout.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
        )
