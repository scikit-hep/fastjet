import warnings

import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402

_default_taus_njettiness = [1, 2, 3, 4]


class _classmultievent:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.data = data
        px, py, pz, E, starts, stops = self.extract_cons(self.data)
        px = self.correct_byteorder(px)
        py = self.correct_byteorder(py)
        pz = self.correct_byteorder(pz)
        E = self.correct_byteorder(E)
        starts = self.correct_byteorder(starts)
        stops = self.correct_byteorder(stops)
        self._results = fastjet._ext.interfacemulti(
            px, py, pz, E, starts, stops, jetdef
        )

    def _check_record(self, data):
        return data.layout.is_record or data.layout.is_numpy

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
        starts = np.asarray(array.layout.starts)
        stops = np.asarray(array.layout.stops)
        return px, py, pz, E, starts, stops

    def single_to_jagged(self, array):
        single = ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(np.array([0, len(array)])),
                ak.contents.RecordArray(
                    [
                        ak.contents.NumpyArray(array.px),
                        ak.contents.NumpyArray(array.py),
                        ak.contents.NumpyArray(array.pz),
                        ak.contents.NumpyArray(array.E),
                    ],
                    ["px", "py", "pz", "E"],
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=array.behavior,
            attrs=array.attrs,
        )
        return single

    def _warn_for_exclusive(self):
        if (
            self.jetdef
            not in [
                fastjet.kt_algorithm,
                fastjet.cambridge_algorithm,
                fastjet.ee_kt_algorithm,
                fastjet.plugin_algorithm,
            ]
        ) and (
            (self.jetdef not in [fastjet.kt_algorithm, fastjet.cambridge_algorithm])
            or self.jetdef.extra_param() < 0
        ):
            warnings.formatwarning = fastjet.formatwarning
            warnings.warn(
                "dcut and exclusive jets for jet-finders other than kt, C/A or genkt with p>=0 should be interpreted with care.",
                stacklevel=2,
            )
        return

    def inclusive_jets(self, min_pt):
        np_results = self._results.to_numpy(min_pt)
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def unclustered_particles(self):
        np_results = self._results.to_numpy_unclustered_particles()
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def exclusive_jets(self, n_jets, dcut):
        self._warn_for_exclusive()
        of = 0
        np_results = 0
        if n_jets == 0:
            raise ValueError("Njets cannot be 0")
        if dcut == -1 and n_jets != -1:
            np_results = self._results.to_numpy_exclusive_njet(n_jets)
            of = np_results[-1]
        if n_jets == -1 and dcut != -1:
            np_results = self._results.to_numpy_exclusive_dcut(dcut)
            of = np_results[-1]
        if np_results == 0 and of == 0:
            raise ValueError("Either NJets or Dcut sould be entered")
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def exclusive_jets_up_to(self, n_jets):
        self._warn_for_exclusive()
        if n_jets == 0:
            raise ValueError("Njets cannot be 0")
        np_results = self._results.to_numpy_exclusive_njet_up_to(n_jets)
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def exclusive_jets_ycut(self, ycut):
        np_results = self._results.to_numpy_exclusive_ycut(ycut)
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def constituent_index(self, min_pt):
        np_results = self._results.to_numpy_with_constituents(min_pt)
        off = np_results[-1]
        out = ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(np_results[0]), ak.contents.NumpyArray(np_results[1])
            )
        )
        out = ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
        return out

    def exclusive_jets_constituent_index(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        np_results = self._results.to_numpy_exclusive_njet_with_constituents(njets)
        off = np_results[-1]
        out = ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(np_results[0]), ak.contents.NumpyArray(np_results[1])
            )
        )
        out = ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
        return out

    def exclusive_jets_softdrop_grooming(
        self,
        njets=1,
        beta=0.0,
        symmetry_cut=0.1,
        symmetry_measure="scalar_z",
        R0=0.8,
        recursion_choice="larger_pt",
        # subtractor = 0,
        mu_cut=float("inf"),
    ):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")
        np_results = self._results.to_numpy_softdrop_grooming(
            njets,
            beta,
            symmetry_cut,
            symmetry_measure,
            R0,
            recursion_choice,  # subtractor,
            mu_cut,
        )

        px = ak.unflatten(
            ak.Array(ak.contents.NumpyArray(np_results[0])),
            ak.Array(ak.contents.NumpyArray(np_results[4])),
            highlevel=False,
        )
        py = ak.unflatten(
            ak.Array(ak.contents.NumpyArray(np_results[1])),
            ak.Array(ak.contents.NumpyArray(np_results[4])),
            highlevel=False,
        )
        pz = ak.unflatten(
            ak.Array(ak.contents.NumpyArray(np_results[2])),
            ak.Array(ak.contents.NumpyArray(np_results[4])),
            highlevel=False,
        )
        E = ak.unflatten(
            ak.Array(ak.contents.NumpyArray(np_results[3])),
            ak.Array(ak.contents.NumpyArray(np_results[4])),
            highlevel=False,
        )
        jetpt = ak.Array(ak.contents.NumpyArray(np_results[5]))
        jeteta = ak.Array(ak.contents.NumpyArray(np_results[6]))
        jetphi = ak.Array(ak.contents.NumpyArray(np_results[7]))
        jetmass = ak.Array(ak.contents.NumpyArray(np_results[8]))
        jetE = ak.Array(ak.contents.NumpyArray(np_results[9]))
        jetpz = ak.Array(ak.contents.NumpyArray(np_results[10]))
        jetdeltaR = ak.Array(ak.contents.NumpyArray(np_results[11]))
        jetsymmetry = ak.Array(ak.contents.NumpyArray(np_results[12]))

        out = ak.zip(
            {
                "constituents": ak.zip(
                    {"px": px, "py": py, "pz": pz, "E": E}, depth_limit=2
                ),
                "msoftdrop": jetmass,
                "ptsoftdrop": jetpt,
                "etasoftdrop": jeteta,
                "phisoftdrop": jetphi,
                "Esoftdrop": jetE,
                "pzsoftdrop": jetpz,
                "deltaRsoftdrop": jetdeltaR,
                "symmetrysoftdrop": jetsymmetry,
            },
            depth_limit=1,
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )
        return out

    def njettiness(
        self,
        measure_definition="NormalizedMeasure",
        axes_definition="OnePass_KT_Axes",
        njets=_default_taus_njettiness,
        beta=1.0,
        R0=0.8,
        Rcutoff=None,
        nPass=None,
        akAxesR0=None,
    ):
        if isinstance(njets, (int, float)):
            njets = [njets]
        if len(njets) == 0:
            raise ValueError("Must provide at least one njets!")
        if any(njet <= 0 for njet in njets):
            raise ValueError("Requested njets must be > 0!")

        double_max = 999.0
        int_max = 999

        np_results = self._results.to_numpy_njettiness(
            measure_definition,
            axes_definition,
            njets,
            beta,
            R0,
            Rcutoff or double_max,
            nPass or int_max,
            akAxesR0 or double_max,
        )
        out = ak.Array(
            ak.contents.NumpyArray(np_results[0]),
        )
        return out

    def exclusive_jets_energy_correlator(
        self,
        njets=1,
        beta=1,
        npoint=0,
        angles=-1,
        alpha=0,
        func="generalized",
        normalized=True,
    ):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")
        np_results = self._results.to_numpy_energy_correlators(
            njets,
            beta,
            npoint,
            angles,
            alpha,
            func,
            normalized,
        )
        out = ak.Array(ak.contents.NumpyArray(np_results))
        return out

    def exclusive_jets_lund_declusterings(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        np_results = self._results.to_numpy_exclusive_njet_lund_declusterings(njets)
        off = np_results[-1]
        out = ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(np_results[0]),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                    ),
                    ("Delta", "kt"),
                ),
            )
        )
        out = ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
        return out

    def unique_history_order(self):
        np_results = self._results.to_numpy_unique_history_order()
        off = np_results[-1]
        out = ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(off), ak.contents.NumpyArray(np_results[0])
            )
        )
        return out

    def exclusive_dmerge(self, njets):
        np_results = self._results.to_numpy_exclusive_dmerge(njets)
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def exclusive_dmerge_max(self, njets):
        np_results = self._results.to_numpy_exclusive_dmerge_max(njets)
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def exclusive_ymerge_max(self, njets):
        np_results = self._results.to_numpy_exclusive_ymerge_max(njets)
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def exclusive_ymerge(self, njets):
        np_results = self._results.to_numpy_exclusive_ymerge(njets)
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def Q(self):
        np_results = self._results.to_numpy_q()
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def Q2(self):
        np_results = self._results.to_numpy_q2()
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def constituents(self, min_pt):
        outputs_to_inputs = self.constituent_index(min_pt)
        shape = ak.num(outputs_to_inputs)
        total = np.sum(shape)
        duplicate = ak.unflatten(np.zeros(total, np.int64), shape)
        prepared = self.data[:, np.newaxis][duplicate]
        return prepared[outputs_to_inputs]

    def exclusive_jets_constituents(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        outputs_to_inputs = self.exclusive_jets_constituent_index(njets)
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
            of = np_results[-1]
        if nsub == -1 and dcut != -1:
            np_results = self._results.to_numpy_exclusive_subjets_dcut(
                px, py, pz, E, dcut
            )
            of = np_results[-1]
        if np_results == 0 and of == 0:
            raise ValueError("Either NJets or Dcut sould be entered") from None

        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
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
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
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
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def n_particles(self):
        np_results = self._results.to_numpy_n_particles()
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def n_exclusive_jets(self, dcut):
        np_results = self._results.to_numpy_n_exclusive_jets(dcut)
        out = ak.Array(ak.contents.NumpyArray(np_results[0]))
        return out

    def childless_pseudojets(self):
        np_results = self._results.to_numpy_childless_pseudojets()
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )

    def jets(self):
        np_results = self._results.to_numpy_jets()
        of = np_results[-1]
        return ak.Array(
            ak.contents.ListOffsetArray(
                ak.index.Index64(of),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
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
            ak.contents.ListOffsetArray(
                ak.index.Index64(np_results[-1]),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
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
            ak.contents.ListOffsetArray(
                ak.index.Index64(np_results[-1]),
                ak.contents.RecordArray(
                    (
                        ak.contents.NumpyArray(np_results[0]),
                        ak.contents.NumpyArray(np_results[1]),
                        ak.contents.NumpyArray(np_results[2]),
                        ak.contents.NumpyArray(np_results[3]),
                    ),
                    ("px", "py", "pz", "E"),
                    parameters={"__record__": "Momentum4D"},
                ),
            ),
            behavior=self.data.behavior,
            attrs=self.data.attrs,
        )
