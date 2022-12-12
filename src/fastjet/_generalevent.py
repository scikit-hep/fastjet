import warnings

import awkward as ak
import numpy as np

import fastjet._ext  # noqa: F401, E402


class _classgeneralevent:
    def __init__(self, data, jetdef):
        self.jetdef = jetdef
        self.data = data
        self._mod_data = data
        self._bread_list = []
        self._clusterable_level = []
        self._results = []

        self._input_mapping = []
        self.multi_layered_listoffset(self.data, ())
        for i in range(len(self._clusterable_level)):
            self._clusterable_level[i] = ak.Array(
                self._clusterable_level[i].layout.to_ListOffsetArray64(True)
            )
            px, py, pz, E, offsets = self.extract_cons(self._clusterable_level[i])
            px = self.correct_byteorder(px)
            py = self.correct_byteorder(py)
            pz = self.correct_byteorder(pz)
            E = self.correct_byteorder(E)
            offsets = self.correct_byteorder(offsets)
            self._results.append(
                fastjet._ext.interfacemulti(px, py, pz, E, offsets, jetdef)
            )

    def _check_listoffset_subtree(self, data):
        return data.layout.is_list

    def _check_listoffset(self, data):
        return isinstance(data.layout, ak.contents.ListOffsetArray)

    def _check_record(self, data):
        return data.layout.is_record or data.layout.is_numpy

    def _check_indexed(self, data):
        return data.layout.is_indexed

    def multi_layered_listoffset_input(self, data, crumb_list):
        if data.layout.is_union:
            for i in range(len(data.layout.contents)):
                temp_crumb = crumb_list + (i,)
                if self._check_subtree_input(ak.Array(data.layout.contents[i])):
                    self._bread_list_input.append(crumb_list)
                    return
                self.multi_layered_listoffset_input(
                    ak.Array(data.layout.contents[i]), temp_crumb
                )
        elif data.layout.is_record:
            for elem in data.layout.fields:
                temp_crumb = crumb_list + (elem,)
                if self._check_subtree_input(ak.Array(data.layout.content(elem))):
                    self._bread_list_input.append(crumb_list)
                    return
                self.multi_layered_listoffset_input(
                    ak.Array(data.layout.content(elem)), temp_crumb
                )
            return
        elif self._check_record(ak.Array(data.layout.content)):
            attributes = dir(data)
            if (
                "px" in attributes
                and "py" in attributes
                and "pz" in attributes
                and "E" in attributes
            ):
                crumb_list = crumb_list + (None,)
                self._bread_list_input.append(crumb_list)
                self._cluster_inputs.append(ak.Array(data.layout.content))
            else:
                crumb_list = crumb_list + (None,)
                self.multi_layered_listoffset_input(
                    ak.Array(data.layout.content), crumb_list
                )
        else:
            crumb_list = crumb_list + (None,)
            self.multi_layered_listoffset_input(
                ak.Array(data.layout.content), crumb_list
            )

    def multi_layered_listoffset(self, data, crumb_list):
        if data.layout.is_union:
            for i in range(len(data.layout.contents)):
                temp_crumb = crumb_list + (i,)
                if self._check_subtree(ak.Array(data.layout.contents[i])):
                    self._bread_list.append(crumb_list)
                    return
                self.multi_layered_listoffset(
                    ak.Array(data.layout.contents[i]), temp_crumb
                )
        elif data.layout.is_record:
            for elem in data.layout.fields:
                temp_crumb = crumb_list + (elem,)
                if self._check_subtree(ak.Array(data.layout.content(elem))):
                    self._bread_list.append(crumb_list)
                    return
                self.multi_layered_listoffset(
                    ak.Array(data.layout.content(elem)), temp_crumb
                )
            return
        elif self._check_listoffset_subtree(ak.Array(data.layout.content)):
            if self._check_record(
                ak.Array(ak.Array(data.layout.content).layout.content),
            ):
                attributes = dir(data)
                if (
                    "px" in attributes
                    and "py" in attributes
                    and "pz" in attributes
                    and "E" in attributes
                ):
                    crumb_list = crumb_list + (None,)
                    self._bread_list.append(crumb_list)
                    self._clusterable_level.append(ak.Array(data.layout.content))
            elif self._check_indexed(
                ak.Array(ak.Array(data.layout.content).layout.content),
            ):
                if self._check_record(
                    ak.Array(
                        ak.Array(
                            ak.Array(data.layout.content).layout.content
                        ).layout.content
                    ),
                ):
                    attributes = dir(data)
                    if (
                        "px" in attributes
                        and "py" in attributes
                        and "pz" in attributes
                        and "E" in attributes
                    ):
                        crumb_list = crumb_list + (None,)
                        self._bread_list.append(crumb_list)
                        self._clusterable_level.append(ak.Array(data.layout.content))
                else:
                    crumb_list = crumb_list + (None,)
                    self.multi_layered_listoffset(
                        ak.Array(data.layout.content), crumb_list
                    )
            else:
                crumb_list = crumb_list + (None,)
                self.multi_layered_listoffset(ak.Array(data.layout.content), crumb_list)
        else:
            crumb_list = crumb_list + (None,)
            self.multi_layered_listoffset(ak.Array(data.layout.content), crumb_list)

    def correct_byteorder(self, data):
        if data.dtype.byteorder == "=":
            pass
        else:
            data = data.dtype.newbyteorder("=")
        return data

    def _check_subtree(self, data):
        if self._check_listoffset_subtree(ak.Array(data.layout)):
            if self._check_record(
                ak.Array(ak.Array(data.layout.content)),
            ):
                attributes = dir(data)
                if (
                    "px" in attributes
                    and "py" in attributes
                    and "pz" in attributes
                    and "E" in attributes
                ):
                    self._clusterable_level.append(data)
                    return True
            elif self._check_indexed(
                ak.Array(data.layout.content),
            ):
                if self._check_record(
                    ak.Array(ak.Array(data.layout.content).layout.content)
                ):
                    attributes = dir(data)
                    if (
                        "px" in attributes
                        and "py" in attributes
                        and "pz" in attributes
                        and "E" in attributes
                    ):
                        self._clusterable_level.append(data)
                        return True
            else:
                return False
        else:
            return False

    def _check_subtree_input(self, data):
        if self._check_record(
            data,
        ):
            attributes = dir(data)
            if (
                "px" in attributes
                and "py" in attributes
                and "pz" in attributes
                and "E" in attributes
            ):
                self._cluster_inputs.append(data)
                return True
            else:
                return False
        else:
            return False

    def extract_cons(self, array):
        px = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).px)
        py = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).py)
        pz = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).pz)
        E = np.asarray(ak.Array(array.layout.content, behavior=array.behavior).E)
        off = np.asarray(array.layout.stops)
        off = np.insert(off, 0, 0)
        return px, py, pz, E, off

    def _replace_multi(self):
        self._mod_data = self.data
        if self._input_flag == 0:
            for i in range(len(self._clusterable_level)):
                self._cur_idx = i
                self._mod_data = ak.Array(self.replace(self._mod_data.layout, i, 0))
            return self._mod_data
        else:
            for i in range(len(self._input_mapping)):
                self._cur_idx = i
                self._mod_data_input = ak.Array(
                    self.replace(self._mod_data_input.layout, self._input_mapping[i], 0)
                )
            return self._mod_data_input.layout

    def replace(self, layout, cluster, level):
        if level == len(self._bread_list[cluster]):
            return self._out[self._cur_idx].layout

        elif self._check_listoffset(ak.Array(layout)):
            if isinstance(layout, ak.contents.ListOffsetArray):
                return layout.copy(
                    content=self.replace(layout.content, cluster, level + 1)
                )
        elif isinstance(
            layout,
            (
                ak.contents.BitMaskedArray,  # all of the layouts with one 'content'
                ak.contents.ByteMaskedArray,
                ak.contents.IndexedArray,
                ak.contents.IndexedOptionArray,
                ak.contents.ListArray,
                ak.contents.RegularArray,
                ak.contents.UnmaskedArray,
            ),
        ):
            return layout.copy(content=self.replace(layout.content, cluster, level + 1))
        elif layout.is_record:
            nextcontents = []
            for elem in layout.fields:
                if elem == self._bread_list[cluster][level]:
                    nextcontents.append(
                        self.replace(
                            layout.content(self._bread_list[cluster][level]),
                            cluster,
                            level + 1,
                        )
                    )
                else:
                    nextcontents.append(
                        layout.content(elem),
                    )
            return layout.copy(contents=nextcontents)
        elif isinstance(layout, ak.record.Record):
            return ak.record.Record(
                self.replace(layout.array, cluster, level + 1),
                layout.at,
            )
        elif layout.is_union:
            nextcontents = []
            for i in range(len(layout.contents)):
                if i == self._bread_list[cluster][level]:
                    nextcontents.append(
                        self.replace(
                            layout.contents[self._bread_list[cluster][level]],
                            cluster,
                            level + 1,
                        )
                    )
                else:
                    nextcontents.append(
                        layout.contents[i],
                    )
            return layout.copy(contents=nextcontents)

        raise AssertionError(layout)

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
                "dcut and exclusive jets for jet-finders other than kt, C/A or genkt with p>=0 should be interpreted with care."
            )
        return

    def inclusive_jets(self, min_pt):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy(min_pt)
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def constituents(self, min_pt):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_with_constituents(min_pt)
            off = np_results[-1]
            out = ak.Array(
                ak.contents.ListOffsetArray(
                    ak.index.Index64(np_results[0]),
                    ak.contents.NumpyArray(np_results[1]),
                ),
                behavior=self.data.behavior,
            )
            outputs_to_inputs = ak.Array(
                ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout)
            )
            shape = ak.num(outputs_to_inputs)
            total = np.sum(shape)
            duplicate = ak.unflatten(np.zeros(total, np.int64), shape)
            prepared = self._clusterable_level[i][:, np.newaxis][duplicate]
            self._out.append(prepared[outputs_to_inputs])
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_jets_constituents(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_njet_with_constituents(
                njets
            )
            off = np_results[-1]
            out = ak.Array(
                ak.contents.ListOffsetArray(
                    ak.index.Index64(np_results[0]),
                    ak.contents.NumpyArray(np_results[1]),
                ),
                behavior=self.data.behavior,
            )
            outputs_to_inputs = ak.Array(
                ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout)
            )
            shape = ak.num(outputs_to_inputs)
            total = np.sum(shape)
            duplicate = ak.unflatten(np.zeros(total, np.int64), shape)
            prepared = self._clusterable_level[i][:, np.newaxis][duplicate]
            self._out.append(prepared[outputs_to_inputs])
        res = ak.Array(self._replace_multi())
        return res

    def constituent_index(self, min_pt):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_with_constituents(min_pt)
            off = np_results[-1]
            out = ak.Array(
                ak.contents.ListOffsetArray(
                    ak.index.Index64(np_results[0]),
                    ak.contents.NumpyArray(np_results[1]),
                ),
                behavior=self.data.behavior,
            )
            self._out.append(
                ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_jets_constituent_index(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_njet_with_constituents(
                njets
            )
            off = np_results[-1]
            out = ak.Array(
                ak.contents.ListOffsetArray(
                    ak.index.Index64(np_results[0]),
                    ak.contents.NumpyArray(np_results[1]),
                ),
                behavior=self.data.behavior,
            )
            self._out.append(
                ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_jets_lund_declusterings(self, njets):
        if njets <= 0:
            raise ValueError("Njets cannot be <= 0")

        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_njet_lund_declusterings(
                njets
            )
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
                ),
                behavior=self.data.behavior,
            )
            self._out.append(
                ak.Array(ak.contents.ListOffsetArray(ak.index.Index64(off), out.layout))
            )
        res = ak.Array(self._replace_multi())
        return res

    def unclustered_particles(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_unclustered_particles()
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def n_particles(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_n_particles()
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def n_exclusive_jets(self, dcut):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_n_exclusive_jets(dcut)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def childless_pseudojets(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_childless_pseudojets()
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def jets(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_jets()
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_jets(self, n_jets, dcut):
        self._warn_for_exclusive()
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            of = 0
            np_results = 0
            if n_jets == 0:
                raise ValueError("Njets cannot be 0")
            if dcut == -1 and n_jets != -1:
                np_results = self._results[i].to_numpy_exclusive_njet(n_jets)
                of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            if n_jets == -1 and dcut != -1:
                np_results = self._results[i].to_numpy_exclusive_dcut(dcut)
                of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            if np_results == 0 and of == 0:
                raise ValueError("Either NJets or Dcut sould be entered")
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_jets_ycut(self, ycut):
        self._warn_for_exclusive()
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_ycut(ycut)
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def unique_history_order(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_unique_history_order()
            off = np_results[-1]
            self._out.append(
                ak.Array(
                    ak.contents.ListOffsetArray(
                        ak.index.Index64(off), ak.contents.NumpyArray(np_results[0])
                    )
                )
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_dmerge(self, njets):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_dmerge(njets)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_dmerge_max(self, njets):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_dmerge_max(njets)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_ymerge_max(self, njets):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_ymerge_max(njets)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_ymerge(self, njets):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_exclusive_ymerge(njets)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def Q(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_q()
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def Q2(self):
        self._out = []
        self._input_flag = 0
        for i in range(len(self._clusterable_level)):
            np_results = self._results[i].to_numpy_q2()
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def get_parents(self, data_inp):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            idx = -1
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_get_parents(px, py, pz, E)
            self._out.append(
                ak.Array(
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
                )
            )

        res = ak.Array(self._replace_multi())
        return res

    def exclusive_subdmerge(self, data_inp, nsub):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):

            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            idx = -1
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_exclusive_subdmerge(
                px, py, pz, E, nsub
            )
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))

        res = ak.Array(self._replace_multi())
        return res

    def exclusive_subjets(self, data_inp, dcut, nsub):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            of = 0
            np_results = 0
            if nsub == 0:
                raise ValueError("Njets cannot be 0")
            if dcut == -1 and nsub != -1:
                np_results = self._results[idx].to_numpy_exclusive_subjets_nsub(
                    px, py, pz, E, nsub
                )
                of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            if nsub == -1 and dcut != -1:
                np_results = self._results[idx].to_numpy_exclusive_subjets_dcut(
                    px, py, pz, E, dcut
                )
                of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            if np_results == 0 and of == 0:
                raise ValueError("Either NJets or Dcut sould be entered")

            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_subjets_up_to(self, data_inp, nsub):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_exclusive_subjets_up_to(
                px, py, pz, E, nsub
            )
            of = np.insert(np_results[-1], len(np_results[-1]), len(np_results[0]))
            self._out.append(
                ak.Array(
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
            )
        res = ak.Array(self._replace_multi())
        return res

    def exclusive_subdmerge_max(self, data_inp, nsub):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_exclusive_subdmerge_max(
                px, py, pz, E, nsub
            )
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def n_exclusive_subjets(self, data_inp, dcut):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_n_exclusive_subjets(
                px, py, pz, E, dcut
            )
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def has_parents(self, data_inp):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_has_parents(px, py, pz, E)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def has_child(self, data_inp):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_has_child(px, py, pz, E)
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def jet_scale_for_algorithm(self, data_inp):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_jet_scale_for_algorithm(
                px, py, pz, E
            )
            self._out.append(ak.Array(ak.contents.NumpyArray(np_results[0])))
        res = ak.Array(self._replace_multi())
        return res

    def get_child(self, data_inp):
        self._cluster_inputs = []
        self._bread_list_input = []
        self._input_mapping = []
        self._out = []
        self._mod_data_input = data_inp
        self._input_flag = 1
        self.multi_layered_listoffset_input(data_inp, ())
        if len(self._cluster_inputs) == 0:
            raise TypeError("The Awkward Array is not valid")
        for i in range(len(self._cluster_inputs)):
            px = self._cluster_inputs[i].px
            py = self._cluster_inputs[i].py
            pz = self._cluster_inputs[i].pz
            E = self._cluster_inputs[i].E
            for j in range(len(self._bread_list)):
                if self._bread_list[j] == self._bread_list_input[i]:
                    idx = j
                    self._input_mapping.append(j)
            if idx == -1:
                continue
            assert len(self._cluster_inputs[i]) == len(self._clusterable_level[idx])
            np_results = self._results[idx].to_numpy_get_child(px, py, pz, E)
            self._out.append(
                ak.Array(
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
                )
            )
        res = ak.Array(self._replace_multi())
        return res
