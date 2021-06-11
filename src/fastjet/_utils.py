import fastjet._swig

import vector
import awkward as ak
import numpy as np


def sorted_by_E(data):
    if isinstance(data, ak.Array):
        try:
            tempE = data.E
        except:
            raise ValueError(
                "Needs either correct coordinates or embedded vector backend"
            )
        tmpsort = ak.argsort(tempE, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_E(data)


def sorted_by_pt(data):
    if isinstance(data, ak.Array):
        try:
            temppt = data.pt
        except:
            raise ValueError(
                "Needs either correct coordinates or embedded vector backend"
            )
        tmpsort = ak.argsort(temppt, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_pt(data)


def sorted_by_pz(data):
    if isinstance(data, ak.Array):
        try:
            temppz = data.pz
        except:
            raise ValueError(
                "Needs either correct coordinates or embedded vector backend"
            )
        tmpsort = ak.argsort(temppz, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_pz(data)


def sorted_by_rapidity(data):
    if isinstance(data, ak.Array):
        try:
            temprap = data.eta
        except:
            raise ValueError(
                "Needs either correct coordinates or embedded vector backend"
            )
        tmpsort = ak.argsort(temprap, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_pz(data)
