import awkward as ak

import fastjet._swig


def sorted_by_E(data):
    if isinstance(data, ak.Array):
        try:
            tempE = data.E
        except AttributeError:
            raise AttributeError(
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
        except AttributeError:
            raise AttributeError(
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
        except AttributeError:
            raise AttributeError(
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
        except AttributeError:
            raise AttributeError(
                "Needs either correct coordinates or embedded vector backend"
            )
        tmpsort = ak.argsort(temprap, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_pz(data)


def join(*argv):
    if isinstance(argv[0], ak.Array):
        for arg in argv:
            if not isinstance(arg, ak.Array):
                raise AttributeError("All arguments need to be of the same type")
        # CODE FOR THIS NEEDS TO BE WRITTEN
    else:
        for arg in argv:
            if isinstance(arg, ak.Array):
                raise AttributeError("All arguments need to be of the same type")
        if len(argv) == 1:  # Calling different constructors
            return fastjet._swig.join(argv[0])
        if len(argv) == 2:
            return fastjet._swig.join(argv[0], argv[1])
        if len(argv) == 3:
            return fastjet._swig.join(argv[0], argv[1], argv[2])
        if len(argv) == 4:
            return fastjet._swig.join(argv[0], argv[1], argv[2], argv[3])
        if len(argv) > 4:
            raise ValueError("Length exceeded")
