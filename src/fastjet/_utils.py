import awkward as ak

import fastjet._swig

# light wrapping for the functions to raise an error if the user inputs awkward arrays into functions meant for swig


def sorted_by_E(data):
    if isinstance(data, ak.Array):
        try:
            tempE = data.E
        except AttributeError:
            raise AttributeError(
                "Needs either correct coordinates or embedded vector backend"
            ) from None
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
            ) from None
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
            ) from None
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
            ) from None
        tmpsort = ak.argsort(temprap, axis=-1)
        return data[tmpsort]
    else:
        return fastjet._swig.sorted_by_pz(data)


def join(*argv):
    if isinstance(argv[0], ak.Array):
        raise TypeError("Use inbuilt methods for awkward arrays")
    else:
        for arg in argv:
            if isinstance(arg, ak.Array):
                raise AttributeError(
                    "All arguments need to be of the same type"
                ) from None
        if len(argv) == 1:  # Calling different constructors
            return fastjet._swig.join(argv[0])
        if len(argv) == 2:
            return fastjet._swig.join(argv[0], argv[1])
        if len(argv) == 3:
            return fastjet._swig.join(argv[0], argv[1], argv[2])
        if len(argv) == 4:
            return fastjet._swig.join(argv[0], argv[1], argv[2], argv[3])
        if len(argv) > 4:
            raise ValueError("Length exceeded") from None


def dot_product(a, b):
    if isinstance(a, ak.Array) or isinstance(b, ak.Array):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return fastjet._swig.dot_product(a, b)


def sort_indices(indices, values):
    if isinstance(indices, ak.Array) or isinstance(values, ak.Array):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return sort_indices(indices, values)


def theta(a, b):
    if isinstance(a, ak.Array) or isinstance(b, ak.Array):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return fastjet._swig.theta(a, b)


def have_same_momentum(a, b):
    if isinstance(a, ak.Array) or isinstance(b, ak.Array):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return fastjet._swig.have_same_momentum(a, b)


def cos_theta(a, b):
    if isinstance(a, ak.Array) or isinstance(b, ak.Array):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return fastjet._swig.cos_theta(a, b)


def PtYPhiM(pt, y, phi, m):
    if (
        isinstance(pt, ak.Array)
        or isinstance(y, ak.Array)
        or isinstance(phi, ak.Array)
        or isinstance(m, ak.Array)
    ):
        raise TypeError("Use inbuilt methods for Awkward Array") from None
    else:
        return fastjet._swig.PtYPhiM(pt, y, phi, m)
