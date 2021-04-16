// BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

#include <pybind11/pybind11.h>

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(_ext, m) {
    m.def("add", &add);
    m.def("subtract", [](int i, int j) { return i - j; });
}
