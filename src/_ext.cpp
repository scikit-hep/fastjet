// BSD 3-Clause License; see
// https://github.com/scikit-hep/fastjet/blob/main/LICENSE

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <unordered_map>
#include <vector>

#include <fastjet/AreaDefinition.hh>
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/LundGenerator.hh>

#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace fj = fastjet;
namespace py = pybind11;
using namespace pybind11::literals;

typedef struct {
  PyObject_HEAD void *ptr;
  void *ty;
  int own;
  PyObject *next;
} SwigPyObject;

template <typename T>
T swigtocpp(py::object obj) { // unwraps python object to get the cpp pointer
                              // from the swig bindings
  auto upointer = obj.attr("this").ptr();
  auto swigpointer = reinterpret_cast<SwigPyObject *>(upointer);
  auto objpointervoid = swigpointer->ptr;
  auto objpointer = reinterpret_cast<T>(objpointervoid);
  return objpointer;
}
class output_wrapper {
public:
  std::vector<std::shared_ptr<fj::ClusterSequence>> cse;
  std::vector<std::shared_ptr<std::vector<fj::PseudoJet>>> parts;

  std::shared_ptr<fj::ClusterSequence> getCluster() {
    auto a = cse[0];
    return a;
  }
  void setCluster() {}
};

output_wrapper interfacemulti(
    py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
    py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
    py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
    py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
    py::array_t<int, py::array::c_style | py::array::forcecast> offsets,
    py::object jetdef) {
  py::buffer_info infooff = offsets.request();
  py::buffer_info infopx = pxi.request();
  py::buffer_info infopy =
      pyi.request(); // requesting buffer information of the input
  py::buffer_info infopz = pzi.request();
  py::buffer_info infoE = Ei.request();

  auto offptr = static_cast<int *>(infooff.ptr);
  auto pxptr = static_cast<double *>(infopx.ptr);
  auto pyptr =
      static_cast<double *>(infopy.ptr); // pointer to the initial value
  auto pzptr = static_cast<double *>(infopz.ptr);
  auto Eptr = static_cast<double *>(infoE.ptr);

  int dimoff = infooff.shape[0];
  output_wrapper ow;
  std::vector<double> nevents;
  std::vector<double> offidx;
  std::vector<double> constphi;
  std::vector<double> idx;
  std::vector<double> idxo;
  for (int i = 0; i < dimoff - 1; i++) {
    std::vector<fj::PseudoJet> particles;
    for (int j = *offptr; j < *(offptr + 1); j++) {
      particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
      pxptr++;
      pyptr++;
      pzptr++;
      Eptr++;
    }

    std::vector<fj::PseudoJet> jets;
    auto jet_def = swigtocpp<fj::JetDefinition *>(jetdef);
    std::shared_ptr<std::vector<fj::PseudoJet>> pj =
        std::make_shared<std::vector<fj::PseudoJet>>(particles);
    std::shared_ptr<fastjet::ClusterSequence> cs =
        std::make_shared<fastjet::ClusterSequence>(*pj, *jet_def);
    auto j = cs->inclusive_jets();
    offptr++;
    ow.cse.push_back(cs);
    ow.parts.push_back(pj);
  }
  return ow;
}

PYBIND11_MODULE(_ext, m) {
  using namespace fastjet;
  m.def("interfacemulti", &interfacemulti,
        py::return_value_policy::take_ownership);
  /// Jet algorithm definitions

  py::class_<output_wrapper>(m, "output_wrapper")
    .def_property("cse", &output_wrapper::getCluster,&output_wrapper::setCluster)
    .def("to_numpy",
      [](const output_wrapper ow, double min_pt = 0) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->inclusive_jets(min_pt).size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->inclusive_jets(min_pt);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets from multievent clustering and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_with_constituents",
      [](const output_wrapper ow, double min_pt = 0) {
        auto css = ow.cse;
        int64_t len = css.size();
        auto jk = 0;
        auto sizepar = 0;

        for(int i = 0; i < len; i++){
        jk += css[i]->inclusive_jets(min_pt).size();
        sizepar += css[i]->n_particles();
        }
        jk++;

        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {sizepar}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;

        ptreventoffsets[eventidx] = 0;
        eventidx++;

        auto jetoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufjetoffsets = jetoffsets.request();
        int *ptrjetoffsets = (int *)bufjetoffsets.ptr;
        size_t jetidx = 0;

        size_t idxh = 0;
        ptrjetoffsets[jetidx] = 0;
        jetidx++;
        auto eventprev = 0;


        for (unsigned int i = 0; i < css.size(); i++){

        auto jets = css[i]->inclusive_jets(min_pt);
        int size = css[i]->inclusive_jets(min_pt).size();
        auto idx = css[i]->particle_jet_indices(jets);
        int64_t sizz = css[i]->n_particles();
        auto prev = ptrjetoffsets[jetidx-1];

        for (unsigned int j = 0; j < jets.size(); j++){
        ptrjetoffsets[jetidx] = jets[j].constituents().size() + prev;
        prev = ptrjetoffsets[jetidx];
        jetidx++;
        }
        for(int k = 0; k < size; k++){
          for(int j = 0; j <sizz; j++){
            if(idx[j] == k){
              ptrid[idxh] = j;
              idxh++;
            }
          }
        }
        ptreventoffsets[eventidx] = jets.size()+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            jetoffsets,
            parid,
            eventoffsets
          );
      }, "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_njet",
      [](const output_wrapper ow, const int n_jets = 0) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->exclusive_jets(n_jets).size();
        }

        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->exclusive_jets(n_jets);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, "n_jets"_a = 0, R"pbdoc(
        Retrieves the exclusive jets upto n jets from multievent clustering and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_njet_with_constituents",
      [](const output_wrapper ow, const int n_jets = 0) {
        auto css = ow.cse;
        int64_t len = css.size();
        auto jk = 0;
        auto sizepar = 0;

        for(int i = 0; i < len; i++){
          jk += css[i]->exclusive_jets(n_jets).size();
          sizepar += css[i]->n_particles();
        }
        jk++;

        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {sizepar}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;

        ptreventoffsets[eventidx] = 0;
        eventidx++;

        auto jetoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufjetoffsets = jetoffsets.request();
        int *ptrjetoffsets = (int *)bufjetoffsets.ptr;
        size_t jetidx = 0;

        size_t idxh = 0;
        ptrjetoffsets[jetidx] = 0;
        jetidx++;
        auto eventprev = 0;


        for (unsigned int i = 0; i < css.size(); i++){  // iterate through events
            auto jets = css[i]->exclusive_jets(n_jets);
            int size = css[i]->exclusive_jets(n_jets).size();
            auto idx = css[i]->particle_jet_indices(jets);
            int64_t sizz = css[i]->n_particles();
            auto prev = ptrjetoffsets[jetidx-1];

            for (unsigned int j = 0; j < jets.size(); j++){
              ptrjetoffsets[jetidx] = jets[j].constituents().size() + prev;
              prev = ptrjetoffsets[jetidx];
              jetidx++;
            }
            for(int k = 0; k < size; k++){  // iterate through jets in event
              for(int j = 0; j <sizz; j++){  // iterate through particles in event
                if(idx[j] == k){  // if particle jet index matches jet index assign it to jet j
                  ptrid[idxh] = j;
                  idxh++;
                }
              }
            }
            ptreventoffsets[eventidx] = jets.size()+eventprev;
            eventprev = ptreventoffsets[eventidx];
            eventidx++;
          }
        return std::make_tuple(
            jetoffsets,
            parid,
            eventoffsets
          );
      }, "n_jets"_a = 0, R"pbdoc(
        Retrieves the constituents of exclusive jets upto n jets from multievent clustering and converts them to numpy arrays.
        Args:
          n_jets: Number of exclusive subjets. Default: 0.
        Returns:
          jet offsets, particle indices, and event offsets
      )pbdoc")
      .def("to_numpy_exclusive_dcut",
      [](const output_wrapper ow, const double dcut = 100) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->exclusive_jets(dcut).size();
        }

        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->exclusive_jets(dcut);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, "dcut"_a = 100, R"pbdoc(
        Retrieves the exclusive jets upto the given dcut from multievent clustering and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_ycut",
      [](const output_wrapper ow, const double ycut = 100) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->exclusive_jets_ycut(ycut).size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->exclusive_jets_ycut(ycut);
        for (unsigned int j = 0; j < jets.size(); j++){
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, "dcut"_a = 100, R"pbdoc(
        Retrieves the exclusive jets upto the given dcut from multievent clustering and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_dmerge",
      [](const output_wrapper ow, int njets = 0) {
        auto css = ow.cse;
        int64_t len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->exclusive_dmerge(njets);
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, "njets"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_dmerge_max",
      [](const output_wrapper ow, int njets = 0) {
        auto css = ow.cse;
        auto len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->exclusive_dmerge_max(njets);
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, "njets"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_ymerge_max",
      [](const output_wrapper ow, int njets = 0) {
        auto css = ow.cse;
        auto len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->exclusive_ymerge_max(njets);
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, "njets"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_ymerge",
      [](const output_wrapper ow, int njets = 0) {
        auto css = ow.cse;
        auto len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->exclusive_ymerge(njets);
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, "njets"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_q",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        auto len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->Q();
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_q2",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        auto len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {len}, {sizeof(double)}));
        auto bufparid = parid.request();
        double *ptrid = (double *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->Q2();
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_subjets_dcut",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          double dcut = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto jk = 0;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->inclusive_jets();
        jk += css[i]->exclusive_subjets(jets[indices[i]],dcut).size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto jets = css[i]->exclusive_subjets(incjets[indices[i]],dcut);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_subjets_nsub",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          int nsub = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto jk = 0;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->inclusive_jets();
        jk += css[i]->exclusive_subjets(jets[indices[i]],nsub).size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto jets = css[i]->exclusive_subjets(incjets[indices[i]],nsub);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_subjets_up_to",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          int nsub = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto jk = 0;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->inclusive_jets();
        jk += css[i]->exclusive_subjets_up_to(jets[indices[i]],nsub).size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto jets = css[i]->exclusive_subjets_up_to(incjets[indices[i]],nsub);
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_subdmerge",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          int nsub = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {dimpy}, {sizeof(double)}));
        auto bufpx = out_value.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto value = css[i]->exclusive_subdmerge(incjets[indices[i]],nsub);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_subdmerge_max",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          int nsub = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {dimpy}, {sizeof(double)}));
        auto bufpx = out_value.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto value = css[i]->exclusive_subdmerge_max(incjets[indices[i]],nsub);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_n_exclusive_subjets",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei,
          double dcut = 0
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {dimpy}, {sizeof(int)}));
        auto bufpx = out_value.request();
        int *ptrpx = (int *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto value = css[i]->n_exclusive_subjets(incjets[indices[i]],dcut);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_has_parents",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(bool), py::format_descriptor<bool>::value, 1, {dimpy}, {sizeof(bool)}));
        auto bufpx = out_value.request();
        bool *ptrpx = (bool *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        fj::PseudoJet pj2(0,0,0,0);
        auto value = css[i]->has_parents(incjets[indices[i]],pj1, pj2);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Tells whether the given jet has parents or not.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_has_child",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(bool), py::format_descriptor<bool>::value, 1, {dimpy}, {sizeof(bool)}));
        auto bufpx = out_value.request();
        bool *ptrpx = (bool *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        auto value = css[i]->has_child(incjets[indices[i]],pj1);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Tells whether the given jet has children or not.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
    .def("to_numpy_jet_scale_for_algorithm",
      [](
          const output_wrapper ow,
          py::array_t<double,
          py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        int dimpy = infopy.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }
        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        auto out_value = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {dimpy}, {sizeof(double)}));
        auto bufpx = out_value.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;

        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        auto value = css[i]->jet_scale_for_algorithm(incjets[indices[i]]);
        ptrpx[idxe] = value;
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            out_value,
            off
          );
      }, R"pbdoc(
        Retrieves the exclusive subjets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_unique_history_order",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        int64_t len = css.size();
        int jk = 0;
        for(unsigned int i = 0; i<len; i++){

          jk += css[i]->unique_history_order().size();
        }
        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;
        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        ptreventoffsets[eventidx] = 0;
        eventidx++;
        size_t idxh = 0;
        auto eventprev = 0;
        for (unsigned int i = 0; i < css.size(); i++){
        auto info= css[i]->unique_history_order();
        for(unsigned int j =0; j < info.size(); j++){
        ptrid[idxh] = info[j];
        idxh++;}
        ptreventoffsets[eventidx] = info.size()+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_n_particles",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        int64_t len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->n_particles();
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, R"pbdoc(
        Gets n_particles.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_n_exclusive_jets",
      [](const output_wrapper ow, double dcut) {
        auto css = ow.cse;
        int64_t len = css.size();

        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
        size_t idxh = 0;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){
        ptrid[idxh] = css[i]->n_exclusive_jets(dcut);
        idxh++;
        ptreventoffsets[eventidx] = 1+eventprev;
        eventprev = ptreventoffsets[eventidx];
        eventidx++;
          }
        return std::make_tuple(
            parid,
            eventoffsets
          );
      }, R"pbdoc(
        Gets n_exclusive_jets.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_exclusive_njet_lund_declusterings",
      [](const output_wrapper ow, const int n_jets = 0) {
        auto css = ow.cse;
        int64_t len = css.size();
        auto jk = 0;

        for(int i = 0; i < len; i++){
          jk += css[i]->exclusive_jets(n_jets).size();
        }
        jk++;

        auto lund_generator = fastjet::contrib::LundGenerator();
        std::vector<double> Delta_vec;
        std::vector<double> kt_vec;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;

        ptreventoffsets[eventidx] = 0;
        eventidx++;

        auto jetoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufjetoffsets = jetoffsets.request();
        int *ptrjetoffsets = (int *)bufjetoffsets.ptr;
        size_t jetidx = 0;

        size_t idxh = 0;
        ptrjetoffsets[jetidx] = 0;
        jetidx++;
        auto eventprev = 0;

        for (unsigned int i = 0; i < css.size(); i++){  // iterate through events
          auto jets = css[i]->exclusive_jets(n_jets);
          int size = css[i]->exclusive_jets(n_jets).size();
          auto prev = ptrjetoffsets[jetidx-1];

          for (unsigned int j = 0; j < jets.size(); j++){
            auto lund_result = lund_generator.result(jets[j]);
            auto splittings = lund_result.size();
            for (unsigned int k = 0; k < splittings; k++){
              Delta_vec.push_back(lund_result[k].Delta());
              kt_vec.push_back(lund_result[k].kt());
            }

            ptrjetoffsets[jetidx] = splittings + prev;
            prev = ptrjetoffsets[jetidx];
            jetidx++;
          }

          ptreventoffsets[eventidx] = jets.size() + eventprev;
          eventprev = ptreventoffsets[eventidx];
          eventidx++;
        }

        auto Deltas = py::array(Delta_vec.size(), Delta_vec.data());
        auto kts = py::array(kt_vec.size(), kt_vec.data());

        return std::make_tuple(
            jetoffsets,
            Deltas,
            kts,
            eventoffsets
          );
      }, "n_jets"_a = 0, R"pbdoc(
        Calculates the Lund declustering Delta and k_T parameters from exclusive n_jets and converts them to numpy arrays.
        Args:
          n_jets: Number of exclusive subjets. Default: 0.
        Returns:
          jet offsets, splitting Deltas, kts, and event offsets.
      )pbdoc")
      .def("to_numpy_unclustered_particles",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->unclustered_particles().size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->unclustered_particles();
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the unclustered particles from multievent clustering and converts them to numpy arrays.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_childless_pseudojets",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->childless_pseudojets().size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->childless_pseudojets();
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the childless pseudojets from multievent clustering and converts them to numpy arrays.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_jets",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        jk += css[i]->jets().size();
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto jets = ow.cse[i]->jets();
        for (unsigned int j = 0; j < jets.size(); j++)
        {
          ptrpx[idxe] = jets[j].px();
          ptrpy[idxe] = jets[j].py();
          ptrpz[idxe] = jets[j].pz();
          ptrE[idxe] = jets[j].E();
          idxe++;
        }
        *ptroff = jets.size()+*(ptroff-1);
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the childless pseudojets from multievent clustering and converts them to numpy arrays.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_get_parents",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }

        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        fj::PseudoJet pj2(0,0,0,0);
        auto value = css[i]->has_parents(incjets[indices[i]],pj1, pj2);
        if(value == true){
        jk += 2;}
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        fj::PseudoJet pj2(0,0,0,0);
        auto value = css[i]->has_parents(incjets[indices[i]],pj1, pj2);
        if(value == true){
        ptrpx[idxe] = pj1.px();
        ptrpy[idxe] = pj1.py();
        ptrpz[idxe] = pj1.pz();
        ptrE[idxe] = pj1.E();
        idxe++;
        ptrpx[idxe] = pj2.px();
        ptrpy[idxe] = pj2.py();
        ptrpz[idxe] = pj2.pz();
        ptrE[idxe] = pj2.E();
        idxe++;
        *ptroff = 2+ *(ptroff-1);
        }
        else{
          *ptroff = *(ptroff-1);
        }
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the unclustered particles from multievent clustering and converts them to numpy arrays.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
    .def("to_numpy_get_child",
      [](
          const output_wrapper ow,
          py::array_t<double, py::array::c_style | py::array::forcecast> pxi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pyi,
          py::array_t<double, py::array::c_style | py::array::forcecast> pzi,
          py::array_t<double, py::array::c_style | py::array::forcecast> Ei
        ) {
        py::buffer_info infopx = pxi.request();
        py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
        py::buffer_info infopz = pzi.request();
        py::buffer_info infoE = Ei.request();

        auto pxptr = static_cast<double *>(infopx.ptr);
        auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
        auto pzptr = static_cast<double *>(infopz.ptr);
        auto Eptr = static_cast<double *>(infoE.ptr);

        int dimpx = infopx.shape[0];
        auto css = ow.cse;
        int64_t len = css.size();
        // Don't specify the size if using push_back.

        std::vector<fj::PseudoJet> particles;
        for(int j = 0; j < dimpx; j++ ){
          particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
          pxptr++;
          pyptr++;
          pzptr++;
          Eptr++;
          }

        std::vector<int> indices;
        for(unsigned int i = 0 ; i < len; i++){
          std::unordered_map<double, int> umap;
          auto jets = ow.cse[i]->inclusive_jets();
          for(unsigned int j = 0 ; j < jets.size(); j++){
            umap.insert({jets[j].rap(),j});
          }
          auto got = umap.find(particles[i].rap());
          if (got == umap.end()){
              throw "Jet Not in this ClusterSequence";
          }
          if(got == umap.end()){
          }
          indices.push_back(got->second);
        }
        // Don't specify the size if using push_back.
        auto jk = 0;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        auto value = css[i]->has_child(incjets[indices[i]],pj1);
        if(value == true){
        jk += 1;}
        }
        auto px = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpx = px.request();
        double *ptrpx = (double *)bufpx.ptr;

        auto py = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpy = py.request();
        double *ptrpy = (double *)bufpy.ptr;

        auto pz = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufpz = pz.request();
        double *ptrpz = (double *)bufpz.ptr;

        auto E = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
        auto bufE = E.request();
        double *ptrE = (double *)bufE.ptr;

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len+1}, {sizeof(int)}));
        auto bufoff = off.request();
        int *ptroff = (int *)bufoff.ptr;
        size_t idxe = 0;
        *ptroff = 0;
        ptroff++;
        for(int i = 0; i < len; i++){
        auto incjets = ow.cse[i]->inclusive_jets();
        fj::PseudoJet pj1(0,0,0,0);
        auto value = css[i]->has_child(incjets[indices[i]],pj1);
        if(value == true){
        ptrpx[idxe] = pj1.px();
        ptrpy[idxe] = pj1.py();
        ptrpz[idxe] = pj1.pz();
        ptrE[idxe] = pj1.E();
        idxe++;
        *ptroff = 1+ *(ptroff-1);
        }
        else{
          *ptroff = *(ptroff-1);
        }
        ptroff++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E,
            off
          );
      }, R"pbdoc(
        Retrieves the unclustered particles from multievent clustering and converts them to numpy arrays.
        Args:
          None.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc");
  py::class_<ClusterSequence>(m, "ClusterSequence")
      .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &,
                    const bool &>(),
           "pseudojets"_a, "jet_definition"_a,
           "write_out_combinations"_a = false,
           "Create a ClusterSequence, starting from the supplied set of "
           "PseudoJets and clustering them with jet definition specified by "
           "jet_definition (which also specifies the clustering strategy)");
}
