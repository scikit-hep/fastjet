// BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <unordered_map>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/AreaDefinition.hh>
#include <fastjet/GhostedAreaSpec.hh>

#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace fj = fastjet;
namespace py = pybind11;
using namespace pybind11::literals;

typedef struct{
  PyObject_HEAD
  void *ptr;
  void *ty;
  int own;
  PyObject *next;
} SwigPyObject;

template <typename T>
T swigtocpp(py::object obj) {  // unwraps python object to get the cpp pointer from the swig bindings
  auto upointer = obj.attr("this").ptr();
  auto swigpointer = reinterpret_cast<SwigPyObject*>(upointer);
  auto objpointervoid = swigpointer->ptr;
  auto objpointer = reinterpret_cast<T>(objpointervoid);
  return objpointer;
}
class output_wrapper{
  public:
  std::vector<std::shared_ptr<fj::ClusterSequence>> cse;
  std::vector<std::shared_ptr<std::vector<fj::PseudoJet>>> parts;

  std::shared_ptr<fj::ClusterSequence> getCluster(){
    auto a = cse[0];
    return a;
  }
  void setCluster(){}
};

fj::ClusterSequence interface(py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, py::object jetdef)
{
  // py::buffer_info infooff = offsets.request();
  py::buffer_info infopx = pxi.request();
  py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
  py::buffer_info infopz = pzi.request();
  py::buffer_info infoE = Ei.request();

  // auto offptr = static_cast<int *>(infooff.ptr);
  auto pxptr = static_cast<double *>(infopx.ptr);
  auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
  auto pzptr = static_cast<double *>(infopz.ptr);
  auto Eptr = static_cast<double *>(infoE.ptr);

  // int dimoff = infooff.shape[0];
  int dimpx = infopx.shape[0];
  int dimpy = infopy.shape[0];
  int dimpz = infopz.shape[0];
  int dimE = infoE.shape[0];

  std::vector<double> nevents;
  std::vector<double> offidx;
  std::vector<double> constphi;
  std::vector<double> idx;
  std::vector<double> idxo;

  std::vector<fj::PseudoJet> particles;
  for (int i = 0; i < dimpx; i++) {
    particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
    pxptr++;
    pyptr++;
    pzptr++;
    Eptr++;
    }
  std::vector<fj::PseudoJet> jets;
  auto jet_def = swigtocpp<fj::JetDefinition*>(jetdef);
  fj::ClusterSequence cs(particles, *jet_def);
  jets = fj::sorted_by_pt(cs.inclusive_jets());
  return cs;
}

output_wrapper interfacemulti(py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei,py::array_t<int, py::array::c_style | py::array::forcecast> offsets, py::object jetdef)
{
  py::buffer_info infooff = offsets.request();
  py::buffer_info infopx = pxi.request();
  py::buffer_info infopy = pyi.request();  // requesting buffer information of the input
  py::buffer_info infopz = pzi.request();
  py::buffer_info infoE = Ei.request();

  auto offptr = static_cast<int *>(infooff.ptr);
  auto pxptr = static_cast<double *>(infopx.ptr);
  auto pyptr = static_cast<double *>(infopy.ptr);  // pointer to the initial value
  auto pzptr = static_cast<double *>(infopz.ptr);
  auto Eptr = static_cast<double *>(infoE.ptr);

  int dimoff = infooff.shape[0];
  int dimpx = infopx.shape[0];
  int dimpy = infopy.shape[0];
  int dimpz = infopz.shape[0];
  int dimE = infoE.shape[0];
  output_wrapper ow;
  std::vector<double> nevents;
  std::vector<double> offidx;
  std::vector<double> constphi;
  std::vector<double> idx;
  std::vector<double> idxo;
  for (int i = 0; i < dimoff-1; i++) {
    std::vector<fj::PseudoJet> particles;
    for(int j = *offptr; j < *(offptr+1); j++ ){
    particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
    pxptr++;
    pyptr++;
    pzptr++;
    Eptr++;
    }

  std::vector<fj::PseudoJet> jets;
  auto jet_def = swigtocpp<fj::JetDefinition*>(jetdef);
  std::shared_ptr<std::vector<fj::PseudoJet>> pj = std::make_shared<std::vector<fj::PseudoJet>>(particles);
  std::shared_ptr<fastjet::ClusterSequence> cs = std::make_shared<fastjet::ClusterSequence>(*pj, *jet_def);
  auto j = cs->inclusive_jets();
  offptr++;
  ow.cse.push_back(cs);
  ow.parts.push_back(pj);
  }
  return ow;
}

template <typename T>
class IterableWrapper
{
private:
  T fIterable;
  size_t fIndex;
  size_t fSize;

public:
  IterableWrapper(T &iterable) : fIterable(iterable), fIndex(0), fSize(iterable.size()){};
  typename T::value_type &operator*() { return fIterable.at(fIndex); }
  IterableWrapper &operator++()
  {
    fIndex++;
    return *this;
  }
  const size_t index() const { return fIndex; }
  const size_t size() const { return fSize; }
};

/*template<typename T>
class IterableWrapper {
 private:
  T fIterable;
  typename T::iterator fIterator;
  typename T::iterator fEnd;
 public:
  IterableWrapper(T& iterable): fIterable(iterable), fIterator(iterable.begin()), fEnd(iterable.end()) {};
  typename T::value_type& operator*() { return *fIterator; }
  IterableWrapper& operator++() {
    ++fIterator; return *this;
  }
  typename T::iterator begin() const { return fIterable.begin(); }
  typename T::iterator end() const { return fEnd; }
};*/

class IterableWrapperSentinel
{
};

/**
 * Determines whether we reached the end of the iterator.
 *
 * @return True if we reached the end of the iterator.
 */
template <typename T>
bool operator==(const IterableWrapper<T> &it, const IterableWrapperSentinel &)
{
  return it.index() == it.size();
}

/*template<typename T>
bool operator==(const IterableWrapper<T> & it, const IterableWrapperSentinel &) {
  return it == it.end();
}*/

/**
 * Create PseudoJet objects from a numpy array of px, py, pz, E. Axis 0 is the number of particles,
 * while axis 1 must be the 4 parameters.
 *
 * Note: The aray is required to be c-style, which ensures that it works with other packages. For example,
 *       pandas caused a problem in some cases without that argument.
 *
 * @param[jets] Numpy input array.
 * @returns Vector of PseudoJets.
 */
std::vector<fastjet::PseudoJet> constructPseudojetsFromNumpy(const py::array_t<double, py::array::c_style | py::array::forcecast> &jets)
{
  // Retrieve array and relevant information
  py::buffer_info info = jets.request();
  // I'm not sure which one of these is better.
  //auto inputJets = static_cast<double *>(info.ptr);
  auto inputJets = jets.data();
  std::vector<fastjet::PseudoJet> outputJets;
  // This defines our numpy array shape.
  int nParticles = info.shape[0];
  int nParams = info.shape[1];

  // Validation.
  if (nParams != 4)
  {
    throw std::runtime_error("Number of params is not correct. Should be four per particle.");
  }
  // Convert the arrays
  for (size_t i = 0; i < nParticles; ++i)
  {

    outputJets.push_back(fastjet::PseudoJet(
        inputJets[i * nParams + 0], inputJets[i * nParams + 1],
        inputJets[i * nParams + 2], inputJets[i * nParams + 3]));
  }

  return outputJets;
}

struct JetFinderSettings
{
  fastjet::JetDefinition _jetDefinition;
  fastjet::AreaDefinition _areaDefinition;

  const fastjet::JetDefinition JetDefinition() { return _jetDefinition; }
  void SetJetDefinition(fastjet::JetDefinition &def) { _jetDefinition = def; }
  const fastjet::AreaDefinition AreaDefinition() { return _areaDefinition; }
  void SetAreaDefinition(fastjet::AreaDefinition &area) { _areaDefinition = area; }
};

PYBIND11_MODULE(_ext, m) {
  using namespace fastjet;
  m.def("interface", &interface, py::return_value_policy::take_ownership);
  m.def("interfacemulti", &interfacemulti, py::return_value_policy::take_ownership);
  /// Jet algorithm definitions
  py::enum_<JetAlgorithm>(m, "JetAlgorithm", py::arithmetic(), "Jet algorithms")
    .value("kt_algorithm", JetAlgorithm::kt_algorithm, "the longitudinally invariant kt algorithm")
    .value("cambridge_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("aachen_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("cambridge_aachen_algorithm", JetAlgorithm::cambridge_algorithm, "the longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).")
    .value("antikt_algorithm", JetAlgorithm::antikt_algorithm, "like the k_t but with distance measures dij = min(1/kti^2,1/ktj^2) Delta R_{ij}^2 / R^2 diB = 1/kti^2")
    .value("genkt_algorithm", JetAlgorithm::genkt_algorithm, "like the k_t but with distance measures dij = min(kti^{2p},ktj^{2p}) Delta R_{ij}^2 / R^2 diB = 1/kti^{2p} where p = extra_param()")
    .value("cambridge_for_passive_algorithm", JetAlgorithm::cambridge_for_passive_algorithm, "a version of cambridge with a special distance measure for particles whose pt is < extra_param(); this is not usually intended for end users, but is instead automatically selected when requesting a passive Cambridge area.")
    .value("genkt_for_passive_algorithm", JetAlgorithm::genkt_for_passive_algorithm, "a version of genkt with a special distance measure for particles whose pt is < extra_param() [relevant for passive areas when p<=0] ***** NB: THERE IS CURRENTLY NO IMPLEMENTATION FOR THIS ALG *******")
    .value("ee_kt_algorithm", JetAlgorithm::ee_kt_algorithm, "the e+e- kt algorithm")
    .value("ee_genkt_algorithm", JetAlgorithm::ee_genkt_algorithm, "the e+e- genkt algorithm (R > 2 and p=1 gives ee_kt)")
    .value("plugin_algorithm", JetAlgorithm::plugin_algorithm, "any plugin algorithm supplied by the user")
    .value("undefined_jet_algorithm", JetAlgorithm::undefined_jet_algorithm, "the value for the jet algorithm in a JetDefinition for which no algorithm has yet been defined")
    .export_values();

  // Recombination scheme definitions
  py::enum_<RecombinationScheme>(m, "RecombinationScheme", py::arithmetic(), "Recombination schemes")
    .value("E_scheme", RecombinationScheme::E_scheme, "summing the 4-momenta")
    .value("pt_scheme", RecombinationScheme::pt_scheme, "pt weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling E=| p|")
    .value("pt2_scheme", RecombinationScheme::pt2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling E=| p|")
    .value("Et_scheme", RecombinationScheme::Et_scheme, "pt weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling | p|->=E")
    .value("Et2_scheme", RecombinationScheme::Et2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) with preprocessing to make things massless by rescaling | p|->=E")
    .value("BIpt_scheme", RecombinationScheme::BIpt_scheme, "pt weighted recombination of y,phi (and summing of pt's), with no preprocessing")
    .value("BIpt2_scheme", RecombinationScheme::BIpt2_scheme, "pt^2 weighted recombination of y,phi (and summing of pt's) no preprocessing")
    .value("WTA_pt_scheme", RecombinationScheme::WTA_pt_scheme, "pt-based Winner-Takes-All (WTA) recombination: the result of the recombination has the rapidity, azimuth and mass of the PseudoJet with the larger pt, and a pt equal to the sum of the two pt's")
    .value("WTA_modp_scheme", RecombinationScheme::WTA_modp_scheme, "mod-p-based Winner-Takes-All (WTA) recombination: the result of the recombination gets the 3-vector direction and mass of the PseudoJet with the larger |3-momentum| (modp), and a |3-momentum| equal to the scalar sum of the two |3-momenta|.")
    .value("external_scheme", RecombinationScheme::external_scheme, "for the user's external scheme")
    .export_values();

  // Jet reconstruction strategy definitions
  py::enum_<Strategy>(m, "Strategy", "Jet reconstruction strategies")
    .value("N2MHTLazy9AntiKtSeparateGhosts", Strategy::N2MHTLazy9AntiKtSeparateGhosts, "Like N2MHTLazy9 in a number of respects, but does not calculate ghost-ghost distances and so does not carry out ghost-ghost recombination.")
    .value("N2MHTLazy9", Strategy::N2MHTLazy9, "Only looks into a neighbouring tile for a particle's nearest neighbour (NN) if that particle's in-tile NN is further than the distance to the edge of the neighbouring tile. Uses tiles of size R and a 3x3 tile grid around the particle.")
    .value("N2MHTLazy25", Strategy::N2MHTLazy25, " Similar to N2MHTLazy9, but uses tiles of size R/2 and a 5x5 tile grid around the particle.")
    .value("N2MinHeapTiled", Strategy::N2MinHeapTiled, "faster that N2Tiled above about 500 particles; differs from it by retainig the di(closest j) distances in a MinHeap (sort of priority queue) rather than a simple vector. ")
    .value("N2Tiled", Strategy::N2Tiled, "fastest from about 50..500")
    .value("N2Plain", Strategy::N2Plain, "fastest below 50")
    .value("Best", Strategy::Best, "automatic selection of the best (based on N), including the LazyTiled strategies that are new to FJ3.1")
    .value("NlnN", Strategy::NlnN, "best of the NlnN variants -- best overall for N>10^4. (Does not work for R>=2pi)")
    .value("NlnN3pi", Strategy::NlnN3pi, "legacy N ln N using 3pi coverage of cylinder. (Does not work for R>=2pi)")
    .value("NlnN4pi", Strategy::NlnN4pi, "legacy N ln N using 4pi coverage of cylinder")
    .value("NlnNCam4pi", Strategy::NlnNCam4pi, "Chan's closest pair method (in a variant with 4pi coverage), for use exclusively with the Cambridge algorithm. (Does not work for R>=2pi)")
    .value("NlnNCam2pi2R", Strategy::NlnNCam2pi2R, "Chan's closest pair method (in a variant with 2pi+2R coverage), for use exclusively with the Cambridge algorithm.  (Does not work for R>=2pi)")
    .value("NlnNCam", Strategy::NlnNCam, "Chan's closest pair method (in a variant with 2pi+minimal extra variant), for use exclusively with the Cambridge algorithm. (Does not work for R>=2pi)")
    .value("BestFJ30", Strategy::BestFJ30, "the automatic strategy choice that was being made in FJ 3.0 (restricted to strategies that were present in FJ 3.0)")
    .value("plugin_strategy", Strategy::plugin_strategy, "the plugin has been used...")
    .export_values();

  py::class_<output_wrapper>(m, "output_wrapper")
    .def_property("cse", &output_wrapper::getCluster,&output_wrapper::setCluster)
    .def("to_numpy",
      [](const output_wrapper ow, double min_pt = 0) {
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
        auto jk = 0;
        auto sizepar = 0;

        for(int i = 0; i < len; i++){
        jk += css[i]->inclusive_jets().size();
        sizepar += css[i]->n_particles();
        }
        jk++;

        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {sizepar}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;

        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;

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
        int size = css[i]->inclusive_jets().size();
        auto idx = css[i]->particle_jet_indices(jets);
        auto sizz = css[i]->n_particles();
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
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      .def("to_numpy_exclusive_dcut",
      [](const output_wrapper ow, const double dcut = 100) {
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
        auto jk = len;

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
        auto jk = len;

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
        auto jk = len;

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
        auto jk = len;

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
        auto jk = len;

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
        auto jk = len;

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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, double dcut = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, int nsub = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, int nsub = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, int nsub = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, int nsub = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, double dcut = 0) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei) {

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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
        int jk = 0;
        for(unsigned int i = 0; i<len; i++){

          jk += css[i]->unique_history_order().size();
        }
        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;
        auto eventoffsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
        auto bufeventoffsets = eventoffsets.request();
        int *ptreventoffsets = (int *)bufeventoffsets.ptr;
        size_t eventidx = 0;
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
        auto len = css.size();
        auto jk = len;

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
        auto len = css.size();
        auto jk = len;

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
      .def("to_numpy_unclustered_particles",
      [](const output_wrapper ow) {
        auto css = ow.cse;
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
        auto len = css.size();
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

        auto off = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {len}, {sizeof(int)}));
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei) {
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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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
      [](const output_wrapper ow, py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei) {
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
        int dimpz = infopz.shape[0];
        int dimE = infoE.shape[0];
        auto css = ow.cse;
        auto len = css.size();
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
  py::class_<JetDefinition>(m, "JetDefinition", "Jet definition")
    .def(py::init<JetAlgorithm, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    .def(py::init<JetAlgorithm, double, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "R"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    .def(py::init<JetAlgorithm, double, double, RecombinationScheme, Strategy>(), "jet_algorithm"_a, "R"_a, "extra_parameter"_a, "recombination_scheme"_a = E_scheme, "strategy"_a = Best)
    // Leave out implementation of an external recombiner because it's not entirely clear how we should handle that in the bindings.
    .def(py::init<JetAlgorithm, double, RecombinationScheme, Strategy, int>(), "jet_algorithm"_a, "R"_a, "recombination_scheme"_a, "strategy"_a, "n_parameters"_a)
    .def_property("recombination_scheme", &JetDefinition::recombination_scheme, &JetDefinition::set_recombination_scheme, "The recombination scheme.")
    .def_property("jet_algorithm", &JetDefinition::jet_algorithm, &JetDefinition::set_jet_algorithm, "The jet algorithm.")
    .def_property_readonly("R", &JetDefinition::R, "The jet resolution parameter.")
    .def_property("extra_parameter", &JetDefinition::extra_param, &JetDefinition::set_extra_param, "A general purpose extra parameter, whose meaning depends on the algorithm, and may often be unused.")
    .def_property_readonly("strategy", &JetDefinition::strategy, "Jet finding strategy")
    .def("description", [](JetDefinition & jet, const bool include_recombiner) { if (include_recombiner == true) { return jet.description(); } return jet.description_no_recombiner(); }, "include_recombiner"_a = true, "A textual description of the current jet definition. The recombiner description is included by default.")
    .def_static("algorithm_description", &JetDefinition::algorithm_description, "jet_algorithm"_a, "A short textual description of the algorithm jet_algorithm")
    .def("__repr__", [](JetDefinition& jetDef){
        std::stringstream s;
        //s << "<JetDefinition jet_algorithm=" << jetDef.jet_algorithm() << ", R=" << jetDef.R() << " at " << &jetDef << ">";
        s << "<JetDefinition R=" << jetDef.R() << " at " << &jetDef << ">";
        return s.str();
      });

  py::class_<PseudoJet>(m, "PseudoJet")
    .def(py::init<double, double, double, double>(), "px"_a, "py"_a, "pz"_a, "E"_a)
    .def_property_readonly("e", &PseudoJet::E)
    // Alias for `e`
    .def_property_readonly("E", &PseudoJet::E)
    .def_property_readonly("px", &PseudoJet::px)
    .def_property_readonly("py", &PseudoJet::py)
    .def_property_readonly("pz", &PseudoJet::pz)
    .def_property_readonly("phi", &PseudoJet::phi, "phi (in the range 0..2pi)")
    .def_property_readonly("phi_std", &PseudoJet::phi_std, "phi in the range -pi..pi")
    .def_property_readonly("phi_02pi", &PseudoJet::phi_02pi, "phi in the range 0..2pi")
    .def_property_readonly("rap", &PseudoJet::rap, "the rapidity or some large value when the rapidity is infinite.")
    // Alias for `rap`
    .def_property_readonly("rapidity", &PseudoJet::rapidity, "the rapidity or some large value when the rapidity is infinite.")
    .def_property_readonly("pseudorapidity", &PseudoJet::pseudorapidity, "the pseudo-rapidity or some large value when the rapidity is infinite.")
    // Alias for `pseudorapidity`
    .def_property_readonly("eta", &PseudoJet::eta, "the pseudo-rapidity or some large value when the rapidity is infinite.")
    .def_property_readonly("pt", &PseudoJet::pt, "the scalar transverse momentum")
    .def_property_readonly("pt2", &PseudoJet::pt2, "the squared transverse momentum")
    .def_property_readonly("perp", &PseudoJet::perp, "the scalar transverse momentum")
    .def_property_readonly("perp2", &PseudoJet::perp2, "the squared transverse momentum")
    .def_property_readonly("m", &PseudoJet::m, "the squared invariant mass")
    .def_property_readonly("m2", &PseudoJet::m2, "the invariant mass")
    .def_property_readonly("mt", &PseudoJet::mt, "the squared transverse mass = kt^2+m^2")
    .def_property_readonly("mt2", &PseudoJet::mt2, "the transverse mass = sqrt(kt^2+m^2)")
    .def_property_readonly("mperp", &PseudoJet::mperp, "the squared transverse mass = kt^2+m^2")
    .def_property_readonly("mperp2", &PseudoJet::mperp2, "the transverse mass = sqrt(kt^2+m^2)")
    .def_property_readonly("modp", &PseudoJet::modp, "the 3-vector modulus = sqrt(px^2+py^2+pz^2)")
    .def_property_readonly("modp2", &PseudoJet::modp, "the squared 3-vector modulus = sqrt(px^2+py^2+pz^2)")
    .def_property_readonly("et", &PseudoJet::Et, "the transverse energy")
    .def_property_readonly("et2", &PseudoJet::Et2, "the squared transverse energy")
    // Intro ducted in a more recent versions of fastjet.
    /*.def_property_readonly("cos_theta", &PseudoJet::cos_theta, "cosine of the polar angle")
    .def_property_readonly("theta", &PseudoJet::theta, "the polar angle")*/
    .def_property_readonly("beam_distance", &PseudoJet::beam_distance, "distance between this jet and the beam")
    // Methods
    .def("delta_R", &PseudoJet::delta_R, "other"_a)
    .def("boost", &PseudoJet::boost, "rest_frame"_a, "transform this jet (given in the rest frame of the given PseudoJet) into a jet in the lab frame")
    .def("unboost", &PseudoJet::boost, "rest_frame"_a, " transform this jet (given in lab) into a jet in the rest frame of the given PseudoJet")
    // Operations
    .def(py::self + py::self, "other"_a, "Add the given PseudoJets.")
    .def(py::self - py::self, "other"_a, "Subtract the given PseudoJets.")
    .def(py::self += py::self, "other"_a, "Add the given PseudoJet.")
    .def(py::self -= py::self, "other"_a, "Subtract the given PseudoJet.")
    .def(py::self * double(), "val"_a, "Scale by the given value.")
    .def(py::self *= double(), "val"_a, "Scale by the given value.")
    .def(py::self / double(), "val"_a, "Divide by the given value.")
    .def(py::self /= double(), "val"_a, "Divide by the given value.")
    .def(py::self == py::self, "other"_a, "Compare the given PseudoJets. Returns true if the 4 momentum components of the two PseudoJets are identical and all the internal indices (user, cluster_history)")
    .def(py::self != py::self, "other"_a, "Compare the given PseudoJets. Returns true if they are not equal.")
    // Reset
    // Cannot use the templated version directly - we have to instantiate types explicitly.
    .def("reset", (void (PseudoJet::*)(double, double, double, double)) &PseudoJet::reset, "px"_a, "py"_a, "pz"_a, "E"_a, "Reset the PseudoJet to the properties of the given jet px, py, pz, E.")
    .def("reset", &PseudoJet::reset<std::vector<double>>, "jet_four_momentum"_a, "Reset the PseudoJet to the properties of the given jet px, py, py, E.")
    .def("reset", &PseudoJet::reset<PseudoJet>, "new_psueodjet"_a, "Reset the PseudoJet to the properties of the given jet.")
    .def("reset_momentum", (void (PseudoJet::*)(double, double, double, double)) &PseudoJet::reset_momentum, "px"_a, "py"_a, "pz"_a, "E"_a, "Reset the momentum PseudoJet to the momentum of the given jet.")
    .def("reset_momentum", (void (PseudoJet::*)(const PseudoJet &)) &PseudoJet::reset_momentum, "new_psueodjet"_a, "Reset the momentum PseudoJet to the momentum of the given jet.")
    // User index
    .def_property("user_index", &PseudoJet::user_index, &PseudoJet::set_user_index, "The user index allows the user to add simple identifying information to a particle/jet")
    // Area
    .def_property_readonly("area", &PseudoJet::area, "The jet (scalar) area.")
    .def_property_readonly("area_error", &PseudoJet::area_error, "The error (uncertainty) associated with the determination of the area of this jet.")
    .def_property_readonly("area_4vector", &PseudoJet::area_4vector, "The jet 4-vector area.")
    .def_property_readonly("is_pure_ghost", &PseudoJet::is_pure_ghost, "True if this jet is made exclusively of ghosts.")
    // Constituents
    .def_property_readonly("constituents", &PseudoJet::constituents, "The constituents.")
    .def("__getitem__", [](const PseudoJet &jet, size_t i) {
      auto constituents = jet.constituents();
      if (i >= constituents.size()) throw py::index_error();
      return constituents[i];
    }, "i"_a, R"pbdoc(
      Retrieve an individual constituent at a given index.
      Args:
        i: Index of the constituent to retrieve.
      Returns:
        Constituent at that index.
    )pbdoc")
    .def("__iter__",
      [](const PseudoJet &jet) {
        auto constituents = jet.constituents();
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(constituents), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), R"pbdoc(
        Iterates over the constituents in the PseudoJet.
      )pbdoc")
    // Cluster sequence related
    .def_property_readonly("has_associated_cluster_sequence", &PseudoJet::has_associated_cluster_sequence, "true if this PseudoJet has an associated ClusterSequence")
    .def_property_readonly("has_associated_cs", &PseudoJet::has_associated_cluster_sequence, "true if this PseudoJet has an associated ClusterSequence")
    .def_property_readonly("has_valid_cluster_sequence", &PseudoJet::has_valid_cluster_sequence, "true if this PseudoJet has an associated and still valid(ated) ClusterSequence.")
    .def_property_readonly("has_valid_cs", &PseudoJet::has_valid_cluster_sequence, "true if this PseudoJet has an associated and still valid(ated) ClusterSequence.")
    .def("has_partner", &PseudoJet::has_partner, "partner"_a, "Check if it has been recombined with another PseudoJet in which case, return its partner through the argument. Otherwise, 'partner' is set to 0.")
    .def("has_child", &PseudoJet::has_child, "child"_a, "Check if it has been recombined with another PseudoJet in which case, return its child through the argument. Otherwise, 'child' is set to 0.")
    .def("parents", [](PseudoJet & jet) {
        PseudoJet parent1, parent2;
        jet.has_parents(parent1, parent2);
        // We'd like to return None if the parents aren't available. However, we can't handle this nicely without c++17.
        // Instead, we just return both the parents. The user has to check...
        //return (parent1.pt() != 0 || parent2.pt() != 0) ? std::optional<std::tuple<PseudoJet, PseudoJet>>{std::make_tuple(parent1, parent2)} : py::cast<py::none>(Py_None);
        return std::make_tuple(parent1, parent2);
      }, "Return the parents of the PseudoJet if it is the product of a recombination. If not, returns two empty PseudoJets.")
    .def("contains", &PseudoJet::contains, "other"_a, " Check if the current PseudoJet contains the one passed as argument.")
    .def("is_inside", &PseudoJet::is_inside, "other"_a, "Check if the current PseudoJet is contained the one passed as argument.")
    .def("__repr__", [](PseudoJet & jet){
        std::stringstream s;
        //s << "<PseudoJet px=" << jet.px() << ", py=" << jet.py() << ", pz=" << jet.pz() << ", E=" << jet.E();
        s << "<PseudoJet pt=" << jet.pt() << ", eta=" << jet.eta() << ", phi=" << jet.phi() << ", m=" << jet.m() << " at " << &jet << ">";
        return s.str();
      });

  // Helper functions
  m.def("dot_product", &dot_product, "jet_1"_a, "jet_2"_a, "Returns the 4-vector dot product of a and b");
  m.def("have_same_momentum", &have_same_momentum, "jet_1"_a, "jet_2"_a, "Returns true if the momenta of the two input jets are identical");
  m.def("sorted_by_pt", &sorted_by_pt, "jets"_a, "Return a vector of jets sorted into decreasing transverse momentum");
  m.def("sorted_by_pz", &sorted_by_pz, "jets"_a, "Return a vector of jets sorted into increasing pz");
  m.def("sorted_by_rapidity", &sorted_by_rapidity, "jets"_a, "Return a vector of jets sorted into increasing rapidity");
  m.def("sorted_by_E", &sorted_by_E, "jets"_a, "Return a vector of jets sorted into decreasing energy");

  py::class_<ClusterSequence>(m, "ClusterSequence")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const bool &>(), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    // numpy constructor.
    .def(py::init([](const py::array_t<double> & pseudojets, const JetDefinition & jetDef, const bool & writeOutCombination){ auto convertedPseudojets = constructPseudojetsFromNumpy(pseudojets); return ClusterSequence(convertedPseudojets, jetDef, writeOutCombination); }), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    .def("inclusive_jets", &ClusterSequence::inclusive_jets, "pt_min"_a = 0., "Return a vector of all jets (in the sense of the inclusive algorithm) with pt >= ptmin. Time taken should be of the order of the number of jets returned.")
    .def("__getitem__", [](const ClusterSequence &cs, size_t i) {
      auto inclusive_jets = cs.inclusive_jets();
      if (i >= inclusive_jets.size()) throw py::index_error();
      return inclusive_jets[i];
    }, "i"_a, R"pbdoc(
      Retrieve an individual jet at a given index.
      Args:
        i: Index of the jet to retrieve.
      Returns:
        Jet at that index.
    )pbdoc")
    .def("__iter__",
      [](const ClusterSequence &cs) {
        auto jets = cs.inclusive_jets();
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), R"pbdoc(
        Finds the include jets and iterates over them.
      )pbdoc")
    .def("__call__",
      [](const ClusterSequence &cs, double min_pt = 0) {
        auto jets = cs.inclusive_jets(min_pt);
        return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
      }, py::keep_alive<0, 1>(), "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          List of inclusive jets.
      )pbdoc")
    .def("to_numpy",
      [](const ClusterSequence &cs, double min_pt = 0) {
        auto jets = cs.inclusive_jets(min_pt);
        // Don't specify the size if using push_back.
        auto jk = jets.size();
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
        size_t idxe = 0;
        for (unsigned int i = 0; i < jets.size(); i++)
        {

          ptrpx[idxe] = jets[i].px();
          ptrpy[idxe] = jets[i].py();
          ptrpz[idxe] = jets[i].pz();
          ptrE[idxe] = jets[i].E();
          idxe++;
        }
        return std::make_tuple(
            px,
            py,
            pz,
            E
          );
      }, "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc")
      .def("to_numpy_with_constituents",
      [](const ClusterSequence &cse, double min_pt = 0) {
        auto jets = cse.inclusive_jets(min_pt);
        // Don't specify the size if using push_back.
        auto jk = jets.size();
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
        size_t idxe = 0;

        auto offsets = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {jk}, {sizeof(int)}));
        auto bufoffsets = offsets.request();
        int *ptroffsets = (int *)bufoffsets.ptr;

        size_t off = 0;
        int prev = 0;
        for (unsigned int i = 0; i < jets.size(); i++)
        {
          ptrpx[idxe] = jets[i].px();
          ptrpy[idxe] = jets[i].py();
          ptrpz[idxe] = jets[i].pz();
          ptrE[idxe] = jets[i].E();
          idxe++;
          std::vector<fj::PseudoJet> constituents = jets[i].constituents();
          ptroffsets[off] = constituents.size() + prev;
          prev = ptroffsets[off];
          off++;
        }
        auto size = jets.size();
        auto sizepar = cse.n_particles();



        auto parid = py::array(py::buffer_info(nullptr, sizeof(int), py::format_descriptor<int>::value, 1, {sizepar}, {sizeof(int)}));
        auto bufparid = parid.request();
        int *ptrid = (int *)bufparid.ptr;
        auto idx = cse.particle_jet_indices(jets);

        size_t idxh = 0;
        for(int i = 0; i < size; i++){
          for(int j = 0; j <sizepar; j++){
            if(idx[j] == i){
              ptrid[idxh] = j;
              idxh++;
            }
          }
        }
        for(int j = 0; j <sizepar; j++){
            if(idx[j] == -1){
              ptrid[idxh] = j;
              idxh++;
            }
          }

        return std::make_tuple(
            px,
            py,
            pz,
            E,
            idx,
            parid,
            offsets
          );
      }, "min_pt"_a = 0, R"pbdoc(
        Retrieves the inclusive jets and converts them to numpy arrays.
        Args:
          min_pt: Minimum jet pt to include. Default: 0.
        Returns:
          pt, eta, phi, m of inclusive jets.
      )pbdoc");

  py::class_<ClusterSequenceArea, ClusterSequence>(m, "ClusterSequenceArea")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const AreaDefinition &>(), "pseudojets"_a, "jet_definition"_a, "area_definition"_a, "Create a ClusterSequenceArea, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
    .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const GhostedAreaSpec &>(), "pseudojets"_a, "jet_definition"_a, "ghost_spec"_a, "Create a ClusterSequenceArea, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)");

  py::enum_<AreaType>(m, "AreaType", py::arithmetic(), "Jet area definitions")
    .value("invalid_area", AreaType::invalid_area, "Invalid jet area")
    .value("active_area", AreaType::active_area, "Active jet area")
    .value("active_area_explicit_ghosts", AreaType::active_area_explicit_ghosts, "Active jet area with explicit ghosts")
    .value("one_ghost_passive_area", AreaType::one_ghost_passive_area, "One ghost with passive area")
    .value("passive_area", AreaType::passive_area, "Passive area")
    .value("voronoi_area", AreaType::voronoi_area, "Voronoi area")
    .export_values();

  py::class_<GhostedAreaSpec>(m, "GhostedAreaSpec", "Ghost area specification")
    .def(py::init<double, int, double, double, double, double>(),
        "max_rapidity"_a,
        "repeat_in"_a = gas::def_repeat,
        "ghost_area_in"_a = gas::def_ghost_area,
        "grid_scatter_in"_a = gas::def_grid_scatter,
        "pt_scatter_in"_a = gas::def_pt_scatter,
        "mean_ghost_pt_in"_a = gas::def_mean_ghost_pt)
    .def_property("max_rapidity", &GhostedAreaSpec::ghost_maxrap, &GhostedAreaSpec::set_ghost_maxrap, "Max ghost rapidity.")
    .def_property("ghost_area", &GhostedAreaSpec::ghost_area, &GhostedAreaSpec::set_ghost_area, "Ghost area.")
    .def_property("grid_scatter", &GhostedAreaSpec::grid_scatter, &GhostedAreaSpec::set_grid_scatter, "Grid scatter.")
    .def_property("pt_scatter", &GhostedAreaSpec::pt_scatter, &GhostedAreaSpec::set_pt_scatter, "pt scatter.")
    .def_property("mean_ghost_pt", &GhostedAreaSpec::mean_ghost_pt, &GhostedAreaSpec::set_mean_ghost_pt, "Mean ghost pt.");

  py::class_<AreaDefinition>(m, "AreaDefinition", "Area definition")
    .def(py::init<AreaType, const GhostedAreaSpec &>(), "area_type"_a, "ghost_spec"_a)
    .def(py::init<AreaType>(), "area_type"_a)
    .def_property_readonly("area_type", &AreaDefinition::area_type)
    .def_property_readonly("ghost_spec", (const GhostedAreaSpec & (AreaDefinition::*)() const) &AreaDefinition::ghost_spec)
    .def("__repr__", [](AreaDefinition& areaDef){
        std::stringstream s;
        //s << "<AreaDefinition area_type=" << areaDef.area_type() << ", ghost_spec=" << areaDef.ghost_spec() << " at " << &areaDef << ">";
        s << "<AreaDefinition at " << &areaDef << ">";
        return s.str();
      })
    ;

  // Awkward array
  // Ensure dependencies are loaded.
  //py::module::import("awkward");

  py::class_<JetFinderSettings>(m, "JetFinderSettings", "Encompasses jet finder settings")
    .def(py::init<const JetDefinition, const AreaDefinition>(), "jet_definition"_a, "area_definition"_a)
    .def_property("jet_definition", &JetFinderSettings::JetDefinition, &JetFinderSettings::SetJetDefinition)
    .def_property("area_definition", &JetFinderSettings::AreaDefinition, &JetFinderSettings::SetAreaDefinition)
    .def("__repr__", [](JetFinderSettings & settings){
        std::stringstream s;
        //s << "<JetFinderSettings jet_definition R=" << settings.JetDefinition() << ", area_definition=" << settings.AreaDefinition()  << " at " << &settings << ">";
        s << "<JetFinderSettings jet_def R=" << settings.JetDefinition().R() << " at " << &settings << ">";
        return s.str();
      })
    ;

}
