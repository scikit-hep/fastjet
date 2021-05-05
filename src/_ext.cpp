// BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

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

py::dict interface(py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, std::map<std::string,std::string> params, std::map<std::string,float> paramf)//, py::object jetdef)
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

  //auto aa = jetdef.attr("description").cast<std::string>();
  //std::cout << aa << std::endl;

  ///auto roo;
  // vector<int> dims;  // the dimensions of the input array
  // int n = 1;
  // vector<vector<float>> grid = {{1,2,3},{1,2,3}};
  // grid.clear();
  // vector <float> too;

  // for (auto r: info.shape) {
  // dims.push_back(r);
  // n *= r;					// total number of elements
  // }

  // for (int i = 0; i < dims[0]; i++) {
  // Vector to store column elements
  // vector <float> too;
  // for (int j = 0; j < dims[1]; j++) {
  // too.push_back(*ptr);
  // ptr++;

  // }
  // Pushing back above 1D vector
  // to create the 2D vector
  // grid.push_back(too);
  // }
  // for(int i = 0; i <dims[0]; i++){  //for debugging
  // for(int j = 0; j <dims[1]; j++){
  // std::cout<<grid[i][j];
  // }
  // }
  auto algo =  (std::string)params["algor"];
  auto R = paramf["R"];
  std::vector<double> nevents;
  std::vector<double> offidx;
  std::vector<double> constphi;
  std::vector<double> idx;
  std::vector<double> idxo;

 // nevents.push_back(dimoff - 1);
  //offidx.push_back(0);
  //idx.push_back(0);
  //for (int k = 0; k < dimoff - 1; k++)
  //{
    std::vector<fj::PseudoJet> particles;
    // an event with three particles:    px    py  pz      E
    //if(*offptr == *(offptr+1)){
      //idx.push_back(idx[idx.size() - 1]);
      //offptr++;
      //continue;
    //}
    for (int i = 0; i < dimpx; i++) {
      particles.push_back(fj::PseudoJet(*pxptr, *pyptr, *pzptr, *Eptr));
      pxptr++;
      pyptr++;
      pzptr++;
      Eptr++;
    }
    //offptr++;


    // for(int i = *(offptr-1); i < *offptr; i++){
    // std::cout<<particles[i].px()<<" "<<particles[i].py()<<" "<<particles[i].pz()<<" "<<particles[i].E()<<endl;
    // }
    std::vector<fj::PseudoJet> jets;

    if (algo.compare("anti-kt") == 0) {
      fj::JetDefinition jet_def(fj::antikt_algorithm, R);

      // run the clustering, extract the jets
      fj::ClusterSequence cs(particles, jet_def);
      //std::vector<int> input = cs.unique_history_order();
      //cout << "----------------------------------------------------------" << //endl;
      // for (int i = 0; i < input.size(); i++) {
      // std::cout << input.at(i) << ' ';
      // }
      //cout << endl;
      //cout << "----------------------------------------------------------" << endl;
      jets = fj::sorted_by_pt(cs.inclusive_jets());
      std::cout << "Clustering with " << jet_def.description() << std::endl;
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
      std::cout << "        pt y phi" << std::endl;
      for (unsigned int i = 0; i < jets.size(); i++) {
        std::cout << "jet " << i << ": " << jets[i].pt() << " "
                  << jets[i].rap() << " " << jets[i].phi() << std::endl;
        ptrpx[idxe] = jets[i].px();
        ptrpy[idxe] = jets[i].py();
        ptrpz[idxe] = jets[i].pz();
        ptrE[idxe] = jets[i].E();
        idxe++;
        //std::vector<fj::PseudoJet> constituents = jets[i].constituents();
        //unsigned int j;
        //std::cout << constituents.size() << endl;
        //for (j = 0; j < constituents.size(); j++) {
        //std::cout << "    constituent " << j << " " << constituents[j].px() << constituents[j].py() << constituents[j].pz() << constituents[j].E() << std::endl;
        //constpt.push_back(constituents[j].pt());
        //constrap.push_back(constituents[j].rap());
        //constphi.push_back(constituents[j].phi());
        //for (auto k = 0; k < particles.size(); k++) {
        //if (constituents[j].px() == particles[k].px() && constituents[j].py() == particles[k].py() && constituents[j].pz() == particles[k].pz() && constituents[j].E() == particles[k].E()) {
        //idxo.push_back(k);
        //std::cout << endl;
        //std::cout << "---------------------------------------------------------------" << endl;
        //std::cout << k;
        //std::cout << endl;
        //std::cout << "---------------------------------------------------------------" << endl;
        //}
        //}
        //}
        //offidx.push_back(idxo.size());
        //if (idx.size() == 0)
        //{
        //idx.push_back(j);
        //}
        //else
        //{
          //idx.push_back(j + idx[idx.size() - 1]);
        //}
      }
      py::dict out;
      out["px"] = px; // throws exception error - ptyes.h Line 546
      out["py"] = py;
      out["pz"] = pz;
      out["E"] = E;

      //std::map<std::string, py::array> out = {{"part0-node1-data", pt},{"part0-node2-data", rap},{"part0-node3-data", phi}};
      return out;
    }


  if (algo.compare("cambridge_algorithm") == 0) {
      fj::JetDefinition jet_def(fj::cambridge_algorithm, R);

      // run the clustering, extract the jets
      fj::ClusterSequence cs(particles, jet_def);

      jets = cs.inclusive_jets();
      std::cout << "Clustering with " << jet_def.description() << std::endl;
      std::cout << "        pt y phi" << std::endl;
      for (unsigned int i = 0; i < jets.size(); i++) {
        std::cout << "jet " << i << ": " << jets[i].pt() << " "
                  << jets[i].rap() << " " << jets[i].phi() << std::endl;
        std::vector<fj::PseudoJet> constituents = jets[i].constituents();
        unsigned int j;
        for (j = 0; j < constituents.size(); j++) {
          std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                    << std::endl;
          //constpt.push_back(constituents[j].pt());
          //constrap.push_back(constituents[j].rap());
          //constphi.push_back(constituents[j].phi());
        }
        if (idx.size() == 0) {
          idx.push_back(j);
        } else {
          idx.push_back(j + idx[idx.size() - 1]);
        }
      }
    }

    if (algo.compare("ee_kt_algorithm") == 0) {
      fj::JetDefinition jet_def(fj::ee_kt_algorithm);

      // run the clustering, extract the jets
      fj::ClusterSequence cs(particles, jet_def);

      jets = fj::sorted_by_pt(cs.inclusive_jets());
      std::cout << "Clustering with " << jet_def.description() << std::endl;
      std::cout << "        pt y phi" << std::endl;
      for (unsigned int i = 0; i < jets.size(); i++) {
        std::cout << "jet " << i << ": " << jets[i].pt() << " "
                  << jets[i].rap() << " " << jets[i].phi() << std::endl;
        std::vector<fj::PseudoJet> constituents = jets[i].constituents();
        unsigned int j;
        for (j = 0; j < constituents.size(); j++) {
          std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                    << std::endl;
          //constpt.push_back(constituents[j].pt());
          //constrap.push_back(constituents[j].rap());
          //constphi.push_back(constituents[j].phi());
        }
        if (idx.size() == 0) {
          idx.push_back(j);
        } else {
          idx.push_back(j + idx[idx.size() - 1]);
        }
      }
    }
  //}

  // print out some infos

  // print the jets

}

PYBIND11_MODULE(_ext, m) {
  m.def("interface", &interface);
  // py::class_<ClusterSequence>(m, "ClusterSequence")
  //     .def(py::init<const std::vector<PseudoJet> &, const JetDefinition &, const bool &>(), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
  //     // numpy constructor.
  //     .def(py::init([](const py::array_t<double> &pseudojets, const JetDefinition &jetDef, const bool &writeOutCombination) { auto convertedPseudojets = constructPseudojetsFromNumpy(pseudojets); return ClusterSequence(convertedPseudojets, jetDef, writeOutCombination); }), "pseudojets"_a, "jet_definition"_a, "write_out_combinations"_a = false, "Create a ClusterSequence, starting from the supplied set of PseudoJets and clustering them with jet definition specified by jet_definition (which also specifies the clustering strategy)")
  //     .def("inclusive_jets", &ClusterSequence::inclusive_jets, "pt_min"_a = 0., "Return a vector of all jets (in the sense of the inclusive algorithm) with pt >= ptmin. Time taken should be of the order of the number of jets returned.")
  //     .def(
  //         "__getitem__", [](const ClusterSequence &cs, size_t i) {
  //           auto inclusive_jets = cs.inclusive_jets();
  //           if (i >= inclusive_jets.size())
  //             throw py::index_error();
  //           return inclusive_jets[i];
  //         },
  //         "i"_a, R"pbdoc(
  //     Retrieve an individual jet at a given index.
  //     Args:
  //       i: Index of the jet to retrieve.
  //     Returns:
  //       Jet at that index.
  //   )pbdoc")
  //     .def(
  //         "__iter__", [](const ClusterSequence &cs) {
  //           auto jets = cs.inclusive_jets();
  //           return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
  //         },
  //         py::keep_alive<0, 1>(), R"pbdoc(
  //       Finds the include jets and iterates over them.
  //     )pbdoc")
  //     .def(
  //         "__call__", [](const ClusterSequence &cs, double min_pt = 0) {
  //           auto jets = cs.inclusive_jets(min_pt);
  //           return py::make_iterator(IterableWrapper<std::vector<PseudoJet>>(jets), IterableWrapperSentinel());
  //         },
  //         py::keep_alive<0, 1>(), "min_pt"_a = 0, R"pbdoc(
  //       Retrieves the inclusive jets.
  //       Args:
  //         min_pt: Minimum jet pt to include. Default: 0.
  //       Returns:
  //         List of inclusive jets.
  //     )pbdoc")
  //     .def(
  //         "to_numpy", [](const fj::ClusterSequence &cs, double min_pt = 0) {
  //           auto jets = cs.inclusive_jets(min_pt);
  //           // Don't specify the size if using push_back.
  //           std::vector<double> pt, eta, phi, m;
  //           for (const auto &jet : jets)
  //           {
  //             pt.push_back(jet.pt());
  //             eta.push_back(jet.eta());
  //             phi.push_back(jet.phi());
  //             m.push_back(jet.m());
  //           }
  //           return std::make_tuple(
  //               py::array(py::cast(pt)),
  //               py::array(py::cast(eta)),
  //               py::array(py::cast(phi)),
  //               py::array(py::cast(m)));
  //         },
  //         "min_pt"_a = 0, R"pbdoc(
  //       Retrieves the inclusive jets and converts them to numpy arrays.
  //       Args:
  //         min_pt: Minimum jet pt to include. Default: 0.
  //       Returns:
  //         pt, eta, phi, m of inclusive jets.
  //     )pbdoc");
}
