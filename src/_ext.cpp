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

py::dict interface(py::array_t<double, py::array::c_style | py::array::forcecast> pxi, py::array_t<double, py::array::c_style | py::array::forcecast> pyi, py::array_t<double, py::array::c_style | py::array::forcecast> pzi, py::array_t<double, py::array::c_style | py::array::forcecast> Ei, std::map<std::string,std::string> params, std::map<std::string,float> paramf)
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
      auto pt = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
      auto bufpt = pt.request();
      double *ptrpt = (double *)bufpt.ptr;

      auto phi = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
      auto bufphi = phi.request();
      double *ptrphi = (double *)bufphi.ptr;

      auto rap = py::array(py::buffer_info(nullptr, sizeof(double), py::format_descriptor<double>::value, 1, {jk}, {sizeof(double)}));
      auto bufrap = rap.request();
      double *ptrrap = (double *)bufrap.ptr;
      size_t idxe = 0;
      std::cout << "        pt y phi" << std::endl;
      for (unsigned int i = 0; i < jets.size(); i++) {
        std::cout << "jet " << i << ": " << jets[i].pt() << " "
                  << jets[i].rap() << " " << jets[i].phi() << std::endl;
        ptrpt[idxe] = jets[i].pt();
        ptrphi[idxe] = jets[i].phi();
        ptrrap[idxe] = jets[i].rap();
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
      out["part0-node1-data"] = pt; // throws exception error - ptyes.h Line 546
      out["part0-node2-data"] = rap;
      out["part0-node3-data"] = phi;

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
}
