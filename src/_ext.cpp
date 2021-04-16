// BSD 3-Clause License; see https://github.com/scikit-hep/fastjet/blob/main/LICENSE

#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include <fastjet/ClusterSequence.hh>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace fj = fastjet;

namespace py = pybind11;

void interface(){
  std::vector<fj::PseudoJet> particles;
  // an event with three particles:    px    py  pz      E
  particles.push_back(fj::PseudoJet( 99.0,  0.1,  0, 100.0));
  particles.push_back(fj::PseudoJet(  4.0, -0.1,  0,   5.0));
  particles.push_back(fj::PseudoJet(-99.0,    0,  0,  99.0));

  // choose a jet definition
  double R = 0.8;
  fj::JetDefinition jet_def(fj::antikt_algorithm, R);

  // run the clustering, extract the jets
  fj::ClusterSequence cs(particles, jet_def);
  std::vector<fj::PseudoJet> jets = fj::sorted_by_pt(cs.inclusive_jets());

  // print out some infos
  std::cout << "Clustering with " << jet_def.description() << std::endl;

  // print the jets
  std::cout <<   "        pt y phi" << std::endl;
  for (unsigned int i = 0;  i < jets.size();  i++) {
    std::cout << "jet " << i << ": " << jets[i].pt() << " "
              << jets[i].rap() << " " << jets[i].phi() << std::endl;
    std::vector<fj::PseudoJet> constituents = jets[i].constituents();
    for (unsigned int j = 0;  j < constituents.size();  j++) {
      std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                << std::endl;
    }
  }
}

PYBIND11_MODULE(_ext, m) {
    m.def("interface", &interface);
}
