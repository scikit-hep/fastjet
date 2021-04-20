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
using namespace std;
namespace py = pybind11;

std::map<string,vector<float>> interface(py::array_t<float> xs, double Rv, string algor ){
  
    py::buffer_info info = xs.request(); 	// requesting buffer information of the input
    auto ptr = static_cast<float *>(info.ptr);	// pointer to the initial value
    vector<int> dims;				// the dimensions of the input array
    int n = 1;
    vector<vector<float>> grid = {{1,2,3},{1,2,3}};
    grid.clear();
    vector <float> too;   

    for (auto r: info.shape) {
      dims.push_back(r);
      n *= r;					// total number of elements
    }
      
    for (int i = 0; i < dims[0]; i++) { 
        // Vector to store column elements 
        vector <float> too; 
  
        for (int j = 0; j < dims[1]; j++) { 
            too.push_back(*ptr); 
            ptr++; 
        } 
  
        // Pushing back above 1D vector 
        // to create the 2D vector 
        grid.push_back(too); 
    } 
   //for(int i = 0; i <dims[0]; i++){  //for debugging
    //for(int j = 0; j <dims[1]; j++){
    //	std::cout<<grid[i][j];
    //}
   //} 
  std::vector<fj::PseudoJet> particles;
  // an event with three particles:    px    py  pz      E
  for(int i = 0; i < dims[0]; i++){
  particles.push_back(fj::PseudoJet( grid[i][0],  grid[i][1],  grid[i][2], grid[i][3]));
  }
  for(int i = 0; i < dims[0]; i++){
  std::cout<<particles[i][0]<<" "<<particles[i][1]<<" "<<particles[i][2]<<" "<<particles[i][3]<<endl;
  }
  // choose a jet definition
  double R = Rv;
  string algo = algor;
  std::vector<fj::PseudoJet> jets;
  //fj::JetDefinition jet_def;
  std::vector<float> pt;
  std::vector<float> rap;
  std::vector<float> phi;
  std::vector<float> constpt;
  std::vector<float> constrap;
  std::vector<float> constphi;
  std::vector<float> idx;
  
  
  
  
  if(algo.compare("antikt_algorithm")==0){
  fj::JetDefinition jet_def(fj::antikt_algorithm, R);

  // run the clustering, extract the jets
  fj::ClusterSequence cs(particles, jet_def);
  
  jets = fj::sorted_by_pt(cs.inclusive_jets());
  std::cout << "Clustering with " << jet_def.description() << std::endl;
    std::cout <<   "        pt y phi" << std::endl;
  for (unsigned int i = 0;  i < jets.size();  i++) {
    std::cout << "jet " << i << ": " << jets[i].pt() << " "
              << jets[i].rap() << " " << jets[i].phi() << std::endl;
    pt.push_back(jets[i].pt());rap.push_back(jets[i].rap());phi.push_back(jets[i].phi());
    std::vector<fj::PseudoJet> constituents = jets[i].constituents();
    unsigned int j;
    for (j = 0;  j < constituents.size();  j++) {
      std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                << std::endl;
    constpt.push_back(constituents[j].pt());constrap.push_back(constituents[j].rap());constphi.push_back(constituents[j].phi());
    }
    if (idx.size() == 0){
    idx.push_back(j);}
    else{
    idx.push_back(j+idx[idx.size()-1]);
    }
  }
  }
  
  
  
  
  if(algo.compare("kt_algorithm")==0){
  
  fj::JetDefinition jet_def(fj::kt_algorithm, R);

  // run the clustering, extract the jets
  fj::ClusterSequence cs(particles, jet_def);
  
  jets = fj::sorted_by_pt(cs.inclusive_jets());
  std::cout << "Clustering with " << jet_def.description() << std::endl;
  }
  if(algo.compare("cambridge_algorithm")==0){
  
  fj::JetDefinition jet_def(fj::cambridge_algorithm, R);

  // run the clustering, extract the jets
  fj::ClusterSequence cs(particles, jet_def);
  
  jets = fj::sorted_by_pt(cs.inclusive_jets());
  std::cout << "Clustering with " << jet_def.description() << std::endl;
    std::cout <<   "        pt y phi" << std::endl;
  for (unsigned int i = 0;  i < jets.size();  i++) {
    std::cout << "jet " << i << ": " << jets[i].pt() << " "
              << jets[i].rap() << " " << jets[i].phi() << std::endl;
    pt.push_back(jets[i].pt());rap.push_back(jets[i].rap());phi.push_back(jets[i].phi());
    std::vector<fj::PseudoJet> constituents = jets[i].constituents();
    unsigned int j;
    for (j = 0;  j < constituents.size();  j++) {
      std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                << std::endl;
    constpt.push_back(constituents[j].pt());constrap.push_back(constituents[j].rap());constphi.push_back(constituents[j].phi());
    }
    if (idx.size() == 0){
    idx.push_back(j);}
    else{
    idx.push_back(j+idx[idx.size()-1]);
    }
  }
  }
  
  
  
  if(algo.compare("ee_kt_algorithm")==0){
  fj::JetDefinition jet_def(fj::ee_kt_algorithm);

  // run the clustering, extract the jets
  fj::ClusterSequence cs(particles, jet_def);
  
  jets = fj::sorted_by_pt(cs.inclusive_jets());
  std::cout << "Clustering with " << jet_def.description() << std::endl;
    std::cout <<   "        pt y phi" << std::endl;
  for (unsigned int i = 0;  i < jets.size();  i++) {
    std::cout << "jet " << i << ": " << jets[i].pt() << " "
              << jets[i].rap() << " " << jets[i].phi() << std::endl;
    pt.push_back(jets[i].pt());rap.push_back(jets[i].rap());phi.push_back(jets[i].phi());
    std::vector<fj::PseudoJet> constituents = jets[i].constituents();
    unsigned int j;
    for (j = 0;  j < constituents.size();  j++) {
      std::cout << "    constituents " << j << "'s pt: " << constituents[j].pt()
                << std::endl;
    constpt.push_back(constituents[j].pt());constrap.push_back(constituents[j].rap());constphi.push_back(constituents[j].phi());
    }
    if (idx.size() == 0){
    idx.push_back(j);}
    else{
    idx.push_back(j+idx[idx.size()-1]);
    }
  }
  }
  
  // print out some infos
  

  // print the jets

  std::map<string,vector<float>> out = {{"pt",pt},{"rap",rap},{"phi",phi},{"constpt",constpt},{"constrap",constrap},{"constphi",constphi},{"idx",idx}};
  return out;
}

PYBIND11_MODULE(_ext, m) {
    m.def("interface", &interface);
}
