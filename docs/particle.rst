The Object Oriented Interface
==============================

Clustering the data
--------------------

The fastjet library provides many options for the user to perform clustering on HEP data. The library has been designed keeping in mind the different requirements of users. The basic clustering process is described below.


Clustering Specification
-------------------------

The fastjet library has some clases specifically made to provide the different parameters for clustering. This includes the following classes :

* `JetDefinition <http://fastjet.fr/repo/doxygen-3.4.0/classfastjet_1_1JetDefinition.html>`__
* `AreaDefinition <http://fastjet.fr/repo/doxygen-3.4.0/classfastjet_1_1AreaDefinition.html>`__
* `RangeDefinition <http://fastjet.fr/repo/doxygen-3.4.0/classfastjet_1_1RangeDefinition.html>`__

For example, the JetDefinition class can be instantiated in the following way: ::

	import fastjet
	jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)

The JetDefinition class takes varied number of arguments, the first argument is always the type of algorithm, the number of rest of the arguments depends on how many parameters the given algorithm requires.

The JetAlgorithms
----------------------
The JetDefinition class takes `JetAlgorithms <http://fastjet.fr/repo/doxygen-3.4.0/namespacefastjet.html#a6377b557cbb936d4046d2aa936170dc0>`__  as arguments. In the above example we have chosen the Anti-kt algorithm. The list of algorithms is as following:

* ``ee_genkt_algorithm`` : The e+e- genkt algorithm (R > 2 and p=1 gives ee_kt)
* ``ee_kt_algorithm`` : The e+e- kt algorithm
* ``genkt_algorithm`` : Like the k_t but with distance measures dij = min(kti^{2p},ktj^{2p}) Delta R_{ij}^2 / R^2 diB = 1/kti^{2p} where p = extra_param()
* ``kt_algorithm`` : The longitudinally invariant kt algorithm
* ``cambridge_for_passive_algorithm`` : A version of cambridge with a special distance measure for particles whose pt is < extra_param(); This is not usually intended for end users, but is instead automatically selected when requesting a passive Cambridge area.
* ``cambridge_algorithm`` : The longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm).
* ``antikt_algorithm`` : Like the k_t but with distance measures dij = min(1/kti^2,1/ktj^2) Delta R_{ij}^2 / R^2 diB = 1/kti^2
* There are other algorithms mentioned in the link that do not work.

The Data
--------

The input for the classic interface is a list of PseudoJets. To use the classic interface here's what the data should look like (This is a single event interface, one function call can only process one event): ::

	>>> array = [fastjet.PseudoJet(1.1,1.2,1.3,1.4),
	... fastjet.PseudoJet(2.1,2.2,2.3,2.4),
	... fastjet.PseudoJet(3.1,3.2,3.3,3.4)]


ClusterSequence Class
----------------------

After defining the JetDefinition class, the user can provide this instance to the ClusterSequence class as an argument, along with the input data to perform the clustering: ::

	fastjet.ClusterSequence(inputs, jetdef)


Extracting Information
-----------------------
Any output that has to be an Array will be a list of PseudoJets if it's particle data. For example: ::

	>>> inc_jets = cluster.inclusive_jets()
	>>> for elem in inc_jets:
        ...     print("px:", elem.px(),"py:", elem.py(),"pz:", elem.pz(),"E:", elem.E(),)
        px: 6.300000000000001 py: 6.6000000000000005 pz: 6.8999999999999995 E: 7.199999999999999
