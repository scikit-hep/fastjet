The Array Oriented Interface
============================

The tutorial on this page describes how the user can use `Awkward Arrays <https://awkward-array.org/quickstart.html>`__  to perform clustering on the particle data.

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
* There are other algorithms mentioned do not work.

The Data
---------
The input data for the Multi-event interface has to be an Awkward Array. One such example is as follows: ::

	>>> import awkward as ak
	>>> array = ak.Array(
        ... [
        ... 	{"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
        ... 	{"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
        ... 	{"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ... ],
    	... )

The Awkward Array here is a Record Array of Lorentz Vectors.

.. note::
   The inputs can be provided in more Awkward Array formats than described here.


ClusterSequence Class
----------------------

After defining the JetDefinition class, the user can provide this instance to the ClusterSequence class as an argument, along with the input data to perform the clustering: ::

	>>> cluster = fastjet.ClusterSequence(array, jetdef)
           <fastjet._pyjet.AwkwardClusterSequence object at 0x7f1413120a90>


Extracting Information
-----------------------
Any output that has to be an Array will be an Awkward Array in the array oriented interface. For example: ::

	>>> cluster.inclusive_jets()
	   <Array [{px: 1.2, py: 3.2, ... E: 48.7}] type='2 * Momentum4D["px": float64, "py...'>

Limitations
-----------
The Awkward Array interface is only available for the fastjet.ClusterSequence class. The Awkward Array functionality is likely to be expanded to other classes in the future.
