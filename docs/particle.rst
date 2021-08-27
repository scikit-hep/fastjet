The Object Oriented Interface
==============================

Clustering the data
--------------------

The fastjet library provides many options for the user to perform clustering on HEP data. The library has been designed keeping in mind the different requirements of users. The basic clustering process is described below.


The interfaces
---------------
The fastjet library contains two interfaces within:

* **The Awkward interface**
* **The Classic interface**

The Awkward interface is the new interface made to handle multi-event data, whereas the classic interface is the same as the C++ library, designed to handle the data in a particle-at-a-time fashion. The tutorials on this page are for classes (and methods) which are **common** to both, the Awkward interface and the Classic interface.

Clustering specification
-------------------------

The fastjet library has some clases specifically made to provide the different parameters for clustering. This includes the following classes :

* JetDefinition
* AreaDefinition
* RangeDefinition

For example, the JetDefinition class can be instantiated in the following way: ::

	import fastjet
	jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)

The JetDefinition class takes varied number of arguments, the first argument is always the type of algorithm, the number of rest of the arguments depends on how many parameters the given algorithm requires.

The algorithm classes
----------------------
The JetDefinition class takes algorithm classes as arguments. In the above example we have chosen the Anti-kt algorithm. The list of algorithms is as following:

* ee_genkt_algorithm
* ee_kt_algorithm
* genkt_algorithm
* genkt_for_passive_algorithm
* kt_algorithm
* cambridge_for_passive_algorithm
* cambridge_algorithm
* cambridge_aachen_algorithm
* antikt_algorithm

ClusterSequence Class
----------------------

After defining the JetDefinition class, the user can provide this instance to the ClusterSequence class as an argument, along with the input data to perform the clustering: ::

	fastjet.ClusterSequence(inputs, jetdef)

.. note::
   The inputs can be provided in different ways depending on the user requirements, to know how to do that, please refer to interface specific pages of this documentation
