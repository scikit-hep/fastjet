The Object Oriented Interface
==============================

Clustering the data
--------------------

The fastjet library provides many options for the user to perform clustering on HEP data. The library has been designed keeping in mind the different requirements of users. The basic clustering process is described below.


Clustering Specification
-------------------------

The fastjet library has some clases specifically made to provide the different parameters for clustering. This includes the following classes :

* `JetDefinition <http://fastjet.fr/repo/doxygen-3.4.1/classfastjet_1_1JetDefinition.html>`__
* `AreaDefinition <http://fastjet.fr/repo/doxygen-3.4.1/classfastjet_1_1AreaDefinition.html>`__
* `RangeDefinition <http://fastjet.fr/repo/doxygen-3.4.1/classfastjet_1_1RangeDefinition.html>`__

For example, the JetDefinition class can be instantiated in the following way: ::

	import fastjet
	jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)

The JetDefinition class takes varied number of arguments, the first argument is always the type of algorithm, the number of rest of the arguments depends on how many parameters the given algorithm requires.

The JetAlgorithms
----------------------
The JetDefinition class takes `JetAlgorithms <http://fastjet.fr/repo/doxygen-3.4.1/namespacefastjet.html#a6377b557cbb936d4046d2aa936170dc0>`__  as arguments. In the above example we have chosen the Anti-kt algorithm. The list of algorithms is as following:

* ``kt_algorithm`` : The longitudinally invariant :math:`k_t` algorithm with distance measures

  .. math::
	\begin{align*}
	d_{ij} &= \min(p_{Ti}^2,p_{Tj}^2) \frac{\Delta R_{ij}^2}{R^2}, \\
	d_{iB} &= p_{Ti}^2.
	\end{align*}
* ``cambridge_algorithm`` : The longitudinally invariant variant of the cambridge algorithm (aka Aachen algoithm) with distance measures

  .. math::
	\begin{align*}
	d_{ij} &= \frac{\Delta R_{ij}^2}{R^2}, \\
	d_{iB} &= 1.
	\end{align*}
* ``antikt_algorithm`` : Like the :math:`k_t` but with distance measures

  .. math::
	\begin{align*}
	d_{ij} &= \min(p_{Ti}^{-2},p_{Tj}^{-2}) \frac{\Delta R_{ij}^2}{R^2}, \\
	d_{iB} &= p_{Ti}^{-2}.
	\end{align*}
* ``genkt_algorithm`` : The generalized :math:`k_t` algorithm but with distance measures

  .. math::
	\begin{align*}
    d_{ij} &= \min(p_{Ti}^{2p},p_{Tj}^{2p}) \frac{\Delta R_{ij}^2}{R^2}, \\
    d_{iB} &= p_{Ti}^{2p}
	\end{align*}

  where ``p = extra_param()``. Special cases are for :math:`p=1` (gives the :math:`k_t` algorithm), :math:`p=0` (gives the Cambridge/Aachen algorithm), and :math:`p=-1` (gives the anti-:math:`k_t` algorithm),.
* ``ee_kt_algorithm`` : The :math:`e^+e^-` :math:`k_t` algorithm (also known as Durham algoithm) with distance measure

  .. math::
    d_{ij} = 2\min(E_{i}^{2},E_{j}^{2}) (1-\cos\theta_{ij}).

* ``ee_genkt_algorithm`` : The :math:`e^+e^-` ``genkt`` algorithm with distance measures

  .. math::
	\begin{align*}
    d_{ij} &= \min(E_{i}^{2p},E_{j}^{2p}) \frac{1-\cos\theta_{ij}}{1-\cos R}, \\
    d_{iB} &= E_{i}^{2p}
	\end{align*}

  For :math:`\pi < R < 3\pi` and :math:`p=1` it gives the ``ee_kt`` algorithm.
* ``cambridge_for_passive_algorithm`` : A version of cambridge with a special distance measure for particles with ``pt < extra_param()``. This is not usually intended for end users, but is instead automatically selected when requesting a passive Cambridge area.
* Other algorithms that are not mentioned are most likely not implemented yet.

For a more complete description of the algorithms, please refer to the `fastjet manual <http://fastjet.fr/repo/fastjet-doc-3.4.1.pdf>`__.

The Data
--------

The input for the classic interface is a list of PseudoJets. To use the classic interface here's what the data should look like (This is a single event interface, one function call can only process one event): ::

	>>> array = [
	... 	fastjet.PseudoJet(1.1,1.2,1.3,11.4),
	... 	fastjet.PseudoJet(2.1,2.2,2.3,12.4),
	... 	fastjet.PseudoJet(3.1,3.2,3.3,13.4)
	... ]


ClusterSequence Class
----------------------

After defining the JetDefinition class, the user can provide this instance to the ClusterSequence class as an argument, along with the input data to perform the clustering: ::

	>>> cluster = fastjet.ClusterSequence(array, jetdef)
	#--------------------------------------------------------------------------
	#                         FastJet release 3.4.1
	#                 M. Cacciari, G.P. Salam and G. Soyez
	#     A software package for jet finding and analysis at colliders
	#                           http://fastjet.fr
	#
	# Please cite EPJC72(2012)1896 [arXiv:1111.6097] if you use this package
	# for scientific work and optionally PLB641(2006)57 [hep-ph/0512210].
	#
	# FastJet is provided without warranty under the GNU GPL v2 or higher.
	# It uses T. Chan's closest pair algorithm, S. Fortune's Voronoi code,
	# CGAL and 3rd party plugin jet algorithms. See COPYING file for details.
	#--------------------------------------------------------------------------
	>>> cluster
	<fastjet._swig.ClusterSequence; proxy of <Swig Object of type 'fastjet::ClusterSequence *' at 0x11b15bc90> >


Extracting Information
-----------------------
Any output that has to be an Array will be a list of PseudoJets if it's particle data. For example: ::

	>>> inc_jets = cluster.inclusive_jets()
	>>> for elem in inc_jets:
	... 	print("px:", elem.px(),"py:", elem.py(),"pz:", elem.pz(),"E:", elem.E(),)
	px: 6.300000000000001 py: 6.6000000000000005 pz: 6.8999999999999995 E: 37.2
