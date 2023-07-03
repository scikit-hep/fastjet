The Array Oriented Interface
============================

The tutorial on this page describes how the user can use `Awkward Arrays <https://awkward-array.org/quickstart.html>`__  to perform clustering on the particle data.

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
---------
The input data for the Multi-event interface has to be an Awkward Array. One such example is as follows: ::

	>>> import awkward as ak
	>>> array1 = ak.Array(
    ...     [
    ...         {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 23.5},
    ...         {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 755.12},
    ...         {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 835.56},
    ...     ],
    ... )


The Awkward Array here is a Record Array of Lorentz Vectors. ``fastjet`` is able to handle record arrays automatically when the fields are named ``("px", "py", "pz", "E")``.
In order to pass arrays with ``("pt", "eta", "phi", "M")`` fields, the user has to use ``vector`` and register the behavior. For example: ::

	>>> import vector
	>>> vector.register_awkward()
	>>> array2 = ak.Array(
	...     [
	...         {"pt": 3.42, "eta": 1.24, "phi": 1.21, "M": 22.6},
	...         {"pt": 71.8, "eta": 2.72, "phi": 1.11, "M": 519},
	...         {"pt": 71.1, "eta": 2.73, "phi": 1.1, "M": 631},
	...     ],
	...     with_name="Momentum4D",
	... )


The ``with_name`` argument is used to specify the name of the vector that will be used to override the default behavior. The ``vector.register_awkward()`` function registers the behavior of ``Momentum4D`` in ``awkward``.
The previous example with ``("px", "py", "pz", "E")`` fields still works: ::

	>>> array3 = ak.Array(
	...     [
	...         {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 23.5},
	...         {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 755.12},
	...         {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 835.56},
	...     ],
	...     with_name="Momentum4D",
	... )


It is good practice to register the behavior of the vector before using it in ``fastjet``.

ClusterSequence Class
----------------------

After defining the JetDefinition class, the user can provide this instance to the ClusterSequence class as an argument, along with the input data to perform the clustering: ::

	>>> cluster = fastjet.ClusterSequence(array1, jetdef)
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
	<fastjet._pyjet.AwkwardClusterSequence object at 0x7f1413120a90>


Extracting Information
-----------------------
Any output that has to be an Array will be an Awkward Array in the array oriented interface. For example: ::

	>>> cluster.inclusive_jets()
	<MomentumArray4D [{px: 1.2, py: 3.2, pz: 5.4, ...}, ...] type='2 * Momentum...'>

Limitations
-----------
The Awkward Array interface is only available for the fastjet.ClusterSequence class. The Awkward Array functionality is likely to be expanded to other classes in the future.
