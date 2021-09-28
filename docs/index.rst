.. toctree::
   :caption: Tutorial
   :hidden:

   Awkward
   particle

.. toctree::
   :caption: References
   :hidden:

   pseudojet
   clustersequence

.. |br| raw:: html

    <br/>

.. image:: logo.svg
    :width: 300px
    :alt: fastjet
    :target: https://github.com/scikit-hep/fastjet


Fastjet is a library for performing Jet-Finding *within* the Scikit-HEP ecosystem.
The library includes the classic interface, and a new interface built to perform clustering on multi-event Awkward Array objects.

.. note::
   Any questions about the C++ library, the jet finding algorithms, etc. should be directed to `fastjet.fr <http://fastjet.fr>`__. This page is meant to document the Python interfaces.

Documentation
---------------

* Python interface - This site.
* `GitHub <https://github.com/scikit-hep/fastjet/>`__

Installation
-------------

fastjet can be installed from `pypi <https://pypi.org/project/fastjet/>`__ using pip: ::

	pip install fastjet

Most users will get a precompiled binary (wheel) for your operating system and Python version. If not, the above attempts to compile from source.

The Interfaces
---------------
The fastjet library provides many options for the user to perform clustering on HEP data. The library has been designed keeping in mind the different requirements of users.

The fastjet library contains two interfaces within:

* **The Awkward interface**
* **The Classic interface**

The Awkward interface is the new interface made to handle multi-event data, whereas the classic interface is the same as the `C++ library <http://fastjet.fr/>`__, designed to handle the data in a particle-at-a-time fashion. The tutorials are divided into two to explain how each of the interfaces work. Please take a look at the tutorial section to get started.



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
