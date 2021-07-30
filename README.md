# fastjet

[![Actions Status][actions-badge]][actions-link]
[![Documentation Status][rtd-badge]][rtd-link]
[![Code style: black][black-badge]][black-link]

[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]
[![PyPI platforms][pypi-platforms]][pypi-link]

[![GitHub Discussion][github-discussions-badge]][github-discussions-link]
[![Gitter][gitter-badge]][gitter-link]
[![Scikit-HEP][sk-badge]](https://scikit-hep.org/)



[actions-badge]:            https://github.com/scikit-hep/fastjet/workflows/CI/badge.svg
[actions-link]:             https://github.com/scikit-hep/fastjet/actions
[black-badge]:              https://img.shields.io/badge/code%20style-black-000000.svg
[black-link]:               https://github.com/psf/black
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/fastjet
[conda-link]:               https://github.com/conda-forge/fastjet-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/aryan26roy/fastjet.git/discussions
[gitter-badge]:             https://badges.gitter.im/https://github.com/aryan26roy/fastjet.git/community.svg
[gitter-link]:              https://gitter.im/https://github.com/aryan26roy/fastjet.git/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
[pypi-link]:                https://pypi.org/project/fastjet/
[pypi-platforms]:           https://img.shields.io/pypi/pyversions/fastjet
[pypi-version]:             https://badge.fury.io/py/fastjet.svg
[rtd-badge]:                https://readthedocs.org/projects/fastjet/badge/?version=latest
[rtd-link]:                 https://fastjet.readthedocs.io/en/latest/?badge=latest
[sk-badge]:                 https://scikit-hep.org/assets/images/Scikit--HEP-Project-blue.svg

Official FastJet bindings to Python and Awkward Array.

Main features of Fastjet:
  * Contains Vectorized as well as Non-Vectorized interface for Fastjet.
  * Compiled against the complete Fastjet library in C++.
  * Has Awkward Array and Vector as dependency.
  * Provides the functionality to cluster multiple events at a time.
  * Input data can be in any coordinate system.

# Installation
The package can be installed from PyPI using the following command:
``` bash
python -m pip install fastjet
```
# Overview

Some of the basic functionalities of Fastjet and how to use them are listed below.

``` python
import fastjet
import awkward as ak
import vector
```
The input data can be either a awkward array or a list of Pseudojets.

## Awkward Array
```python
input_data = ak.Array(
    [
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
        [
            {"px": 1.2, "py": 3.2, "pz": 5.4, "E": 2.5, "ex": 0.78},
            {"px": 32.2, "py": 64.21, "pz": 543.34, "E": 24.12, "ex": 0.35},
            {"px": 32.45, "py": 63.21, "pz": 543.14, "E": 24.56, "ex": 0.0},
        ],
    ],
    with_name="Momentum4D",
)
```
## List of PseudoJets
```Python
input_data = [fastjet.PseudoJet(1,1,1,1)
             ,fastjet.PseudoJet(1.2,1.2,1.2,1.2)
             ,fastjet.PseudoJet(3,3,3,3)
             ,fastjet.PseudoJet(-1,-12,2,1)
             ,fastjet.PseudoJet(-1,-12,2.1,0.9)]
```
## Clustering
The classes (clustering and clustering specification) are the same for all the input types:

```python
jetdef = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.6)
cluster = fastjet.ClusterSequence(input_data, jetdef)
```
## Outputs
The output can be extracted using function calls (output can be an Awkward Array or a list of PseudoJets depending on the input).

```python
inclusive_jets = cluster.inclusive_jets()
```
# Installation For Developers
Clone this repository recursively to get the dependencies.

```
git clone --recursive https://github.com/scikit-hep/fastjet.git
```

## Build dependencies

There are still external build-time dependencies for the C++ components of `fastjet` that need to be installed on the build machine.
To install the build-time dependencies run the following installation commands for your respective operating system:

### Debian/Ubuntu

``` bash
sudo apt-get update && sudo apt-get install -y libboost-dev gfortran swig autoconf libtool
```

## Build and install

Then you can build it using the following command:
``` bash
python -m pip install .[test]
```
