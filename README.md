<img src="https://raw.githubusercontent.com/scikit-hep/fastjet/main/docs/logo.svg">

[![Actions Status][actions-badge]][actions-link]
[![cirrus-ci Status][cirrus-ci-badge]][cirrus-ci-link]
[![PyPI version][pypi-version]][pypi-link]
[![Conda-Forge][conda-badge]][conda-link]

[![PyPI platforms][pypi-platforms]][pypi-link]
[![GitHub Discussion][github-discussions-badge]][github-discussions-link]
[![Scikit-HEP][sk-badge]](https://scikit-hep.org/)



[actions-badge]:            https://github.com/scikit-hep/fastjet/workflows/CI/badge.svg
[actions-link]:             https://github.com/scikit-hep/fastjet/actions
[cirrus-ci-badge]:          https://api.cirrus-ci.com/github/scikit-hep/fastjet.svg?branch=main
[cirrus-ci-link]:           https://cirrus-ci.com/github/scikit-hep/fastjet
[conda-badge]:              https://img.shields.io/conda/vn/conda-forge/fastjet
[conda-link]:               https://github.com/conda-forge/fastjet-feedstock
[github-discussions-badge]: https://img.shields.io/static/v1?label=Discussions&message=Ask&color=blue&logo=github
[github-discussions-link]:  https://github.com/scikit-hep/fastjet/discussions
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
# Tutorial

For a tutorial please look at the tutorial section of readthedocs page of this package.

<br>
<p align = "center">
<a href = "https://fastjet.readthedocs.io">
<img src = "https://img.shields.io/badge/read-documentation-blue" width="300" height="50">
</a>
</p>
<br>

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
sudo apt-get update && sudo apt-get install -y libboost-dev libmpfr-dev libgmp-dev swig autoconf libtool
```

## Build and install

Then you can build it using the following command:
``` bash
python -m pip install '.[test]'
```
