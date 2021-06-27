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



[actions-badge]:            https://github.com/aryan26roy/fastjet.git/workflows/CI/badge.svg
[actions-link]:             https://github.com/aryan26roy/fastjet.git/actions
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

Main features of Fastjet :
  * Contains Vectorized as well as Non-Vectorized interface for Fastjet.
  * Compiled against the complete Fastjet library in C++.
  * Has Awkward Array and Vector as dependency.
  * Provides the functionality to cluster multiple events at a time.
  * Input data can be in any coordinate system.

# Installation
The package can be installed from pypi using the following command :
``` bash
pip install fastjet
```
# Overview

Some of the basic functionalities of Fastjet and how to use them are listed below.

``` python
import fastjet
import awkward as ak
import vector
```


