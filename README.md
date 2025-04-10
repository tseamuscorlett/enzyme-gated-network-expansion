# Enzyme-gated network expansion

Jupyter notebooks for running and analyzing data from enzyme-gated network expansion.

Simulations were run using 'folds_multiprocessing.py' and 'RUN_foldgated.ipynb'

Analysis files (.ipynb) include dependencies on assets from 'data' directory which can be downloaded at XYZ.

Phyletic distribution of ECOD X-groups calculated using code here: https://github.com/tseamuscorlett/ECOD_domain_phyletic_distribution

Enzyme-gated network expansion was run using the python package below.


# networkExpansionPy

Python package to construct biosphere-level metabolic networks, and run network expansion algorithms.  This package contains functions to prune biochemical reactions based on thermodynamic constraints.

The package is built from algorithms described in the following papers:

Goldford J.E. et al, Remnants of an ancient metabolism without phosphate. Cell 168, 1–9, March 9, 2017

Goldford J.E. et al, Environmental boundary conditions for the origin of life converge to an organo-sulfur metabolism. Nature Ecol Evo 3,12 1715-1724, November 11, 2019

## Installation

In a conda or virtual environment, clone git repo and install using pip.

```sh
git clone https://github.com/jgoldford/networkExpansionPy.git
cd networkExpansionPy
pip install -e .
```
