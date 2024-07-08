[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1003173.svg)](https://doi.org/10.5281/zenodo.1003173)

# StructureAnalysis

This package provides tools for statistical analysis of the structure of a plant.

> warning::
> This package in under development, even the master branch is broken

## Install

### Dependencies
We recommand to use python virtual environments to install the package (conda is a good choice, especialy since conda-only packages are required as dependencies).

Required dependencies available from `conda-forge` channel are :
- boost
- scons
- openalea.deploy
- openalea.sconsx
- nose (for tests, only via pip install)

### Doing the development install
To install the package for developments, you can use the `multisetup.py` script. It will use the `Multisetup` from `openalea.deploy` to compile `stat_tools` and `sequence_analysis` packages.

```bash
python multisetup.py develop
```

However, compilation can be a bit long and you may want to install using parallel compilation. To do so, you can use the following sequence (compilation with 4 cores, but you can adapt to your ressources):

```bash
cd stat_tools
scons -j4
pip install -e .
cd ../sequence_analysis
scons -j4
pip install -e .
cd ..
```

### After-install issues
Observed only on OSX (so far), you may need to extend your `LD_LIBRARY_PATH` environment variable to include the path of your newly built libraries:
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/your/conda/env/:/path/to/StructureAnalysis/stat_tool/:/path/to/StructureAnalysis/sequence_analysis/
```