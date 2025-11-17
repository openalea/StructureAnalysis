# Installation

You must use conda environment : <https://docs.conda.io/en/latest/index.html>

## Users

### Create a new environment with stat_tool installed in there

```bash

mamba create -n stat_tool -c openalea3 -c conda-forge  openalea.stat_tool
mamba activate stat_tool
```

Install stat_tool in a existing environment

```bash
mamba install -c openalea3 -c conda-forge openalea.stat_tool
```

### (Optional) Test your installation

```bash
mamba install -c conda-forge pytest
git clone https://github.com/openalea/stat_tool.git
cd stat_tool/test; pytest
```

## Developers

### Install From source

```bash
# Install dependency with conda
mamba env create -n phm -f conda/environment.yml
mamba activate stat_tool

# Clone stat_tool and install
git clone https://github.com/openalea/stat_tool.git
cd stat_tool
pip install .

# (Optional) Test your installation
cd test; pytest
```
