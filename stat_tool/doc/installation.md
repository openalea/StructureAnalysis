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
mamba env create -n stat -f conda/environment.yml
mamba activate stat_tool

# Clone stat_tool and install
git clone https://github.com/openalea/stat_tool.git
cd stat_tool
pip install --no-build-isolation -e .
# (Optional) Test your installation
cd test; pytest
```

Compilation options for developers such as WITH_TEST, WITH_EFENCE are defined in pyproject.toml. They can be used with

```bash
pip install --no-build-isolation --config-settings=cmake.define.WITH_TEST=TRUE -e .
pip install --no-build-isolation --config-settings=cmake.build-type="Debug" -e .
```
and also combined:
```bash
pip install --no-build-isolation --config-settings=cmake.define.WITH_TEST=TRUE --config-settings=cmake.define.WITH_EFENCE=TRUE --config-settings=cmake.build-type="Debug" -e .
```
