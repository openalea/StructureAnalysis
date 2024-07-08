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

To install the package for developments, you can use the `multisetup.py` script. It will use the `Multisetup` from `openalea.deploy` to compile `stat_tool` and `sequence_analysis` packages.

```bash
python multisetup.py develop
```

However, compilation can be a bit long and you may want to install using parallel compilation. To do so, you can use the following sequence (compilation with 4 cores, but you can adapt to your ressources):

```bash
cd stat_tool
scons -j4
pip install -e .
cd ../sequence_analysis
scons -j4
pip install -e .
cd ..
```

### After-install issues

Observed only on OSX (so far), you may need to extend your `LD_LIBRARY_PATH` environment variable to include the path of your newly built libraries:

```bash
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:/path/to/your/conda/env/lib:/path/to/StructureAnalysis/stat_tool/build-scons/lib/:/path/to/StructureAnalysis/sequence_analysis/build-scons/lib/
```

## Wrapping C++ code

The C++ code is wrapped using `boost::python`. The wrapping is done using `Scons` and `openalea.sconsx` (thin additional layer on top of `Scons` for `OpenAlea` projects).

> warning::
> We still need to make sure that all the C++ code is properly wrapped, and no requirements are missing. Is running all the tests enough ?

## Code layout

As a first step, we focus on two packages: `stat_tool` and `sequence_analysis`. The other packages are either becoming irrelevant (based on original `AML` code) or are not yet implemented.

### stat_tool

You can find below the structure of the `stat_tool` package. The `Scosntruct` file is used to compile the C++ code. The `src/cpp` directory contains the C++ code, and the `src/wrapper` directory contains the code to wrap the C++ code into Python. `src/openalea/stat_tool` contains the Python code that uses the wrapped C++ code.

```bash
.
├── AUTHORS.txt
├── ChangeLog.txt
├── LICENSE.txt
├── README.txt
├── SConstruct
├── conda
│   └── meta.yaml
├── debian
│   ├── changelog
│   ├── compat
│   ├── control
│   ├── copyright
│   ├── rules
│   └── source
│       └── format
├── doc
│   ├── Doxyfile
│   ├── Makefile
│   ├── README
│   ├── _static
│   │   ├── fig7_1.png
│   │   └── fig7_2.png
│   ├── conf.py
│   ├── contents.rst
│   ├── cpp
│   │   └── doxygen.conf
│   ├── make.bat
│   ├── pyplots
│   │   └── example1.py
│   └── user
│       ├── admin.rst
│       ├── ...
│       └── visualea_oak_demo.rst
├── examples
│   ├── Model
│   │   ├── Markov
│   │   │   ├── belren1.hsc
│   │   │   ├── ...
│   │   │   └── wij1.hsc
│   │   └── align1.a
│   ├── Sample
│   │   ├── Histogram
│   │   │   ├── meri1.his
│   │   │   ├── ...
│   │   │   └── peup6.his
│   │   └── Sequences
│   │       ├── belren1.seq
│   │       ├── ...
│   │       └── wij1.seq
│   ├── exploratory.aml
│   └── exploratory.py
├── install_stool.sh
├── merge_info.txt
├── methodo.txt
├── setup.py
├── share
│   └── data
│       ├── chene_sessile.vec
│       ├── ...
│       └── peup6.his
├── src
│   ├── cpp
│   │   ├── SConscript
│   │   ├── SConscriptWIG
│   │   ├── categorical_process.cpp
│   │   ├── chain.cpp
│   │   ├── chain_algorithms.cpp
|   |   ├── ...
│   │   ├── chain_reestimation.h
│   │   ├── chain_reestimation.hpp
│   │   ├── vectors.cpp
│   │   └── vectors.h
│   ├── openalea
│   │   └── stat_tool
│   │       ├── __init__.py
│   │       ├── cluster.py
│   │       ├── ...
│   │       └── vectors.py
│   └── wrapper
│       ├── SConscript
│       ├── boost_python_aliases.h
│       ├── export_base.cpp
│       ├── export_base.h
|       ├── ...
│       ├── export_vectors.cpp
│       ├── export_vectors.h
│       ├── stat_tool_wrap.cpp
│       └── wrapper_util.h
└── test
    ├── README.md
    ├── aml
    │   ├── stat_tool_test.aml
    │   ├── stat_tool_test.ipynb
    │   └── stat_tool_test.py
    ├── cpp
    │   ├── SConscript
    │   ├── np_model.mix
    │   └── test_multivariate_mixture.cpp
    ├── data
    │   ├── angles_ahp6_10.seq
    |   ├── ...
    │   └── vectors2.vec
    ├── stat_tool_examples.py
    ├── stat_tool_sequence_analysis_class.txt
    ├── test_cluster.py
    ├── ...
    ├── test_vectors_functional.py
    └── tools.py
```

### sequence_analysis

You can find below the structure of the `sequence_analysis` package, which is similar to `stat_tool`. The `Scosntruct` file is used to compile the C++ code. The `src/cpp` directory contains the C++ code, and the `src/wrapper` directory contains the code to wrap the C++ code into Python. `src/openalea/sequint`and `src/openalea/sequence_analysis` contains the Python code that uses the wrapped C++ code.

```bash
.
├── AUTHORS.txt
├── ChangeLog.txt
├── LICENSE.txt
├── README.txt
├── SConstruct
├── conda
│   └── meta.yaml
├── debian
│   ├── changelog
│   ├── compat
│   ├── control
│   ├── copyright
│   ├── rules
│   └── source
│       └── format
├── doc
│   ├── Doxyfile
│   ├── Makefile
│   ├── conf.py
│   ├── contents.rst
│   ├── cpp
│   │   └── doxygen.conf
│   ├── make.bat
│   ├── pyplots
│   │   ├── example_oak_1.py
│   │   ├── ...
│   │   └── sequence_plot_values.py
│   └── user
│       ├── admin.rst
│       ├── autosum.rst
│       ├── index.rst
│       ├── overview.txt
│       ├── scripts
│       │   ├── example_oak.py
│       │   └── example_oak.rst
│       └── tutorial.rst
├── sequence_analysis.py
├── setup.py
├── share
│   ├── README.txt
│   └── data
│       ├── abri13.ren
│       ├── abricotier_suivi_11.seq
│       ├── belren1.hsc
│       ├── belren1.seq
│       ├── ...
│       ├── wij1.hsc
│       └── wij1.seq
├── src
│   ├── cpp
│   │   ├── SConscript
│   │   ├── alignment.cpp
│   │   ├── ...
│   │   └── vomc_distributions2.cpp
│   ├── openalea
│   │   ├── seqint
│   │   │   ├── README.txt
│   │   │   ├── __init__.py
│   │   │   ├── config.py
│   │   │   ├── ...
│   │   │   └── xl_io.py
│   │   └── sequence_analysis
│   │       ├── __init__.py
│   │       ├── _sequence_analysis.so
│   │       ├── compare.py
│   │       ├── ...
│   │       └── variable_order_markov.py
│   ├── sequence_analysis_wralea
│   │   ├── __init__.py
│   │   ├── __wralea__.py
│   │   ├── demo
│   │   │   ├── __init__.py
│   │   │   ├── change_point
│   │   │   │   ├── Demo_ChangePoint_stat_tool_wralea.py
│   │   │   │   ├── __init__.py
│   │   │   │   ├── angelique_internode_length.seq
│   │   │   │   ├── change_point_demo.aml
│   │   │   │   ├── icon.png
│   │   │   │   ├── pin_laricio_7x.seq
│   │   │   │   └── reinet1.hsc
│   │   │   └── stat_tool_tutorial
│   │   │       ├── __init__.py
│   │   │       ├── __wralea__.py
│   │   │       ├── angelique_internode_length.seq
│   │   │       ├── icon.png
│   │   │       ├── pin_laricio_7x.seq
│   │   │       └── reinet1.hsc
│   │   ├── icon.png
│   │   └── stat.py
│   └── wrapper
│       ├── SConscript
│       ├── boost_python_aliases.h
│       ├── csequence.cpp
│       ├── export_base.cpp
│       ├── export_base.h
│       ├── ...
│       ├── sequence_analysis_wrap.cpp
│       └── wrapper_util.h
├── test
│   ├── _test_variable_order_markov.py
│   ├── functional1.py
│   ├── functional2.py
│   ├── functional3.py
│   ├── test_add_absorbing_run.py
│   ├── ...
│   ├── test_transcode.py
│   └── tools.py
└── tutorials
    ├── Code
    │   ├── __init__.py
    │   ├── amlseq2R.py
    │   ├── matrix_plot.py
    │   └── python_dics2R.py
    ├── Utils
    │   ├── __init__.py
    │   ├── __pycache__
    │   │   └── __init__.cpython-310.pyc
    │   ├── dos2unix.py
    │   └── unix2dos.py
    ├── scratch
    ├── seq1v_5s_LR_init.hsmc
    ├── sequences.ipynb
    └── sim_v_5s_LR.hsmc
```
