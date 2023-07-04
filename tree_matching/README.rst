=======
TreeMatching
=======


**TreeMatching** is an open-source toolkit for the comparison of tree graph structure.


=============
Installation
=============


``TreeMatching`` distribution is based on the ``conda`` software environment management system.
To install conda, you may refer to its installation page: https://docs.conda.io/projects/conda/en/latest/user-guide/install/

To install TreeMatching, you need to create an environment (named for instance pgl) :

.. code-block:: bash

        conda create -n treematching openalea.treematching -c openalea3 -c conda-forge

The package ``openalea.treematching`` is retrieved from the ``openalea3`` channel (developement) and its dependencies will be taken from ``conda-forge`` channel.

Then, you need to activate the treematching environment

.. code-block:: bash

        conda activate treematching


Or use the treematching modules in Python

.. code-block:: bash

        ipython


.. code-block:: python

        >>> from openalea.tree_matching import *



=============
Compiling
=============


The simplest way to build TreeMatching is to use conda (see below).

Then, setup your Conda environment with all required dependencies :

.. code:: bash

    # Linux or macOS
    conda env create -f environment.yaml

Now, you can build, then install TreeMatching :

.. code:: bash

    cd tree_matching
    mkdir build
    cd build

    # Linux
    cmake .. -DCMAKE_INSTALL_PREFIX=${CONDA_PREFIX}

    # Windows -> Visual Studio 2015 is required
    cmake .. -G "NMake Makefiles" -DCMAKE_INSTALL_PREFIX=%LIBRARY_PREFIX%  ..

    cmake --build . --target install --config Release

    cd ..
    python setup.py install --prefix=${CONDA_PREFIX}

You're done !


Help and Support
----------------

Please open an **Issue** if you need support or that you run into any error (Installation, Runtime, etc.).
We'll try to resolve it as soon as possible.

==============
Authors
==============

TreeMatching was developed by Pascal Ferraro, Aida Ouangraoua, Christophe Godin and Frédéric Boudon.

