============================
Delft3D Flexible Mesh Python
============================


.. image:: https://img.shields.io/pypi/v/delft3dfmpy.svg
        :target: https://pypi.python.org/pypi/delft3dfmpy

.. image:: https://img.shields.io/travis/grongen/delft3dfmpy.svg
        :target: https://travis-ci.org/grongen/delft3dfmpy

.. image:: https://readthedocs.org/projects/delft3dfmpy/badge/?version=latest
        :target: https://delft3dfmpy.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status




Python package to generate delft3dfm models from standardized data models or other flow models.


* Free software: MIT license
* Documentation: https://dhydamo.readthedocs.io.
* WIKI (in Dutch): https://hkvconfluence.atlassian.net/wiki/spaces/DHYD/pages/222396421/Achtergrond

Features
--------

* Read Hydamo or shape files to Python data structure
* Build Delft3D FM model from Hydamo or seperately defined elements
* Generate rectangular meshes with refinement
* Create 1D 2D connections
* Create RR model schematisation
* Write to Delft3D FM model

Installation
------------
By the following procedure you can install the delft3dfmpy module:

1.  Install a Anaconda or Minoconda Python distribution:
        https://www.anaconda.com/products/individual
        https://docs.conda.io/en/latest/miniconda.html

    During installation, tick the box “Add Anaconda to PATH”, even though it colors a suggestive red.

2.  Clone or download the repository. This can be done with the "Code" button at the upper right of this page. Unpack the downloaded file if you've chosen to download the repository as .zip.

3.  Open a command prompt and navigate (with the "cd" command) to the directory were you've cloned or unpacked delft3dfmpy. This directory should contain amongst other the files environment.yml and setup.py. Create the delft3dfmpy environment with the required modules, by executing the following commands in an anaconda command prompt:

        conda env create -f environment.yml

3.  Activate the created environment by the following command in command prompt:

       conda activate delft3dfmpy

 Install the delft3dfmpy in the active environment by the following command in command prompt:

      pip install .

 Alternatively, you can link the directory with delft3dfmpy to your anaconda environment in develop-mode:
 
      pip install -e .


Activate the created environment in an anaconda command prompt (conda activate delft3dfmpy) before running your notebook or script. A jupyter notebook or command prompt for the environment can also be launched from the Anaconda Navigator. 
For more information on how to use environments, see: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

Usage
-----
The usage is best described by the example notebook: https://github.com/openearth/delft3dfmpy/blob/master/notebooks/Usage_introduction_coupled_RRFM.ipynb
