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

Environment Preparation
^^^^^^^^^^^^^^^^^^^^^^^
Prepare an environment with the correct dependencies for delft3dfmpy
1.  Install a Anaconda or Minoconda Python distribution:

     https://www.anaconda.com/products/individual
     https://docs.conda.io/en/latest/miniconda.html

2.  Save the content of https://raw.githubusercontent.com/openearth/delft3dfmpy/master/environment.yml and store this in a local file named ``environment.yml``

3.  Open Ananconda prompt and enter the directory where the ``envrionment.yml`` from step 2 is stored.

4.  Install the ``delft3dfmpy`` environment with the required modules, by executing the following command in the opened command prompt.

        ``conda env create -f environment.yml``

    This should create a ``delft3dfmpy`` environment with the required dependencies.

5.  Activate the created environment by the following command in command prompt:

        ``conda activate delft3dfmpy``

    You now have a correct and activated environment for installation of the ``delft3dfmpy`` Python package

Python package ``delft3dfmpy``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1.  Install the Python package ``delft3dfmpy`` directly, by executing the following command in the opened command prompt.

        ``python -m pip install git+https://github.com/openearth/delft3dfmpy.git``

2.  If the command prompt states ``Successfully built delft3dfmpy`` then installation is succesful.

Usage
-----

Activate the created environment in an (anaconda) command prompt (conda activate delft3dfmpy) before running your notebook or script. A jupyter notebook or command prompt for the environment can also be launched from the Anaconda Navigator. 
For more information on how to use environments, see: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

The usage is best described by the example notebook: https://github.com/openearth/delft3dfmpy/blob/master/notebooks/Usage_introduction_coupled_RRFM.ipynb

Contribution
------------

Contributions are much welcome for the documentation, missing but desirable features and bugs. The code and development happens in this GitHub repository: https://github.com/openearth/delft3dfmpy.

For contributions, use the following guidelines:

1.  Fork the project on GitHub, clone the fork to your Operating System and open the repository as folder/workspace in your favorite IDE (eg. VSCode + Python extension).

2.  Make your contributions and test the changes locally.

3.  Once satisfied, push your changes as a new branch to your fork and create a Pull Request to the original repository.

4.  A maintainer on the main GitHub repository will review your PR and guide the merging process. 
