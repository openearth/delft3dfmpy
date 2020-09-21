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

2.  Save the content of https://raw.githubusercontent.com/openearth/delft3dfmpy/master/environment.yml and store this in a local file named environment.yml

3.  from this repository and create a delft3dfmpy environment with the required modules, by executing the following command in the opened command prompt.

     conda env create -f environment.yml

     This should create a delft3dfmpy environment with the required dependencies.

3.  Activate the created environment by the following command in command prompt:

     conda activate delft3dfmpy

4.  Clone this github directory by the following command in the command prompt:

     git clone https://github.com/openearth/delft3dfmpy.git

5.  Make sure your enviornment is still active and move to the delft3dfm folder and install the delft3dfmpy in the active environment by the following commands:

     cd delft3dfmpy

     pip install .

    if you installed Anaconda for the current user only, you can install delft3dfmpy by:

    python -m pip install --user .

    Or in case you want to further develop replace the last command by:

     pip install -e .

    or:

     python -m pip install --user -e .

    Including the keyword `-e` makes Python recognize the module during import in your script: `import delft3dfmpy` without moving it to the site-packages directory. This means that you can edit the code from the directory where you cloned or unzipped the repository, and the changes will automatically be recognized on the next import. If you do not plan to modify the code you can leave the keyword `-e` out. In this case pip will install the package in the site-packages directory, as any other pip or conda installed package. You can still make changes to the code, but in order for them to be picked up they need to be done in the site-packages/delft3dfmpy directory

Activate the created environment in an (anaconda) command prompt (conda activate delft3dfmpy) before running your notebook or script. A jupyter notebook or command prompt for the environment can also be launched from the Anaconda Navigator. 
For more information on how to use environments, see: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

Usage
-----
The usage is best described by the example notebook: https://github.com/openearth/delft3dfmpy/blob/master/notebooks/Usage_introduction_coupled_RRFM.ipynb
