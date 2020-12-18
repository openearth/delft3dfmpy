============================
Delft3D Flexible Mesh Python
============================


.. image:: https://img.shields.io/pypi/v/delft3dfmpy.svg
        :target: https://pypi.python.org/pypi/delft3dfmpy

.. image:: https://img.shields.io/pypi/l/delft3dfmpy.svg
        :target: https://img.shields.io/pypi/l/delft3dfmpy

Python package to generate delft3dfm models from standardized data models or other flow models.


* Free software: MIT license
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

The package ``delft3dfmpy`` requires you to have (a) an environment with required dependencies and (b) an Integrated Development Envrionment (IDE) that can access this envrionment. 

Please meet these two conditions first with instructions below before installing the Python package ``delft3dfmpy``.


Python package ``delft3dfmpy``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
1.  Install the Python package ``delft3dfmpy`` directly, by executing the following command in an Anaconda Prompt 

        ``python -m pip install delft3dfmpy``

2.  If the command prompt states ``Successfully built delft3dfmpy`` then installation is succesful.


Environment Preparation
^^^^^^^^^^^^^^^^^^^^^^^
Prepare an environment with the correct dependencies for ``delft3dfmpy``.

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


Envrionment Integration in your IDE 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Integration of the prepared envrionment depends on the IDE of usage. Here we mention briefly some options for the user (pick one!): 

1.  For a new instance of Jupyter within the activated environment:

        Using Notebook

        ``conda install -c conda-forge notebook``

        or using JupyterLab

        ``conda install -c conda-forge jupyterlab`` 

2.  To register the newly created environment as a new kernel for Jupyter (Notebook or JupyterLab):

        ``python -m ipykernel install --user --name=delft3dfmpy``

3.  No extra actions are required for PyCharm, Spyder or VSCode + Python extension.


Usage
-----

Activate the created environment in an (Anaconda) command prompt (``conda activate delft3dfmpy``) before running your notebook or script. A Jupyter notebook or command prompt for the environment can also be launched from the Anaconda Navigator. 
For more information on how to use environments, see: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html

The usage is best described by the example notebook: https://github.com/openearth/delft3dfmpy/blob/master/notebooks/Usage_introduction_coupled_RRFM.ipynb

Contribution
------------

Contributions are much welcome for documentation, desirable features and bugs. The code and development happens in this GitHub repository: https://github.com/openearth/delft3dfmpy.

For contributions, use the following guidelines:

1.  Fork the project on GitHub and clone the fork to your Operating System.

2.  Make sure you have installed and activated the environment as is described above.

3.  Delft3dfmpy uses ``flit`` to to build, package and publish the project. To install the development dependencies and register the cloned fork as a Python package for development purposes do the following:

        From an elevated Anaconda Prompt (run as Administrator) within the activated ``delft3dfmpy`` environment:

        ``conda install -c conda-forge flit``

        ``flit install --deps develop --symlink``

        This installs the development dependencies and creates a symbolic link in the Python site-packages folder of the activated environment.

4.  Open the repository as folder/workspace in your favorite IDE (eg. VSCode + Python extension)

5.  Make your contributions and test the changes locally.

6.  Once satisfied, push your changes as a new branch to your fork and create a Pull Request to the original repository.

7.  A maintainer on the main GitHub repository will review your PR and guide the merging process. 
