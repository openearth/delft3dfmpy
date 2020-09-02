#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

from setuptools import setup, find_packages

here = os.path.abspath(os.path.dirname(__file__))

with open("README.rst") as readme_file:
    readme = readme_file.read()

setup(
    name="delft3dfmpy",
    description="Delft3dFM Py, automated 1D-2D and RR model setup for the Delft-3D FM suite",
    long_description=readme + "\n\n",
    url="https://github.com/openearth/delft3dfmpy/",
    author="Deltares - HKV",
    author_email="",
    packages=find_packages(),
    package_dir={"delft3dfmpy": "delft3dfmpy"},
    test_suite="tests",
    use_scm_version=True,
    setup_requires=["setuptools_scm"],
    python_requires=">=3.6",
    install_requires=[
        "pip==18.1",
        "bumpversion==0.5.3",
        "wheel==0.32.1",
        "watchdog==0.9.0",
        "flake8==3.5.0",
        "tox==3.5.2",
        "coverage==4.5.1",
        "Sphinx==1.8.1",
        "twine==1.12.1",
        "GDAL>=2.4.1",
        "geopandas>=0.5.1",
        "Fiona>=1.8.4",
        "Pillow>=6.1.0",
        "pandas>=1.0.3",
        "numpy>=1.16.4",
        "netCDF4>=1.5.1.2",
        "imod>=0.8.0",
        "rasterio>=1.1.2",
        "rasterstats>=0.14.0",
        "descartes>=1.1.0"
    ],
    extras_require={
        "dev": ["pytest", "pytest-cov", "black"],
        "optional": [],
    },
    entry_points="""
    """,
    include_package_data=True,
    license="MIT",
    zip_safe=False,
    classifiers=[
        # https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Photogrammetry",
        "License :: OSI Approved :: ???",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
    ],
    keywords="hydraulics hydrology data-science delft3dfm deltares hkv hydamo",
)
