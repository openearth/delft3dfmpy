# -*- coding: utf-8 -*-

"""Top-level package for delft3dfmpy."""

__author__ = """Guus Rongen"""
__email__ = 'rongen@hkv.nl'
__version__ = '1.1.3'

from delft3dfmpy.core.dfm import DFlowFMModel
from delft3dfmpy.core.drr import DFlowRRModel
from delft3dfmpy.core.mesh2d import Rectangular
from delft3dfmpy.datamodels.hydamo import HyDAMO
from delft3dfmpy.datamodels.osm import OSM
from delft3dfmpy.io.dflowfmwriter import DFlowFMWriter
from delft3dfmpy.io.dflowrrwriter import DFlowRRWriter

from delft3dfmpy.core.logging import initialize_logger
initialize_logger()
