import itertools
import logging
import os
   
import geopandas as gpd
import numpy as np
import pandas as pd
import tqdm
from scipy.spatial import KDTree
from shapely.geometry import LineString, Point, Polygon

from delft3dfmpy.converters import hydamo_to_dflowrr
from delft3dfmpy.core import checks, geometry
from delft3dfmpy.datamodels.common import ExtendedGeoDataFrame
from delft3dfmpy.datamodels.cstructures import meshgeom, meshgeomdim
from delft3dfmpy.io import drrreader

logger = logging.getLogger(__name__)

class DFlowRRModel:
    """Main data structure for RR-model in DflowFM. Contains subclasses
    for unpaved, paved,greehouse and open water nodes and external forcings (seepage, precipitation, evaporation)
      """

    def __init__(self):

        
        self.d3b_parameters = {}       

        self.unpaved = Unpaved(self)
        
        self.paved = Paved(self)
        
        self.greenhouse = Greenhouse(self)               
        
        self.openwater = Openwater(self)               
                
        self.external_forcings = ExternalForcings(self)
        
        self.dimr_path = ''
    
class ExternalForcings:
    """
    Class for external forcings, which contains the boundary
    conditions and the initial conditions.
    """

    def __init__(self, dflowrrmodel):
        # Point to relevant attributes from parent
        self.dflowrrmodel = dflowrrmodel
        self.io = drrreader.ExternalForcingsIO(self)
        
        self.boundary_nodes = {}
        self.seepage = {}
        self.precip = {}
        self.evap = {}
        
    def add_precip(self, id, series):
        self.precip[id] = {
            'precip' : series
        }
        
    def add_evap(self, id, series):
        self.evap[id] = {
            'evap' : series
        }
        
    def add_seepage(self, id, series):
        self.seepage[id] = {
            'seepage' : series
        }
    
    def add_boundary_node(self, id, px, py):
        self.boundary_nodes[id] = {
            'id' : id,  
            'px' : px,
            'py' : py            
        }
    
class Unpaved:
    """
    Class for unpaved nodes
    """
        
    def __init__(self, dflowrrmodel):
        # Point to relevant attributes from parent
        self.dflowrrmodel = dflowrrmodel
        
        # initialize a dataframe for every type of nodes related to 'unpaved'
        self.unp_nodes = {}    
        self.ernst_defs = {}
        
        self.io = drrreader.UnpavedIO(self)

    def add_unpaved(self,id, total_area, lu_areas, surface_level, soiltype, surface_storage, infiltration_capacity, initial_gwd, meteo_area, px, py, boundary_node):
        self.unp_nodes[id] = {
            'id' : 'unp_'+id,  
            'na' : '16',
            'ar' : lu_areas,
            'ga' : total_area,
            'lv' : surface_level,
            'co' : '3',
            'su' : '0',
            'sd' : surface_storage,
            'sp' : 'sep_'+id,
            'ic' : infiltration_capacity,
            'ed' : 'ernst_'+id,
            'bt' : soiltype,
            'ig' : initial_gwd,
            'mg' : surface_level,
            'gl' : '1.5',
            'is' : '0',
            'ms' : 'ms_'+meteo_area,
            'px': px,
            'py': py,
            'boundary_node': boundary_node            
        }
        
    def add_ernst_def(self,id, cvo, lv, cvi, cvs):
        self.ernst_defs[id] = {
            'id' : 'ernst_'+id,            
            'cvi' : cvi,
            'cvs' : cvs,
            'cvo' : cvo,
            'lv'  : lv            
        }        
        
class Paved:
    """
    Class for paved nodes.
    """

    def __init__(self, dflowrrmodel):
        # Point to relevant attributes from parent
        self.dflowrrmodel = dflowrrmodel        
        self.pav_nodes = {}
        self.io = drrreader.PavedIO(self)
        
        self.node_geom = {}
        self.link_geom = {}
        
    #PAVE id 'pav_Nde_n003' ar 16200 lv 1 sd '1' ss 0 qc 0 1.94E-05 0 qo 2 2 ms 'Station1' aaf 1 is 0 np 0 dw '1' ro 0 ru 0 qh '' pave#
    def add_paved(self,id,  area, surface_level, street_storage, sewer_storage, pump_capacity, meteo_area, px, py, boundary_node):
        self.pav_nodes[id] = {
            'id' : 'pav_'+id,        
            'ar' : area,
            'lv' : surface_level,            
            'qc' : pump_capacity,            
            'strs' : street_storage,
            'sews' : sewer_storage,            
            'ms' : 'ms_'+meteo_area,
            'is' : '0',
            'np' : '0',
            'ro' : '0',
            'ru' : '0',
            'px': px,
            'py': py,
            'boundary_node': boundary_node
        }                        
                
class Greenhouse:
    """
    Class for greenhouse nodes
    """
    def __init__(self, dflowrrmodel):        
                
        self.dflowrrmodel = dflowrrmodel
        self.gh_nodes = {}
        # Create the io class
        self.io = drrreader.GreenhouseIO(self)

#    GRHS id ’1’ na 10 ar 1000. 0. 0. 3000. 0. 0. 0. 0. 0. 0. sl 1.0 as 0. sd ’roofstor 1mm’ si
#    ’silo typ1’ ms ’meteostat1’ is 50.0 grhs

    def add_greenhouse(self, id, area, surface_level, roof_storage, meteo_area, px, py, boundary_node):
        self.gh_nodes[id] = {            
            'id': 'gh_'+id,
            'ar' : area,
            'sl': surface_level,            
            'sd': roof_storage,
            'ms' : 'ms_'+meteo_area,
            'is' : '0',
            'px': px,
            'py': py,
            'boundary_node': boundary_node
        }

class Openwater:
    """
    Class for open water nodes
    """
    def __init__(self, dflowrrmodel):        
                
        self.dflowrrmodel = dflowrrmodel
        self.ow_nodes = {}
        # Create the io class
        self.io = drrreader.OpenwaterIO(self)

    def add_openwater(self, id, area, meteo_area, px, py, boundary_node):
        self.ow_nodes[id] = {            
            'id': 'ow_'+id,
            'ar' : area,
            'ms' : 'ms_'+meteo_area,
            'px': px,
            'py': py,
            'boundary_node': boundary_node
        }
         