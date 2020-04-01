import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from delft3dfmpy.converters import hydamo_to_dflowrr
from delft3dfmpy.core import checks
from delft3dfmpy.datamodels.common import ExtendedGeoDataFrame
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

logger = logging.getLogger(__name__)

class UnpavedIO:
    def __init__(self, unpaved):
        self.unpaved = unpaved

    def unpaved_from_input(self, catchments, landuse, surface_level, soiltype,  surface_storage, infiltration_capacity , initial_gwd , meteo_areas):
        """
        Method to fill an unpaved RR-node from input.        

        Parameters
        ----------
        catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        landuse : raster
            Land use information
        surface_level : raster
            Elevation [cm+NAP]
        soiltype : raster
            Soil type information
        surface_storage : float or raster
            Storage on the surface. If a float, it is spatially uniform.
        infiltration_capacity : float or raster
            Infiltration capacity. If a float, it is spatially uniform.
        initial_gwd : float or raster
            Initial groundwater depth. If a float, it is spatially uniform.
        meteo_areas : shapely polygon object
            Representative area for a given meteo-station. 

        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_unpaved(catchments, landuse, surface_level, soiltype, surface_storage, infiltration_capacity, initial_gwd, meteo_areas)
        for unpav in geconverteerd.itertuples():
            self.unpaved.add_unpaved(  
                id = unpav.code,
                total_area = unpav.total_area,  
                lu_areas = unpav.lu_areas, 
                surface_level = unpav.mvlevel,
                soiltype = unpav.soiltype,
                surface_storage = unpav.surstor,
                infiltration_capacity = unpav.infcap,
                initial_gwd = unpav.initial_gwd,
                meteo_area = unpav.meteostat,                
                px = unpav.px, 
                py = unpav.py,
                boundary_node = unpav.boundary
            )
            
    def ernst_from_input(self, catchments, depths=None, resistance=None):
        """
        Method to create a dictionary with Ernst layer depths and resistances.

        Parameters
        ----------
         catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        depths : list optional
            Depth boundaries for the Ernst resistances.
        resistance : list, optional
            Resistances [d-1] for every Ernst layer

        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_ernst(catchments, depths, resistance)
        for ernst in geconverteerd.itertuples():
            self.unpaved.add_ernst_def(  
                id = ernst.code,
                cvo = ernst.reslist,
                lv = ernst.lvs
            )

class PavedIO:
    def __init__(self, paved):
        self.paved = paved

    def paved_from_input(self, catchments, landuse, surface_level, street_storage, sewer_storage, pumpcapacity, meteo_areas):
        """
        Method to create an RR-paved node from input.
        
        Parameters
        ----------
         catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        landuse : raster
            Land use information
        surface_level : raster
            Elevation [cm+NAP]
        street_storage : float or raster
            Storage on the street. If a float, it is spatially uniform.
        sewer_storage : float or raster
            Storage on the surface. If a float, it is spatially uniform.
        pumpcapaicty : TYPE
            Sewage pump capacity. If a float, it is spatially uniform.
        meteo_areas : shapely polygon object
            Representative area for a given meteo-station. 

        Returns
        -------
        None.

        """
               
        geconverteerd = hydamo_to_dflowrr.generate_paved(catchments, landuse, surface_level, street_storage, sewer_storage, pumpcapacity, meteo_areas)
        #PAVE id 'pav_Nde_n003' ar 16200 lv 1 sd '1' ss 0 qc 0 1.94E-05 0 qo 2 2 ms 'Station1' aaf 1 is 0 np 0 dw '1' ro 0 ru 0 qh '' pave#
        for pav in geconverteerd.itertuples():
            self.paved.add_paved(
                id = pav.code,           
                area = pav.area,
                surface_level = pav.mvlevel,
                street_storage = pav.streetstor,
                sewer_storage = pav.sewstor,                
                pump_capacity = pav.pumpcap,
                meteo_area = pav.meteostat,                               
                px = pav.px,
                py = pav.py,
                boundary_node = pav.boundary
            )
          

class GreenhouseIO:
    def __init__(self, greenhouse):
        self.greenhouse = greenhouse

    def greenhouse_from_input(self, catchments, landuse, surface_level, roof_storage, meteo_areas ):
        """
        Method to create an RR greenhouse-node from input.
        
        Parameters
        ----------
        catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        landuse : raster
            Land use information
        surface_level : raster
            Elevation [cm+NAP]        
        roof_storage : float or raster
            Storage on the greenhouse roof. If a float, it is spatially uniform.
        meteo_areas : shapely polygon object
            Representative area for a given meteo-station. 

        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_greenhouse(catchments, landuse, surface_level, roof_storage, meteo_areas)

        for gh in geconverteerd.itertuples():
             self.greenhouse.add_greenhouse(
                 id = gh.code,
                 area = gh.area,
                 surface_level = gh.mvlevel,
                 roof_storage = gh.roofstor,                 
                 meteo_area = gh.meteostat,
                 px = gh.px,
                 py = gh.py,
                 boundary_node = gh.boundary
         )

class OpenwaterIO:
    def __init__(self, openwater):
        self.openwater = openwater

    def openwater_from_input(self, catchments, landuse, meteo_areas):
        """
         Method to create an RR openwater-node from input.
        
        Parameters
        ----------
        catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        landuse : raster
            Land use information
        meteo_areas : shapely polygon object
            Representative area for a given meteo-station. 

        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_openwater(catchments, landuse, meteo_areas)

        for ow in geconverteerd.itertuples():
             self.openwater.add_openwater(
                 id = ow.code,
                 area = ow.area,
                 meteo_area = ow.meteostat,
                 px = ow.px,
                 py = ow.py,
                 boundary_node = ow.boundary
         )
          
            
class ExternalForcingsIO:

    def __init__(self, external_forcings):
        self.external_forcings = external_forcings

    def seepage_from_input(self, catchments, seepage_folder):  
        """
        Method to get seepage fluxes from raster input.
        
        Parameters
        ----------
        catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
        seepage_folder : location where rasters are stored.

        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_seepage( catchments, seepage_folder)            
        for cat in geconverteerd.iteritems():
            self.external_forcings.add_seepage(
                 id = cat[0],                
                 series = cat[1] 
            )        
        
    def precip_from_input(self, areas, precip_folder): 
         """
        Method to get precipitation fluxes from raster input.
        
        Parameters
        ----------
        areas : shapely polygon object.
            Representative areas for meteo-stations
        precip_folder : location where rasters are stored.

        Returns
        -------
        None.

        """
         geconverteerd = hydamo_to_dflowrr.generate_precip( areas, precip_folder)
        
         for cat in geconverteerd.iteritems():
            self.external_forcings.add_precip(
                 id = cat[0],                
                 series = cat[1] 
            )
            
    def evap_from_input(self, areas, evap_folder, dissolve_field=None):        
         """
        Method to get evaporation fluxes from raster input.
        
        Parameters
        ----------
        areas : shapely polygon object.
            Representative areas for meteo-stations
        evaporation_folder : location where rasters are stored.
        dissolve_field: field name of areas on wchich to dissolve to obtain an area (Sobek supports only one location for evaporation).
        Returns
        -------
        None.

        """
         geconverteerd = hydamo_to_dflowrr.generate_evap( areas, evap_folder, dissolve_field=dissolve_field)
         for cat in geconverteerd.iteritems():
            self.external_forcings.add_evap(
                 id = cat[0],                
                 series = cat[1] 
            )
            
    def boundary_from_input(self, boundary_nodes, catchments):    
        """
        Method to generate RR boundary nodes from input.

        Parameters
        ----------
        boundary_nodes : shapely point object
            Lateral nodes in DFM that are associated with RR catchments 
        catchments : shapely polygon object.
            Catchment areas; every cachtment gets a set of RR-nodes.
            
        Returns
        -------
        None.

        """
        geconverteerd = hydamo_to_dflowrr.generate_boundary( boundary_nodes, catchments)
        for bn in geconverteerd.itertuples():
             self.external_forcings.add_boundary_node(
                 id = bn.code,
                 px = bn.px,
                 py = bn.py
             )