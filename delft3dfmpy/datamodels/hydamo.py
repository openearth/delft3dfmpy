import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point, Polygon

class HyDAMO:
    """
    HyDAMO model
    """

    def __init__(self, extent_file=None):

        # Read geometry to clip data
        if extent_file is not None:
            self.clipgeo = gpd.read_file(extent_file).unary_union
        else:
            self.clipgeo = None
            
        # Create standard dataframe for network, crosssections, orifices, weirs
        self.branches = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry',
            'typeruwheid',            
        ])
        
        self.profile = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry',
            'globalid',
            'profiellijnid'            
            
        ])
        self.profile_roughness = ExtendedDataFrame(required_columns=[
            'code',
            'globalid',
            'profielpuntid'        
        ])
        
        
        # self.crosssections = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
        #     'code',
        #     'geometry'            
        # ])
        

        self.param_profile = ExtendedDataFrame(required_columns=[
            'globalid',
            'normgeparamprofielid',
            'hydroobjectid'
        ])
        
        self.param_profile_values = ExtendedDataFrame(required_columns=[
            'normgeparamprofielid',
            'soortparameter',
            'waarde',            
            'ruwheidlaag',
            'ruwheidhoog',
            'typeruwheid'
        ])
        
        # Weirs
        self.weirs = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',
            'geometry',
            'globalid',
            'soortstuw',                          
            'afvoercoefficient'
        ])         
        
        # opening
        self.opening = ExtendedDataFrame(required_columns=[            
            'vormopening',            
            'globalid',
            'hoogstedoorstroombreedte',
            'hoogstedoorstroomhoogte',
            'laagstedoorstroombreedte',
            'laagstedoorstroomhoogte',            
            'vormopening',
            'afvoercoefficient'
         ])
        
        # opening
        self.closing_device = ExtendedDataFrame(required_columns=[
            'code'    
            
        ])
        
        # opening
        self.management_device = ExtendedDataFrame(required_columns=[
            'code',
            'soortregelbaarheid',            
            'maximalehoogtebovenkant',
            'maximalehoogtebovenkant'          
         ])
        
        # Bridges
        self.bridges = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',                        
            'globalid',
            'geometry',
            'lengte',
            'intreeverlies',
            'uittreeverlies'
            
        ])
        
        # # Orifices
        # self.orifices = ExtendedGeoDataFrame(geotype=Point, required_columns=[
        #      'code',
        #      'geometry',
        #      'laagstedoorstroomhoogte',
        #      'laagstedoorstroombreedte',
        #      'schuifhoogte',
        #      'afvoercoefficient'
        #  ])

        # Culverts
        self.culverts = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',            
            'geometry',
            'lengte',
            'hoogteopening',
            'breedteopening',
            'hoogtebinnenonderkantbene',
            'hoogtebinnenonderkantbov',
            'vormkoker',
            'intreeverlies',
            'uittreeverlies',
            'typeruwheid',
            'ruwheid'
        ])
        
        # Laterals
        self.laterals = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',            
            'geometry'
        ])

        # Gemalen
        self.pumpstations = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',            
            'globalid',
            'geometry',
        ])
        self.pumps = ExtendedDataFrame(required_columns=[
            'code',
            'globalid',
            'gemaalid',
            'maximalecapaciteit'            
        ])
        self.management = ExtendedDataFrame(required_columns=[
            'code',
            'globalid'            
        ])
        # self.afsluitmiddel = ExtendedDataFrame(required_columns=[
        #     'code',
        #     'soortafsluitmiddelcode',
        #     'codegerelateerdobject'
        # ])

        # Hydraulische randvoorwaarden
        self.boundary_conditions = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',
            'typerandvoorwaarde',
            'geometry'
        ])
        
        # RR catchments
        self.catchments = ExtendedGeoDataFrame(geotype=Polygon, required_columns=[
            'code',
            'geometry',
            'globalid',
            'lateraleknoopid'            
        ])
        
        # Laterals
        self.laterals = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',
            'geometry',
            'globalid'            
        ])
        
        # RR overflows
        self.overflows = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',
            'geometry',
            'codegerelateerdobject',
            'fractie'
            
        ])

        # RR sewer areas
        self.sewer_areas = ExtendedGeoDataFrame(geotype=Polygon, required_columns=[
            'code',
            'geometry'            
        ])
        
        
    def roughness_int2str(self, df):
        roughness = {
            1: "Chezy",
            2: "Manning",
            3: "StricklerNikuradse",
            4: "Strickler",
            5: "WhiteColebrook",
            6: "deBosBijkerk"
        }
        for val in roughness.keys():
            if isinstance(df.typeruwheid.iloc[0],str):            
                idx = df[df['typeruwheid'] == str(val)].index
            else:
                idx = df[df['typeruwheid'] == val].index
            df.loc[idx, 'typeruwheid'] = roughness[val]
        return df
    
    def roughness_str2int(self, df):
        roughness = {
            1: "Chezy",
            2: "Manning",
            3: "StricklerNikuradse",
            4: "Strickler",
            5: "WhiteColebrook",
            6: "deBosBijkerk"
        }
        for strval in roughness.values():
            idx = df[df['typeruwheid'] == strval].index
            df.loc[idx, 'typeruwheid'] = [key for key, val in roughness.items() if val==strval][0]
        return df
    
        
    def to_pickle(self, filename, overwrite=False):
        # Check if path exists
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(f'File "{filename}" alraedy exists.')

        # Dump object
        with open(filename, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def from_pickle(cls, filename):
        # Read object
        with open(filename, 'rb') as handle:
            loaded_cls = pickle.load(handle)

        return loaded_cls



