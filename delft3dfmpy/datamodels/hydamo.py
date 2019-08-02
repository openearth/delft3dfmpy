import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point

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
            'ruwheidswaarde',
            'ruwheidstypecode'
        ])
        
        self.crosssections = ExtendedGeoDataFrame(geotype=LineString, required_columns=['code', 'geometry', 'ruwheidswaarde'])

        self.parametrised_profiles = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'bodemhoogtebenedenstrooms',
            'bodemhoogtebovenstrooms',
            'bodembreedte',
            'taludhellinglinkerzijde',
            'taludhellingrechterzijde',
            'hoogteinsteeklinkerzijde',
            'hoogteinsteekrechterzijde'
        ])
        
        # Weirs
        self.weirs = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code',
            'geometry',
            'laagstedoorstroomhoogte',
            'hoogstedoorstroomhoogte',
            'laagstedoorstroombreedte',
            'hoogstedoorstroombreedte',
            'kruinbreedte',
            'afvoercoefficient'
        ])

        # Orifices
        self.orifices = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry',
            'hoogteopening',
            'breedteopening',
            'indpeilregulpeilscheidend',
            'hoogtebinnenonderkantbenedenstrooms',
            'hoogtebinnenonderkantbovenstrooms'
        ])

        # Culverts
        self.culverts = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry',
            'lengte',
            'hoogteopening',
            'breedteopening',
            'hoogtebinnenonderkantbenedenstrooms',
            'hoogtebinnenonderkantbovenstrooms',
            'vormcode',
            'intreeverlies',
            'uittreeverlies',
            'ruwheidstypecode',
            'ruwheidswaarde'
        ])
        
        # Laterals
        self.laterals = ExtendedGeoDataFrame(geotype=Point, required_columns=['code', 'name', 'geometry'])

        # Gemalen
        self.gemalen = ExtendedGeoDataFrame(geotype=Point, required_columns=['code', 'naam', 'geometry'])
        self.pumps = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            'code', 'maximalecapaciteit', 'geometry', 'codegerelateerdobject'])
        self.sturing = ExtendedDataFrame(required_columns=[
            'code', 'streefwaarde', 'bovenmarge', 'ondermarge', 'codegerelateerdobject'])

        # Hydraulische randvoorwaarden
        self.boundary_conditions = ExtendedGeoDataFrame(
            geotype=Point, required_columns=['code', 'typerandvoorwaardecode', 'geometry'])

    def to_pickle(filename, overwrite=False):
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


