import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point, Polygon

class 3Di2FM:
    """
    3Di to D-Flow FM datamodel
    """

    def __init__(self, extent_file=None):

        # Read geometry to clip data
        if extent_file is not None:
            self.clipgeo = gpd.read_file(extent_file).unary_union
        else:
            self.clipgeo = None
            
        # Create standard dataframe for network, crosssections, structures
        self.branches = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry'
        ])
        
        self.crslocs = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "shift",
            "crosssectiondefinitionid",
            "geometry"
        ])
        
        self.crsdefs = pd.DataFrame(columns=[
            "id",
            "type",
            "thalweg",
            "numlevels",
            "levels",
            "flowwidths",
            "totalwidths"
        ])

        self.generalstructures = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            # Only require id. Other columns can be specified if needed, but requires knowledge of general structure
            # keys
            "id",
            "geometry"
        ])

        self.pumps = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "controlside",
            "maximumcapacity",
            "startlevelsuctionside",
            "stoplevelsuctionside",
            "startleveldeliveryside",
            "stopleveldeliveryside",
            "geometry"
        ])

        self.culverts = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "leftlevel",
            "rightlevel",
            "crosssectiondefinitionid",
            "length",
            "inletlosscoeff",
            "outletlosscoeff",
            "frictiontype",
            "frictionvalue",
            "geometry"
        ])

        self.culvert_crsdefs = pd.DataFrame(columns=[
            "id",
            "type",
            "thalweg",
            "numlevels",
            "levels",
            "flowwidths",
            "totalwidths"
        ])

        self.orifices = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "crestlevel",
            "crestwidth",
            "gateloweredgelevel",
            "corrcoef",
            "uselimitflowpos",
            "limitflowpos",
            "uselimitflowneg",
            "limitflowneg",
            "geometry"
        ])

        self.weirs = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "crestlevel",
            "crestwidth",
            "corrcoeff",
            "geometry"
        ])

        self.uweirs = ExtendedGeoDataFrame(geotype=Point, required_columns=[
            "id",
            "yvalues",
            "zvalues",
            "crestlevel",
            "dischargecoeff",
            "geometry"
        ])

        # Just regular dataframe, because storage nodes are connected to nodes and thus have no geo information
        self.storagenodes = pd.DataFrame(columns=[
            "id",
            "usestreetstorage",
            #"nodeid",  # will be set after generating FM 1D network
            "bedlevel",
            "area",
            "storagetype",
        ])


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


