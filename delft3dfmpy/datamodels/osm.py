import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point, Polygon

class OSM:
    """
    OpenStreetMap model
    """

    def __init__(self, extent_file=None, data_columns=None):

        # Read geometry to clip data
        if extent_file is not None:
            self.clipgeo = gpd.read_file(extent_file).unary_union
        else:
            self.clipgeo = None

        # Get required columns of OSM data
        self.data_columns = data_columns

        # Create standard dataframe for network, cross sections, orifices, weirs
        # FIXME: check available columns and required columns for the OSM data, and apply these here
        self.branches = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('branches'))

        # FIXME: in openstreetmap, cross sections are not linestrings, perpendicular to stream, but profile types and dimensions of a channel
        # It may be that this is "parameterised cross sections, and we simply don't need the property below.
        self.crosssections = ExtendedGeoDataFrame(geotype=LineString, required_columns=[
            'code',
            'geometry',
            'ruwheidswaarde',
            'ruwheidstypecode'
        ])

        # FIXME: ensure that all required parameterised properties are provided. I can imagine this is a matter of making
        # several parameterised profiles for different profile types (e.g. trapezoidal, rectangular, circular, etc.)
        self.parametrised_profiles = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('crosssections'))


        # FIXME: ensure that all culvert types and properties can be handled. We probably have circular and box-shaped culverts, sometimes with multiple openings
        self.culverts = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('structures'))

        # # FIXME: not sure what laterals in this context mean, but I don't think we need it at this stage.
        # self.laterals = ExtendedGeoDataFrame(geotype=Point, required_columns=[
        #     'code',
        #     'geometry'
        # ])
        #

    def get_columns(self, key):
        cols = [x.strip() for x in self.data_columns[key].split('#')[0].strip().split(',')]
        return cols

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

