import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point, Polygon
import logging

class OSM:
    """
    OpenStreetMap model
    """

    def __init__(self, extent_file=None, data_columns=None, proj_crs = None, logger=logging):

        # Read geometry to clip data
        self.logger = logger
        self.crs_out = proj_crs

        if extent_file is not None:
            self.logger.debug(f'extent file found, reading from {extent_file}')
            self.clipgdf = gpd.read_file(extent_file)
            self.clipgeo = self.clipgdf.unary_union

            if self.crs_out is not None:
                if self.crs_out!=self.clipgdf.crs:
                    self.clipgdf.to_crs(self.crs_out, inplace=True)
        else:
            self.logger.debug(f'No extent file found, assuming the entire file is needed.')
            self.clipgdf = None
            self.clipgeo = None



        # Get required columns of OSM data
        self.data_columns = data_columns

        # Create standard dataframe for network, cross sections, orifices, weirs
        # FIXME: check available columns and required columns for the OSM data, and apply these here
        self.branches = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('branches'))

        # FIXME: ensure that all required parameterised properties are provided. I can imagine this is a matter of making
        # several parameterised profiles for different profile types (e.g. trapezoidal, rectangular, circular, etc.)
        self.profiles = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('crosssections'))

        # FIXME: ensure that all culvert types and properties can be handled. We probably have circular and box-shaped culverts, sometimes with multiple openings
        self.culverts = ExtendedGeoDataFrame(geotype=LineString, required_columns=self.get_columns('culverts'))

        # # FIXME: not sure what laterals in this context mean, but I don't think we need it at this stage.
        # self.laterals = ExtendedGeoDataFrame(geotype=Point, required_columns=[
        #     'code',
        #     'geometry'
        # ])
        #

    def get_columns(self, key):
        cols = [x.strip() for x in self.data_columns[key].replace(' ', '').split(',')]
        return cols

    def to_pickle(self, filename, overwrite=False):
        # Check if path exists
        if os.path.exists(filename) and not overwrite:
            raise FileExistsError(f'File "{filename}" already exists.')

        # Dump object
        with open(filename, 'wb') as handle:
            pickle.dump(self, handle)

    @classmethod
    def from_pickle(cls, filename):
        # Read object
        with open(filename, 'rb') as handle:
            loaded_cls = pickle.load(handle)

        return loaded_cls

