import os
import pickle

import geopandas as gpd
import pandas as pd

from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame
from shapely.geometry import LineString, Point, Polygon
import logging

class datamodel:
    """
    OpenStreetMap model
    """

    #FIXME BMA: we do not need to specify all columns on forehand.
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

    def get_columns(self, key):
        if self.data_columns[key] is not None:
            cols = [x.strip() for x in self.data_columns[key].replace(' ', '').split(',')]
        else:
            cols = None
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

