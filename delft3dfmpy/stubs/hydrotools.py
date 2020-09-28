# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 23:18:13 2020

@author: danie
"""

import geopandas as gpd

def _filter(gdf, attribute_filter):
    if isinstance(attribute_filter, dict):
        for key, value in attribute_filter.items():
            if not isinstance(value, list):
                value = [value]
                gdf = gdf[gdf[key].isin(value)]
        return gdf
    else:
       raise IOError('attribute_filter should be dictionary') 


def read_file(path,
              hydamo_property,
              attribute_filter=None,
              column_mapping=None,
              ):
        """
        Read any OGR supported feature-file to match hydamo-property.

        A mask file can be specified to clip the selection.

        Parameters
        ----------
        path : str or Path
            Path to the feature file
        hydamo_property : HyDAMO property 
            property to map to (HyDAMO.branches, HyDAMO.crosssections)
        attribute_filter: dict
            dict with lists or strings of the format {'column_name': [values to keep]}
        column_mapping: dict
            dict for renaming input colunns to required columns
            
        Result: GeoDataFrame matching the HyDAMO property
        """
      
        gdf = gpd.read_file(path)
        
        #filter by attribute
        if attribute_filter:
            gdf = _filter(gdf, attribute_filter)
            
        #map to hydamo columns
        if column_mapping:
            gdf.rename(columns=column_mapping, inplace=True)
            
        # drop all columns not needed   
        drop_cols = [col for col in gdf.columns if not col in hydamo_property.required_columns]
        if len(drop_cols) > 0:
            gdf = gdf.drop(drop_cols, axis=1)
        
        return gdf