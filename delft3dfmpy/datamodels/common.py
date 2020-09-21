import os
import re

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import ogr
from shapely import wkb
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

from copy import deepcopy

from delft3dfmpy.core import geometry

class ExtendedGeoDataFrame(gpd.GeoDataFrame):

    # normal properties
    _metadata = ['required_columns', 'geotype'] + gpd.GeoDataFrame._metadata

    def __init__(self, geotype, required_columns=None, *args, **kwargs):
        # Check type
        if required_columns is None:
            required_columns = []
        elif not isinstance(required_columns, list):
            required_columns = [required_columns]
        
        # Add required columns to column list
        if 'columns' in kwargs.keys():
            kwargs['columns'] += required_columns
        else:
            kwargs['columns'] = required_columns
        
        super(ExtendedGeoDataFrame, self).__init__(*args, **kwargs)

        self.required_columns = required_columns[:]
        self.geotype = geotype
        
    # def drop(self,edf, item, index_col=None,axis=None):
    #     #edf = ExtendedGeoDataFrame(geotype=LineString)
    #     temp = gpd.GeoDataFrame()        
    #     for field in self.iteritems(): 
    #         temp[field[0]] = field[1]        
    #     temp.drop(item,axis=axis, inplace=True)  
    #     edf.set_data(temp, index_col=index_col, check_columns=True)              
    #     return edf
                
    def copy(self, deep=True):
        """
        Create a copy
        """

        index = self.index.tolist() if not deep else deepcopy(self.index.tolist())
        columns = self.columns.tolist() if not deep else deepcopy(self.columns.tolist())

        edf = ExtendedGeoDataFrame(
            geotype=self.geotype,
            required_columns=[],
            index=index,
            columns=columns,
        )
        edf.required_columns.extend(self.required_columns[:])
        edf.loc[:, :] = self.values if not deep else deepcopy(self.values)
        
        return edf
    
    def delete_all(self):
        """
        Empty the dataframe
        """
        if not self.empty:
            self.iloc[:, 0] = np.nan
            self.dropna(inplace=True)

    def read_shp(self, path, index_col=None, column_mapping=None, check_columns=True, clip=None, check_geotype=True):
        """
        Import function, extended with type checks. Does not destroy reference to object.
        """
        # Read GeoDataFrame
        gdf = gpd.read_file(path)
        
        if 'MultiPolygon' in str(gdf.geometry.type):
            #gdf = gdf[gdf.geometry.type != 'MultiPolygon']  
            sfx = ['_'+str(i) for i in range(100)]
            gdf = gdf.explode()
            for ftc in gdf.code.unique():                                
                if len(gdf[gdf.code==ftc])>1:
                    gdf.loc[gdf.code==ftc,'code'] = [i+sfx[ii] for ii,i in enumerate(gdf[gdf.code==ftc].code)]
                    print(f'{ftc} is MultiPolygon; split into single parts.')               
                
            #print('Features of type \"Multipolygon\" encountered: they are skipped.')
        
        # Check number of entries
        if gdf.empty:
            raise IOError('Imported shapefile contains no rows.')

        # Rename columns:
        if column_mapping is not None:
            gdf.rename(columns=column_mapping, inplace=True)

        # Add data to class GeoDataFrame
        self.set_data(gdf, index_col=index_col, check_columns=check_columns, check_geotype=check_geotype)

        # Clip if extent is provided
        if clip is not None:
            self.clip(clip)


    def set_data(self, gdf, index_col=None, check_columns=True, check_geotype=True):

        if not self.empty:
            self.delete_all()

        # Check columns
        if check_columns:
            self._check_columns(gdf)

        # Copy content
        for col, values in gdf.iteritems():
            self[col] = values.values

        if index_col is None:
            self.index = gdf.index
            self.index.name = gdf.index.name

        else:
            self.index = gdf[index_col]
            self.index.name = index_col
        
        # Check geometry types
        if check_geotype:
            self._check_geotype()

    def _check_columns(self, gdf):
        """
        Check presence of columns in geodataframe
        """
        present_columns = gdf.columns.tolist()
        for column in self.required_columns:
            if column not in present_columns:
                raise KeyError('Column "{}" not found. Got {}, Expected at least {}'.format(
                    column, ', '.join(present_columns), ', '.join(self.required_columns)))

    def _check_geotype(self):
        """
        Check geometry type
        """
        if not all(isinstance(geo, self.geotype) for geo in self.geometry):
            raise TypeError('Geometrietype "{}" vereist. De ingevoerde shapefile heeft geometrietype(n) {}.'.format(
                re.findall('([A-Z].*)\'', repr(self.geotype))[0],
                self.geometry.type.unique().tolist()
            ))

    def read_gml(self, gml_path, index_col=None, groupby_column=None, order_column=None, column_mapping={}, check_columns=True, check_geotype=True, clip=None):
        """
        Read GML file to GeoDataFrame.

        This function has the option to group Points into LineStrings. To do so,
        specify the groupby column (which the set has in common) and the order_column,
        the column which indicates the order of the grouping.

        A mask file can be specified to clip the selection.

        Parameters
        ----------
        gml_path : str
            Path to the GML file
        groupby_column : str
            Optional, column to group points by
        order_column : str
            Optional, columns to specify the order of the grouped points
        mask_file : str
            File containing the mask to clip points.
        """

        if not os.path.exists(gml_path):
            raise OSError(f'File not found: "{gml_path}"')

        ogr.UseExceptions()
        gml = ogr.Open(gml_path)
        layer = gml.GetLayer()
        layerDefinition = layer.GetLayerDefn()
        nfields = layerDefinition.GetFieldCount()

        # Get column names for features
        columns = [layerDefinition.GetFieldDefn(i).GetName() for i in range(nfields)]
              
        # Collect features
        features = [f for f in layer]

        # Get geometries
        georefs = [f.GetGeometryRef() for f in features]
        for i, geo in enumerate(georefs):
            if geo is None:
                print('Skipping invalid geometry.')
                del georefs[i]
                del features[i]

        geometries = []
        new_feats = []
        for i,f in enumerate(features):
            geometry = wkb.loads(georefs[i].ExportToWkb())            
            if(geometry.type=='MultiPolygon'):
                new_geoms = list(geometry)                                
                geometries.extend(new_geoms)
                new_features = [f]*len(new_geoms)
                new_feats.extend(new_features)                
            else:
                geometries.append(geometry)
                new_feats.append(f)
        features = new_feats        
        #geometries = [wkb.loads(geo.ExportToWkb()) for geo in georefs]
                      
                
        # Get group by columns
        if groupby_column is not None:
            # Check if the group by column is found
            if groupby_column not in columns:
                raise ValueError('Groupby column not found in feature list.')

            # Check if the geometry is as expected
            if not isinstance(geometries[0], Point):
                raise ValueError('Can only group Points to LineString')

            # Get the values from the column that is grouped
            columnid = columns.index(groupby_column)
            groupbyvalues = [f.GetField(columnid) for f in features]
            volgid = columns.index(order_column)
            order = [f.GetField(volgid) for f in features]

            # Create empty dict for lines
            branches, counts = np.unique(groupbyvalues, return_counts=True)
            lines = {branch: [0] * count for branch, count in zip(branches, counts)}

            # Since the order does not always start at 1, find the starting number per group
            startnr = {branch: 999 for branch in branches}
            for branch, volgnr in zip(groupbyvalues, order):
                startnr[branch] = min(volgnr, startnr[branch])

            # Determine relative order of points in profile (required if the point numbering is not subsequent)
            order_rel = []
            for branch, volgnr in zip(groupbyvalues, order):
                lst_volgnr = [x[1] for x in zip(groupbyvalues, order) if x[0] == branch]
                lst_volgnr.sort()
                for i, x in enumerate(lst_volgnr):
                    if volgnr == x:
                        order_rel.append(i)

            # Filter branches with too few points
            singlepoint = (counts < 2)

            # Assign points
            for point, volgnr, branch, volgnr_rel in zip(geometries, order, groupbyvalues, order_rel):
                #lines[branch][volgnr - startnr[branch]] = point
                lines[branch][volgnr_rel] = point
                
            # Group geometries to lines
            for branch in branches[~singlepoint]:
                if any(isinstance(pt, int) for pt in lines[branch]):
                    print(f'Points are not properly assigned for branch "{branch}". Check the GML.')
                    lines[branch] = [pt for pt in lines[branch] if not isinstance(pt, int)]
                lines[branch] = LineString(lines[branch])

            # Set order for branches with single point to 0, so features are not loaded
            for branch in branches[singlepoint]:
                order[groupbyvalues.index(branch)] = 0

            # Read fields at first occurence
            startnrs = [startnr[branch] for branch in groupbyvalues]
            fields = [list(map(f.GetField, range(nfields))) for i, volgnr, f in zip(order, startnrs, features) if i == volgnr]

            # Get geometries in correct order for featuresF
            geometries = [lines[row[columnid]] for row in fields]

        else:
            fields = [list(map(f.GetField, range(nfields))) for f in features]

        # Create geodataframe
        gdf = gpd.GeoDataFrame(fields, columns=columns, geometry=geometries)
        gdf.rename(columns=column_mapping, inplace=True)

        # add a letter to 'exploded' multipolygons
        sfx = ['_'+str(i) for i in range(100)]
        for ftc in gdf.code.unique():                                
            if len(gdf[gdf.code==ftc])>1:
                gdf.loc[gdf.code==ftc,'code'] = [i+sfx[ii] for ii,i in enumerate(gdf[gdf.code==ftc].code)]
                print(f'{ftc} is MultiPolygon; split into single parts.')   
                    
        # Add data to class GeoDataFrame
        self.set_data(gdf, index_col=index_col, check_columns=check_columns, check_geotype=check_geotype)

        if clip is not None:
            self.clip(clip)

    def clip(self, geometry):
        """
        Clip geometry
        """
        if not isinstance(geometry, (Polygon, MultiPolygon)):
            raise TypeError('Expected geometry of type Polygon or MultiPolygon')

        # Clip if needed
        gdf = self.loc[self.intersects(geometry).values]
        if gdf.empty:
            raise ValueError('Found no features within extent geometry.')

        self.set_data(gdf)


    def snap_to_branch(self, branches, snap_method, maxdist=5):
        """Snap the geometries to the branch"""
        geometry.find_nearest_branch(branches=branches, geometries=self, method=snap_method, maxdist=maxdist)

class ExtendedDataFrame(pd.DataFrame):

    _metadata = ['required_columns'] + pd.DataFrame._metadata

    def __init__(self, required_columns=None, *args, **kwargs):
        super(ExtendedDataFrame, self).__init__(*args, **kwargs)

        if required_columns is None:
            required_columns = []

        self.required_columns = required_columns[:] if isinstance(required_columns, list) else [required_columns]

    def delete_all(self):
        """
        Empty the dataframe
        """
        if not self.empty:
            self.iloc[:, 0] = np.nan
            self.dropna(inplace=True)

    def set_data(self, df, index_col):

        if not self.empty:
            self.delete_all()

        # Copy content
        for col, values in df.iteritems():
            self[col] = values.values

        if index_col is None:
            self.index = df.index
            self.index.name = df.index.name

        else:
            self.index = df[index_col]
            self.index.name = index_col

        # Check columns and types
        self._check_columns()

    def add_data(self, df):

        if not np.in1d(df.columns, self.columns).all():
            raise KeyError('The new df contains columns that are not present in the current df.')

        # Concatenate data
        current = pd.DataFrame(self.values, index=self.index, columns=self.columns)
        newdf = pd.concat([current, df], ignore_index=False, sort=False)

        # Empty current df
        self.delete_all()

        # Add values again
        self.set_data(newdf, index_col=self.index.name)

    def _check_columns(self):
        """
        Check presence of columns in geodataframe
        """
        present_columns = self.columns.tolist()
        for i, column in enumerate(self.required_columns):
            if column not in present_columns:
                raise KeyError('Column "{}" not found. Got {}, Expected at least {}'.format(
                    column, ', '.join(present_columns), ', '.join(self.required_columns)))


    def read_gml(self, gml_path, index_col=None):
        """
        Read GML file to GeoDataFrame.

        This function has the option to group Points into LineStrings. To do so,
        specify the groupby column (which the set has in common) and the order_column,
        the column which indicates the order of the grouping.

        A mask file can be specified to clip the selection.

        Parameters
        ----------
        gml_path : str
            Path to the GML file
        index_col : str
            Optional, column to be set as index
        """

        if not os.path.exists(gml_path):
            raise OSError(f'File not found: "{gml_path}"')

        ogr.UseExceptions()
        gml = ogr.Open(gml_path)
        layer = gml.GetLayer()
        layer_definition = layer.GetLayerDefn()
        nfields = layer_definition.GetFieldCount()

        # Get column names for features
        columns = [layer_definition.GetFieldDefn(i).GetName() for i in range(nfields)]

        # Collect features
        features = [f for f in layer]

        # Read fields
        fields = [list(map(f.GetField, range(nfields))) for f in features]

        # Create geodataframe
        gmldf = pd.DataFrame(fields, columns=columns)

        # Add data to class GeoDataFrame
        self.set_data(gmldf, index_col=index_col)

