# test workflow for OpenStreetMap data

import os
import configparser, json
from delft3dfmpy import OSM, DFlowFMModel
from delft3dfmpy.datamodels.common import ExtendedGeoDataFrame
from delft3dfmpy.core.logging import initialize_logger
import matplotlib.pyplot as plt
import pandas as pd
import geopandas as gpd

import logging

root = os.path.abspath('../data/osm')
fn_ini = os.path.join(root, 'osm_settings.ini')

logger = initialize_logger('osm2fm.log', log_level=10)

# Read ini file
logger.info(f'Read config from {fn_ini}')
config = configparser.ConfigParser(inline_comment_prefixes=[";", "#"])
config.read(fn_ini)

# Path to data
path = config.get('input', 'DataPath')

# Extend of study area
fn_pilot_area = os.path.join(path, config.get('input', 'studyareafile'))

logger.info(f'All data is expected to be in {path}')

# Get required columns
required_columns_data = config._sections['datacolumns']

# Get parameters
parameters = config._sections['parameter']

# Initialise osm data
osm = OSM(fn_pilot_area, required_columns_data,parameters['projectedcrs'], logger=logger)

# print(type(osm))

# Id column
id = config.get('datacolumns','idcolumn')


# Read branches and store in OSM data model
osm.branches.read_shp(os.path.join(path,config.get('input','datafile')),
                      index_col=id,
                      proj_crs=osm.crs_out,
                      clip=osm.clipgeo,
                      id_col=id,
                      filter_cols=True,
                      logger=logger)


# Read cross-sections and store in OSM data model
osm.profiles.read_shp(os.path.join(path,config.get('input','datafile')),
                      index_col=id,
                      proj_crs=osm.crs_out,
                      clip=osm.clipgeo,
                      id_col=id,
                      filter_cols=True,
                      logger=logger
                      )

# retrieve profiles at start of each line segment
profiles_start = osm.profiles.branch_to_prof(offset=0.5, prefix = 'Prof_', suffix='_A', rename_col='id')
# retrieve profiles at end of each line segment
profiles_end = osm.profiles.branch_to_prof(offset=0.5, prefix = 'Prof_', suffix='_B', rename_col='id', vertex_end=True)
# combine into one and add set to OSM
osm.profiles.set_data(pd.concat([profiles_start, profiles_end]), check_columns=False, check_geotype=False)
# Merge profile_cl and profile_open for profiles
osm.profiles.merge_columns(col1='profile_cl', col2='profile_op', rename_col='profile')

# # Read culverts into OSM
osm.culverts.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs= osm.crs_out, clip = osm.clipgeo,
                      id_col=id, filter_cols=True, filter_rows={'drain_type': 'culvert'}, logger=logger)
# Merge profile_cl and profile_open for culverts
osm.culverts.merge_columns(col1='profile_cl', col2='profile_op', rename_col='profile')
# Snap culvert to branches and determine centroid.
osm.culverts.snap_to_branch(osm.branches, snap_method='ends')


# Plot branches, cross-sections and culverts
fig1, ax1 = plt.subplots(figsize=(10, 10))

ax1.xaxis.set_visible(False)
ax1.yaxis.set_visible(False)

background = plt.imread(path+'/background_projected.png')
ax1.imshow(background, extent=(524564.3221, 529442.7747, 9246725.9975, 9249557.8336), interpolation='lanczos')
osm.clipgdf.plot(ax=ax1, color='w', alpha=0.5)
osm.branches.plot(ax=ax1, label='Channel')
osm.profiles.geometry.plot(ax=ax1, marker='.', color='r' , markersize=5, label='Cross section')
osm.culverts.centroid.plot(ax=ax1, color='yellow', label='Culvert', markersize=5, zorder=10)
#plt.show()

# TODO: CROSS SECTION LOCATION - add cross sections at start and end of branch. Take the longitudinal slope with SHIFT parameter into account
# TODO: replace dummy dem values of dem at cross-sections, dem at leftlevel, dem at rightlevel for profiles and culverts
# TODO: compute depth_to_bottom for culverts and profiles
# TODO: CROSS SECTION DEFINITION -  assign elevation value to cross sections and depth. this needs to be retrieved from a DEM (which we have!)
#osm.profiles.sample_raster(rasterio,offset=None,geometry)
# Temporary DEM_crsloc and depth_to_bottom
osm.profiles['DEM_crsloc'] = 12
osm.profiles['depth_to_bottom'] = osm.profiles[['depth', 'diameter']].max(axis=1)+0.5
osm.profiles['bottom_level'] = osm.profiles.DEM_crsloc - osm.profiles.depth_to_bottom

# Temporary DEM_leftlevel, DEM_rightlevel and depth_to_bottom
osm.culverts['DEM_leftlevel'] = 11
osm.culverts['DEM_rightlevel'] = 10
osm.culverts['depth_to_bottom'] = osm.culverts[['depth', 'diameter']].max(axis=1)+0.3
osm.culverts['leftlevel'] = osm.culverts.DEM_leftlevel - osm.culverts.depth_to_bottom
osm.culverts['rightlevel'] = osm.culverts.DEM_rightlevel - osm.culverts.depth_to_bottom

# needed variables: id, branchid, chainage, leftlevel, rightlevel, crosssection, outletlosscoef, allowedflowdir
#  valveonoff, numlosscoeff, valveopeningheight, relopening, losscoeff,
#                     frictiontype, frictionvalue

# Start dfmmodel
dfmmodel = DFlowFMModel()


# Collect culverts
# TODO: id  --> c_+ id of branchid
osm.culverts['id'] = 'C_' + osm.culverts['id']

# TODO: type  --> 'culvert' --> drain_type
osm.culverts.drain_type #However in d-hydamo this is fixed

# TODO: branchid --> branch_id
osm.culverts.branch_id

# TODO: chainage --> branch_offset
osm.culverts.branch_offset

# TODO: left level/ rightlevel

# TODO: allowed flowdir --> none assumption everyhting is both direction
osm.culverts['allowedflowdir'] = 'both'

# TODO: inletloss - set standaard in inifile
osm.culverts['inletlosscoeff'] = parameters['inletlosscoefculvert']
# TODO: outletloss - set standaard in inifile
osm.culverts['outletlosscoeff'] = parameters['outletlosscoefculvert']

# TODO: csDefid - get branchid --> search profile
osm.culverts['crosssection'] = [{} for _ in range(len(osm.culverts))]

for culvert in osm.culverts.itertuples():
    # Generate cross section definition name
    if culvert.profile == 'round':
        crosssection = {'shape': 'circle', 'diameter': culvert.diameter}

    elif culvert.profile == 'boxed_rectangular':
        crosssection = {'shape': 'rectangle', 'height': culvert.depth, 'width': culvert.width,
                        'closed': 'yes'}
    else:
        crosssection = {'shape': 'circle', 'diameter': 0.50}
        print(f'Culvert {culvert.id} has an unknown shape: {culvert.profile}. Applying a default profile (round - 50cm)')

# TODO: valve onoff - no valve = 0
osm.culverts['valveOnOff'] = 0

# TODO: numlosscoeff - standard on 0
osm.culverts['numlosscoef'] = 0

# TODO: bedfrictiontype - inifile
friction_type = parameters['frictiontype']

# TODO: bedfrictionvalue - material and smoothness, inifile
friction_values = dict(zip(parameters['frictionmaterials'].split(','), list(map(float,parameters['frictionvalues'].split(',')))))

osm.culverts.columns
# TODO: CROSS SECTIONS DEFINTION - specify roughness dependent on material add this

# TODO: CROSS SECTION DEFINITION - create circular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create rectangular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create trapezoid profiles --> prof_idofbranch
# TODO: CROSS SECTION LOCATION - select rows with drain type other than culvert
# TODO: STRUCTURE - determine length
# TODO: STRUCTURE - snapping of open drains over culverts. May not be needed as wel only use parameterized profiles
# TODO: STRUCTURE - where cross sections are meeting a culvert, it may make sense to straighten the elevation value up and downstream of culvert
# TODO: plot branches + cross sections locations + structures

# TODO: create DFM parsers for drains, culverts and cross sections based on dfmmodel.structures.io.xxxx_from_hydamo as example

# TODO: create 1D network

# TODO: create 2D grid with refinement

# TODO: create 1D2D links

print("Hello world")
