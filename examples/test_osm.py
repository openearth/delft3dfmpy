# test workflow for OpenStreetMap data

import os
import configparser, json
from delft3dfmpy import OSM
from delft3dfmpy.core.logging import initialize_logger
import matplotlib.pyplot as plt
import geopandas as gpd

import logging

root = os.path.abspath('../data/osm')
fn_ini = os.path.join(root, 'osm_settings.ini')

logger = initialize_logger('osm2fm.log', log_level=10)
# Read ini file
logger.info(f'Read config from {fn_ini}')
config = configparser.ConfigParser()
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

# TODO: BRANCHES - read id column from json.  Do not deviate between drain type
# Read branches and store in OSM data model
osm.branches.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs = osm.crs_out, clip = osm.clipgeo
                      , id_col=id, filter_cols=True, logger=logger)

# TODO: CROSS SECTIONS DEFINTION - read id, drain_type, material, width, depth, top_width, diameter, profile_op, profile_cl, bottom_width columns from json
# read cross-sections
osm.parametrised_profiles.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs = osm.crs_out, clip = osm.clipgeo
                      , id_col=id, filter_cols=True, logger=logger)

# TODO: CROSS SECTIONS DEFINTION - specify roughness dependent on material add this
# TODO: CROSS SECTION DEFINITION -  assign elevation value to cross sections. this needs to be retrieved from a DEM (which we have!)
# TODO: CROSS SECTION DEFINITION - create circular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create rectangular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create trapezoid profiles --> prof_idofbranch
# TODO: CROSS SECTION LOCATION - select rows with drain type other than culvert
# TODO: CROSS SECTION LOCATION - add cross sections at start and end of branch. Take the longitudinal slope with SHIFT parameter into account
osm.parametrised_profiles.snap_to_branch(osm.branches, snap_method='ends')

# FIXME: cross-sections are plotted as branches. This means that snapping should be changed.
# FIXME: when adding the culvert. I encounter crs problem.
# Plot branches and cross-sections
plt.rcParams['axes.edgecolor'] = 'w'

fig, ax = plt.subplots(figsize=(10, 10))

#ax.fill(*osm.clipgeo.exterior.xy, color='w', alpha=0.5)
ax.xaxis.set_visible(False)
ax.yaxis.set_visible(False)

background = plt.imread(path+'/background_projected.png')
ax.imshow(background, extent=(524564.3221, 529442.7747, 9246725.9975, 9249557.8336), interpolation='lanczos')
osm.clipgdf.plot(ax=ax, color='w', alpha=0.5)
osm.branches.plot(ax=ax, label='Channel')
osm.parametrised_profiles.plot(ax=ax, color='C3', label='Cross section')
plt.show()

# # TODO: STRUCTURE - read id, draintype
# # Read culverts
# osm.culverts.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs= crs_proj, clip = osm.clipgeo,
#                       id_col=id, filter_cols=True, draintype_col=config.get('datacolumns','draintypecolumn')
#                       , filter_culverts=True, logger=logger)

# TODO: STRUCTURE - select rows with draintype culvert
# TODO: STRUCTURE - determine length and midpoint location culvert
# TODO: STRUCTURE - snapping of open drains over culverts. May not be needed as wel only use parameterized profiles
# TODO: STRUCTURE - where cross sections are meeting a culvert, it may make sense to straighten the elevation value up and downstream of culvert
# TODO: plot branches + cross sections locations + structures

# TODO: create DFM parsers for drains, culverts and cross sections based on dfmmodel.structures.io.xxxx_from_hydamo as example

# TODO: create 2D grid with refinement

# TODO: create 1D2D links

print("Hello world")
