# test workflow for OpenStreetMap data

import os
import configparser, json
from delft3dfmpy import OSM
from delft3dfmpy.core.logging import initialize_logger
import matplotlib.pyplot as plt
import pandas as pd

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
osm.branches.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs = osm.crs_out, clip = osm.clipgeo
                      , id_col=id, filter_cols=True, logger=logger)

# Read cross-sections and store in OSM data model
osm.profiles.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs = osm.crs_out,
                      clip = osm.clipgeo, id_col=id, filter_cols=True, geometry_duplicater=1, logger=logger)

# Snap profiles to branches: create two profile locations along branches.
# osm.profiles.snap_to_branch(osm.branches, snap_method='one_data')
# retrieve profiles at start of each line segment
profiles_start = osm.profiles.branch_to_prof(offset=0.5, prefix='A_', rename_col='id')
# retrieve profiles at end of each line segment
profiles_end = osm.profiles.branch_to_prof(offset=0.5, prefix='B_', rename_col='id', vertex_end=True)
# concat into a new profiles object
osm.profiles = pd.concat([profiles_start, profiles_end])

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
profiles_start.plot(ax=ax, marker='.', alpha=0.3, color='g')
profiles_end.plot(ax=ax, marker='.', alpha=0.3, color='r')

plt.show()

# # Read culverts into OSM
osm.culverts.read_shp(os.path.join(path,config.get('input','datafile')),index_col=id, proj_crs= osm.crs_out, clip = osm.clipgeo,
                      id_col=id, filter_cols=True, draintype_col=config.get('datacolumns','draintypecolumn')
                      , filter_culverts=True, logger=logger)

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
osm.profiles.geometry.interpolate(osm.profiles.branch_offset).plot(ax=ax1, marker='*', markersize=5, color='C3', label='Cross section', zorder=5)
osm.culverts.centroid.plot(ax=ax1, color='yellow', label='Culvert', markersize=5, zorder=10)
plt.show()


# TODO: CROSS SECTION LOCATION - add cross sections at start and end of branch. Take the longitudinal slope with SHIFT parameter into account
# TODO: CROSS SECTIONS DEFINTION - specify roughness dependent on material add this
# TODO: CROSS SECTION DEFINITION -  assign elevation value to cross sections. this needs to be retrieved from a DEM (which we have!)
# TODO: CROSS SECTION DEFINITION - create circular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create rectangular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create trapezoid profiles --> prof_idofbranch
# TODO: CROSS SECTION LOCATION - select rows with drain type other than culvert
# TODO: STRUCTURE - select rows with draintype culvert
# TODO: STRUCTURE - determine length and midpoint location culvert
# TODO: STRUCTURE - snapping of open drains over culverts. May not be needed as wel only use parameterized profiles
# TODO: STRUCTURE - where cross sections are meeting a culvert, it may make sense to straighten the elevation value up and downstream of culvert
# TODO: plot branches + cross sections locations + structures

# TODO: create DFM parsers for drains, culverts and cross sections based on dfmmodel.structures.io.xxxx_from_hydamo as example

# TODO: create 2D grid with refinement

# TODO: create 1D2D links

print("Hello world")
