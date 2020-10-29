# test workflow for OpenStreetMap data

import os
import configparser, json
from delft3dfmpy import OSM

# Read ini file
config = configparser.ConfigParser()
config.read(os.path.abspath('../data/osm')+'/osm_settings.ini')

# Path to data
path = config.get('input', 'DataPath')

# Extend of study area
fn_pilot_area = os.path.join(path, config.get('input', 'studyareafile'))

print(f'All data is expected to be in {path}')

# Get required columns
required_columns_data = config._sections['datacolumns']

osm = OSM(fn_pilot_area, required_columns_data)

print(type(osm))

# TODO: read branches, culvert, cross sections and properties from file and plot network to check
# TODO: BRANCHES - read id column from json.  Do not deviate between drain type
osm.branches.read_shp(os.path.join(path,config.get('input','datafile')),index_col='id',clip = osm.clipgeo, id_col='id')

print('hello world')
# TODO: BRANCHES - connect branches on the right locations
# TODO: plot branches

# TODO: CROSS SECTIONS DEFINTION - read id, drain_type, material, width, depth, top_width, diameter, profile_op, profile_cl, bottom_width columns from json
# TODO: CROSS SECTIONS DEFINTION - specify roughness dependent on material add this
# TODO: CROSS SECTION DEFINITION -  assign elevation value to cross sections. this needs to be retrieved from a DEM (which we have!)
# TODO: CROSS SECTION DEFINITION - create circular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create rectangular profiles --> prof_idofbranch
# TODO: CROSS SECTION DEFINITION - create trapezoid profiles --> prof_idofbranch
# TODO: CROSS SECTION LOCATION - select rows with drain type other than culvert
# TODO: CROSS SECTION LOCATION - add cross sections at start and end of branch. Take the longitudinal slope with SHIFT parameter into account
# TODO: plot branches + cross sections locations

# TODO: STRUCTURE - read id, draintype
# TODO: STRUCTURE - select rows with draintype culvert
# TODO: STRUCTURE - determine length and midpoint location culvert
# TODO: STRUCTURE - snapping of open drains over culverts. May not be needed as wel only use parameterized profiles
# TODO: STRUCTURE - where cross sections are meeting a culvert, it may make sense to straighten the elevation value up and downstream of culvert
# TODO: plot branches + cross sections locations + structures

# TODO: create DFM parsers for drains, culverts and cross sections based on dfmmodel.structures.io.xxxx_from_hydamo as example

# TODO: create 2D grid with refinement

# TODO: create 1D2D links

