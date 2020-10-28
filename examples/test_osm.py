# test workflow for OpenStreetMap data

import os

from delft3dfmpy import OSM

# path to the package containing the dummy-data
path = os.path.abspath('../data/osm')
fn_pilot_area = os.path.join(path, 'study_area.geojson')

print(f'All data is expected to be in {path}')

# initialize the class
osm = OSM(fn_pilot_area)
print(type(osm))

# TODO: read branches, culvert, cross sections and properties from file and plot network to check
# TODO: BRANCHES - read id column from json.  Do not deviate between drain type
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

