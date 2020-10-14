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

# TODO: snapping. May not be needed as wel only use parameterized profiles

# TODO assign elevation value to cross sections. this needs to be retrieved from a DEM (which we have!)

# TODO: where cross sections are meeting a culvert, it may make sense to straighten the elevation value up and downstream of culvert

# TODO: create DFM parsers for drains, culverts and cross sections based on dfmmodel.structures.io.xxxx_from_hydamo as example

