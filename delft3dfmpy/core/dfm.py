import itertools
import logging
import os
from typing import Union, List

import geopandas as gpd
import numpy as np
import pandas as pd
from tqdm.auto import tqdm
from scipy.spatial import KDTree
from shapely.geometry import LineString, Point, Polygon
import shapely

from delft3dfmpy.converters import hydamo_to_dflowfm
from delft3dfmpy.core import checks, geometry
from delft3dfmpy.datamodels.common import ExtendedGeoDataFrame
from delft3dfmpy.datamodels.cstructures import meshgeom, meshgeomdim
from delft3dfmpy.io import dfmreader

logger = logging.getLogger(__name__)


class DFlowFMModel:
    """Main data structure for dflowfm model. Contains subclasses
    for network, structures, cross sections, observation points
    and external forcings.
    """

    def __init__(self):

        self.mdu_parameters = {}

        self.network = Network(self)

        self.structures = Structures(self)

        self.crosssections = CrossSections(self)

        self.observation_points = ObservationPoints(self)

        self.external_forcings = ExternalForcings(self)

        self.storagenodes = StorageNodes(self)

        self.roughness_mapping = {
            "Chezy": "Chezy",
            "Manning": "Manning",
            "StricklerKn": "StricklerNikuradse",
            "StricklerKs": "Strickler",
            "White Colebrook": "WhiteColebrook",
            "Bos en Bijkerk": "deBosBijkerk",
            "Onbekend": "Strickler",
            "Overig": "Strickler",
        }

    def export_network(self, output_dir, overwrite=False):
        """
        Expert network to shapefiles

        Three files are exported. One for the 2d mesh, one for the 1d mesh
        and one for the links between the 1d and 2d mesh.

        Parameters
        ----------
        output_dir : str
            Path to output directory
        overwrite : bool, optional
            Whether to overwrite existing files, by default False

        Raises
        ------
        FileExistsError
            If one of the shape file names already exists.
        """
        # Check if files already exist
        files = ["mesh1d.shp", "mesh2d.shp", "links1d2d.shp"]
        paths = [os.path.join(output_dir, file) for file in files]
        if not overwrite:
            for path in paths:
                if os.path.exists(path):
                    raise FileExistsError(
                        f'Path "{path}" already exists. Choose another output folder or specify overwrite=True.'
                    )

        # Links
        links = self.network.links1d2d.get_1d2dlinks(as_gdf=True)
        links.crs = {"init": "epsg:28992"}
        links.to_file(paths[2])

        # Mesh2d
        mesh2d = gpd.GeoDataFrame(
            geometry=[Polygon(poly) for poly in self.network.mesh2d.get_faces()],
            crs="epsg:28992",
        )
        # Add properties
        # Save
        mesh2d.to_file(paths[1])

        # Mesh1d
        mesh1d = gpd.GeoDataFrame(
            geometry=[LineString(line) for line in self.network.mesh1d.get_segments()],
            crs="epsg:28992",
        )
        # Add properties
        mesh1d["node1"], mesh1d["node2"] = self.network.mesh1d.get_values(
            "edge_nodes", as_array=True
        ).T
        # Save
        mesh1d.to_file(paths[0])


class ExternalForcings:
    """
    Class for external forcings, which contains the boundary
    conditions and the initial conditions.
    """

    def __init__(self, dflowfmmodel):
        # Point to relevant attributes from parent
        self.dflowfmmodel = dflowfmmodel
        self.initial_waterlevel_polygons = gpd.GeoDataFrame(
            columns=["waterlevel", "geometry", "locationType"]
        )
        self.initial_waterdepth_polygons = gpd.GeoDataFrame(
            columns=["waterdepth", "geometry", "locationType"]
        )
        self.missing = None
        self.mdu_parameters = dflowfmmodel.mdu_parameters

        # GeoDataFrame for saving boundary conditions
        self.boundaries = {}  # gpd.GeoDataFrame()

        # Dataframe for saving time series for structure
        self.structures = pd.DataFrame(
            columns=[
                "id",
                "type",
                "parameter",
                "time",
                "time_unit",
                "value",
                "value_unit",
            ]
        )

        # Dictionary for saving laterals
        self.laterals = {}

        self.io = dfmreader.ExternalForcingsIO(self)

    def set_initial_waterlevel(self, level, polygon=None, name=None, locationtype="1d"):
        """
        Method to set initial water level. A polygon can be given to
        limit the initial water level to a certain extent.

        """
        # Get name is not given as input
        if name is None:
            name = "wlevpoly{:04d}".format(len(self.initial_waterlevel_polygons) + 1)

        # Add to geodataframe
        if polygon == None:
            new_df = pd.DataFrame(
                {
                    "waterlevel": level,
                    "geometry": polygon,
                    "locationType": locationtype,
                },
                index=[name],
            )
            self.initial_waterlevel_polygons = new_df
        else:
            self.initial_waterlevel_polygons.loc[name] = {
                "waterlevel": level,
                "geometry": polygon,
                "locationType": locationtype,
            }

    def set_missing_waterlevel(self, missing):
        """
        Method to set the missing value for the water level.
        this overwrites the water level at missing value in the mdu file.

        Parameters
        ----------
        missing : float
            Water depth
        """
        self.mdu_parameters["WaterLevIni"] = missing

    def set_initial_waterdepth(self, depth, polygon=None, name=None, locationtype="1d"):
        """
        Method to set the initial water depth in the 1d model. The water depth is
        set by determining the water level at the locations of the cross sections.

        Parameters
        ----------
        depth : float
            Water depth
        """
        # Get name is not given as input
        if name is None:
            name = "wlevpoly{:04d}".format(len(self.initial_waterdepth_polygons) + 1)
        # Add to geodataframe
        if polygon == None:

            new_df = pd.DataFrame(
                {
                    "waterdepth": depth,
                    "geometry": polygon,
                    "locationType": locationtype,
                },
                index=[name],
            )

            self.initial_waterdepth_polygons = new_df
        else:
            self.initial_waterdepth_polygons.loc[name] = {
                "waterdepth": depth,
                "geometry": polygon,
                "locationType": locationtype,
            }

    def add_rainfall_2D(self, fName, bctype="rainfall"):
        """
        Parameters
        ----------
        fName : str
            Location of netcdf file containing rainfall rasters
        bctype : str
            Type of boundary condition. Currently only rainfall is supported
        """

        assert bctype in ["rainfall"]

        # Add boundary condition
        self.boundaries["rainfall_2D"] = {
            "file_name": fName,
            "bctype": bctype + "bnd",
        }

    def add_boundary_condition(self, name, pt, bctype, series, branchid=None):
        """
        Add boundary conditions to model:
        - The boundary condition can be discharge or waterlevel
        - Is specified by a geographical location (pt) and a branchid
        - If no branchid is given, the nearest is searched
        - The boundary condition is added to the end of the given or nearest branch.
        - To specify a time dependendend boundary: a timeseries with values should be given
        - To specify a constant boundary a float should be given

        Parameters
        ----------
        name : str
            ID of the boundary condition
        pt : tuple or shapely.geometry.Point
            Location of the boundary condition
        bctype : str
            Type of boundary condition. Currently only discharge and waterlevel are supported
        series : pd.Series or float
            If a float, a constant in time boundary condition is used. If a pandas series,
            the values per time step are used. Index should be in datetime format
        branchid : str, optional
            ID of the branch. If None, the branch nearest to the given location (pt) is
            searched, by default None
        """

        assert bctype in ["discharge", "waterlevel"]

        unit = "m3/s" if bctype.lower() == "discharge" else "m"
        if name in self.boundaries.keys():
            raise KeyError(
                f'A boundary condition with name "{name}" is already present.'
            )

        if isinstance(pt, tuple):
            pt = Point(*pt)

        # Find the nearest node
        nodes1d = np.asarray(
            [
                n.split("_")
                for n in self.dflowfmmodel.network.mesh1d.description1d[
                    "network_node_ids"
                ]
            ],
            dtype="float",
        )
        get_nearest = KDTree(nodes1d)
        distance, idx_nearest = get_nearest.query(pt)
        nodeid = f"{nodes1d[idx_nearest][0]:6f}_{nodes1d[idx_nearest][1]:6f}"

        # # If   branch is not given explicitly, find the nearest
        # if branchid is None:
        #     branchid = self.dflowfmmodel.network.branches.distance(pt).idxmin()

        # # Find nearest branch for geometry
        # branch = self.dflowfmmodel.network.schematised.at[branchid, 'geometry']
        # extended_line = geometry.extend_linestring(line=branch, near_pt=pt, length=1.0)

        # # Create intersection line for boundary condition
        # bcline = LineString(geometry.orthogonal_line(line=extended_line, offset=0.1, width=0.1))

        # Convert time to minutes
        if isinstance(series, pd.Series):
            times = ((series.index - series.index[0]).total_seconds() / 60.0).tolist()
            values = series.values.tolist()
            startdate = pd.datetime.strftime(series.index[0], "%Y-%m-%d %H:%M:%S")
        else:
            times = None
            values = series
            startdate = "0000-00-00 00:00:00"

        # Add boundary condition
        self.boundaries[name] = {
            "code": name,
            "bctype": bctype + "bnd",
            "value": values,
            "time": times,
            "filetype": 9,
            "operand": "O",
            "method": 3,
            "time_unit": f"minutes since {startdate}",
            "value_unit": unit,
            "nodeid": nodeid,
        }

        # Check if a 1d2d link should be removed
        self.dflowfmmodel.network.links1d2d.check_boundary_link(self.boundaries[name])

    def add_rain_series(self, name, values, times):
        """
        Adds a rain series a boundary condition.
        Specify name, values, and times

        Parameters
        ----------
        name : str
            ID of the condition
        values : list of floats
            Values of the rain intensity
        times : list of datetime
            Times for the values
        """
        # Add boundary condition
        self.boundaries[name] = {
            "code": name,
            "bctype": "rainfall",
            "filetype": 1,
            "method": 1,
            "operand": "O",
            "value": values,
            "time": times,
            "geometry": None,
            "branchid": None,
        }

    def set_structure_series(self, structure_id, structure_type, parameter, series):
        if structure_type.lower() == "weir":
            weirnames = list(self.dflowfmmodel.structures.weirs.keys())
            if structure_id not in weirnames:
                raise IndexError(
                    f'"{structure_id}" not in index: "{",".join(weirnames)}"'
                )
            self.dflowfmmodel.structures.weirs[structure_id][
                parameter
            ] = "boundaries.bc"
            if "crest" in parameter.lower():
                parameter = "weir_crestLevel"
                unit = "m+NAP"
        elif structure_type.lower() == "orifice":
            orificenames = list(self.dflowfmmodel.structures.orifices.keys())
            if structure_id not in orificenames:
                raise IndexError(
                    f'"{structure_id}" not in index: "{",".join(orificenames)}"'
                )
            if "gate" in parameter.lower():
                parameter = "orifice_gateLowerEdgeLevel"
                unit = "m"
            self.dflowfmmodel.structures.orifices[structure_id][
                parameter
            ] = "boundaries.bc"
        elif structure_type.lower() == "culvert":
            culvertnames = list(self.dflowfmmodel.structures.culverts.keys())
            if structure_id not in culvertnames:
                raise IndexError(
                    f'"{structure_id}" not in index: "{",".join(culvertnames)}"'
                )
            if "valve" in parameter.lower():
                parameter = "culvert_valveOpeningHeight"
                unit = "m"
            self.dflowfmmodel.structures.culverts[structure_id][
                parameter
            ] = "boundaries.bc"
        else:
            raise NotImplementedError(
                "Only implemented for weirs, culverts and orifices."
            )

        #  Convert time to minutes
        if isinstance(series, pd.DataFrame):
            series = series.iloc[:, 0]
        if isinstance(series, pd.Series):
            times = ((series.index - series.index[0]).total_seconds() / 60.0).tolist()
            values = series.values.tolist()
            startdate = series.index[0]
        else:
            times = None
            values = series
            startdate = "0000000000"

        # Add to dataframe
        self.structures.loc[structure_id] = {
            "id": structure_id,
            "type": structure_type,
            "parameter": parameter,
            "time_unit": f'minutes since {pd.datetime.strftime(startdate,"%Y-%m-%d %H:00:00")}',
            "value_unit": unit,
            "time": times,
            "value": values,
        }


class CrossSections:
    """
    Cross section class.
    Contains dictionaries for the cross section locations and
    cross section definitions. Also has an subclass for input
    and output.
    """

    def __init__(self, dflowfmmodel):
        """
        Constructor

        Parameters
        ----------
        DFlowFMMOdel
            Model schematisation of which the cross sections are part.
        """
        self.io = dfmreader.CrossSectionsIO(self)

        self.crosssection_loc = {}
        self.crosssection_def = {}

        self.dflowfmmodel = dflowfmmodel

        self.default_definition = None
        self.default_locations = None
        self.default_definition_shift = 0.0

        self.get_roughnessname = self.dflowfmmodel.network.get_roughness_description

    def set_default_definition(self, definition, shift=0.0):
        """
        Add default profile
        """
        if definition not in self.crosssection_def.keys():
            raise KeyError(f'Cross section definition "{definition}" not found."')

        self.default_definition = definition
        self.default_definition_shift = shift

    def set_default_locations(self, locations):
        """
        Add default profile locations
        """

        self.default_locations = locations

    def add_yz_definition(
        self, yz=None, thalweg=None, roughnesstype=None, roughnessvalue=None, name=None
    ):
        """
        Add xyz crosssection

        Parameters
        ----------
        code : str
            Id of cross section
        branch : str
            Name of branch
        offset : float
            Position of cross section along branch. If not given, the position is determined
            from the branches in the network. These should thus be given in this case.
        crds : np.array
            Nx2 array with y, z coordinates
        """

        # get coordinates
        length, z = yz.T
        if name is None:
            name = f"yz_{len(crosssection_def):08d}"

        # Get roughnessname
        roughnessname = self.get_roughnessname(roughnesstype, roughnessvalue)

        # Add to dictionary
        self.crosssection_def[name] = {
            "id": name,
            "type": "yz",
            "thalweg": np.round(thalweg, decimals=3),
            "yzCount": len(z),
            "yCoordinates": list_to_str(length),
            "zCoordinates": list_to_str(z),
            "sectionCount": 1,
            "frictionIds": roughnessname,
            "frictionPositions": list_to_str([length[0], length[-1]]),
        }

        return name

    def add_circle_definition(self, diameter, roughnesstype, roughnessvalue, name=None):
        """
        Add circle cross section. The cross section name is derived from the shape and roughness,
        so similar cross sections will result in a single definition.
        """
        # Get name if not given
        if name is None:
            name = f"circ_d{diameter:.3f}"

        # Get roughnessname
        roughnessname = self.get_roughnessname(roughnesstype, roughnessvalue)

        # Add to dictionary
        self.crosssection_def[name] = {
            "id": name,
            "type": "circle",
            "thalweg": 0.0,
            "diameter": diameter,
            "frictionId": roughnessname,
        }

        return name

    def add_rectangle_definition(
        self, height, width, closed, roughnesstype, roughnessvalue, name=None
    ):
        """
        Add rectangle cross section. The cross section name is derived from the shape and roughness,
        so similar cross sections will result in a single definition.
        """
        # Get name if not given
        if name is None:
            name = f"rect_h{height:.3f}_w{width:.3f}"

        # Get roughnessname
        roughnessname = self.get_roughnessname(roughnesstype, roughnessvalue)

        # Add to dictionary
        self.crosssection_def[name] = {
            "id": name,
            "type": "rectangle",
            "thalweg": 0.0,
            "height": height,
            "width": width,
            "closed": int(closed),
            "frictionId": roughnessname,
        }

        return name

    def add_trapezium_definition(
        self,
        slope,
        maximumflowwidth,
        bottomwidth,
        closed,
        roughnesstype,
        roughnessvalue,
        bottomlevel=None,
        name=None,
    ):
        """
        Add rectangle cross section. The cross section name is derived from the shape and roughness,
        so similar cross sections will result in a single definition.
        """
        # Get name if not given
        if name is None:
            name = f"trapz_s{slope:.1f}_bw{bottomwidth:.1f}_bw{maximumflowwidth:.1f}"

        # Get roughnessname
        roughnessname = self.get_roughnessname(roughnesstype, roughnessvalue)

        if bottomlevel is None:
            bottomlevel = 0.0

        if not closed:
            levels = f"{bottomlevel} 100"
            flowwidths = (
                f"{bottomwidth:.2f} {bottomwidth + 2.*((100.0-bottomlevel)*slope):.2f}"
            )
        else:
            levels = f"0 {((maximumflowwidth - bottomwidth)/2.0) / slope:.2f}"
            flowwidths = f"{bottomwidth:.2f} {maximumflowwidth:.2f}"

        # Add to dictionary
        self.crosssection_def[name] = {
            "id": name,
            "type": "zw",
            "thalweg": 0.0,
            "numLevels": 2,
            "levels": levels,
            "flowWidths": flowwidths,
            "totalWidths": flowwidths,
            "frictionId": roughnessname,
        }

        return name

    def add_zw_definition(
        self,
        numLevels,
        levels,
        flowWidths,
        totalWidths,
        roughnesstype,
        roughnessvalue,
        name=None,
    ):
        """
        Add zw cross section. The cross section name is derived from the shape and roughness,
        so similar cross sections will result in a single definition.
        """
        # Get name if not given
        if name is None:
            name = (
                f'zw_h{levels.replace(" ","_"):.1f}_w{flowWidths.replace(" ","_"):.1f}'
            )

        # Get roughnessname
        roughnessname = self.get_roughnessname(roughnesstype, roughnessvalue)

        # Add to dictionary
        self.crosssection_def[name] = {
            "id": name,
            "type": "zw",
            "thalweg": 0.0,
            "numLevels": int(numLevels),
            "levels": levels,
            "flowWidths": flowWidths,
            "totalWidths": totalWidths,
            "frictionId": roughnessname,
        }

        return name

    def add_crosssection_location(
        self, branchid, chainage, definition, minz=np.nan, shift=0.0
    ):

        descr = f"{branchid}_{chainage:.1f}"
        # Add cross section location
        self.crosssection_loc[descr] = {
            "id": descr,
            "branchid": branchid,
            "chainage": chainage,
            "shift": shift,
            "definitionId": definition,
        }

    def get_branches_without_crosssection(self):
        # First find all branches that match a cross section
        branch_ids = {dct["branchid"] for _, dct in self.crosssection_loc.items()}
        # Select the branch-ids that do nog have a matching cross section
        branches = self.dflowfmmodel.network.branches
        no_crosssection = branches.index[~np.isin(branches.index, list(branch_ids))]

        return no_crosssection.tolist()

    def get_structures_without_crosssection(self):
        struc_ids = [dct["id"] for _, dct in self.crosssection_def.items()]
        bridge_ids = [
            dct["csDefId"] for _, dct in self.dflowfmmodel.structures.bridges.items()
        ]
        no_cross_bridge = np.asarray(bridge_ids)[
            ~np.isin(bridge_ids, struc_ids)
        ].tolist()
        no_crosssection = no_cross_bridge
        # uweir_ids = [dct['numlevels'] for _, dct in self.dflowfmmodel.structures.uweirs.items()]
        # no_cross_uweir = uweirs.index[~np.isin(uweirs.index, list(uweir_ids))]
        # no_crosssection = no_crosssection.append(no_cross_uweir)
        return no_crosssection

    def get_bottom_levels(self):
        """Method to determine bottom levels from cross sections"""

        # Initialize lists
        data = []
        geometry = []

        for key, css in self.crosssection_loc.items():
            # Get location
            geometry.append(
                self.dflowfmmodel.network.schematised.at[
                    css["branchid"], "geometry"
                ].interpolate(css["chainage"])
            )
            shift = css["shift"]

            # Get depth from definition if yz and shift
            definition = self.crosssection_def[css["definitionId"]]
            minz = shift
            if definition["type"].lower() == "yz":
                minz += min(float(z) for z in definition["zCoordinates"].split())

            data.append([css["branchid"], css["chainage"], minz])

        # Add to geodataframe
        gdf = gpd.GeoDataFrame(
            data=data, columns=["branchid", "chainage", "minz"], geometry=geometry
        )
        return gdf


class Links1d2d:
    def __init__(self, network):
        self.mesh1d = network.mesh1d
        self.mesh2d = network.mesh2d
        self.network = network

        # List for 1d 2d links
        self.nodes1d = []
        self.faces2d = []

    def generate_1d_to_2d(self, max_distance=np.inf, branchid=None):
        """
        Generate 1d2d links from 1d nodes. Each 1d node is connected to
        the nearest 2d cell. A maximum distance can be specified to remove links
        that are too long. Also the branchid can be specified, if you only want
        to generate links from certain 1d branches.

        Parameters
        ----------
        max_distance : int, float
            The maximum length of a link. All longer links are removed.
        branchid : str or list
            ID's of branches for which the connection from 1d to 2d is made.
        """
        logger.info(f"Generating links from 1d to 2d based on distance.")

        # Create KDTree for faces
        faces2d = np.c_[
            self.mesh2d.get_values("facex"), self.mesh2d.get_values("facey")
        ]
        get_nearest = KDTree(faces2d)

        # Get network geometry
        nodes1d = self.mesh1d.get_nodes()
        idx = self.mesh1d.get_nodes_for_branch(branchid)

        # Get nearest 2d nodes
        distance, idx_nearest = get_nearest.query(nodes1d[idx])
        close = distance < max_distance

        # Add link data
        nodes1didx = np.arange(len(nodes1d))[idx]
        self.nodes1d.extend(nodes1didx[close] + 1)
        self.faces2d.extend(idx_nearest[close] + 1)

        # Remove conflicting 1d2d links
        for bc in self.network.dflowfmmodel.external_forcings.boundaries.values():
            if "geometry" not in bc:
                continue
            self.check_boundary_link(bc)

    def generate_2d_to_1d(
        self,
        max_distance=np.inf,
        intersecting=True,
        branchid=None,
        shift_to_centroid=True,
    ):
        """
        Generate 1d2d links from 2d cells. A maximum distance can be specified
        to remove links that are too long. Also a branchid can be specified to only
        generate links to certain branches.

        In case of a 1D and 2D grid that is on top of each other the user might want
        to generate links only for intersecting cells, where in case of non-overlapping meshes
        the use might want to use the shortest distance. This behaviour can be specified with
        the option intersecting:
        1. intersecting = True: each 2d cell crossing a 1d branch segment is connected to
            the nearest 1d cell.
        2. intersecting = False: each 2d cell is connected to the nearest 1d cell,
            If the link crosses another cell it is removed.
        In case of option 2. setting a max distance will speed up the the process a bit.

        Parameters
        ----------
        max_distance : int, float
            Maximum allowed length for a link
        intersecting : bool
            Make connections for intersecting 1d and 2d cells or based on
            nearest neighbours
        branchid : str or list of str
            Generate only to specified 1d branches
        """
        logger.info(
            f'Generating links from 2d to 1d based on {"intersection" if intersecting else "distance"}.'
        )

        # Collect polygons for cells
        centers2d = self.mesh2d.get_faces(geometry="center")
        idx = np.arange(len(centers2d), dtype="int")
        # Create KDTree for 1d cells
        nodes1d = self.mesh1d.get_nodes()
        nodes1didx = self.mesh1d.get_nodes_for_branch(branchid)
        get_nearest = KDTree(nodes1d[nodes1didx])

        # Make a pre-selection
        if max_distance < np.inf:
            # Determine distance from 2d to nearest 1d
            distance, _ = get_nearest.query(centers2d)
            idx = idx[distance < max_distance]

        # Create GeoDataFrame
        logger.info(f"Creating GeoDataFrame of ({len(idx)}) 2D cells.")
        cells = gpd.GeoDataFrame(
            data=centers2d[idx],
            columns=["x", "y"],
            index=idx + 1,
            geometry=[
                Polygon(cell)
                for i, cell in enumerate(self.mesh2d.get_faces())
                if i in idx
            ],
        )

        # Find intersecting cells with branches
        logger.info("Determine intersecting or nearest branches.")
        if branchid is None:
            branches = self.network.branches
        elif isinstance(branchid, str):
            branches = self.network.branches.loc[[branchid]]
        else:
            branches = self.network.branches.loc[branchid]

        if intersecting:
            geometry.find_nearest_branch(branches, cells, method="intersecting")
        else:
            geometry.find_nearest_branch(
                branches, cells, method="overal", maxdist=max_distance
            )

        # Drop the cells without intersection
        cells.dropna(subset=["branch_offset"], inplace=True)
        faces2d = np.c_[cells.x, cells.y]

        # Get nearest 1d nodes
        distance, idx_nearest = get_nearest.query(faces2d)
        close = distance < max_distance

        # Add link data
        nodes1didx = np.where(nodes1didx)[0][idx_nearest]
        self.nodes1d.extend(nodes1didx[close] + 1)
        self.faces2d.extend(cells.index.values[close])

        if not intersecting:
            logger.info("Remove links that cross another 2D cell.")
            # Make sure only the nearest cells are accounted by removing all links that also cross another cell
            links = self.get_1d2dlinks(as_gdf=True)
            todrop = []

            # Remove links that intersect multiple cells
            cellbounds = cells.bounds.values.T
            for link in tqdm(
                links.itertuples(),
                total=len(links),
                desc="Removing links crossing mult. cells",
            ):
                selectie = cells.loc[
                    geometry.possibly_intersecting(cellbounds, link.geometry)
                ].copy()
                if selectie.intersects(link.geometry).sum() > 1:
                    todrop.append(link.Index)
            links.drop(todrop, inplace=True)

            # Re-assign
            del self.nodes1d[:]
            del self.faces2d[:]

            self.nodes1d.extend(links["node1did"].values.tolist())
            self.faces2d.extend(links["face2did"].values.tolist())

        # Shift centers of 2d faces to centroid if they are part of a 1d-2d link
        if shift_to_centroid:
            # Get current centers
            cx, cy = self.mesh2d.get_faces(geometry="center").T
            # Calculate centroids for cells with link
            idx = np.array(self.faces2d) - 1
            centroids = np.vstack(
                [cell.mean(axis=0) for cell in np.array(self.mesh2d.get_faces())[idx]]
            ).T
            cx[idx] = centroids[0]
            cy[idx] = centroids[1]
            # Set values back to geometry
            self.mesh2d.set_values("facex", cx)
            self.mesh2d.set_values("facey", cy)

        # Remove conflicting 1d2d links
        for bc in self.network.dflowfmmodel.external_forcings.boundaries.values():
            if bc["geometry"] is None:
                continue
            self.check_boundary_link(bc)

    def check_boundary_link(self, bc):
        """
        Since a boundary conditions is not picked up when there is a bifurcation
        in the first branch segment, potential 1d2d links should be removed.

        This function should be called whenever a boundary conditions is added,
        or the 1d2d links are generated.
        """

        # Can only be done after links have been generated
        if not self.nodes1d or not self.faces2d:
            return None

        # Find the nearest node with the KDTree
        nodes1d = self.mesh1d.get_nodes()
        get_nearest = KDTree(nodes1d)
        distance, idx_nearest = get_nearest.query(
            [float(pt) for pt in bc["nodeid"].split("_")]
        )
        node_id = idx_nearest + 1

        # Check 1. Determine if the nearest node itself is not a bifurcation
        edge_nodes = self.mesh1d.get_values("edge_nodes", as_array=True)
        counts = {u: c for u, c in zip(*np.unique(edge_nodes, return_counts=True))}
        if counts[node_id] > 1:
            logger.warning(
                f"The boundary condition at {node_id} is not a branch end. Check if it is picked up by dflowfm."
            )

        # Check 2. Check if any 1d2d links are connected to the node or next node. If so, remove.
        # Find the node(s) connect to 'node_id'
        to_remove = np.unique(edge_nodes[(edge_nodes == node_id).any(axis=1)])
        for item in to_remove:
            while item in self.nodes1d:
                loc = self.nodes1d.index(item)
                self.nodes1d.pop(loc)
                self.faces2d.pop(loc)
                nx, ny = nodes1d[item - 1]
                # bcx, bcy = bc['geometry'].centroid.coords[0]
                logger.info(
                    f"Removed link(s) from 1d node: ({nx:.2f}, {ny:.2f}) because it is too close to boundary condition at node {node_id:.0f}."
                )

    def get_1d2dlinks(self, as_gdf=False):
        """
        Method to get 1d2d links as array with coordinates or geodataframe.

        Parameters
        ----------
        as_gdf : bool
            Whether to export as geodataframe (True) or numpy array (False)
        """

        if not any(self.nodes1d):
            return None

        # Get 1d nodes and 2d faces
        nodes1d = self.mesh1d.get_nodes()
        faces2d = self.mesh2d.get_faces(geometry="center")

        # Get links
        links = np.dstack(
            [nodes1d[np.array(self.nodes1d) - 1], faces2d[np.array(self.faces2d) - 1]]
        )

        if not as_gdf:
            return np.array([line.T for line in links])
        else:
            return gpd.GeoDataFrame(
                data=np.c_[self.nodes1d, self.faces2d],
                columns=["node1did", "face2did"],
                geometry=[LineString(line.T) for line in links],
            )

    def remove_1d2d_from_numlimdt(self, file, threshold, node="2d"):
        """
        Remove 1d2d links based on numlimdt file
        """
        if node == "1d":
            links = self.get_1d2dlinks(as_gdf=True)

        with open(file) as f:
            for line in f.readlines():
                x, y, n = line.split()
                if int(n) >= threshold:
                    if node == "2d":
                        self.remove_1d2d_link(
                            float(x), float(y), mesh=node, max_distance=2.0
                        )
                    else:
                        # Find the 1d node connected to the link
                        idx = links.distance(Point(float(x), float(y))).idxmin()
                        x, y = links.at[idx, "geometry"].coords[0]
                        self.remove_1d2d_link(x, y, mesh=node, max_distance=2.0)

    def remove_1d2d_link(self, x, y, mesh, max_distance):
        """
        Remove 1d 2d link based on x y coordinate.
        Mesh can specified, 1d or 2d.
        """
        if mesh == "1d":
            pts = self.mesh1d.get_nodes()
            ilink = 0
        elif mesh == "2d":
            pts = np.c_[self.mesh2d.get_faces(geometry="center")]
            ilink = 1
        else:
            raise ValueError()

        # Find nearest link
        dists = np.hypot(pts[:, 0] - x, pts[:, 1] - y)
        if dists.min() > max_distance:
            return None
        imin = np.argmin(dists)

        # Determine what rows to remove (if any)
        linkdim = self.nodes1d if mesh == "1d" else self.faces2d
        to_remove = [link for link in (linkdim) if link == (imin + 1)]
        for item in to_remove:
            while item in linkdim:
                loc = linkdim.index(item)
                self.nodes1d.pop(loc)
                self.faces2d.pop(loc)

    def remove_1d_endpoints(self):
        """Method to remove 1d2d links from end points of the 1d mesh. The GUI
        will interpret every endpoint as a boundary conditions, which does not
        allow a 1d 2d link at the same node. To avoid problems with this, use
        this method.
        """
        # Can only be done after links have been generated
        if not self.nodes1d or not self.faces2d:
            return None

        nodes1d = self.mesh1d.get_nodes()
        edge_nodes = self.mesh1d.get_values("edge_nodes", as_array=True)

        # Select 1d nodes that are only present in a single edge
        edgeid, counts = np.unique(edge_nodes, return_counts=True)
        to_remove = edgeid[counts == 1]

        for item in to_remove:
            while item in self.nodes1d:
                loc = self.nodes1d.index(item)
                self.nodes1d.pop(loc)
                self.faces2d.pop(loc)
                nx, ny = nodes1d[item - 1]
                logger.info(
                    f"Removed link(s) from 1d node: ({nx:.2f}, {ny:.2f}) because it is connected to an end-point."
                )

    def remove_links1d2d_within_polygon(self, polygon) -> None:
        """Remove 1d2d links within a given polygon or multipolygon

        Args:
            network (Network): The network from which the links are removed
            polygon (Union[Polygon, MultiPolygon]): The polygon that indicates which to remove
        """

        # Create an array with 2d facecenters and 1d nodes, that form the links
        nodes1d = self.mesh1d.get_nodes()[np.array(self.nodes1d) - 1]
        faces2d = self.mesh2d.get_faces(geometry="center")[np.array(self.faces2d) - 1]

        # Check which links intersect the provided area
        index = np.zeros(len(nodes1d), dtype=bool)
        for part in geometry.as_polygon_list(polygon):
            index |= geometry.points_in_polygon(nodes1d, part)
            index |= geometry.points_in_polygon(faces2d, part)

        # Remove these links
        for i in reversed(np.where(index)[0]):
            self.nodes1d.pop(i)
            self.faces2d.pop(i)


class Network:
    def __init__(self, dflowfmmodel):
        # Link dflowmodel
        self.dflowfmmodel = dflowfmmodel

        # Mesh 1d offsets
        self.offsets = {}

        # Branches and schematised branches
        self.branches = ExtendedGeoDataFrame(
            geotype=LineString, required_columns=["code", "geometry"]
        )
        self.schematised = ExtendedGeoDataFrame(
            geotype=LineString, required_columns=["geometry"]
        )

        # Create mesh for the 1d network
        self.mesh1d = meshgeom(meshgeomdim())
        self.mesh1d.meshgeomdim.dim = 1

        self.mesh1d.edge_branchoffsets = []
        self.mesh1d.edge_branchidx = []
        self.mesh1d.edge_x = []
        self.mesh1d.edge_y = []

        # Create mesh for the 1d network
        self.mesh2d = meshgeom(meshgeomdim())
        self.mesh2d.meshgeomdim.dim = 2

        # Create 1d2dlinks
        self.links1d2d = Links1d2d(self)

        # Dictionary for roughness definitions
        self.roughness_definitions = {}

        # Link mdu parameters
        self.mdu_parameters = dflowfmmodel.mdu_parameters

    def set_branches(self, branches):
        """
        Set branches from geodataframe
        """
        # Check input
        checks.check_argument(
            branches, "branches", (ExtendedGeoDataFrame, gpd.GeoDataFrame)
        )
        # Add data to branches
        self.branches.set_data(branches[self.branches.required_columns])
        # Copy branches to schematised
        self.schematised.set_data(self.branches[self.schematised.required_columns])

    def snap_branch_ends(self, offset, output=None):
        """
        Method to snap branch ends to other branch ends within a given offset.

        Parameters
        offset : float
            Maximum distance between end points. If the distance is larger, they are not snapped.
        """
        # Collect endpoints
        endpoints = []
        for branch in self.branches.itertuples():
            endpoints.append((branch.geometry.coords[0], branch.Index, 0))
            endpoints.append((branch.geometry.coords[-1], branch.Index, -1))

        # Create KDTree of endpoints
        snapped = 0
        df = []

        # For every endpoint determine the distance to the nearest other endpoint
        for i, (endpoint, branchid, side) in tqdm(
            enumerate(endpoints), total=len(endpoints), desc="Snapping branch ends"
        ):
            # Find distance and nearest coordinate of other points
            other_pts = [pt[0] for j, pt in enumerate(endpoints) if j != i]
            mindist, minidx = KDTree(other_pts).query(endpoint)

            # Change endpoint if dist is not 0.0 but smaller than offset.
            if mindist != 0.0 and mindist <= offset:

                # Change coordinates of branch
                crds = self.branches.at[branchid, "geometry"].coords[:]
                crds[side] = other_pts[minidx]
                self.branches.at[branchid, "geometry"] = LineString(crds)
                snapped += 1

                df.append([endpoint, other_pts[minidx]])
                endpoints[i] = (other_pts[minidx], branchid, side)

        if output is not None:
            df = pd.DataFrame(df)
            df.to_csv(os.path.join(output, f"snappedpoints_{offset}.csv"))
            print(f"Snapped {snapped} points.")

    def set_branch_order(self, branchids, idx=None):
        """
        Group branch ids so that the cross sections are
        interpolated along the branch.

        Parameters
        ----------
        branchids : list
            List of branches to group
        """
        # Get the ids (integers) of the branch names given by the user
        branchidx = np.isin(self.mesh1d.description1d["network_branch_ids"], branchids)
        # Get current order
        branchorder = self.mesh1d.get_values("nbranchorder", as_array=True)
        # Update
        if idx is None:
            branchorder[branchidx] = branchorder.max() + 1
        else:
            if not isinstance(idx, int):
                raise TypeError("Expected integer.")
            branchorder[branchidx] = idx
        # Save
        self.mesh1d.set_values("nbranchorder", branchorder)

    def set_branch_interpolation_modelwide(self):
        """
        Set cross-section interpolation over nodes on all branches model-wide. I

        Note:
            - Interpolation will be set using branchorder property.
            - Branches are grouped between bifurcations where multiple branches meet.
            - No interpolation is applied over these type of bifurcations.
            - No branch order is set on branch groups consisting of 1 branch.
        """
        self.get_grouped_branches()
        for group in self.branch_groups.values():
            if len(group) > 1:
                self.set_branch_order(group)

    def make_nodes_to_branch_map(self):
        # Note: first node is upstream, second node is downstream
        self.nodes_to_branch_map = {
            b: [self.mesh1d.description1d["network_node_ids"][_idx - 1] for _idx in idx]
            for b, idx in zip(
                self.mesh1d.description1d["network_branch_ids"],
                self.mesh1d.get_values("nedge_nodes", as_array=True),
            )
        }

    def make_branches_to_node_map(self):
        self.make_nodes_to_branch_map()
        self.branches_to_node_map = {
            n: [k for k, v in self.nodes_to_branch_map.items() if n in v]
            for n in self.mesh1d.description1d["network_node_ids"]
        }

    def generate_nodes_with_bedlevels(
        self, resolve_at_bifurcation_method="min", return_reversed_branches=False
    ):
        """
        Generate nodes with upstream and downstream bedlevels derived from set cross-sections on branch. It takes into
        account whether or not branch order is specified (so interpolation over nodes is set).

        Nodes in between cross-sections on same branch/branch-group are linearly interpolated. Outside are extrapolated
        constant (e.g. end points of branch).

        Branch groups which include branch with reversed direction compared to whole group will be taken into account.
        Use return_reversed_branches=True to return that information

        Specify with resolve_at_bifurcation_method how to resolve bedlevel at bifurcation of more than 2 branches using
        options 'min' (minimum), 'max' (maximum), 'mean' (average).
        """
        assert resolve_at_bifurcation_method in ["min", "max", "mean"], (
            f"Incorrect value for "
            f"'resolve_at_bifurcation_method' supplied. "
            f"Either use 'min', 'max' or 'mean'"
        )
        bedlevels_crs_branches = self.dflowfmmodel.crosssections.get_bottom_levels()
        branch_order = self.mesh1d.get_values("nbranchorder", as_array=True)
        self.make_branches_to_node_map(), self.make_nodes_to_branch_map()
        nodes_dict = {
            n: {"up": [], "down": []} for n in self.branches_to_node_map.keys()
        }
        reserved_branches = []
        for order, (branch, nodes) in tqdm(
            zip(branch_order, self.nodes_to_branch_map.items()),
            total=len(branch_order),
            desc="Getting bedlevels",
        ):
            if order == -1:
                # No branch order so just get upstream and downstream levels
                branch_length = self.branches.loc[branch, "geometry"].length
                subset = bedlevels_crs_branches.loc[
                    bedlevels_crs_branches["branchid"] == branch
                ]
                if subset.empty:
                    continue  # if this happens, crs is not defined. This can be a problem.
                nodes_dict[nodes[0]]["up"].append(
                    np.interp(0.0, subset["chainage"], subset["minz"])
                )
                nodes_dict[nodes[1]]["down"].append(
                    np.interp(branch_length, subset["chainage"], subset["minz"])
                )
            else:
                # In case of branch order, first collect all branches and set in them in order of up- to downstream
                all_branches = [
                    self.mesh1d.description1d["network_branch_ids"][i]
                    for i in np.argwhere(order == branch_order).ravel()
                ]
                all_nodes = [self.nodes_to_branch_map[b] for b in all_branches]
                # First check if any of the branches has a bedlevel from a cross-section profile otherwise skip
                check = all(
                    [
                        bedlevels_crs_branches.loc[
                            bedlevels_crs_branches["branchid"] == b
                        ].empty
                        for b in all_branches
                    ]
                )
                if check:
                    continue  # if this happens, cross-section is not defined. This can be a problem.

                # Check if every branch is from up to down direction. Otherwise fix by reversing
                n = 0
                n_length = len(all_nodes)
                direction = list(np.ones(n_length))
                max_tries = 0
                while n < n_length:
                    up = (
                        np.count_nonzero(
                            all_nodes[n][0] == np.array([x[0] for x in all_nodes])
                        )
                        == 1
                    )
                    down = (
                        np.count_nonzero(
                            all_nodes[n][1] == np.array([x[1] for x in all_nodes])
                        )
                        == 1
                    )
                    if (not up) or (not down):
                        # Reverse
                        all_nodes[n] = list(np.flip(all_nodes[n]))
                        direction[n] = direction[n] * -1
                    n += 1
                    # Check if indeed everything is now in proper direction. Otherwise try again
                    if n == n_length:
                        up = all(
                            [
                                np.count_nonzero(
                                    node[0] == np.array([x[0] for x in all_nodes])
                                )
                                == 1
                                for node in all_nodes
                            ]
                        )
                        down = all(
                            [
                                np.count_nonzero(
                                    node[1] == np.array([x[1] for x in all_nodes])
                                )
                                == 1
                                for node in all_nodes
                            ]
                        )
                        if (not up) or (not down):
                            n = 0
                            max_tries += 1
                    if max_tries > 500:
                        print(
                            f"Can't fix correct directions branch groups {all_branches}"
                        )
                        break

                # Add reserved branches to return
                reserved_branches.extend(
                    [b for b, d in zip(all_branches, direction) if d == -1]
                )

                # Get most upstream node. Otherwise just pick i_upstream = 0 as starting point
                i_upstream = [
                    i
                    for i, n in enumerate([x[0] for x in all_nodes])
                    if n not in [x[1] for x in all_nodes]
                ]
                if len(i_upstream) == 1:
                    i_upstream = i_upstream[0]
                else:
                    # It could be that branch order group forms a ring. In this case check first which node has more
                    # than 2 branches (bifurcation) or just 1 branch (boundary) connected.
                    i_upstream = [
                        i
                        for i, n in enumerate([x[0] for x in all_nodes])
                        if (len(self.branches_to_node_map[n]) > 2)
                        or (len(self.branches_to_node_map[n]) == 1)
                    ]
                    if len(i_upstream) == 1:
                        i_upstream = i_upstream[0]
                    else:
                        raise ValueError(
                            f"Something is not right with the branch order group {all_branches}"
                        )

                # Now put branch list in correct order
                all_branches_sorted = []
                all_nodes_sorted = []
                direction_sorted = []
                for ii in range(len(all_branches)):
                    all_branches_sorted.append(all_branches[i_upstream])
                    all_nodes_sorted.append(all_nodes[i_upstream])
                    direction_sorted.append(direction[i_upstream])
                    try:
                        i_upstream = [
                            i
                            for i, n in enumerate([x[0] for x in all_nodes])
                            if [x[1] for x in all_nodes][i_upstream] == n
                        ][0]
                    except IndexError:
                        break
                # Stitch chainages and bedlevels together
                chainage, bedlevel = [], []
                branch_length = 0
                for b, d in zip(all_branches_sorted, direction_sorted):
                    subset = bedlevels_crs_branches.loc[
                        bedlevels_crs_branches["branchid"] == b
                    ]
                    chain, bed = subset["chainage"], subset["minz"]
                    # Reverse chainage and bedlevel arrays
                    if d == -1:
                        chain = np.flip(chain)
                        bed = np.flip(bed)
                    chainage.extend(chain + branch_length)
                    bedlevel.extend(bed)
                    branch_length = self.branches.loc[b, "geometry"].length
                # Get chainage of up- and downstream node of loop branch within the overall branch
                if len(all_branches_sorted) == 1:
                    up_node_chainage = 0
                    down_node_chainage = self.branches.loc[
                        all_branches_sorted[0], "geometry"
                    ].length
                else:
                    i = np.argmax(
                        [
                            1 if ((nodes == n) or (list(np.flip(nodes)) == n)) else 0
                            for n in all_nodes_sorted
                        ]
                    )
                    up_node_chainage = sum(
                        [0]
                        + [
                            self.branches.loc[b, "geometry"].length
                            for b, n in zip(
                                all_branches_sorted[:-1], all_nodes_sorted[:-1]
                            )
                        ][: i + 1]
                    )
                    down_node_chainage = sum(
                        [
                            self.branches.loc[b, "geometry"].length
                            for b, n in zip(all_branches_sorted, all_nodes_sorted)
                        ][: i + 1]
                    )
                # Finally interpolate
                nodes_dict[nodes[0]]["up"].append(
                    np.interp(up_node_chainage, chainage, bedlevel)
                )
                nodes_dict[nodes[1]]["down"].append(
                    np.interp(down_node_chainage, chainage, bedlevel)
                )

        # Summarize everything and save
        nodes = list(nodes_dict.keys())
        node_geom = [
            Point(x, y)
            for x, y in zip(
                self.mesh1d.get_values("nnodex"), self.mesh1d.get_values("nnodey")
            )
        ]
        if resolve_at_bifurcation_method == "min":
            upstream_bedlevel = [
                np.min(v["up"]) if len(v["up"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
            downstream_bedlevel = [
                np.min(v["down"]) if len(v["down"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
        elif resolve_at_bifurcation_method == "max":
            upstream_bedlevel = [
                np.max(v["up"]) if len(v["up"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
            downstream_bedlevel = [
                np.max(v["down"]) if len(v["down"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
        elif resolve_at_bifurcation_method == "mean":
            upstream_bedlevel = [
                np.average(["up"]) if len(v["up"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
            downstream_bedlevel = [
                np.average(v["down"]) if len(v["down"]) > 0 else np.nan
                for v in nodes_dict.values()
            ]
        else:
            raise NotImplementedError

        self.nodes = gpd.GeoDataFrame(
            index=nodes_dict.keys(),
            data={
                "code": nodes_dict.keys(),
                "upstream_bedlevel": upstream_bedlevel,
                "downstream_bedlevel": downstream_bedlevel,
            },
            geometry=node_geom,
        )

        if return_reversed_branches:
            return list(np.unique(reserved_branches))

    def get_grouped_branches(self):
        """
        Get grouped branch ids to use in set_branch_order function
        """
        # Get all network data
        branch_ids = self.mesh1d.description1d["network_branch_ids"]
        node_ids = self.mesh1d.description1d["network_node_ids"]
        branch_edge_nodes_idx = self.mesh1d.get_values("nedge_nodes", as_array=True)
        # Collect all node ids per branch and all branches per node id
        self.make_nodes_to_branch_map()
        self.make_branches_to_node_map()

        branch_ids_checked = []
        groups = {0: []}
        for branch_id in branch_ids:
            if branch_id in branch_ids_checked:
                continue

            connected_nodes = self.nodes_to_branch_map[branch_id]  # get connected nodes
            for n in connected_nodes:
                b = self.branches_to_node_map[n]  # get connected branches
                if len(b) > 2:
                    continue  # in this case there's a bifurcation so skip
                elif len(b) == 1:
                    groups[list(groups.keys())[-1] + 1] = b  # b is already a list
                    branch_ids_checked.extend(b)
                    continue  # in this case the branch is not connected to other branches. make separate group
                else:
                    # remove branch from connected_branches because we are not interested in it
                    b = [_b for _b in b if _b != branch_id][0]
                    if b in branch_ids_checked:
                        # connected branch is already added to a group, so this means that branch should be added to
                        # that group
                        groups[[k for k, v in groups.items() if b in v][0]].append(
                            branch_id
                        )
                        branch_ids_checked.extend([branch_id])
                    else:
                        groups[list(groups.keys())[-1] + 1] = [
                            branch_id
                        ]  # otherwise add to group
                        branch_ids_checked.extend([branch_id])
                branch_ids_checked = list(np.unique(branch_ids_checked))
        groups.pop(0)  # remove the 0th group because empty

        # The branches are grouped but not fully connected (although routine above should do the trick in theory).
        # Try merging groups if a branch is in multiple groups.
        # In that case we know that branch groups are connected to eachother and is in fact a bigger group
        for b in branch_ids:
            _groups = groups.copy()  # make copy to apply changes on
            _k = -1  # default index
            # Loop over groups
            for k, v in groups.items():
                if b in v:  # if branch in grouped branches
                    if _k == -1:
                        _k = k  # set index to first group found
                    else:
                        _groups[_k].extend(
                            v
                        )  # otherwise add group to first found group
                        _groups[_k] = list(
                            np.unique(_groups[_k])
                        )  # remove duplicates due to add
                        _groups.pop(k)  # and remove group from groups
            groups = _groups.copy()  # copy changed dict over original
        # One pass over all branches should be sufficient to group everything together. Otherwise raise error
        if (
            max(
                [
                    sum([1 if b in v else 0 for k, v in groups.items()])
                    for b in branch_ids
                ]
            )
            > 1
        ):
            raise ValueError(
                f"Still branches contained in multiple groups. Maximum number of groups where this "
                f"happens: {max([sum([1 if b in v else 0 for k, v in groups.items()]) for b in branch_ids])}"
            )

        # save
        self.branch_groups = groups.copy()

    # generate network and 1d mesh
    def generate_1dnetwork(
        self,
        one_d_mesh_distance: float = 40.0,
        seperate_structures: bool = True,
        max_dist_to_struc: Union[float, None] = None,
        urban_branches: Union[List[str], None] = None,
    ) -> None:
        """Generate a 1D mesh based on the preferred node distance and the loaded branches.
        The node positions are determined based on three criteria:
        - First, such that the distance in between nodes is closest to "one_d_mesh_distance"
        - Second, at least one node should be added in between two subsequent structures. If
          you don't want this, set "seperate_structures" to False
        - Third, a maximum distance between a node and a structure can be specified with
          "max_dist_to_struc"

        Parameters
        ----------
        one_d_mesh_distance : float, optional
            Preferred distance in between nodes, by default 40.0
        seperate_structures : bool, optional
            Whether to add a node in between subsequent structures, by default True
        max_dist_to_struc : Union[float, None], optional
            A maximum distance between a node and a structure, by default None
        urban_branches : Union[List[str], None], optional
            A list of branch IDs for which only nodes are placed at the end points, by default None
        """

        if self.branches.empty:
            raise ValueError(
                "Branches should be added before 1d network can be generated."
            )

        checks.check_argument(one_d_mesh_distance, "one_d_mesh_distance", (float, int))

        # Temporary dictionary to store the id number of the nodes and branches
        nodes = []
        edge_nodes_dict = {}

        # Check if any structures present (if not, structures will be None)
        structures = self.dflowfmmodel.structures.as_dataframe(
            generalstructures=True,
            weirs=True,
            bridges=True,
            culverts=True,
            pumps=True,
            uweirs=True,
            orifices=True,
            compounds=True,
        )

        # If offsets are not predefined, generate them base on one_d_mesh_distance
        if not self.offsets:
            self.generate_offsets(
                one_d_mesh_distance,
                structures=structures if seperate_structures else None,
                max_dist_to_struc=max_dist_to_struc,
                urban_branches=urban_branches,
            )

        # Add the network data to the 1d mesh structure
        sorted_branches = self.branches.iloc[self.branches.length.argsort().values]

        # Add network branch data
        dimensions = self.mesh1d.meshgeomdim
        dimensions.nbranches = len(sorted_branches)
        self.mesh1d.set_values(
            "nbranchorder", (np.ones(dimensions.nbranches, dtype=int) * -1).tolist()
        )
        self.mesh1d.set_values(
            "nbranchlengths", sorted_branches.geometry.length + 1e-12
        )
        self.mesh1d.description1d["network_branch_ids"] = sorted_branches.index.astype(
            str
        ).tolist()
        self.mesh1d.description1d[
            "network_branch_long_names"
        ] = sorted_branches.index.astype(str).tolist()

        # Add network branch geometry
        coords = [line.coords[:] for line in sorted_branches.geometry]
        geomx, geomy = list(zip(*list(itertools.chain(*coords))))[:2]
        dimensions.ngeometry = len(geomx)
        self.mesh1d.set_values("nbranchgeometrynodes", [len(lst) for lst in coords])
        self.mesh1d.set_values("ngeopointx", geomx)
        self.mesh1d.set_values("ngeopointy", geomy)

        branch_names = sorted_branches.index.astype(str).tolist()
        branch_longnames = ["long_" + s for s in branch_names]

        network_edge_nodes = []
        mesh1d_edge_nodes = []
        mesh1d_node_branchidx = []
        mesh1d_node_branchoffset = []
        mesh1d_edge_branchidx = []
        mesh1d_edge_branchoffset = []
        mesh1d_edge_x = []
        mesh1d_edge_y = []
        mesh1d_node_names = []

        # For each branch
        for i_branch, branch in enumerate(sorted_branches.itertuples()):

            # Get branch coordinates. Round to 6 decimals to make sure that, due to precision errors in coordinates,
            # already existing first and/or end nodes are not recognized
            points = shapely.wkt.loads(
                shapely.wkt.dumps(branch.geometry, rounding_precision=6)
            ).coords[:]

            # Network edge node administration
            # -------------------------------
            first_point = points[0]
            last_point = points[-1]

            # Get offsets from dictionary
            offsets = self.offsets[branch.Index]
            # The number of links on the branch
            nlinks = len(offsets) - 1

            # also get the offsets of the edge nodes, halfway the segments
            edge_offsets = [
                (offsets[i] + offsets[i + 1]) / 2.0
                for i in range(np.max([1, len(offsets) - 1]))
            ]

            # Check if the first and last point of the branch are already in the set
            if first_point not in nodes:
                first_present = False
                nodes.append(first_point)
            else:
                first_present = True
                offsets = offsets[1:]

            if last_point not in nodes:
                last_present = False
                nodes.append(last_point)
            else:
                last_present = True
                offsets = offsets[:-1]

            # If no points remain, add an extra halfway: each branch should have at least 1 node
            if len(offsets) == 0:
                if urban_branches is not None and branch.Index not in urban_branches:
                    offsets = np.array([branch.geometry.length / 2.0])
                    edge_offsets = np.array(
                        [i * branch.geometry.length for i in [0.25, 0.75]]
                    )
                    nlinks += 1
                else:
                    offsets = self.offsets[branch.Index]

            # Get the index of the first and last node in the dictionary (1 based, so +1)
            i_from = nodes.index(first_point) + 1
            i_to = nodes.index(last_point) + 1
            if i_from == i_to:
                raise ValueError(
                    f"For {branch.Index} a ring geometry was found: start and end node are the same. Ring geometries are not accepted."
                )

            network_edge_nodes.append([i_from, i_to])

            # Mesh1d edge node administration
            # -------------------------------
            # First determine the start index. This is equal to the number of already present points (+1, since 1 based)
            start_index = len(mesh1d_node_branchidx) + 1
            # For each link, create a new edge node connection
            if first_present:
                start_index -= 1
            new_edge_nodes = [
                [start_index + i, start_index + i + 1] for i in range(nlinks)
            ]
            # If the first node is present, change the first point of the first edge to the existing point
            if first_present:
                new_edge_nodes[0][0] = edge_nodes_dict[first_point]
            else:
                edge_nodes_dict[first_point] = new_edge_nodes[0][0]
            # If the last node is present, change the last point of the last edge too
            if last_present:
                new_edge_nodes[-1][1] = edge_nodes_dict[last_point]
            else:
                edge_nodes_dict[last_point] = new_edge_nodes[-1][1]
            # Add to edge_nodes
            mesh1d_edge_nodes.extend(new_edge_nodes)

            # Update number of nodes
            mesh_point_names = [f"{branch.Index}_{offset:.2f}" for offset in offsets]

            # Append ids, longnames, branch and offset
            self.mesh1d.description1d["mesh1d_node_ids"].extend(mesh_point_names)
            self.mesh1d.description1d["mesh1d_node_long_names"].extend(mesh_point_names)
            mesh1d_node_branchidx.extend([i_branch + 1] * len(offsets))
            mesh1d_node_branchoffset.extend(offsets.tolist())

            lengths = np.r_[
                0.0,
                np.cumsum(
                    np.hypot(
                        np.diff(coords[i_branch], axis=0)[:, 0],
                        np.diff(coords[i_branch], axis=0)[:, 1],
                    )
                ),
            ]
            edge_x = []
            edge_y = []
            for i_edge, edge in enumerate(edge_offsets):
                closest = (np.abs(lengths - edge)).argmin()
                if lengths[closest] > edge:
                    closest = closest - 1
                edge_x.append(
                    edge
                    / (lengths[closest + 1] - lengths[closest])
                    * (coords[i_branch][closest + 1][0] - coords[i_branch][closest][0])
                    + coords[i_branch][closest][0]
                )
                edge_y.append(
                    edge
                    / (lengths[closest + 1] - lengths[closest])
                    * (coords[i_branch][closest + 1][1] - coords[i_branch][closest][1])
                    + coords[i_branch][closest][1]
                )

            mesh1d_edge_branchidx.extend([i_branch + 1] * len(edge_offsets))
            mesh1d_edge_branchoffset.extend(edge_offsets)
            mesh1d_edge_x.extend(edge_x)
            mesh1d_edge_y.extend(edge_y)

        # Parse nodes
        dimensions.nnodes = len(nodes)
        nodex, nodey = list(zip(*nodes))[:2]
        self.mesh1d.set_values("nnodex", nodex)
        self.mesh1d.set_values("nnodey", nodey)
        self.mesh1d.description1d["network_node_ids"].extend(
            [f"{xy[0]:.6f}_{xy[1]:.6f}" for xy in nodes]
        )
        self.mesh1d.description1d["network_node_long_names"].extend(
            [f"x={xy[0]:.6f}_y={xy[1]:.6f}" for xy in nodes]
        )

        # Add edge node data to mesh
        self.mesh1d.set_values("nedge_nodes", np.ravel(network_edge_nodes))
        self.mesh1d.meshgeomdim.numedge = len(mesh1d_edge_nodes)
        self.mesh1d.set_values("edge_nodes", np.ravel(mesh1d_edge_nodes))

        self.mesh1d.edge_branchidx = mesh1d_edge_branchidx
        self.mesh1d.edge_branchoffset = mesh1d_edge_branchoffset
        self.mesh1d.edge_x = mesh1d_edge_x
        self.mesh1d.edge_y = mesh1d_edge_y

        # Add mesh branchidx and offset to mesh
        dimensions.numnode = len(mesh1d_node_branchidx)
        self.mesh1d.set_values("branchidx", mesh1d_node_branchidx)
        self.mesh1d.set_values("branchoffsets", mesh1d_node_branchoffset)

        # Process the 1d network (determine x and y locations) and determine schematised branches
        schematised, _ = self.mesh1d.process_1d_network()
        for idx, geometry in schematised.items():
            self.schematised.at[idx, "geometry"] = geometry

    def _generate_1d_spacing(self, anchor_pts, one_d_mesh_distance):
        """
        Generates 1d distances, called by function generate offsets
        """
        offsets = []
        for i in range(len(anchor_pts) - 1):
            section_length = anchor_pts[i + 1] - anchor_pts[i]
            if section_length <= 0.0:
                raise ValueError("Section length must be larger than 0.0")
            nnodes = max(2, int(round(section_length / one_d_mesh_distance) + 1))
            offsets.extend(
                np.linspace(
                    anchor_pts[i], anchor_pts[i + 1], nnodes - 1, endpoint=False
                ).tolist()
            )
        offsets.append(anchor_pts[-1])

        return np.asarray(offsets)

    def generate_offsets(
        self,
        one_d_mesh_distance: float,
        structures: Union[ExtendedGeoDataFrame, None] = None,
        max_dist_to_struc: Union[float, None] = None,
        urban_branches: Union[List[str], None] = None,
    ):
        """
        Method to generate 1d network grid point locations. The distances are generated
        based on the 1d mesh distance and anchor points. The anchor points can for
        example be structures; every structure should be seperated by a gridpoint.

        Parameters
        ----------
        one_d_mesh_distance : float
            Preferred distance in between nodes
        structures : Union[ExtendedGeoDataFrame, None], optional
            Whether to add a node in between subsequent structures, by default None
        max_dist_to_struc : Union[float, None], optional
            A maximum distance between a node and a structure, by default None
        urban_branches : Union[List[str], None], optional
            A list of branch IDs for which only nodes are placed at the end points, by default None

        """
        # For each branch
        for branch in self.branches.itertuples():
            # Distribute points along network [1d mesh]
            if urban_branches is not None:
                if branch.Index in urban_branches:
                    offsets = self._generate_1d_spacing(
                        [0.0, branch.geometry.length], branch.geometry.length
                    )
                else:
                    offsets = self._generate_1d_spacing(
                        [0.0, branch.geometry.length], one_d_mesh_distance
                    )
            else:
                offsets = self._generate_1d_spacing(
                    [0.0, branch.geometry.length], one_d_mesh_distance
                )
            self.offsets[branch.Index] = offsets

        if structures is not None:
            # Check argument
            checks.check_argument(
                structures,
                "structures",
                (pd.DataFrame, gpd.GeoDataFrame),
                columns=["branchid", "chainage"],
            )

            # Get structure data from dfs
            ids_offsets = structures[["branchid", "chainage"]]
            idx = structures["branchid"] != ""
            if idx.any():
                logger.warning("Some structures are not linked to a branch.")
            ids_offsets = ids_offsets.loc[idx, :]

            # For each branch
            for branch_id, group in ids_offsets.groupby("branchid"):

                # Check if structures are located at the same offset
                u, c = np.unique(group["chainage"], return_counts=True)
                if any(c > 1):
                    logger.warning(
                        "Structures {} have the same location.".format(
                            ", ".join(
                                group.loc[
                                    np.isin(group["chainage"], u[c > 1])
                                ].index.tolist()
                            )
                        )
                    )

                branch = self.branches.at[branch_id, "geometry"]
                # Limits are the lengths at which the structures are located
                limits = sorted(group["chainage"].unique())

                anchor_pts = [0.0, branch.length]
                offsets = self._generate_1d_spacing(anchor_pts, one_d_mesh_distance)

                # Merge limits with start and end of branch
                limits = [-1e-3] + limits + [branch.length + 1e-3]

                # If any structures
                if len(limits) > 2:

                    # also check if the calculation point are close enough to the structures
                    if max_dist_to_struc is not None:
                        additional = []

                        # Skip the first and the last, these are no structures
                        for i in range(1, len(limits) - 1):
                            # if the distance between two limits is large than twice the max distance to structure,
                            # the mesh point will be too far away. Add a limit on the minimum of half the length and
                            # two times the max distance
                            dist_to_prev_limit = limits[i] - (
                                max(additional[-1], limits[i - 1])
                                if any(additional)
                                else limits[i - 1]
                            )
                            if dist_to_prev_limit > 2 * max_dist_to_struc:
                                additional.append(
                                    limits[i]
                                    - min(2 * max_dist_to_struc, dist_to_prev_limit / 2)
                                )

                            dist_to_next_limit = limits[i + 1] - limits[i]
                            if dist_to_next_limit > 2 * max_dist_to_struc:
                                additional.append(
                                    limits[i]
                                    + min(2 * max_dist_to_struc, dist_to_next_limit / 2)
                                )

                        # Join the limits
                        limits = sorted(limits + additional)

                    # Get upper and lower limits
                    upper_limits = limits[1:]
                    lower_limits = limits[:-1]

                    # Determine the segments that are missing a grid point
                    in_range = [
                        ((offsets > lower) & (offsets < upper)).any()
                        for lower, upper in zip(lower_limits, upper_limits)
                    ]

                    while not all(in_range):
                        # Get the index of the first segment without grid point
                        i = in_range.index(False)

                        # Add it to the anchor pts
                        anchor_pts.append((lower_limits[i] + upper_limits[i]) / 2.0)
                        anchor_pts = sorted(anchor_pts)

                        # Generate new offsets
                        offsets = self._generate_1d_spacing(
                            anchor_pts, one_d_mesh_distance
                        )

                        # Determine the segments that are missing a grid point
                        in_range = [
                            ((offsets > lower) & (offsets < upper)).any()
                            for lower, upper in zip(lower_limits, upper_limits)
                        ]

                    if len(anchor_pts) > 2:
                        logger.info(
                            f"Added 1d mesh nodes on branch {branch_id} at: {anchor_pts}, due to the structures at {limits}."
                        )

                # Set offsets for branch id
                self.offsets[branch_id] = offsets

    def get_roughness_description(self, roughnesstype, value):

        if np.isnan(float(value)):
            raise ValueError("Roughness value should not be NaN.")

        # Check input
        # checks.check_argument(roughnesstype, 'roughness type', (str, int))
        # checks.check_argument(value, 'roughness value', (float, int, np.float, np.integer))

        # Convert integer to string
        # map HyDAMO definition to D-Hydro definition
        roughnesstype = self.dflowfmmodel.roughness_mapping[roughnesstype]

        # Get name
        name = f"{roughnesstype}_{float(value)}"

        # Check if the description is already known
        if name.lower() in map(str.lower, self.roughness_definitions.keys()):
            return name

        # Convert roughness type string to integer for dflowfm
        delft3dfmtype = roughnesstype

        # Add to dict
        self.roughness_definitions[name] = {
            "name": name,
            "code": delft3dfmtype,
            "value": value,
        }

        return name

    def add_mesh2d(self, twodmesh):
        """
        Add 2d mesh to network object.
        """
        if not hasattr(twodmesh, "meshgeom"):
            checks.check_argument(twodmesh, "twodmesh", meshgeom)
            geometries = twodmesh
        else:
            if not isinstance(twodmesh.meshgeom, meshgeom):
                raise TypeError('The given mesh should have an attribute "meshgeom".')
            geometries = twodmesh.meshgeom

        # Add the meshgeom
        self.mesh2d.add_from_other(geometries)

        # Add mdu parameters
        if hasattr(twodmesh, "missing_z_value"):
            if twodmesh.missing_z_value is not None:
                self.mdu_parameters["BedlevUni"] = twodmesh.missing_z_value

    def get_node_idx_offset(self, branch_id, pt, nnodes=1):
        """
        Get the index and offset of a node on a 1d branch.
        The nearest node is looked for.
        """

        # Project the point on the branch
        dist = self.schematised[branch_id].project(pt)

        # Get the branch data from the networkdata
        branchidx = (
            self.mesh1d.description1d["network_branch_ids"].index(
                self.str2chars(branch_id, self.idstrlength)
            )
            + 1
        )
        pt_branch_id = self.mesh1d.get_values("branchidx", as_array=True)
        idx = np.where(pt_branch_id == branchidx)

        # Find nearest offset
        offsets = self.mesh1d.get_values("branchoffset", as_array=True)[idx]
        isorted = np.argsort(np.absolute(offsets - dist))
        isorted = isorted[: min(nnodes, len(isorted))]

        # Get the offset
        offset = [offsets[imin] for imin in isorted]
        # Get the id of the node
        node_id = [idx[0][imin] + 1 for imin in isorted]

        return node_id, offset


class Structures:
    def __init__(self, dflowfmmodel):
        self.generalstructures = {}
        self.pumps = {}
        self.weirs = {}
        self.bridges = {}
        self.uweirs = {}
        self.culverts = {}
        self.orifices = {}
        self.compounds = {}

        self.dflowfmmodel = dflowfmmodel

        # Create the io class
        self.io = dfmreader.StructuresIO(self)

    def add_generalstructure(
        self,
        id,
        branchid,
        chainage,
        name=np.nan,
        allowedflowdir="both",
        upstream1width=np.nan,
        upstream1level=np.nan,
        upstream2width=np.nan,
        upstream2level=np.nan,
        crestwidth=np.nan,
        crestlevel=np.nan,
        crestlength=np.nan,
        downstream1width=np.nan,
        downstream1level=np.nan,
        downstream2width=np.nan,
        downstream2level=np.nan,
        gateloweredgelevel=10e10,
        posfreegateflowcoeff=1.0,
        posdrowngateflowcoeff=1.0,
        posfreeweirflowcoeff=1.0,
        posdrownweirflowcoeff=1.0,
        poscontrcoeffreegate=1.0,
        negfreegateflowcoeff=1.0,
        negdrowngateflowcoeff=1.0,
        negfreeweirflowcoeff=1.0,
        negdrownweirflowcoeff=1.0,
        negcontrcoeffreegate=1.0,
        extraresistance=0.0,
        gateheight=10e10,
        gateopeningwidth=0,
        gateopeninghorizontaldirection="symmetric",
        usevelocityheight="true",
    ):
        self.generalstructures[id] = remove_nan_values(
            {
                "type": "generalStructure",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "upstream1Width": upstream1width,
                "upstream1Level": upstream1level,
                "upstream2Width": upstream2width,
                "upstream2Level": upstream2level,
                "crestWidth": crestwidth,
                "crestLevel": crestlevel,
                "crestLength": crestlength,
                "downstream1Width": downstream1width,
                "downstream1Level": downstream1level,
                "downstream2Width": downstream2width,
                "downstream2Level": downstream2level,
                "gateLowerEdgeLevel": gateloweredgelevel,
                "posFreeGateFlowCoeff": posfreegateflowcoeff,
                "posDrownGateFlowCoeff": posdrowngateflowcoeff,
                "posFreeWeirFlowCoeff": posfreeweirflowcoeff,
                "posDrownWeirFlowCoeff": posdrownweirflowcoeff,
                "posContrCoefFreeGate": poscontrcoeffreegate,
                "negFreeGateFlowCoeff": negfreegateflowcoeff,
                "negDrownGateFlowCoeff": negdrowngateflowcoeff,
                "negFreeWeirFlowCoeff": negfreeweirflowcoeff,
                "negDrownWeirFlowCoeff": negdrownweirflowcoeff,
                "negContrCoefFreeGate": negcontrcoeffreegate,
                "extraResistance": extraresistance,
                "gateHeight": gateheight,
                "gateOpeningWidth": gateopeningwidth,
                "gateOpeningHorizontalDirection": gateopeninghorizontaldirection,
                "useVelocityHeight": usevelocityheight,
            }
        )

    def add_pump(
        self,
        id,
        branchid,
        chainage,
        orientation,
        numstages,
        controlside,
        capacity,
        name=np.nan,
        startlevelsuctionside=np.nan,
        stoplevelsuctionside=np.nan,
        startleveldeliveryside=np.nan,
        stopleveldeliveryside=np.nan,
    ):
        if controlside.lower() == "suctionside":
            pass
        elif controlside.lower() == "deliveryside":
            pass
        elif controlside.lower() != "both":
            raise ValueError(
                "Incorrect controlSide value specified. Either use 'suctionSide', 'deliverySide' or 'both'."
            )
        self.pumps[id] = remove_nan_values(
            {
                "type": "pump",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "orientation": orientation,
                "numstages": numstages,
                "controlSide": controlside,
                "capacity": capacity,
                "startlevelSuctionSide": startlevelsuctionside,
                "stoplevelSuctionSide": stoplevelsuctionside,
                "startlevelDeliverySide": startleveldeliveryside,
                "stoplevelDeliverySide": stopleveldeliveryside,
            }
        )

    def add_weir(
        self,
        id,
        branchid,
        chainage,
        crestlevel,
        crestwidth,
        name=np.nan,
        allowedflowdir="both",
        corrcoeff=1.0,
        usevelocityheight="true",
    ):
        self.weirs[id] = remove_nan_values(
            {
                "type": "weir",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "allowedFlowDir": allowedflowdir,
                "crestLevel": crestlevel,
                "crestWidth": crestwidth,
                "corrCoeff": corrcoeff,
                "useVelocityHeight": usevelocityheight,
            }
        )

    def add_orifice(
        self,
        id,
        branchid,
        chainage,
        crestlevel,
        crestwidth,
        gateloweredgelevel,
        name=np.nan,
        allowedflowdir="both",
        uselimitflowpos=False,
        limitflowpos=0.0,
        uselimitflowneg=False,
        limitflowneg=0.0,
        corrcoeff=1.0,
        usevelocityheight="true",
    ):
        self.orifices[id] = remove_nan_values(
            {
                "type": "orifice",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "allowedFlowDir": allowedflowdir,
                "crestLevel": crestlevel,
                "crestWidth": crestwidth,
                "gateLowerEdgeLevel": gateloweredgelevel,
                "useLimitFlowPos": uselimitflowpos,
                "limitFlowPos": limitflowpos,
                "useLimitFlowNeg": uselimitflowneg,
                "limitFlowNeg": limitflowneg,
                "corrCoeff": corrcoeff,
                "useVelocityHeight": usevelocityheight,
            }
        )

    def add_bridge(
        self,
        id,
        branchid,
        chainage,
        length,
        shift,
        crosssection,
        inletlosscoeff,
        outletlosscoeff,
        name=np.nan,
        allowedflowdir="both",
        frictiontype="StricklerKs",
        frictionvalue=75.0,
    ):
        """
        Add a bridge to the schematisation.

        Note that the cross section should be handed as dictionary. This should contain the
        shape (circle, rectangle) and the required arguments.

        For now, we use the smae cross sections as culverts. This needs to be refined.

        """

        # Check the content of the cross section dictionary
        # checks.check_dictionary(crosssection, required='shape', choice=['diameter', ['width', 'height', 'closed']])
        # map frictiontype from HyDAMO to D-Hydro
        frictiontype = self.dflowfmmodel.roughness_mapping[frictiontype]

        self.bridges[id] = remove_nan_values(
            {
                "type": "bridge",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "allowedFlowDir": allowedflowdir,
                "csDefId": crosssection,
                "shift": shift,
                "inletLossCoeff": inletlosscoeff,
                "outletLossCoeff": outletlosscoeff,
                "frictionType": frictiontype,
                "friction": frictionvalue,
                "length": length,
            }
        )

    def add_uweir(
        self,
        id,
        branchid,
        chainage,
        crestlevel,
        yvalues,
        zvalues,
        name=np.nan,
        allowedflowdir="both",
        dischargecoeff=1.0,
    ):

        """
        Add a universal weir to the schematisation.
        """
        assert len(yvalues.split()) == len(
            zvalues.split()
        ), "Number of y values should equal number of z values"
        self.uweirs[id] = remove_nan_values(
            {
                "type": "universalWeir",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "allowedFlowDir": allowedflowdir,
                "numLevels": len(yvalues.split()),
                "yValues": yvalues,
                "zValues": zvalues,
                "crestLevel": crestlevel,
                "dischargeCoeff": dischargecoeff,
            }
        )

    def add_culvert(
        self,
        id,
        branchid,
        chainage,
        leftlevel,
        rightlevel,
        crosssection,
        length,
        inletlosscoeff,
        outletlosscoeff,
        name=np.nan,
        allowedflowdir="both",
        valveonoff=0,
        numlosscoeff=0,
        valveopeningheight=None,
        relopening=None,
        losscoeff=None,
        frictiontype="StricklerKs",
        frictionvalue=75.0,
    ):
        """
        Add a culvert to the schematisation.

        Note that the cross section should be handed as dictionary. This should contain the
        shape (circle, rectangle) and the required arguments.

        """

        if isinstance(crosssection, dict):
            # Check the content of the cross section dictionary
            checks.check_dictionary(
                crosssection,
                required="shape",
                choice=["diameter", ["width", "height", "closed"]],
            )

            # Get roughnessname
            # map HyDAMO roughness type to D-Hydro roughness type
            frictiontype_mapped = self.dflowfmmodel.roughness_mapping[frictiontype]

            # Add cross section definition
            # WORKAROUND: for the GUI, a profile definition has to be created for every culvert. We do this temporary.
            if crosssection["shape"].lower() == "circle":
                definition = self.dflowfmmodel.crosssections.add_circle_definition(
                    crosssection["diameter"], frictiontype, frictionvalue, name=id
                )
            elif crosssection["shape"].lower() == "rectangle":
                definition = self.dflowfmmodel.crosssections.add_rectangle_definition(
                    crosssection["height"],
                    crosssection["width"],
                    crosssection["closed"],
                    frictiontype,
                    frictionvalue,
                    name=id,
                )
            else:
                NotImplementedError(
                    f'Cross section with shape "{crosssection["shape"]}" not implemented.'
                )
        elif isinstance(crosssection, str):
            # Id of cross-section definition that's already added
            definition = crosssection
        else:
            raise ValueError(
                f"Type of 'crosssection' key is not supported ({type(crosssection)})"
            )

        # Add the culvert to the dictionary
        self.culverts[id] = remove_nan_values(
            {
                "type": "culvert",
                "id": id,
                "name": id if np.isnan(name) else name,
                "branchid": branchid,
                "chainage": chainage,
                "allowedFlowDir": allowedflowdir,
                "leftLevel": leftlevel,
                "rightLevel": rightlevel,
                "csDefId": definition,
                "length": round(length, 3),
                "inletLossCoeff": inletlosscoeff,
                "outletLossCoeff": outletlosscoeff,
                "valveOnOff": int(valveonoff),
                "numLossCoeff": int(numlosscoeff),
                "bedFrictionType": frictiontype_mapped,
                "bedFriction": frictionvalue,
            }
        )
        if valveonoff > 0:
            self.culverts[id]["valveOpeningHeight"] = valveopeningheight
            self.culverts[id]["relOpening"] = relopening
            self.culverts[id]["lossCoeff"] = losscoeff

    def add_compound(self, id, numstructures, structurelist, name=np.nan):
        self.compounds[id] = remove_nan_values(
            {
                "type": "compound",
                "id": id,
                "name": id if np.isnan(name) else name,
                "numStructures": numstructures,
                "structureIds": structurelist,
            }
        )

    def as_dataframe(
        self,
        generalstructures=False,
        pumps=False,
        weirs=False,
        bridges=False,
        culverts=False,
        uweirs=False,
        orifices=False,
        compounds=False,
    ):
        """
        Returns a dataframe with the structures. Specify with the keyword arguments what structure types need to be returned.
        """
        dfs = []
        for df, descr, add in zip(
            [
                self.generalstructures,
                self.culverts,
                self.weirs,
                self.bridges,
                self.pumps,
                self.uweirs,
                self.orifices,
                self.compounds,
            ],
            [
                "generalstructure",
                "culvert",
                "weir",
                "bridge",
                "pump",
                "uweir",
                "orifice",
                "compound",
            ],
            [
                generalstructures,
                culverts,
                weirs,
                bridges,
                pumps,
                uweirs,
                orifices,
                compounds,
            ],
        ):
            if any(df) and add:
                df = pd.DataFrame.from_dict(df, orient="index")
                df.insert(loc=0, column="structype", value=descr, allow_duplicates=True)
                dfs.append(df)

        if len(dfs) > 0:
            return pd.concat(dfs, sort=False)


class ObservationPoints(ExtendedGeoDataFrame):
    def __init__(self, dflowfmmodel):
        super(ObservationPoints, self).__init__(
            geotype=Point,
            required_columns=[
                "name",
                "branchId",
                "chainage",
                "geometry",
                "locationType",
            ],
        )

        self._metadata.append("dflowfmmodel")
        self.dflowfmmodel = dflowfmmodel

    def add_points(self, crds, names, locationTypes=None, snap_distance=None):
        """
        Method to add observation points to schematisation. Observation points can be of type '1d' or '2d'. 1d-points are snapped to the branch.

        Parameters
        ----------
        crds : Nx2 list or array
            x and y coordinates of observation points
        names : str or list
            names of the observation points
        locationTypes:  str or list
            type of the observationpoints: 1d or 2d
        """
        if snap_distance is None:
            snap_distance = 5
        if isinstance(names, str):
            names = [names]
            crds = [crds]
        if locationTypes is not None:
            if isinstance(names, str):
                locationTypes = [locationTypes]
            # split 1d and 2d points, as the first ones need to be snapped to branches
            obs2d = gpd.GeoDataFrame()
            obs2d["name"] = [
                n for nn, n in enumerate(names) if locationTypes[nn] == "2d"
            ]
            obs2d["locationType"] = "2d"
            # obs2d['geometry'] = [Point(*pt) for ipt,pt in enumerate(crds) if (locationTypes[ipt]=='2d')&(not isinstance(pt, Point))]
            obs2d["geometry"] = [
                Point(*pt) if not isinstance(pt, Point) else pt
                for ipt, pt in enumerate(crds)
                if (locationTypes[ipt] == "2d")
            ]
            obs2d["x"] = [pt.coords[0][0] for pt in obs2d["geometry"]]
            obs2d["y"] = [pt.coords[0][1] for pt in obs2d["geometry"]]
            names1d = [n for n_i, n in enumerate(names) if locationTypes[n_i] == "1d"]
            crds1d = [c for c_i, c in enumerate(crds) if locationTypes[c_i] == "1d"]
        else:
            names1d = names
            crds1d = crds

        # Check if data for snapping is available
        network = self.dflowfmmodel.network
        if not network.mesh1d.meshgeomdim.numnode:
            raise ValueError(
                "The network geometry should be generated before the observation points can be snapped to 1d."
            )

        obs1d = gpd.GeoDataFrame()
        obs1d["name"] = names1d
        obs1d["geometry"] = [
            Point(*pt) if not isinstance(pt, Point) else pt for pt in crds1d
        ]
        obs1d["locationType"] = "1d"
        geometry.find_nearest_branch(
            network.branches, obs1d, method="overal", maxdist=snap_distance
        )
        obs1d.rename(
            columns={"branch_id": "branchId", "branch_offset": "chainage"}, inplace=True
        )
        obs = obs1d.append(obs2d, sort=True) if locationTypes is not None else obs1d

        # Add to dataframe
        self.set_data(obs, index_col="name", check_columns=True)


class StorageNodes:
    def __init__(self, dflowfmmodel):
        self.storagenodes = {}

        self.dflowfmmodel = dflowfmmodel

        # Create the io class
        self.io = dfmreader.StorageNodesIO(self)

    def add_storagenode(
        self,
        id,
        nodeid,
        usestreetstorage="true",
        nodetype="unspecified",
        name=np.nan,
        usetable="false",
        bedlevel=np.nan,
        area=np.nan,
        streetlevel=np.nan,
        streetstoragearea=np.nan,
        storagetype="reservoir",
        levels=np.nan,
        storagearea=np.nan,
        interpolate="linear",
    ):
        base = {
            "type": "storageNode",
            "id": id,
            "name": id if np.isnan(name) else name,
            "useStreetStorage": usestreetstorage,
            "nodeType": nodetype,
            "nodeId": nodeid,
            "useTable": usetable,
        }
        if usetable == "false":
            out = {
                **base,
                "bedLevel": bedlevel,
                "area": area,
                "streetLevel": streetlevel,
                "streetStorageArea": streetstoragearea,
                "storageType": storagetype,
            }
        elif usetable == "true":
            assert len(levels.split()) == len(
                storagearea.split()
            ), "Number of levels does not equal number of storagearea"
            out = {
                **base,
                "numLevels": len(levels.split()),
                "levels": area,
                "storageArea": streetlevel,
                "interpolate": interpolate,
            }
        else:
            raise ValueError(
                "Value of key 'usetable' is not supported. Either use 'true' or 'false"
            )
        self.storagenodes[id] = remove_nan_values(out)


def list_to_str(lst):
    string = " ".join([f"{number:6.3f}" for number in lst])
    return string


def remove_nan_values(base):
    base_copy = base.copy()
    for k, v in base.items():
        if isinstance(v, float):
            if np.isnan(v):
                base_copy.pop(k)
    return base_copy
