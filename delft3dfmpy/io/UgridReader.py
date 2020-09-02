import netCDF4
import numpy as np
from delft3dfmpy.datamodels.cstructures import meshgeom, meshgeomdim

import logging

logger = logging.getLogger(__name__)

ugrid_dim_dict = {
    "network1d_nEdges": ('1d', 'nbranches'),
    "network1d_nNodes": ('1d', 'nnodes'),
    "network1d_nGeometryNodes": ('1d', 'ngeometry'),    
    "mesh1d_nEdges": ('1d', 'numedge'),
    "mesh1d_nNodes": ('1d', 'numnode'),
    "max_nmesh2d_face_nodes": ('2d', 'maxnumfacenodes'),
    "mesh2d_nEdges": ('2d', 'numedge'),
    "mesh2d_nFaces": ('2d', 'numface'),
    "mesh2d_nNodes": ('2d', 'numnode')
}

ugrid_var_dict = {
    'network1d_node_id': ('1d', 'network_node_ids'),
    'network1d_node_long_name': ('1d', 'network_node_long_names'),
    'network1d_node_x': ('1d', 'nnodex'),
    'network1d_node_y': ('1d', 'nnodey'),
    'network1d_branch_id': ('1d', 'network_branch_ids'),
    'network1d_branch_long_name': ('1d', 'network_branch_long_names'),
    'network1d_edge_length': ('1d', 'nbranchlengths'),
    'network1d_branch_order': ('1d', 'nbranchorder'),
    'network1d_edge_nodes': ('1d', 'nedge_nodes'),
    'network1d_geom_x': ('1d', 'ngeopointx'),
    'network1d_geom_y': ('1d', 'ngeopointy'),
    'network1d_geom_node_count': ('1d', 'nbranchgeometrynodes'),
    'mesh1d_node_id': ('1d', 'mesh1d_node_ids'),
    'mesh1d_node_long_name': ('1d', 'mesh1d_node_long_names'),
    'mesh1d_edge_nodes': ('1d', 'edge_nodes'),
    'mesh1d_node_branch': ('1d', 'branchidx'),
    'mesh1d_node_offset': ('1d', 'branchoffsets'),
    'mesh2d_node_x': ('2d', 'nodex'),
    'mesh2d_node_y': ('2d', 'nodey'),
    'mesh2d_node_z': ('2d', 'nodez'),
    'mesh2d_edge_nodes': ('2d', 'edge_nodes'),
    'mesh2d_face_x': ('2d', 'facex'),
    'mesh2d_face_y': ('2d', 'facey'),
    'mesh2d_face_z': ('2d', 'facez'),
    'mesh2d_face_nodes': ('2d', 'face_nodes')
}

class UgridReader:

    def __init__(self, network):

        self.network = network

    def read_ugrid(self, path):

        """
        Read Ugrid from netcdf and return dflowfm cstructure with grid
        """
        ncfile = netCDF4.Dataset(path)

        # Read mesh1d    
        mesh1d = meshgeom(meshgeomdim())
        read_dimensions(mesh1d.meshgeomdim, '1d', ncfile)
        read_values(mesh1d, '1d', ncfile)
        schematised, branches = mesh1d.process_1d_network()
        self.network.mesh1d.add_from_other(mesh1d)

        # Add branches
        for idx, geometry in branches.items():
            self.network.branches.at[idx, 'geometry'] = geometry
        for idx, geometry in schematised.items():
            self.network.schematised.at[idx, 'geometry'] = geometry

        # Read mesh2d
        mesh2d = meshgeom(meshgeomdim())
        read_dimensions(mesh2d.meshgeomdim, '2d', ncfile)
        read_values(mesh2d, '2d', ncfile)
        self.network.mesh2d.add_from_other(mesh2d)

        # Read links1d2d
        links_1dnode, links_2dface = read_links(ncfile)
        self.network.links1d2d.nodes1d.extend(links_1dnode)
        self.network.links1d2d.faces2d.extend(links_2dface)
        
        ncfile.close()

def read_dimensions(meshgeomdim, readdim, ncfile):
    """
    Function to read dimensions from netcdf file
    """

    assert readdim in ['1d', '2d']

    meshgeomdim.dim = int(readdim[0])

    # Read dimensions
    for ncname, (dim, cname) in ugrid_dim_dict.items():
        if readdim != dim:
            continue

        # Check if variable is in nc file
        if ncname not in ncfile.dimensions.keys():
            logger.error(f'Failed to read dimension "{ncname}" from ncfile for {readdim} mesh.')
        
        value = ncfile.dimensions[ncname].size
        logger.info(f'Read dimension "{ncname}" from ncfile for {readdim} mesh.')
        setattr(meshgeomdim, cname, value)

def read_values(meshgeom, readdim, ncfile):
    """
    Function to read values from netcdf file
    """

    assert readdim in ['1d', '2d']

    # Read variables
    for ncname, (dim, cname) in ugrid_var_dict.items():
        if readdim != dim:
            continue

        # Check if variable is in nc file
        if ncname not in ncfile.variables.keys():
            logger.error(f'Failed to read variable "{ncname}" from ncfile for {readdim} mesh.')
        
        # Read values
        values = ncfile.variables[ncname][:]
        logger.info(f'Read variable "{ncname}" with shape {values.shape} from ncfile for {readdim} mesh.')
        
        # Read description variables (strings)
        if cname in meshgeom.description1d.keys():
            meshgeom.description1d[cname] = list(map(str.strip, netCDF4.chartostring(values)))
        # Read numerical values
        else:
            if values.mask.any():
                values[values.mask] = ncfile.variables[ncname]._FillValue
            meshgeom.set_values(cname, np.ravel(values))
            
def read_links(ncfile):
    """
    Function to read 1d 2d links
    """
    var = 'link1d2d'

    if var not in ncfile.variables.keys():
        return [[], []]
    else:
        return ncfile.variables[var][:, :].T.tolist()