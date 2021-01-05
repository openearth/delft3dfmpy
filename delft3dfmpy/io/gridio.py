import netCDF4
import numpy as np
import sys
sys.path.append('D:/Documents/GitHub/delft3dfmpy')
from delft3dfmpy.datamodels.cstructures import meshgeom, meshgeomdim

import logging

logger = logging.getLogger(__name__)


def to_netcdf_old(meshgeom, path):

    outformat = "NETCDF3_CLASSIC"
    ncfile = netCDF4.Dataset(path, 'w', format=outformat)

    ncfile.createDimension("nNetNode", meshgeom.meshgeomdim.numnode)
    ncfile.createDimension("nNetLink", meshgeom.meshgeomdim.numedge)
    ncfile.createDimension("nNetLinkPts", 2)

    # Mesh2D
    mesh2d = ncfile.createVariable("Mesh2D", "i4", ())
    mesh2d.cf_role = 'mesh_topology'
    mesh2d.node_coordinates = 'NetNode_x NetNode_y'
    mesh2d.node_dimension = 'nNetNode'
    mesh2d.edge_node_connectivity = 'NetLink'
    mesh2d.edge_dimension = 'nNetLink'
    mesh2d.topology_dimension = 1

    # Nodes:
    for dim in list('xyz'):
        ncvar = ncfile.createVariable(f"NetNode_{dim}", "f8", mesh2d.node_dimension)
        ncvar.units = 'm'

        if dim == 'z':
            ncvar.mesh = 'Mesh2D'
            ncvar.coordinates = 'NetNode_x NetNode_y'
            ncvar[:] = np.zeros(meshgeom.meshgeomdim.numnode)

        else:
            ncvar[:] = meshgeom.get_values(f'node{dim}')

    # Links
    ncvar = ncfile.createVariable("NetLink", "i4", ( mesh2d.edge_dimension, "nNetLinkPts"))
    links = meshgeom.get_values('edge_nodes', as_array=True)
    ncvar[:] = links.tolist()
    
    # NetLinkType
    ncvar = ncfile.createVariable("NetLinkType", "i4", mesh2d.edge_dimension)
    ncvar[:] = np.ones(len(links), dtype=int) * 2

    ncfile.close()

def from_netcdf_old(meshgeom, path, only2d=False):
    """
    Method to read mesh from 'old' netcdf 0.9 format
    Function only suitable for 

    Parameters
    ----------
    meshgeom : [type]
        [description]
    path : str
        Path to netcdf with mesh
    """

    ds = netCDF4.Dataset(path, 'r')

    # Get netlinktype
    netlinktype = ds.variables['NetLinkType'][:].data
    links = ds.variables['NetLink'][:, :] - 1

    if only2d and (netlinktype != 2).any():
        linkidx = netlinktype == 2
        nodeidx = np.unique(links[linkidx, :])
        links = links[linkidx, :]

        meshgeom.meshgeomdim.numnode = len(nodeidx)
        meshgeom.meshgeomdim.numedge = sum(linkidx)

        # The links should be renumbered, to compensate the missing ones
        id_mapping = {old_id: new_id for new_id, old_id in enumerate(nodeidx)}
        links = np.reshape([id_mapping[old_id] for old_id in links.ravel()], links.shape)

    else:
        # Get dimensions
        meshgeom.meshgeomdim.numnode = ds.dimensions['nNetNode'].size
        meshgeom.meshgeomdim.numedge = ds.dimensions['nNetLink'].size
        
        nodeidx = slice(None)

    for dim in list('xyz'):
        # Get values
        data = ds.variables[f'NetNode_{dim}'][nodeidx]
        # Allocate
        meshgeom.allocate(f'node{dim}')
        # Set values
        meshgeom.set_values(f'node{dim}', data)

    # Set links
    meshgeom.allocate(f'edge_nodes')
    meshgeom.set_values('edge_nodes', links.ravel() + 1)

    ds.close()
