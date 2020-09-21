import netCDF4
import numpy as np

def to_netcdf_old(meshgeom, path):

    outformat = "NETCDF3_CLASSIC"
    ncfile = netCDF4.Dataset(path, 'w', format=outformat)

    ncfile.createDimension("nNetNode", meshgeom.meshgeomdim.numnode)
    ncfile.createDimension("nNetLink", meshgeom.meshgeomdim.numedge)
    ncfile.createDimension("Two", 2)

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
    ncvar = ncfile.createVariable("NetLink", "i4", ( mesh2d.edge_dimension, "Two"))
    links = meshgeom.get_values('edge_nodes', as_array=True)
    ncvar[:] = links.tolist()
    # ncvar

    ncvar = ncfile.createVariable("NetLinkType", "i4", mesh2d.edge_dimension)
    ncvar[:] = np.ones(len(links), dtype=int) * 2

    ncfile.close()

def from_netcdf_old(meshgeom, path):
    ds = netCDF4.Dataset(path, 'r')

    meshgeom.meshgeomdim.numnode = ds.dimensions['nNetNode'].size
    meshgeom.meshgeomdim.numedge = ds.dimensions['nNetLink'].size

    for dim in list('xyz'):
        # Get values
        data = ds.variables[f'NetNode_{dim}'][:]
        # Allocate
        meshgeom.allocate(f'node{dim}')
        # Set values
        meshgeom.set_values(f'node{dim}', data)

    # Set links
    meshgeom.allocate(f'edge_nodes')
    meshgeom.set_values('edge_nodes', ds.variables['NetLink'][:, :].ravel())

    ds.close()