# coding: utf-8
import os, sys, math, logging
from netCDF4 import Dataset
from collections import OrderedDict
from datetime import *
import numpy as np


class UgridWriter:
    """Writer for FM files"""
    
    def __init__(self):
        self.idstrlength = 40
        self.longstrlength = 80

    def write(self, dflowfmmodel, path, version):  # write ugrid file from GWSW model

        ncfile = self.create_netcdf(path, version)

        # Write the 1d mesh
        if not dflowfmmodel.network.mesh1d.empty():
            self.init_1dnetwork(ncfile, dflowfmmodel.network.mesh1d)
            self.set_1dmesh(ncfile, dflowfmmodel.network.mesh1d)
            self.set_1dnetwork(ncfile, dflowfmmodel.network.mesh1d)

        # Write the 2d mesh
        if not dflowfmmodel.network.mesh2d.empty():
            self.init_2dmesh(ncfile, dflowfmmodel.network.mesh2d)
            self.set_2dmesh(ncfile, dflowfmmodel.network.mesh2d)

        # Write the 1d 2d links
        if dflowfmmodel.network.links1d2d.nodes1d:
            self.init_1d2dlinks(ncfile, dflowfmmodel.network.links1d2d)
            self.set_1d2dlinks(ncfile, dflowfmmodel.network.links1d2d)

        ncfile.close()

    @staticmethod
    def to_char_list(lst, size):
        """Convert list of strings to list of stings with a fixed number of characters"""
        return [item.ljust(size)[:size] for item in lst]

    def create_netcdf(self, path, version):

        # File format:
        outformat = "NETCDF3_CLASSIC" #"NETCDF4"
        # File where we going to write
        ncfile = Dataset(path, 'w', format=outformat)

        # global attributes
        ncfile.Conventions = "CF-1.8 UGRID-1.0"
        ncfile.title = 'Delft3D-FM 1D2D network for model '+os.path.split(path)[-1].rstrip('_net.nc')
        ncfile.source = f"delft3dfmpy v.{version['number']}, D-HyDAMO, model {os.path.split(path)[-1].rstrip('_net.nc')}"
        ncfile.history = f"Created on {version['date']} by {os.path.split(__file__)[-1]}."        
        ncfile.institution = "Deltares/HKV"
        ncfile.references = "https://github.com/openearth/delft3dfmpy/; https://www.deltares.nl; https://www.hkv.nl"      
        ncfile.comment = f"Tested and compatible with D-Flow FM {version['dfm_version']}, DIMRset {version['dimr_version']} and D-HYDRO suite 1D2D {version['suite_version']}"

        return ncfile

    def init_1dnetwork(self, ncfile, mesh1d):

        # dimensions of the network
        ncfile.createDimension("time", None)
        ncfile.createDimension("network1d_nEdges", mesh1d.meshgeomdim.nbranches)
        ncfile.createDimension("network1d_nNodes", mesh1d.meshgeomdim.nnodes)
        ncfile.createDimension("network1d_nGeometryNodes", mesh1d.meshgeomdim.ngeometry)        
        ncfile.createDimension("idstrlength", self.idstrlength)
        ncfile.createDimension("longstrlength", self.longstrlength)
        ncfile.createDimension("mesh1d_nEdges", mesh1d.meshgeomdim.numedge)
        ncfile.createDimension("mesh1d_nNodes", mesh1d.meshgeomdim.numnode)
        ncfile.createDimension("Two", 2)

    def init_2dmesh(self, ncfile, cmesh2d):

        # Create dimensions
        ncfile.createDimension("max_nmesh2d_face_nodes", cmesh2d.meshgeomdim.maxnumfacenodes)
        ncfile.createDimension("mesh2d_nEdges", cmesh2d.meshgeomdim.numedge)
        ncfile.createDimension("mesh2d_nFaces", cmesh2d.meshgeomdim.numface)
        ncfile.createDimension("mesh2d_nNodes", cmesh2d.meshgeomdim.numnode)

    def init_1d2dlinks(self, ncfile, links1d2d):
        # Dimensions
        ncfile.createDimension("nLink1D2D_edge", len(links1d2d.nodes1d))
        
        cm = ncfile.createVariable("composite_mesh", "i4", ())
        cm.cf_role = 'parent_mesh_topology'
        cm.meshes= 'mesh1D mesh2D'
        cm.mesh_contact = 'link1d2d'

    def set_1dnetwork(self, ncfile, cmesh1d):

        #  network topology
        ntw = ncfile.createVariable("network1d", "i4", ())
        ntw.cf_role = 'mesh_topology'
        ntw.edge_dimension = 'network1d_nEdges'
        ntw.edge_geometry = 'network1d_geometry'
        ntw.edge_node_connectivity = 'network1d_edge_nodes'
        ntw.long_name = "Topology data of 1D network"
        ntw.node_coordinates = 'network1d_node_x network1d_node_y'
        ntw.node_dimension = 'network1d_nNodes'
        ntw.topology_dimension = 1        
        ntw.node_id = "network1d_node_id"
        ntw.node_long_name = "network1d_node_long_name"
        ntw.branch_id = "network1d_branch_id"
        ntw.branch_long_name = "network1d_branch_long_name"
        ntw.edge_length = "network1d_edge_length"
        ntw.branch_order = "network1d_branch_order"

        ntw_node_id = ncfile.createVariable("network1d_node_id", "c", ("network1d_nNodes", "idstrlength"))
        ntw_node_id.long_name = "ID of network nodes"        
        ntw_node_id[:] = self.to_char_list(cmesh1d.description1d["network_node_ids"], self.idstrlength)

        ntw_node_longname = ncfile.createVariable("network1d_node_long_name", "c", ("network1d_nNodes", "longstrlength"))
        ntw_node_longname.long_name = "Long name of network nodes"
        ntw_node_longname[:] = self.to_char_list(cmesh1d.description1d["network_node_long_names"], self.longstrlength)

        # network nodes
        ntw_node_x = ncfile.createVariable("network1d_node_x", np.float64, "network1d_nNodes")
        ntw_node_x.standard_name = 'projection_x_coordinate'
        ntw_node_x.long_name = "x coordinates of network nodes"
        ntw_node_x.units = 'm'
        ntw_node_x[:] = cmesh1d.get_values('nnodex', as_array=True)

        ntw_node_y = ncfile.createVariable("network1d_node_y", np.float64, "network1d_nNodes")
        ntw_node_y.standard_name = 'projection_y_coordinate'
        ntw_node_y.long_name = "y coordinates of network nodes"
        ntw_node_y.units = 'm'
        ntw_node_y[:] = cmesh1d.get_values('nnodey', as_array=True)

        ntw_branch_id_name = ncfile.createVariable("network1d_branch_id", "c", ("network1d_nEdges", "idstrlength"))
        ntw_branch_id_name.long_name = "ID of branch geometries"
        ntw_branch_id_name[:] = self.to_char_list(cmesh1d.description1d["network_branch_ids"], self.idstrlength)

        ntw_branch_id_longname = ncfile.createVariable("network1d_branch_long_name", "c", ("network1d_nEdges", "longstrlength"))
        ntw_branch_id_longname.long_name = "Long name of branch geometries"
        ntw_branch_id_longname[:] = self.to_char_list(cmesh1d.description1d["network_branch_long_names"], self.longstrlength)

        ntw_edge_length = ncfile.createVariable("network1d_edge_length", np.float64, "network1d_nEdges")        
        ntw_edge_length.long_name = "Real length of branch geometries"
        ntw_edge_length.units = "m"
        ntw_edge_length[:] = cmesh1d.get_values('nbranchlengths', as_array=True)

        ntw_branch_order = ncfile.createVariable("network1d_branch_order", "i4", "network1d_nEdges")        
        ntw_branch_order.long_name = "Order of branches for interpolation"
        ntw_branch_order.mesh = "network1d"
        ntw_branch_order.location = "edge"
        ntw_branch_order[:] = cmesh1d.get_values('nbranchorder', as_array=True)

        # network edges
        ntw_edge_node = ncfile.createVariable("network1d_edge_nodes", "i4", ("network1d_nEdges", "Two"))
        ntw_edge_node.cf_role = 'edge_node_connectivity'
        ntw_edge_node.long_name = 'start and end nodes of network edges'
        ntw_edge_node.start_index = 1
        ntw_edge_node[:] = cmesh1d.get_values('nedge_nodes', as_array=True)

        # network geometry
        ntw_geom = ncfile.createVariable("network1d_geometry", "i4", ())
        ntw_geom.geometry_type = 'line'
        ntw_geom.long_name = "1D Geometry"
        ntw_geom.node_count = "network1d_geom_node_count"     
        ntw_geom.node_coordinates = 'network1d_geom_x network1d_geom_y'
        
        ntw_geom_node_count = ncfile.createVariable("network1d_geom_node_count", "i4", "network1d_nEdges")
        ntw_geom_node_count.long_name = "Number of geometry nodes per branch"
        ntw_geom_node_count[:] = cmesh1d.get_values('nbranchgeometrynodes', as_array=True)

        ntw_geom_x = ncfile.createVariable("network1d_geom_x", np.float64, ("network1d_nGeometryNodes"))
        ntw_geom_x.standard_name = 'projection_x_coordinate'
        ntw_geom_x.long_name = 'x-coordinate of branch geometry nodes'
        ntw_geom_x.units = 'm' 
        ntw_geom_x[:] = cmesh1d.get_values('ngeopointx', as_array=True)        
        
        ntw_geom_y = ncfile.createVariable("network1d_geom_y", np.float64, ("network1d_nGeometryNodes"))
        ntw_geom_y.standard_name = 'projection_y_coordinate'
        ntw_geom_y.long_name = 'y-coordinate of branch geometry nodes'
        ntw_geom_y.units = 'm'           
        ntw_geom_y[:] = cmesh1d.get_values('ngeopointy', as_array=True)

    def set_1dmesh(self, ncfile, cmesh1d):

        mesh1d = ncfile.createVariable("mesh1d", "i4", ())
        mesh1d.cf_role = 'mesh_topology'
        mesh1d.long_name = "Topology data of 1D Mesh"
        mesh1d.coordinate_space = 'network1d'
        mesh1d.edge_dimension = 'mesh1d_nEdges'
        mesh1d.edge_node_connectivity = 'mesh1d_edge_nodes'
        mesh1d.edge_coordinates = 'mesh1d_edge_branch mesh1d_edge_offset mesh1d_edge_x mesh1d_edge_y'        
        mesh1d.node_coordinates = 'mesh1d_node_branch mesh1d_node_offset'
        mesh1d.node_dimension = 'mesh1d_nNodes'
        mesh1d.node_id = "mesh1d_node_id"
        mesh1d.node_long_name = "mesh1d_node_long_name"
        mesh1d.topology_dimension = 1

        mesh1d_node_id = ncfile.createVariable("mesh1d_node_id", "c", ("mesh1d_nNodes", "idstrlength"))        
        mesh1d_node_id.long_name = "ID of mesh nodes"        
        mesh1d_node_id[:] = self.to_char_list(cmesh1d.description1d["mesh1d_node_ids"], self.idstrlength)

        mesh1d_node_longname = ncfile.createVariable("mesh1d_node_long_name", "c", ("mesh1d_nNodes", "longstrlength"))
        mesh1d_node_longname.long_name = 'Long name of mesh nodes'                
        mesh1d_node_longname[:] = self.to_char_list(cmesh1d.description1d["mesh1d_node_long_names"], self.longstrlength)

        mesh1d_edge_node = ncfile.createVariable("mesh1d_edge_nodes", "i4", ("mesh1d_nEdges", "Two"))
        mesh1d_edge_node.cf_role  = "edge_node_connectivity"
        mesh1d_edge_node.long_name = 'Start and end nodes of mesh edges'
        mesh1d_edge_node.start_index = 1       
        mesh1d_edge_node[:] = cmesh1d.get_values("edge_nodes", as_array=True) 

        mesh1d_edge_branch = ncfile.createVariable("mesh1d_edge_branch", "i4", "mesh1d_nEdges")
        mesh1d_edge_branch.long_name = "Index of branch on which mesh edges are located"
        mesh1d_edge_branch.start_index = 1
        mesh1d_edge_branch[:] = np.ravel(cmesh1d.edge_branchidx)           
       
        mesh1d_edge_offset = ncfile.createVariable("mesh1d_edge_offset", np.float64, "mesh1d_nEdges")
        mesh1d_edge_offset.long_name = "Offset along branch of mesh edges"
        mesh1d_edge_offset.units = "m"        
        mesh1d_edge_offset[:] = np.ravel(cmesh1d.edge_branchoffset)
        
        mesh1d_edge_x = ncfile.createVariable("mesh1d_edge_x", np.float64, "mesh1d_nEdges")
        mesh1d_edge_x.long_name = "x-coordinate of mesh edges"
        mesh1d_edge_x.units = "m"        
        mesh1d_edge_x[:] = np.ravel(cmesh1d.edge_x)
        
        mesh1d_edge_y = ncfile.createVariable("mesh1d_edge_y", np.float64, "mesh1d_nEdges")
        mesh1d_edge_y.long_name = "y-coordinate of mesh edges"
        mesh1d_edge_y.units = "m"        
        mesh1d_edge_y[:] = np.ravel(cmesh1d.edge_y)
        
        mesh1d_node_branch = ncfile.createVariable("mesh1d_node_branch", "i4", "mesh1d_nNodes")
        mesh1d_node_branch.long_name = "Index of branch on which mesh nodes are located"
        mesh1d_node_branch.start_index = 1
        mesh1d_node_branch[:] = cmesh1d.get_values("branchidx", as_array=True)

        mesh1d_node_offset = ncfile.createVariable("mesh1d_node_offset", np.float64, "mesh1d_nNodes", fill_value= 0.)        
        mesh1d_node_offset.long_name = "Offset along branch of mesh nodes"
        mesh1d_node_offset.units = "m"
        mesh1d_node_offset[:] = cmesh1d.get_values("branchoffsets", as_array=True)
       
    # set 2d mesh data to netcdf file
    def set_2dmesh(self, ncfile, cmesh2d):

        mesh2d = ncfile.createVariable("mesh2d", "i4", ())
        mesh2d.long_name = "Topology data of 2D network"
        mesh2d.topology_dimension = 2
        mesh2d.cf_role = 'mesh_topology'
        mesh2d.node_coordinates = 'mesh2d_node_x mesh2d_node_y'
        mesh2d.node_dimension = 'mesh2d_nNodes'
        mesh2d.edge_coordinates = 'mesh2d_edge_x mesh2d_edge_y'
        mesh2d.edge_dimension = 'mesh2d_nEdges'
        mesh2d.edge_node_connectivity = 'mesh2d_edge_nodes'
        mesh2d.face_node_connectivity = 'mesh2d_face_nodes'
        mesh2d.max_face_nodes_dimension = 'max_nmesh2d_face_nodes'
        mesh2d.face_dimension = "mesh2d_nFaces"
        #mesh2d.edge_face_connectivity = "mesh2d_edge_faces"
        mesh2d.face_coordinates = "mesh2d_face_x mesh2d_face_y"

        # Nodes:
        mesh2d_x = ncfile.createVariable("mesh2d_node_x", np.float64, mesh2d.node_dimension)
        mesh2d_y = ncfile.createVariable("mesh2d_node_y", np.float64, mesh2d.node_dimension)
        mesh2d_z = ncfile.createVariable("mesh2d_node_z", np.float64, mesh2d.node_dimension, fill_value=-999.0)

        mesh2d_x.standard_name = 'projection_x_coordinate'
        mesh2d_y.standard_name = 'projection_y_coordinate'
        mesh2d_z.standard_name = 'altitude'

        for var, dim in zip([mesh2d_x, mesh2d_y, mesh2d_z],  list('xyz')):
            setattr(var, 'units', 'm')
            setattr(var, 'mesh', 'mesh2d')
            setattr(var, 'location', 'node')
            setattr(var, 'long_name', f'{dim}-coordinate of mesh nodes')

        mesh2d_z.coordinates = 'mesh2d_node_x mesh2d_node_y'
        mesh2d_z.grid_mapping = ''
        
        mesh2d_x[:] = cmesh2d.get_values("nodex")
        mesh2d_y[:] = cmesh2d.get_values("nodey")
        
        # Edges:
        # mesh2d_xu = ncfile.createVariable("mesh2d_edge_x", np.float64,  "nmesh2d_edges")
        # mesh2d_yu = ncfile.createVariable("mesh2d_edge_y", np.float64,  "nmesh2d_edges")
        # mesh2d_xu[:] = cmesh2d.get_values("edgex")
        # mesh2d_yu[:] = cmesh2d.get_values("edgey")
        # mesh2d_xu.long_name = 'x-coordinate of mesh edges'
        # mesh2d_yu.long_name = 'y-coordinate of mesh edges'

        # for var, dim in zip([mesh2d_xu, mesh2d_yu],  list('xy')):
        #     setattr(var, 'units', 'm')
        #     setattr(var, 'mesh', 'mesh2d')
        #     setattr(var, 'location', 'edge')
        #     setattr(var, 'standard_name', f'projection_{dim}_coordinate')

        mesh2d_en = ncfile.createVariable("mesh2d_edge_nodes", "i4", ( mesh2d.edge_dimension, "Two"), fill_value=-999)
        mesh2d_en.cf_role = 'edge_node_connectivity'
        mesh2d_en.long_name = 'maps every edge to the two nodes that it connects'
        mesh2d_en.start_index = 1
        mesh2d_en.location = 'edge'
        mesh2d_en.mesh = 'mesh2d'
        mesh2d_en[:] = cmesh2d.get_values('edge_nodes', as_array=True)

        # mesh2d_et = ncfile.createVariable("mesh2d_edge_types", "i4", mesh2d.edge_dimension)
        # mesh2d_et.long_name = 'edge type (relation between edge and flow geometry)'
        # mesh2d_et.coordinates = 'mesh2d_edge_x mesh2d_edge_y'
        # mesh2d_et.location = 'edge'
        # mesh2d_et.mesh = 'mesh2d'
        # mesh2d_et.standard_name = ''
        # mesh2d_et.units = ''
        # mesh2d_et[:] = 2

        mesh2d_fn = ncfile.createVariable("mesh2d_face_nodes", "i4", (mesh2d.face_dimension, mesh2d.max_face_nodes_dimension), fill_value=-999)
        mesh2d_fn.cf_role = 'face_node_connectivity'
        mesh2d_fn.mesh = 'mesh2d'
        mesh2d_fn.location = 'face'
        mesh2d_fn.long_name = 'maps every face to the nodes that it defines'
        mesh2d_fn.start_index = 1
        mesh2d_fn[:] = cmesh2d.get_values('face_nodes', as_array=True)

        mesh2d_face_x = ncfile.createVariable("mesh2d_face_x", np.float64, mesh2d.face_dimension)
        mesh2d_face_y = ncfile.createVariable("mesh2d_face_y", np.float64, mesh2d.face_dimension)
        mesh2d_face_z = ncfile.createVariable("mesh2d_face_z", np.float64, mesh2d.face_dimension, fill_value=-999.0)

        for var, dim in zip([mesh2d_face_x, mesh2d_face_y, mesh2d_face_z],  list('xyz')):
            setattr(var, 'units', 'm')
            setattr(var, 'mesh', 'mesh2d')
            setattr(var, 'location', 'face')
            setattr(var, 'standard_name', f'projection_{dim}_coordinate' if dim != 'z' else 'altitude')
            setattr(var, 'long_name', f'{dim}-coordinate of face nodes')

        mesh2d_face_z.coordinates = 'mesh2d_face_x mesh2d_face_y'
        mesh2d_face_z.grid_mapping = ''
    
        mesh2d_face_x[:] = cmesh2d.get_values("facex")
        mesh2d_face_y[:] = cmesh2d.get_values("facey")

        # Assign altitude data
        # To faces
        if cmesh2d.is_allocated('facez'):
            mesh2d_face_z[:] = cmesh2d.get_values("facez")
        # Assign to nodes
        if cmesh2d.is_allocated('nodez'):
            mesh2d_z[:] = cmesh2d.get_values("nodez")
        # Raise error if none of both is allocated
        if not cmesh2d.is_allocated('nodez') and not cmesh2d.is_allocated('facez'):
            raise ValueError('Assign altitude values either to nodes or faces.')


    def set_1d2dlinks(self, ncfile, links1d2d):
        
        nlinks = len(links1d2d.nodes1d)

        link1d2d = ncfile.createVariable("link1d2d", "i4", ("nLink1D2D_edge", "Two"))
        link1d2d.cf_role = 'mesh_topology_contact'
        link1d2d.contact= 'mesh1D:node mesh2D:face'
        link1d2d.contact_type = 'link1d2d_contact_type'
        link1d2d.contact_ids = 'link1d2d_ids'
        link1d2d.contact_long_names = 'link1d2d_long_names'
        link1d2d.start_index = 1
        link1d2d[:, :] = np.c_[links1d2d.nodes1d, links1d2d.faces2d]

        link1d2d_ids = ncfile.createVariable("link1d2d_ids", "c", ("nLink1D2D_edge", "idstrlength"))
        link1d2d_ids.long_name = 'ids of the contact'
        link1d2d_ids[:] = [self.str2chars(f'linkid{i+1}', self.idstrlength) for i in range(nlinks)]

        link1d2d_long_names = ncfile.createVariable("link1d2d_long_names", "c", ("nLink1D2D_edge", "longstrlength"))
        link1d2d_long_names.long_name = 'long names of the contact'
        link1d2d_long_names[:] = [self.str2chars(f'longnames{i+1}', self.longstrlength) for i in range(nlinks)]

        link1d2d_contact_type = ncfile.createVariable("link1d2d_contact_type", "i4", "nLink1D2D_edge", fill_value=-1)
        link1d2d_contact_type[:] = [3] * nlinks

    def str2chars(self, string, size):
        return string.ljust(size)[:size]
    
    # def add_bounds(self, ncfile):
