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

    def write(self, dflowfmmodel, path):  # write ugrid file from GWSW model

        ncfile = self.create_netcdf(path)

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

    def create_netcdf(self, path):

        # File format:
        outformat = "NETCDF3_CLASSIC" #"NETCDF4"
        # File where we going to write
        ncfile = Dataset(path, 'w', format=outformat)

        # global attributes
        ncfile.Conventions = "CF-1.8 UGRID-1.0"
        ncfile.history = "Created on {} D-Flow 1D, D-Flow FM".format(datetime.now())
        ncfile.institution = "Deltares"
        ncfile.reference = "http://www.deltares.nl"
        ncfile.source = "Python script to prepare HyDAMO import"

        return ncfile

    def init_1dnetwork(self, ncfile, mesh1d):

        # dimensions of the network
        ncfile.createDimension("time", None)
        ncfile.createDimension("nnetwork_branches", mesh1d.meshgeomdim.nbranches)
        ncfile.createDimension("nnetwork_nodes", mesh1d.meshgeomdim.nnodes)
        ncfile.createDimension("nnetwork_geometry", mesh1d.meshgeomdim.ngeometry)
        ncfile.createDimension("nnetwork_edges", mesh1d.meshgeomdim.nbranches)
        ncfile.createDimension("idstrlength", self.idstrlength)
        ncfile.createDimension("longstrlength", self.longstrlength)
        ncfile.createDimension("nmesh1d_edges", mesh1d.meshgeomdim.numedge)
        ncfile.createDimension("nmesh1d_nodes", mesh1d.meshgeomdim.numnode)
        ncfile.createDimension("Two", 2)

    def init_2dmesh(self, ncfile, cmesh2d):

        # Create dimensions
        ncfile.createDimension("max_nmesh2d_face_nodes", cmesh2d.meshgeomdim.maxnumfacenodes)
        ncfile.createDimension("nmesh2d_edges", cmesh2d.meshgeomdim.numedge)
        ncfile.createDimension("nmesh2d_faces", cmesh2d.meshgeomdim.numface)
        ncfile.createDimension("nmesh2d_nodes", cmesh2d.meshgeomdim.numnode)

    def init_1d2dlinks(self, ncfile, links1d2d):
        # Dimensions
        ncfile.createDimension("nlinks_1d2d", len(links1d2d.nodes1d))

        cm = ncfile.createVariable("composite_mesh", "i4", ())
        cm.cf_role = 'mesh_topology_parent'
        cm.meshes= 'mesh1D mesh2D'
        cm.mesh_contact = 'link1d2d'

    def set_1dnetwork(self, ncfile, cmesh1d):

        # geometry
        ntw = ncfile.createVariable("network", "i4", ())
        ntw.cf_role = 'mesh_topology'
        ntw.edge_dimension = 'nnetwork_branches'
        ntw.edge_geometry = 'network_geometry'
        ntw.edge_node_connectivity = 'network_edge_nodes'
        ntw.long_name = "Network topology"
        ntw.node_coordinates = 'network_node_x network_node_y'
        ntw.node_dimension = 'nnetwork_nodes'
        ntw.topology_dimension = 1
        ntw.node_ids = "network_node_ids"
        ntw.node_long_names = "network_node_long_names"
        ntw.branch_ids = "network_branch_ids"
        ntw.branch_long_names = "network_branch_long_names"
        ntw.branch_lengths = "network_branch_lengths"
        ntw.branch_order = "network_branch_order"

        ntw_node_id = ncfile.createVariable("network_node_ids", "c", ("nnetwork_nodes", "idstrlength"))
        ntw_node_id.standard_name = 'network_node_ids'
        ntw_node_id.long_name = "The identification name of the node"
        ntw_node_id.mesh = 'network1D'
        ntw_node_id[:] = self.to_char_list(cmesh1d.description1d["network_node_ids"], self.idstrlength)

        ntw_node_longname = ncfile.createVariable("network_node_long_names", "c", ("nnetwork_nodes", "longstrlength"))
        ntw_node_longname.standard_name = 'network_node_longname'
        ntw_node_longname.long_name = "The long name of the node"
        ntw_node_longname.mesh = 'network1D'
        ntw_node_longname[:] = self.to_char_list(cmesh1d.description1d["network_node_long_names"], self.longstrlength)

        ntw_node_x = ncfile.createVariable("network_node_x", np.float64, "nnetwork_nodes")
        ntw_node_x.standard_name = 'projection_x_coordinate'
        ntw_node_x.long_name = "x coordinates of the network connection nodes"
        ntw_node_x.units = 'm'
        ntw_node_x[:] = cmesh1d.get_values('nnodex', as_array=True)

        ntw_node_y = ncfile.createVariable("network_node_y", np.float64, "nnetwork_nodes")
        ntw_node_y.standard_name = 'projection_y_coordinate'
        ntw_node_y.long_name = "y coordinates of the network connection nodes"
        ntw_node_y.units = 'm'
        ntw_node_y[:] = ntw_node_x[:] = cmesh1d.get_values('nnodey', as_array=True)

        ntw_branch_id_name = ncfile.createVariable("network_branch_ids", "c", ("nnetwork_branches", "idstrlength"))
        ntw_branch_id_name.standard_name = 'network_branch_id_name'
        ntw_branch_id_name.long_name = "The identification name of the branch"
        ntw_branch_id_name[:] = self.to_char_list(cmesh1d.description1d["network_branch_ids"], self.idstrlength)

        ntw_branch_id_longname = ncfile.createVariable("network_branch_long_names", "c", ("nnetwork_branches", "longstrlength"))
        ntw_branch_id_longname.standard_name = 'network_branch_longname'
        ntw_branch_id_longname.long_name = "The long name of the branch"
        ntw_branch_id_longname[:] = self.to_char_list(cmesh1d.description1d["network_branch_long_names"], self.longstrlength)

        ntw_branch_length = ncfile.createVariable("network_branch_lengths", np.float64, "nnetwork_branches")
        ntw_branch_length.standard_name = 'network_branch_length'
        ntw_branch_length.long_name = "The calculation length of the branch"
        ntw_branch_length[:] = cmesh1d.get_values('nbranchlengths', as_array=True)

        ntw_branch_order = ncfile.createVariable("network_branch_order", "i4", "nnetwork_branches")
        ntw_branch_order.standard_name = 'network branch order'
        ntw_branch_order.long_name = "The order of the branches for interpolation"
        ntw_branch_order[:] = cmesh1d.get_values('nbranchorder', as_array=True)

        ntw_edge_node = ncfile.createVariable("network_edge_nodes", "i4", ("nnetwork_edges", "Two"))
        ntw_edge_node.cf_role = 'edge_node_connectivity'
        ntw_edge_node.long_name = 'start and end nodes of each branch in the network'
        ntw_edge_node.start_index = 1
        ntw_edge_node[:] = cmesh1d.get_values('nedge_nodes', as_array=True)

        ntw_geom = ncfile.createVariable("network_geometry", "i4", ())
        ntw_geom.geometry_type = 'multiline'
        ntw_geom.long_name = "1D Geometry"
        ntw_geom.node_count = "nnetwork_geometry"
        ntw_geom.part_node_count = 'network_part_node_count'
        ntw_geom.node_coordinates = 'network_geom_x network_geom_y'

        ntw_geom_x = ncfile.createVariable("network_geom_x", np.float64, ("nnetwork_geometry"))
        ntw_geom_x.standard_name = 'projection_x_coordinate'
        ntw_geom_x.units = 'm'
        ntw_geom_x.cf_role = "geometry_x_node"
        ntw_geom_x.long_name = 'x coordinates of the branch geometries'

        ntw_geom_y = ncfile.createVariable("network_geom_y", np.float64, ("nnetwork_geometry"))
        ntw_geom_y.standard_name = 'projection_y_coordinate'
        ntw_geom_y.units = 'm'
        ntw_geom_y.cf_role = "geometry_y_node"
        ntw_geom_y.long_name = 'y coordinates of the branch geometries'

        ntw_geom_x[:] = cmesh1d.get_values('ngeopointx', as_array=True)
        ntw_geom_y[:] = cmesh1d.get_values('ngeopointy', as_array=True)

    def set_1dmesh(self, ncfile, cmesh1d):

        mesh1d = ncfile.createVariable("mesh1d", "i4", ())
        mesh1d.cf_role = 'mesh_topology'
        mesh1d.coordinate_space = 'network'
        mesh1d.edge_dimension = 'nmesh1d_edges'
        mesh1d.edge_node_connectivity = 'mesh1d_edge_nodes'
        mesh1d.long_name = "1D Mesh"
        mesh1d.node_coordinates = 'mesh1d_nodes_branch_id mesh1d_nodes_branch_offset'
        mesh1d.node_dimension = 'nmesh1d_nodes'
        mesh1d.node_ids = "mesh1d_node_ids"
        mesh1d.node_long_names = "mesh1d_node_long_names"
        mesh1d.topology_dimension = 1

        mesh1d_node_count = ncfile.createVariable("network_part_node_count", "i4", "nnetwork_branches")
        mesh1d_node_count.standard_name = 'network part node count'
        mesh1d_node_count.long_name = "The number of nodes in per branch"
        mesh1d_node_count[:] = cmesh1d.get_values('nbranchgeometrynodes', as_array=True)

        mesh1d_node_id = ncfile.createVariable("mesh1d_node_ids", "c", ("nmesh1d_nodes", "idstrlength"))
        mesh1d_node_id.standard_name = 'mesh1d_node_ids'
        mesh1d_node_id.long_name = "The name of the calculation points"
        mesh1d_node_id.mesh = 'mesh1d'
        mesh1d_node_id[:] = self.to_char_list(cmesh1d.description1d["mesh1d_node_ids"], self.idstrlength)

        mesh1d_node_longname = ncfile.createVariable("mesh1d_node_long_names", "c", ("nmesh1d_nodes", "longstrlength"))
        mesh1d_node_longname.standard_name = 'mesh1d_node_longname'
        mesh1d_node_longname.long_name = "The long name of calculation points"
        mesh1d_node_longname.mesh = 'mesh1d'
        mesh1d_node_longname[:] = self.to_char_list(cmesh1d.description1d["mesh1d_node_long_names"], self.longstrlength)

        mesh1d_edge_node = ncfile.createVariable("mesh1d_edge_nodes", "i4", ("nmesh1d_edges", "Two"))
        mesh1d_edge_node.cf_role = 'edge_node_connectivity'
        mesh1d_edge_node.long_name = 'start and end nodes of each branch in the 1d mesh'
        mesh1d_edge_node.start_index = 1
        mesh1d_edge_node[:] = cmesh1d.get_values("edge_nodes", as_array=True)

        mesh1d_point_branch_id = ncfile.createVariable("mesh1d_nodes_branch_id", "i4", "nmesh1d_nodes")
        mesh1d_point_branch_id.standard_name = 'network calculation point branch id'
        mesh1d_point_branch_id.long_name = "The identification the branch of the calculation point"
        mesh1d_point_branch_id.start_index = 1
        mesh1d_point_branch_id[:] = cmesh1d.get_values("branchidx", as_array=True)

        mesh1d_point_branch_offset = ncfile.createVariable("mesh1d_nodes_branch_offset", np.float64, "nmesh1d_nodes")
        mesh1d_point_branch_offset.standard_name = 'network calculation point branch offset'
        mesh1d_point_branch_offset.long_name = "The offset of the calculation point on the branch"
        mesh1d_point_branch_offset.start_index = 1
        mesh1d_point_branch_offset[:] = cmesh1d.get_values("branchoffsets", as_array=True)

    # set 2d mesh data to netcdf file
    def set_2dmesh(self, ncfile, cmesh2d):

        mesh2d = ncfile.createVariable("mesh2d", "i4", ())
        mesh2d.long_name = "Topology data of 2D network"
        mesh2d.topology_dimension = 2
        mesh2d.cf_role = 'mesh_topology'
        mesh2d.node_coordinates = 'mesh2d_node_x mesh2d_node_y'
        mesh2d.node_dimension = 'nmesh2d_nodes'
        mesh2d.edge_coordinates = 'mesh2d_edge_x mesh2d_edge_y'
        mesh2d.edge_dimension = 'nmesh2d_edges'
        mesh2d.edge_node_connectivity = 'mesh2d_edge_nodes'
        mesh2d.face_node_connectivity = 'mesh2d_face_nodes'
        mesh2d.max_face_nodes_dimension = 'max_nmesh2d_face_nodes'
        mesh2d.face_dimension = "nmesh2d_faces"
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

        link1d2d = ncfile.createVariable("link1d2d", "i4", ("nlinks_1d2d", "Two"))
        link1d2d.cf_role = 'mesh_topology_contact'
        link1d2d.contact= 'mesh1D:node mesh2D:face'
        link1d2d.contact_type = 'link1d2d_contact_type'
        link1d2d.contact_ids = 'link1d2d_ids'
        link1d2d.contact_long_names = 'link1d2d_long_names'
        link1d2d.start_index = 1
        link1d2d[:, :] = np.c_[links1d2d.nodes1d, links1d2d.faces2d]

        link1d2d_ids = ncfile.createVariable("link1d2d_ids", "c", ("nlinks_1d2d", "idstrlength"))
        link1d2d_ids.long_name = 'ids of the contact'
        link1d2d_ids[:] = [self.str2chars(f'linkid{i+1}', self.idstrlength) for i in range(nlinks)]

        link1d2d_long_names = ncfile.createVariable("link1d2d_long_names", "c", ("nlinks_1d2d", "longstrlength"))
        link1d2d_long_names.long_name = 'long names of the contact'
        link1d2d_long_names[:] = [self.str2chars(f'longnames{i+1}', self.longstrlength) for i in range(nlinks)]

        link1d2d_contact_type = ncfile.createVariable("link1d2d_contact_type", "i4", "nlinks_1d2d", fill_value=-1)
        link1d2d_contact_type[:] = [3] * nlinks

    def str2chars(self, string, size):
        return string.ljust(size)[:size]
    
    # def add_bounds(self, ncfile):
