from ctypes import POINTER, Structure, c_char, c_double, c_int
import os

import numpy as np
from shapely.geometry import LineString


# Definieer load structure
class meshgeomdim(Structure):

    _fields_ = [
        ('meshname', POINTER(c_char)),
        ('dim', c_int),
        ('numnode', c_int),
        ('numedge', c_int),
        ('numface', c_int),
        ('maxnumfacenodes', c_int),
        ('numlayer', c_int),
        ('layertype', c_int),
        ('nnodes', c_int),
        ('nbranches', c_int),
        ('ngeometry', c_int),
        ('epgs', c_int)
    ]

    def __repr__(self):
        strings = []
        for field, _ in self._fields_:
            value = getattr(self, field)
            strings.append(f'{field:15s}: {value}')

        return '\n'.join(strings)

class meshgeom(Structure):
    _fields_ = [
        ('edge_nodes', POINTER(c_int)),
        ('face_nodes', POINTER(c_int)),
        ('edge_faces', POINTER(c_int)),
        ('face_edges', POINTER(c_int)),
        ('face_links', POINTER(c_int)),
        ('nnodex', POINTER(c_double)),
        ('nnodey', POINTER(c_double)),
        ('nedge_nodes', POINTER(c_int)),
        ('nbranchlengths', POINTER(c_double)),
        ('nbranchgeometrynodes', POINTER(c_int)),
        ('ngeopointx', POINTER(c_double)),
        ('ngeopointy', POINTER(c_double)),
        ('nbranchorder', POINTER(c_double)),
        ('branchidx', POINTER(c_int)),
        ('branchoffsets', POINTER(c_double)),
        #('edge_branchidx', POINTER(c_int)),
        #('edge_branchoffsets', POINTER(c_double)),        
        ('nodex', POINTER(c_double)),
        ('nodey', POINTER(c_double)),
        ('nodez', POINTER(c_double)),
        ('edgex', POINTER(c_double)),
        ('edgey', POINTER(c_double)),
        ('edgez', POINTER(c_double)),
        ('facex', POINTER(c_double)),
        ('facey', POINTER(c_double)),
        ('facez', POINTER(c_double)),
        ('layer_zs', POINTER(c_int)),
        ('interface_zs', POINTER(c_int))
    ]

    def __init__(self, geometries):
        """
        Constructor
        """
        
        self.meshgeomdim = geometries

        # This dictionary contains some extra variables for 1d meshes
        self.description1d = {
            'mesh1d_node_ids': [],
            'mesh1d_node_long_names': [],
            'network_node_ids': [],
            'network_node_long_names': [],
            'network_branch_ids': [],
            'network_branch_long_names': []
        }

        self._meta_ = {
            'nodex': {'ctype': c_double, 'size': ['numnode'], 'allocated': False},
            'nodey': {'ctype': c_double, 'size': ['numnode'], 'allocated': False},
            'nodez': {'ctype': c_double, 'size': ['numnode'], 'allocated': False},
            'branchoffsets': {'ctype': c_double, 'size': ['numnode'], 'allocated': False},
            'branchidx': {'ctype': c_int, 'size': ['numnode'], 'allocated': False},            
            #'edge_branchoffsets': {'ctype': c_double, 'size': ['numedge'], 'allocated': False},
            #'edge_branchidx': {'ctype': c_int, 'size': ['numedge'], 'allocated': False},
            'nnodex': {'ctype': c_double, 'size': ['nnodes'], 'allocated': False},
            'nnodey': {'ctype': c_double, 'size': ['nnodes'], 'allocated': False},
            'ngeopointx': {'ctype': c_double, 'size': ['ngeometry'], 'allocated': False},
            'ngeopointy': {'ctype': c_double, 'size': ['ngeometry'], 'allocated': False},
            'nbranchlengths': {'ctype': c_double, 'size': ['nbranches'], 'allocated': False},
            'nbranchorder': {'ctype': c_double, 'size': ['nbranches'], 'allocated': False},
            'edgex': {'ctype': c_double, 'size': ['numedge'], 'allocated': False},
            'edgey': {'ctype': c_double, 'size': ['numedge'], 'allocated': False},
            'edgez': {'ctype': c_double, 'size': ['numedge'], 'allocated': False},
            'nedge_nodes': {'ctype': c_int, 'size': ['nbranches', 2], 'allocated': False},
            'nbranchgeometrynodes': {'ctype': c_int, 'size': ['nbranches'], 'allocated': False},
            'facex': {'ctype': c_double, 'size': ['numface'], 'allocated': False},
            'facey': {'ctype': c_double, 'size': ['numface'], 'allocated': False},
            'facez': {'ctype': c_double, 'size': ['numface'], 'allocated': False},
            'face_nodes': {'ctype': c_int, 'size': ['numface', 'maxnumfacenodes'], 'allocated': False},
            'edge_nodes': {'ctype': c_int, 'size': ['numedge', 2], 'allocated': False},
            'edge_faces': {'ctype': c_int, 'size': ['numedge', 2], 'allocated': False},
            'face_edges': {'ctype': c_int, 'size': ['numface', 'maxnumfacenodes'], 'allocated': False}
        }

    def get_dimensions(self, var):
        return tuple(getattr(self.meshgeomdim, fac) if isinstance(fac, str) else fac for fac in self._meta_[var]['size'])

    def get_size(self, var):
        return np.prod(self.get_dimensions(var))

    def allocate(self, var):
        ctype = self._meta_[var]['ctype']
        size = self.get_size(var)
        setattr(self, var, (ctype * size)())
        self._meta_[var]['allocated'] = True

    def is_allocated(self, var):
        return self._meta_[var]['allocated']

    def get_values(self, var, as_array=False, size=None):
        if not self.is_allocated(var):
            return None
        size = self.get_size(var) if size is None else size
        values = [getattr(self, var)[i] for i in range(size)]
        if as_array:
            values = np.reshape(values, self.get_dimensions(var))
        return values

    def set_values(self, var, values):
        # First allocate
        self.allocate(var)

        size = self.get_size(var)
        if size != len(values):
            raise ValueError(f'Size of values ({len(values)}) does not match allocated size ({size}) for "{var}".')
        for i, value in enumerate(values):
            getattr(self, var)[i] = value


    def add_values(self, var, values):
        
        # If not yet allocated, don't add but set
        if not self.is_allocated(var):
            self.set_values(var, values)
        
        # Get old values
        old_values = self.get_values(var, size=(self.get_size(var)-len(values)))
        
        # First allocate
        self.allocate(var)
        size = self.get_size(var)
        
        if size != (len(values) + len(old_values)):
            raise ValueError(f'Size of values ({len(values) + len(old_values)}) does not match allocated size ({size})')

        for i, value in enumerate(old_values + values):
            getattr(self, var)[i] = value

    def add_from_other(self, geometries):
        """
        Method to merge mesh with another mesh
        """

        # If the current mesh is empty, copy the maxfacenumnodes
        if self.empty():
            self.meshgeomdim.maxnumfacenodes = geometries.meshgeomdim.maxnumfacenodes
        # If the maxnumfacenodes is not equal, raise an error. Merging two diffently shaped meshes not implemented
        if self.meshgeomdim.maxnumfacenodes != geometries.meshgeomdim.maxnumfacenodes:
            raise NotImplementedError('The maximum number of face nodes differs between the meshes.')

        # Get the counter offset for the nodes
        startnode = self.meshgeomdim.numnode
        startbranch = self.meshgeomdim.nbranches
        
        # Add dimensions. Sum the new dimensions with the old ones
        for dimension in ['numnode', 'numedge', 'numface', 'nnodes', 'ngeometry', 'nbranches']:
            new = getattr(self.meshgeomdim, dimension) + getattr(geometries.meshgeomdim, dimension)
            setattr(self.meshgeomdim, dimension, new)
        
        # Add variables. Add all new data
        for var in ['nodex', 'nodey', 'nodez', 'facex', 'facey', 'facez', 'nnodex', 'nnodey',
                    'nbranchlengths', 'nbranchorder', 'ngeopointx', 'ngeopointy', 'nbranchgeometrynodes']:
            if geometries.is_allocated(var):
                self.add_values(var, geometries.get_values(var))

        # For variables with indexes. For the indexes, add the old start node
        for indexvar in ['edge_nodes', 'face_nodes', 'branchidx', 'branchoffsets']:
            if geometries.is_allocated(indexvar):
                self.add_values(indexvar, [i + startnode for i in geometries.get_values(indexvar)]) 

        # Network indexvar
        indexvar = 'nedge_nodes'
        if geometries.is_allocated(indexvar):
            self.add_values(indexvar, [i + startbranch for i in geometries.get_values(indexvar)]) 

        # Network descriptive variables
        for var in ['network_node_ids', 'network_node_long_names', 'network_branch_ids',
                    'network_branch_long_names', 'mesh1d_node_ids', 'mesh1d_node_long_names']:
            self.description1d[var] += geometries.description1d[var]
        
        
    def process_1d_network(self):
        """
        Determine x, y locations of 1d network
        """
        assert self.meshgeomdim.dim == 1

        # Read network geometry nodes. Create an inverted list to pop
        ngeom = list(zip(self.get_values('ngeopointx'), self.get_values('ngeopointy')))

        # Get branch data
        nbranchnames = self.description1d['network_branch_ids']
        ngeometrynodes = self.get_values('nbranchgeometrynodes')
        branchids = self.get_values('branchidx', as_array=True)
        offsets = self.get_values('branchoffsets', as_array=True)
        
        # Collect branches
        branches = {}
        schematised = {}
        crds = []
        for i, (name, nnodes) in enumerate(zip(nbranchnames, ngeometrynodes)):
            # Create linestring for network branch
            linestring = LineString([ngeom.pop(0) for _ in range(nnodes)])
            branches[name.strip()] = linestring
            # Determine mesh node position on network branch
            branchoffsets = offsets[branchids == (i + 1)]
            meshcrds = [linestring.interpolate(offset).coords[0] for offset in branchoffsets]
            crds.extend(meshcrds[:])
            # Determine if a start or end coordinate needs to be added for constructing a complete LineString
            if not np.isclose(branchoffsets[0], 0.0):
                meshcrds = [linestring.coords[0]] + meshcrds
            if not np.isclose(branchoffsets[-1], linestring.length):
                meshcrds = meshcrds + [linestring.coords[-1]]
            schematised[name.strip()] = LineString(meshcrds)

        # Add values to mesh
        nodex, nodey = list(zip(*crds))
        self.set_values('nodex', nodex)
        self.set_values('nodey', nodey)

        return schematised, branches

    def empty(self):
        """Determine whether mesh is empty, based on number of nodes"""
        return not bool(self.meshgeomdim.numnode)

    def get_nodes(self):
        """
        Return nodes.
        
        Returns
        -------
        np.ndarray
            Numpy array with x and y coordinates of nodes.
        """

        # Get nodes
        nodes = np.c_[
            np.array(self.get_values('nodex')),
            np.array(self.get_values('nodey'))
        ]

        return nodes

    def get_nodes_for_branch(self, branchid):
        """
        Return index of nodes on branch(es).

        Note that each nodes belongs to a single branch in the network.
        This can provide unexpected results on intersections.
        
        Parameters
        ----------
        branchid : str or list
            Branchid for which to return the nodes.
        
        Returns
        -------
        np.ndarray
            boolean array with True for nodes that are on the branch(es)
        """        
        if branchid is None:
            return np.ones(len(self.get_values('nodex')), dtype=bool)

        # Convert to list if needed
        if isinstance(branchid, str):
            branchid = [branchid]
        # Get the ids (integers) of the branch names given by the user
        branchidx = np.where(np.isin(self.description1d['network_branch_ids'], branchid))[0] + 1
        # Select which of the nodes are in the branches
        idx = np.isin(self.get_values('branchidx'), branchidx)
        
        return idx

    def get_segments(self):

        # Read nodes and links from src
        nodes = self.get_nodes()
        # Get links
        links = self.get_values('edge_nodes', as_array=True)

        return nodes[links - 1]


    def get_faces(self, geometry='exterior'):
        """
        Get cells from 2d mesh

        There are three options for returning coordinates:
        - "center": can be the circumcenter (default) or the centroid if the face centers are moved by the user
        - "centroid": the centroid, which is calculated (again) in this function
        - "exterior" (default): the coordinates of the face edges
        """

        assert self.meshgeomdim.dim == 2
        
        # Read nodes and links from src
        nodes = self.get_nodes()

        if geometry == 'center':
            # Can be both the centroid and the circumcenter
            return np.c_[self.get_values('facex'), self.get_values('facey')]

        elif geometry in ['exterior', 'centroid']:

            # Get face nodes
            face_nodes = self.get_values('face_nodes', as_array=True)
            nanvalues = (face_nodes != -999).sum(axis=1).astype(int)
            unique = np.unique(nanvalues)
            
            # If all cells have the same number of nodes
            if len(unique) == 1:
                faces = nodes[face_nodes - 1]
                

            # Else, combine a list for all nodes
            else:
                faces = [None] * len(face_nodes)
                for n in unique:
                    ncells = nodes[face_nodes[(nanvalues == n), :n] - 1]
                    where = np.where(nanvalues == n)[0]
                    for i, cell in zip(where, ncells):
                        faces[i] = cell

            if geometry == 'exterior':
                return faces
            
            # Determine the centroids by averaging the exteriors
            if geometry == 'centroid':
                return np.vstack([f.mean(axis=0) for f in faces])
        
        else:
            raise ValueError(f'Geometry "{geometry}" not recognized. Pick "center", "centroid" or "exterior"')
            

        
