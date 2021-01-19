import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString

from delft3dfmpy.core import checks, geometry
from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame

def generate_culverts(culverts, id_col='id', roughness_values=None, logger=logging):
    """
    Generate culverts from OSM to D-Flow FM
    """

    culverts_dfm = culverts.copy()
    culverts_dfm['crosssection'] = [{} for _ in range(len(culverts_dfm))]

    # Set culvert id
    culverts_dfm[id_col] = 'C_' + culverts_dfm[id_col]

    # Rename to id
    culverts_dfm.set_index(id_col, inplace=True, drop=False)

    # Set cross-section culvert
    for culvert in culverts_dfm.itertuples():
        # Generate cross section definition name
        if culvert.profile == 'round':
            crosssection = {'shape': 'circle', 'diameter': culvert.diameter}

        elif culvert.profile == 'boxed_rectangular':
            crosssection = {'shape': 'rectangle', 'height': culvert.depth, 'width': culvert.width,
                            'closed': 1}
        else:
            logger.debug(f'Culvert {culvert.Index} on branch {culvert.branch_id} has an unknown shape: {culvert.profile}. Applying a default profile (round - 50cm)')
            crosssection = {'shape': 'circle', 'diameter': 0.50}

        # Set culvert parameters
        culverts_dfm.at[culvert.Index, 'allowedflowdir'] = 'both'
        culverts_dfm.at[culvert.Index, 'valveonoff'] = int(0)
        culverts_dfm.at[culvert.Index, 'numlosscoeff'] = int(0)
        culverts_dfm.at[culvert.Index, 'valveopeningheight'] = 0
        culverts_dfm.at[culvert.Index, 'relopening'] = 0
        culverts_dfm.at[culvert.Index, 'losscoeff'] = 0
        culverts_dfm.at[culvert.Index, 'crosssection'] = crosssection

        # Set roughness values of culvert
        if culvert.material in roughness_values.keys():
            culverts_dfm.at[culvert.Index,'friction_value'] = roughness_values[culvert.material]
        elif culvert.material=='concret':
            culverts_dfm.at[culvert.Index, 'friction_value'] = roughness_values['concrete']
        else:
            culverts_dfm.at[culvert.Index, 'friction_value'] = roughness_values['default']
            logger.debug(f'Material is not known for {culvert.id}, therefore the default friction values is used')
    return culverts_dfm


def move_structure(struc, struc_dict, branch, offset):
    """
    Function the move a structure if needed for a compound event.

    Parameters
    ----------
    struc : string
        current sub-structure id
    struc_dict : dict
        dict with all structures of a certain type
    branch : string
        branch id of the first structure in the compound
    offset : float
        chainage of the first structure in the compound

    Returns
    -------
    Dict with shifted coordinates.

    """
    branch2 = struc_dict[struc]['branchid']
    if branch2!=branch:
       logger.warning(f'Structures of not on the same branche. Moving structure {struc} to branch {branch}.')
    struc_dict[struc]['branchid'] = branch
    struc_dict[struc]['chainage'] = offset
    return struc_dict

def generate_compounds(idlist, structurelist, structures, id_col):
    # probably the coordinates should all be set to those of the first structure (still to do)
    compounds_dfm = ExtendedDataFrame(required_columns=[id_col,'structurelist'])
    compounds_dfm.set_data(pd.DataFrame(np.zeros((len(idlist),3)), columns=[id_col,'numstructures','structurelist'], dtype='str'),index_col=id_col)
    compounds_dfm.index = idlist
    for ii,compound in enumerate(compounds_dfm.itertuples()):
        compounds_dfm.at[compound.Index, id_col] = idlist[ii]
        compounds_dfm.at[compound.Index, 'numstructures'] = len(structurelist[ii])

        # check the substructure coordinates. If they do not coincide, move subsequent structures to the coordinates of the first
        for s_i, struc in enumerate(structurelist[ii]):
            if s_i == 0:
                # find out what type the first structure it is and get its coordinates
                if struc in structures.culverts.keys():
                    branch = structures.culverts[struc]['branchid']
                    offset = structures.culverts[struc]['chainage']
                else:
                    raise IndexError('Structure id not found. Make sure all other structures have been added to the model.')
            else:
                # move a subsequent structure to the location of the first
                 if struc in structures.culverts.keys():
                     structures.culverts = move_structure(struc, structures.culverts, branch, offset)

        compounds_dfm.at[compound.Index, 'structurelist'] = ';'.join([f'{s}' for s in structurelist[ii]])

    return compounds_dfm

