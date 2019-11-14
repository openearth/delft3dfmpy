import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString

from delft3dfmpy.core import checks, geometry

logger = logging.getLogger(__name__)

roughness_gml = {
    1: "Chezy",
    2: "Manning",
    3: "StricklerNikuradse",
    4: "Strickler",
    5: "WhiteColebrook",
    6: "deBosBijkerk"
}


def generate_pumps(pompen, sturing, gemalen):
    """
    Generate pumps from hydamo data. The function combines the pumps
    with its steering and pumping stations to pumps that can be imported
    by dflowfm.

    Note that HYDAMO provides to ways of controlling pumps.
    1. By linking the control to the pumping station, and the pumps to the pumping station
    2. By linking the control directly to the pumps itself.

    The pumping capacity is given in cubic meters per minute in Hydamo. Dflowfm requires
    cubic meters per second, so the capacity is divided by 60.
    
    Parameters
    ----------
    pompen : gpd.GeoDataFrame
        Geometry and attributes of pumps
    sturing : gpd.GeoDataFrame
        Attributes of steering
    gemalen : gpd.GeoDataFrame
        Geometry and attributes of pumping stations
    
    Returns
    -------
    pd.DataFrame
        DataFrame with pump attributes suitable for dflowfm
    """
    # Copy dataframe
    pumps_dfm = pompen.copy()

    # HyDAMO contains m3/min, while D-Hydro needs m3/s
    pumps_dfm['maximalecapaciteit'] /= 60

    # Add sturing to pumps
    for idx, pump in pumps_dfm.iterrows():

        # Find sturing for pump
        sturingidx = (sturing.codegerelateerdobject == idx).values

        # The pump and 'sturing' might be linked to the pumping station,
        # so first check if there are multiple pumps with one 'sturing'
        if not sturingidx.sum() == 1:
            gemaalidx = (gemalen.code == pump.codegerelateerdobject).values
            # If there als multiple pumping stations connected to one pump, raise an error
            if sum(gemaalidx) != 1:
                raise IndexError('Multiple pumping stations (gemalen) found for pump.')

            # Find the idx if the pumping station connected to the pump
            gemaalidx = gemalen.iloc[np.where(gemaalidx)[0][0]]['code']
            # Find the control for the pumping station (and thus for the pump)
            sturingidx = (sturing.codegerelateerdobject == gemaalidx).values

            assert sum(sturingidx) == 1

        # Get the control by index
        pump_control = sturing.iloc[np.where(sturingidx)[0][0]]

        if pump_control.doelvariabelecode != 1 and pump_control.doelvariabelecode != 'waterstand':
            raise NotImplementedError('Sturing not implemented for anything else than water level (1).')

        # Add levels for suction side
        pumps_dfm.at[idx, 'startlevelsuctionside'] = pump_control['streefwaarde'] + pump_control['bovenmarge']
        pumps_dfm.at[idx, 'stoplevelsuctionside'] = pump_control['streefwaarde'] - pump_control['ondermarge']

    # Add name of .pli file for the pump
    pumps_dfm['locationfile'] = 'pump_'+pumps_dfm['code']+'.pli'

    return pumps_dfm

def generate_weirs(weirs):
    """
    Generate weirs from Hydamo input

    Currently only simple weirs can be applied. From Hydamo the attributes 
    'laagstedoorstroomhoogte' and 'kruinbreedte' are used to define the weir dimensions.

    The function already contains code for handling more complex weirs,
    but this code is not reached for now.
    
    Parameters
    ----------
    weirs : gpd.GeoDataFrame
        GeoDataFrame with geometry and attributes for weirs.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with weir attributes suitable for dflowfm
    """

    weirs_dfm = weirs.copy().astype('object')

    logger.info('Currently only simple weirs can be applied. From Hydamo the attributes \'laagstedoorstroomhoogte\' and \'kruinbreedte\' are used to define the weir dimensions.')

    return weirs_dfm

    ### THE CODE BELOW IS NOT REACHED. IT CAN BE USED WHEN UNIVERSAL WEIRS BECOME AVAILABLE IN DFLOWFM

    # # Create copy of geometry object with properties
    # self.weirs.set_data(
    #     self.geometries.weirs.reindex(columns=self.weirs.columns, fill_value='').astype('object'), index_col='code')

    for weir in weirs.itertuples():

        # Check levels
        if weir.laagstedoorstroomhoogte >= weir.hoogstedoorstroomhoogte:
            weirs.at[weir.Index, 'weirtype'] = 'weir'

        else:
            weirs.at[weir.Index, 'weirtype'] = 'weir'
            # The universal weir is not supported yet in D-Hydro!
            # self.weirs.at[weir.Index, 'weirtype'] = 'universal weir'

            # Create y,z-values
            yzvalues = [
                (-0.5 * weir.kruinbreedte, weir.hoogstedoorstroomhoogte),
                (-0.5 * weir.hoogstedoorstroombreedte, weir.hoogstedoorstroomhoogte),
                (-0.5 * weir.laagstedoorstroombreedte, weir.laagstedoorstroomhoogte),
                (0.5 * weir.laagstedoorstroombreedte, weir.laagstedoorstroomhoogte),
                (0.5 * weir.hoogstedoorstroombreedte, weir.hoogstedoorstroomhoogte),
                (0.5 * weir.kruinbreedte, weir.hoogstedoorstroomhoogte)
            ]

            # Remove duplicate values
            counts = [yzvalues.count(yz) for yz in yzvalues]
            while any([c > 1 for c in counts]):
                yzvalues.remove(yzvalues[counts.index(max(counts))])
                counts = [yzvalues.count(yz) for yz in yzvalues]

            # Add values
            weirs.at[weir.Index, 'ycoordinates'] = ' '.join([f'{yz[0]:7.3f}' for yz in yzvalues])
            weirs.at[weir.Index, 'zcoordinates'] = ' '.join([f'{yz[1]:7.3f}' for yz in yzvalues])
            weirs.at[weir.Index, 'levelscount'] = len(yzvalues)

def generate_culverts(culverts):

    culverts_dfm = culverts.copy()
    culverts_dfm['crosssection'] = [{} for _ in range(len(culverts_dfm))]

    for culvert in culverts.itertuples():

        # Generate cross section definition name
        if culvert.vormcode == 1 or culvert.vormcode == 'rond' or culvert.vormcode == 5 or culvert.vormcode == 'ellipsvormig':
            crosssection = {'shape': 'circle', 'diameter': culvert.hoogteopening}
            # definition = f'circ_d{culvert.hoogteopening:.3f}'
        elif culvert.vormcode == 3 or culvert.vormcode == 'rechthoekig' or culvert.vormcode == 99 or culvert.vormcode == 'onbekend':
            crosssection = {'shape': 'rectangle', 'height': culvert.hoogteopening, 'width': culvert.breedteopening, 'closed': 1}
            # definition = f'rect_h{culvert.hoogteopening:.3f}_w{culvert.breedteopening:.3f}'

        # Set cross section definition
        culverts_dfm.at[culvert.Index, 'crosssection'] = crosssection

    return culverts_dfm

def dwarsprofiel_to_yzprofiles(crosssections):
    """
    Function to convert hydamo cross sections 'dwarsprofiel' to
    dflowfm input.
    
    Parameters
    ----------
    crosssections : gpd.GeoDataFrame
        GeoDataFrame with x,y,z-coordinates of cross sections
    
    Returns
    -------
    dictionary
        Dictionary with attributes of cross sections, usable for dflowfm
    """
    cssdct = {}

    for css in crosssections.itertuples():
        # The cross sections from hydamo are all yz profiles
        # Determine yz_values
        xyz = np.vstack(css.geometry.coords[:])
        length = np.r_[0, np.cumsum(np.hypot(np.diff(xyz[:, 0]), np.diff(xyz[:, 1])))]
        yz = np.c_[length, xyz[:, -1]]

        # Add to dictionary
        cssdct[css.code] = {
            'branchid': css.branch_id,
            'chainage': css.branch_offset,
            'yz': yz,
            'ruwheidstypecode': css.ruwheidstypecode,
            'ruwheidswaarde': css.ruwheidswaarde
        }
    
    return cssdct

def parameterised_to_profiles(parameterised, branches):
    """
    Generate parametrised cross sections for all branches,
    or the branches missing a cross section.

    Parameters
    ----------
    parameterised : pd.DataFrame
        GeoDataFrame with geometries and attributes of parameterised profiles.
    branches : list
        List of branches for which the parameterised profiles are derived

    Returns
    -------
    dictionary
        Dictionary with attributes of cross sections, usable for dflowfm
    """

    checks.check_argument(parameterised, 'parameterised', (pd.DataFrame, gpd.GeoDataFrame))
    checks.check_argument(branches, 'branches', (list, tuple))

    # Find
    if branches is not None:
        parambranches = parameterised.reindex(index=branches, columns=parameterised.required_columns + ['css_type'])
    else:
        parambranches = parameterised.reindex(columns=parameterised.required_columns + ['css_type'])

    # Drop profiles for which not enough data is available to write (as rectangle)
    nulls = pd.isnull(parambranches[['bodembreedte', 'bodemhoogtebenedenstrooms', 'bodemhoogtebovenstrooms']]).any(axis=1).values
    parambranches.drop(parambranches.index[nulls], inplace=True)

    # Determine characteristics
    botlev = (parambranches['bodemhoogtebenedenstrooms'] + parambranches['bodemhoogtebovenstrooms']) / 2.0
    dh1 = parambranches['hoogteinsteeklinkerzijde'] - botlev
    dh2 = parambranches['hoogteinsteekrechterzijde'] - botlev
    parambranches['height'] = (dh1 + dh2) / 2.0
    parambranches['bottomlevel'] = botlev

    # Determine maximum flow width and slope (both needed for output)
    parambranches['maxflowwidth'] = parambranches['bodembreedte'] + parambranches['taludhellinglinkerzijde'] * dh1 + parambranches['taludhellingrechterzijde'] * dh2
    parambranches['slope'] = (parambranches['taludhellinglinkerzijde'] + parambranches['taludhellingrechterzijde']) / 2.0

    # Determine profile type
    parambranches.loc[:, 'css_type'] = 'trapezium'
    nulls = pd.isnull(parambranches[parameterised.required_columns]).any(axis=1).values
    parambranches.loc[nulls, 'css_type'] = 'rectangle'

    cssdct = {}
    for branch in parambranches.itertuples():
        # Determine name for cross section
        if branch.css_type == 'trapezium':
            cssdct[branch.Index] = {
                'type': branch.css_type,
                'slope': round(branch.slope, 1),
                'maximumflowwidth': round(branch.maxflowwidth, 1),
                'bottomwidth': round(branch.bodembreedte, 3),
                'closed': 0,
                'ruwheidstypecode': int(branch.ruwheidstypecode) if isinstance(branch.ruwheidstypecode, float) else branch.ruwheidstypecode,
                'ruwheidswaarde': branch.ruwheidswaarde,
                'bottomlevel': branch.bottomlevel
            }
        elif branch.css_type == 'rectangle':
            cssdct[branch.Index] = {
                'type': branch.css_type,
                'height': 5.0,
                'width': round(branch.bodembreedte, 3),
                'closed': 0,
                'ruwheidstypecode': int(branch.ruwheidstypecode) if isinstance(branch.ruwheidstypecode, float) else branch.ruwheidstypecode,
                'ruwheidswaarde': branch.ruwheidswaarde,
                'bottomlevel': branch.bottomlevel
            }

    return cssdct

def generate_boundary_conditions(boundary_conditions, schematised):
    """
    Generate boundary conditions from hydamo 'randvoorwaarden' file.

    Parameters
    ----------
    boundary_conditions: gpd.GeoDataFrame
        geodataframe with the locations and properties of the boundary conditions
    schematised : gpd.GeoDataFrame
        geodataframe with the schematised branches
    
    Returns
    -------
    dictionary
        Dictionary with attributes of boundary conditions, usable for dflowfm
    """
    bcdct = {}

    for bndcnd in boundary_conditions.itertuples():

        # Find nearest branch for geometry
        extended_line = geometry.extend_linestring(
            line=schematised.at[bndcnd.branch_id, 'geometry'], near_pt=bndcnd.geometry, length=1.0)

        # Create intersection line for boundary condition
        bcline = LineString(geometry.orthogonal_line(line=extended_line, offset=0.1, width=0.1))

        if bndcnd.typerandvoorwaardecode in [0, 'waterstand']:
            bctype = 'waterlevel'
        elif bndcnd.typerandvoorwaardecode in [1, 'debiet']:
            bctype = 'discharge'

        # Add boundary condition
        bcdct[bndcnd.code] = {
            'code': bndcnd.code,
            'bctype': bctype+'bnd',
            'value': bndcnd.waterstand if not np.isnan(bndcnd.waterstand) else bndcnd.debiet,
            'time': None,
            'geometry': bcline,
            'filetype': 9,
            'operand': 'O',
            'method': 3,
            'branchid': bndcnd.branch_id
        }

    return bcdct
