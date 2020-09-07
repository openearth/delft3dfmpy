import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString

from delft3dfmpy.core import checks, geometry
from delft3dfmpy.datamodels.common import ExtendedDataFrame, ExtendedGeoDataFrame

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
        pumps_dfm.at[idx, 'startlevelsuctionside'] = pump_control['bovenmarge'] 
        pumps_dfm.at[idx, 'stoplevelsuctionside'] = pump_control['ondermarge'] 
    
    return pumps_dfm

def generate_weirs(weirs, afsluitmiddel=None, sturing=None):
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

def generate_orifices(orifices, afsluitmiddel=None, sturing=None):
    """
    Generate orifices from Hydamo input

    Parameters
    ----------
    orifices : gpd.GeoDataFrame
        GeoDataFrame with geometry and attributes for weirs.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with orifice attributes suitable for dflowfm
    """

    orifices_dfm = orifices.copy().astype('object')
    if 'maximaaldebiet' not in orifices_dfm:  
        orifices_dfm['uselimitflow'] = 'false'
        orifices_dfm['limitflow'] = 0.0
    else:
        orifices_dfm['uselimitflow'] = 'true'
        orifices_dfm['limitflow'] = orifices_dfm['maximaaldebiet']
    return orifices_dfm


def generate_uweirs(uweirs, yz_profiles=None, parametrised_profiles=None):
   """
   Generate universal weirs from Hydamo input

   Crest level is determined from the Hydamo laagstedoorstroomhoogte attribute. The (relative) profile is determined from the crossection with codegerelateerdobject pointing to the universal weir.
    
   Parameters
   ---------
   uweirs : gpd.GeoDataFrame
       GeoDataFrame with geometry and attributes for universal weirs.
    
   Returns
   -------
   pd.DataFrame
       DataFrame with universal weir attributes suitable for dflowfm

   """
   
   uweirs_dfm = uweirs.copy().astype('object')
   uweirs_dfm['crosssection'] = [{} for _ in range(len(uweirs_dfm))]

   for uweir in uweirs.itertuples():    
       
       # first search in yz-profiles
       prof=np.empty(0)
       if yz_profiles is not None:
            if 'codegerelateerdobject' in yz_profiles:  
                prof = yz_profiles[yz_profiles['codegerelateerdobject']==uweir.code]   
                if not prof.empty:                    
                    counts = len(prof.geometry[0].coords[:])
                    xyz = np.vstack(prof.geometry[0].coords[:])
                    length = np.r_[0, np.cumsum(np.hypot(np.diff(xyz[:, 0]), np.diff(xyz[:, 1])))]
                    yzvalues = np.c_[length, xyz[:, -1]-np.min(xyz[:,-1])]            
                
       if (prof.empty) & (parametrised_profiles is not None):
            if 'codegerelateerdobject' in parametrised_profiles:  
                # if not found, check the parametrised profiles
                prof = parametrised_profiles[parametrised_profiles['codegerelateerdobject']==uweir.code]           
                nulls = pd.isnull(prof[['bodembreedte', 'bodemhoogtebenedenstrooms', 'bodemhoogtebovenstrooms','taludhellinglinkerzijde','taludhellingrechterzijde','hoogteinsteeklinkerzijde','hoogteinsteekrechterzijde']]).any(axis=1).values
                if nulls:
                   raise ValueError(f'Insufficient fields avaialable in parametrised profile defintion {prof.code}.')
                bodemhoogte = float((prof.bodemhoogtebenedenstrooms + prof.bodemhoogtebovenstrooms)/2.)
                zvalues = [ float(prof.hoogteinsteeklinkerzijde-bodemhoogte), 0., 0., float( prof.hoogteinsteekrechterzijde-bodemhoogte )]
                yvalues = list( np.cumsum( [0.,float((prof.hoogteinsteeklinkerzijde-bodemhoogte)/prof.taludhellinglinkerzijde),float(prof.bodembreedte), float((prof.hoogteinsteekrechterzijde-bodemhoogte)/prof.taludhellingrechterzijde) ] ) )           
                yzvalues = list(zip(yvalues,zvalues))
                counts = len(zvalues)
           
       if len(prof)==0:
           # return an error it is still not found
           raise ValueError(f'{uweir.code} is not found in any cross-section.')
    
       uweirs_dfm.at[uweir.Index, 'numlevels'] = counts
       uweirs_dfm.at[uweir.Index, 'yvalues'] = ' '.join([f'{yz[0]:7.3f}' for yz in yzvalues])
       uweirs_dfm.at[uweir.Index, 'zvalues'] = ' '.join([f'{yz[1]:7.3f}' for yz in yzvalues])        
             
   return uweirs_dfm
            
        # # Check levels
        # if weir.laagstedoorstroomhoogte >= weir.hoogstedoorstroomhoogte:
        #     weirs.at[weir.Index, 'weirtype'] = 'weir'

        # else:
        #     weirs.at[weir.Index, 'weirtype'] = 'weir'
        #     # The universal weir is not supported yet in D-Hydro!
        #     # self.weirs.at[weir.Index, 'weirtype'] = 'universal weir'

        #     # Create y,z-values
        #     yzvalues = [
        #         (-0.5 * weir.kruinbreedte, weir.hoogstedoorstroomhoogte),
        #         (-0.5 * weir.hoogstedoorstroombreedte, weir.hoogstedoorstroomhoogte),
        #         (-0.5 * weir.laagstedoorstroombreedte, weir.laagstedoorstroomhoogte),
        #         (0.5 * weir.laagstedoorstroombreedte, weir.laagstedoorstroomhoogte),
        #         (0.5 * weir.hoogstedoorstroombreedte, weir.hoogstedoorstroomhoogte),
        #         (0.5 * weir.kruinbreedte, weir.hoogstedoorstroomhoogte)
        #     ]

        #     # Remove duplicate values
        #     counts = [yzvalues.count(yz) for yz in yzvalues]
        #     while any([c > 1 for c in counts]):
        #         yzvalues.remove(yzvalues[counts.index(max(counts))])
        #         counts = [yzvalues.count(yz) for yz in yzvalues]

        #     # Add values
        #     weirs.at[weir.Index, 'ycoordinates'] = ' '.join([f'{yz[0]:7.3f}' for yz in yzvalues])
        #     weirs.at[weir.Index, 'zcoordinates'] = ' '.join([f'{yz[1]:7.3f}' for yz in yzvalues])
        #     weirs.at[weir.Index, 'levelscount'] = len(yzvalues)

def generate_bridges(bridges, yz_profiles=None, parametrised_profiles=None):
    """
    Generate bridges from Hydamo input

    Parameters
    
    ----------
    bridges : gpd.GeoDataFrame
        GeoDataFrame with geometry and attributes for bridges.
    
    Returns
    -------
    pd.DataFrame
        DataFrame with weir attributes suitable for dflowfm
    """

    bridges_dfm = bridges.copy().astype('object')
    bridges_dfm['crosssection'] = [{} for _ in range(len(bridges_dfm))]
    bridges_dfm['bedlevel'] = [{} for _ in range(len(bridges_dfm))]
    
    for bridge in bridges.itertuples():        
        # first search in yz-profiles
        prof = yz_profiles[yz_profiles['codegerelateerdobject']==bridge.code]
        if len(prof) > 0:
            bedlevel = np.min([c[2] for c in prof.geometry[0].coords[:]])  
            profile_id=prof.code.values[0]
        else:
            # if not found, check the parametrised profiles
            prof = parametrised_profiles[parametrised_profiles['codegerelateerdobject']==bridge.code]
            bedlevel = (prof['bodemhoogtebovenstrooms'] + prof['bodemhoogtebenedenstrooms'])/2.
           
        if len(prof)==0:
            # return an error it is still not found
            raise ValueError(f'{bridge.code} is not found in any cross-section.')
        
        profile_id=prof.code.values[0]
        bridges_dfm.at[bridge.Index, 'crosssection'] = profile_id
        bridges_dfm.at[bridge.Index, 'bedlevel'] = float(bedlevel)
             
    return bridges_dfm
          
def generate_culverts(culverts,afsluitmiddel):

    culverts_dfm = culverts.copy()
    culverts_dfm['crosssection'] = [{} for _ in range(len(culverts_dfm))]
        
    for culvert in culverts.itertuples():

        # Generate cross section definition name
        if culvert.vormcode == 1 or culvert.vormcode == 'rond' or culvert.vormcode == 5 or culvert.vormcode == 'ellipsvormig':
            crosssection = {'shape': 'circle', 'diameter': culvert.hoogteopening}
            
        elif culvert.vormcode == 3 or culvert.vormcode == 'rechthoekig' or culvert.vormcode == 99 or culvert.vormcode == 'onbekend':
            crosssection = {'shape': 'rectangle', 'height': culvert.hoogteopening, 'width': culvert.breedteopening, 'closed': 1}
        
        else:
            crosssection = {'shape': 'circle', 'diameter': 0.40}
            print(f'Culvert {culvert.code} has an unknown shape: {culvert.vormcode}. Applying a default profile (round - 40cm)')
        
        # Set cross section definition
        culverts_dfm.at[culvert.Index, 'allowedflowdir'] = 'both'
        culverts_dfm.at[culvert.Index, 'valveonoff'] = 0
        culverts_dfm.at[culvert.Index, 'numlosscoeff'] = 0
        culverts_dfm.at[culvert.Index, 'valveopeningheight'] = 0
        culverts_dfm.at[culvert.Index, 'relopening'] = 0
        culverts_dfm.at[culvert.Index, 'losscoeff'] = 0
        # check whether an afsluitmiddel is present and take action dependent on its settings
        if afsluitmiddel is not None:
            if not afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code].empty:
                if len(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code])!=1:
                    raise IndexError(f'Multiple (or no) instances of afsluitmiddel associated with culvert {culvert.code}')
                if int(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code]['soortafsluitmiddelcode'])==5:
                    culverts_dfm.at[culvert.Index, 'allowedflowdir'] = 'positive'
                if int(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code]['soortafsluitmiddelcode'])==4:
                    culverts_dfm.at[culvert.Index, 'valveonoff'] = 1
                    culverts_dfm.at[culvert.Index, 'valveopeningheight'] = float(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code]['hoogte'])
                    culverts_dfm.at[culvert.Index, 'numlosscoeff'] = 1
                    culverts_dfm.at[culvert.Index, 'relopening'] = float(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code]['hoogte'])/culvert.hoogteopening
                    culverts_dfm.at[culvert.Index, 'losscoeff'] = float(afsluitmiddel[afsluitmiddel.codegerelateerdobject==culvert.code]['afvoercoefficient'])
        culverts_dfm.at[culvert.Index, 'crosssection'] = crosssection

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

def generate_compounds(idlist, structurelist, structures):
    # probably the coordinates should all be set to those of the first structure (still to do)        
    compounds_dfm = ExtendedDataFrame(required_columns=['code','structurelist'])
    compounds_dfm.set_data(pd.DataFrame(np.zeros((len(idlist),3)), columns=['code','numstructures','structurelist'], dtype='str'),index_col='code')
    compounds_dfm.index = idlist    
    for ii,compound in enumerate(compounds_dfm.itertuples()):
        compounds_dfm.at[compound.Index, 'code'] = idlist[ii]
        compounds_dfm.at[compound.Index, 'numstructures'] = len(structurelist[ii])
        
        # check the substructure coordinates. If they do not coincide, move subsequent structures to the coordinates of the first
        for s_i, struc in enumerate(structurelist[ii]):            
            if s_i == 0:                
                # find out what type the first structure it is and get its coordinates
                if struc in structures.pumps.keys():
                    branch = structures.pumps[struc]['branchid']
                    offset = structures.pumps[struc]['chainage']
                elif struc in structures.weirs.keys():
                    branch = structures.weirs[struc]['branchid']
                    offset = structures.weirs[struc]['chainage']
                elif struc in structures.uweirs.keys():
                    branch = structures.uweirs[struc]['branchid']
                    offset = structures.uweirs[struc]['chainage']
                elif struc in structures.culverts.keys():
                    branch = structures.culverts[struc]['branchid']
                    offset = structures.culverts[struc]['chainage']                    
                elif struc in structures.bridges.keys():
                    branch = structures.bridges[struc]['branchid']
                    offset = structures.bridges[struc]['chainage']                    
                elif struc in structures.orifices.keys():
                    branch = structures.orifices[struc]['branchid']
                    offset = structures.orifices[struc]['chainage']                                        
                else:
                    raise IndexError('Structure id not found. Make sure all other structures have been added to the model.')
            else:
                # move a subsequent structure to the location of the first
                 if struc in structures.pumps.keys():                    
                     structures.pumps = move_structure(struc, structures.pumps, branch, offset)
                 if struc in structures.weirs.keys():                    
                     structures.weirs = move_structure(struc, structures.weirs, branch, offset)
                 if struc in structures.uweirs.keys():                    
                     structures.uweirs = move_structure(struc, structures.uweirs, branch, offset)
                 if struc in structures.culverts.keys():                    
                     structures.culverts = move_structure(struc, structures.culverts, branch, offset)
                 if struc in structures.bridges.keys():                    
                     structures.bridges = move_structure(struc, structures.bridges, branch, offset)
                 if struc in structures.orifices.keys():                    
                     structures.orifices = move_structure(struc, structures.orifices, branch, offset)                 
                                     
        compounds_dfm.at[compound.Index, 'structurelist'] = ';'.join([f'{s}' for s in structurelist[ii]])
    
    return compounds_dfm

def dwarsprofiel_to_yzprofiles(crosssections, branches):
    """
    Function to convert hydamo cross sections 'dwarsprofiel' to
    dflowfm input.
    d
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
        
       # if css.Index == 'prof_RS1-DP-27341':
       #     print('stop')
        # Determine yz_values
        xyz = np.vstack(css.geometry.coords[:])
        length = np.r_[0, np.cumsum(np.hypot(np.diff(xyz[:, 0]), np.diff(xyz[:, 1])))]
        yz = np.c_[length, xyz[:, -1]]
        # the GUI cannot cope with identical y-coordinates. Add 1 cm to a 2nd duplicate.
        yz[:,0] = np.round(yz[:,0],3)
        for i in range(1,yz.shape[0]):
            if yz[i,0]==yz[i-1,0]:
                yz[i,0] +=0.01
                
        # determine thalweg
        if branches is not None:
            branche_geom = branches[branches.code==css.branch_id].geometry.values        
            thalweg_xyz = css.geometry.intersection(branche_geom[0]).coords[:][0]                
            # and the Y-coordinate of the thalweg
            thalweg = np.hypot( thalweg_xyz[0]-xyz[0,0], thalweg_xyz[1]-xyz[0,1])
        else: 
            thalweg = 0.0
        
        # Add to dictionary
        cssdct[css.code] = {
            'branchid': css.branch_id,
            'chainage': css.branch_offset,
            'yz': yz,
            'thalweg':thalweg,
            'ruwheidstypecode': css.ruwheidstypecode,
            'ruwheidswaarde': css.ruwheidswaarde
        }
    
    return cssdct

def parametrised_to_profiles(parametrised, branches):
    """
    Generate parametrised cross sections for all branches,
    or the branches missing a cross section.

    Parameters
    ----------
    parametrised : pd.DataFrame
        GeoDataFrame with geometries and attributes of parametrised profiles.
    branches : list
        List of branches for which the parametrised profiles are derived

    Returns
    -------
    dictionary
        Dictionary with attributes of cross sections, usable for dflowfm
    """

    checks.check_argument(parametrised, 'parametrised', (pd.DataFrame, gpd.GeoDataFrame))
    checks.check_argument(branches, 'branches', (list, tuple))

    # Find
    if len(branches) != 0:
         parambranches = ExtendedGeoDataFrame(geotype=LineString, columns = parametrised.columns.tolist()+['css_type'])
         parambranches.set_data(parametrised[np.isin(parametrised.code, branches)], index_col='code', check_columns=True)
    else:
        parambranches = parametrised
     #   
    #else:
        #parambranches = parametrised.reindex(columns=parametrised.requiredcolumns + ['css_type'])
    #    parambranches = parametrised[np.isin(parametrised.code,branches)]# ['css_type'])
   
   
    # Drop profiles for which not enough data is available to write (as rectangle)
    #nulls = pd.isnull(parambranches[['bodembreedte', 'bodemhoogtebenedenstrooms', 'bodemhoogtebovenstrooms']]).any(axis=1).values
    #parambranches = parambranches.drop(ExtendedGeoDataFrame(geotype=LineString), parambranches.index[nulls], index_col='code',axis=0)
    #parambranches.drop(parambranches.index[nulls], inplace=True)

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
    nulls = pd.isnull(parambranches[parametrised.required_columns]).any(axis=1).values
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
                'thalweg': 0.0,
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
                'thalweg': 0.0,
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
        # extended_line = geometry.extend_linestring(
        #     line=schematised.at[bndcnd.branch_id, 'geometry'], near_pt=bndcnd.geometry, length=1.0)

        # # Create intersection line for boundary condition
        # bcline = LineString(geometry.orthogonal_line(line=extended_line, offset=0.1, width=0.1))

        if bndcnd.typerandvoorwaardecode in [0, 'waterstand']:
            bctype = 'waterlevel'
        elif bndcnd.typerandvoorwaardecode in [1, 'debiet']:
            bctype = 'discharge'

        # Add boundary condition
        bcdct[bndcnd.code] = {
            'code': bndcnd.code,
            'bctype': bctype,#+'bnd',
            'value': bndcnd.waterstand if not np.isnan(bndcnd.waterstand) else bndcnd.debiet,
            'time': None,
            'geometry': bndcnd.geometry,#bcline,
            'filetype': 9,
            'operand': 'O',
            'method': 3,
            'branchid': bndcnd.branch_id
        }

    return bcdct
