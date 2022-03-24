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
        sturingidx = (sturing.pompid == pump.globalid).values
        
        # find gemaal for pump
        gemaalidx = (gemalen.globalid == pump.gemaalid).values
        
        # so first check if there are multiple pumps with one 'sturing'
        if sum(sturingidx) != 1:            
            raise IndexError('Multiple or no sturingen found for pump.')
            
        # If there als multiple pumping stations connected to one pump, raise an error
        if sum(gemaalidx) != 1:
            raise IndexError('Multiple or no pumping stations (gemalen) found for pump.')

        # Find the idx if the pumping station connected to the pump
            # gemaalidx = gemalen.iloc[np.where(gemaalidx)[0][0]]['code']
            # Find the control for the pumping station (and thus for the pump)
            #@sturingidx = (sturing.codegerelateerdobject == gemaalidx).values

            #assert sum(sturingidx) == 1

        pumps_dfm.at[idx, 'branch_id'] = gemalen.iloc[np.where(gemaalidx)[0][0]]['branch_id']
        pumps_dfm.at[idx, 'branch_offset'] = gemalen.iloc[np.where(gemaalidx)[0][0]]['branch_offset']
        # Get the control by index
        pump_control = sturing.iloc[np.where(sturingidx)[0][0]]

        if pump_control.doelvariable != 1 and pump_control.doelvariable != 'waterstand':
            raise NotImplementedError('Sturing not implemented for anything else than water level (1).')

        # Add levels for suction side
        pumps_dfm.at[idx, 'startlevelsuctionside'] = pump_control['bovengrens'] 
        pumps_dfm.at[idx, 'stoplevelsuctionside'] = pump_control['ondergrens'] 
    
    return pumps_dfm

def generate_weirs(weirs, opening=None, management_device=None, management=None):
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
    weirs_dfm = gpd.GeoDataFrame()
    orifices_dfm = gpd.GeoDataFrame()
    for idx, weir in weirs.iterrows():
        weir_opening = opening[opening.stuwid == weir.globalid]
        weir_mandev = management_device[management_device.kunstwerkopeningid == weir_opening.globalid.to_string(index=False)]
        if weir_mandev.overlaatonderlaat.to_string(index=False) == 'Overlaat':
            weirs_dfm.at[idx,'code'] = weir.code
            weirs_dfm.at[idx,'branch_id'] = weir.branch_id
            weirs_dfm.at[idx,'branch_offset'] = weir.branch_offset
            weirs_dfm.at[idx,'laagstedoorstroomhoogte'] = weir_opening.laagstedoorstroomhoogte.to_string(index=False)
            weirs_dfm.at[idx,'laagstedoorstroombreedte'] = weir_opening.laagstedoorstroombreedte.to_string(index=False)
            weirs_dfm.at[idx,'afvoercoefficient'] = weir.afvoercoefficient
        elif weir_mandev.overlaatonderlaat.to_string(index=False) == 'Onderlaat':
            orifices_dfm.at[idx,'code'] = weir.code
            orifices_dfm.at[idx,'branch_id'] = weir.branch_id
            orifices_dfm.at[idx,'branch_offset'] = weir.branch_offset
            orifices_dfm.at[idx,'laagstedoorstroomhoogte'] = weir_opening.laagstedoorstroomhoogte.to_string(index=False)
            orifices_dfm.at[idx,'laagstedoorstroombreedte'] = weir_opening.laagstedoorstroombreedte.to_string(index=False)
            orifices_dfm.at[idx,'afvoercoefficient'] = weir.afvoercoefficient
            orifices_dfm.at[idx,'schuifhoogte'] = weir_mandev.hoogteopening.to_string(index=False)
            if 'maximaaldebiet' not in weir_mandev:  
                orifices_dfm.at[idx,'uselimitflow'] = 'false'
                orifices_dfm.at[idx,'limitflow'] = 0.0
            else:
                orifices_dfm.at[idx,'uselimitflow'] = 'true'
                orifices_dfm.at[idx,'limitflow'] = weir_mandev.maximaaldebiet.to_string(index=False)
            
    return [weirs_dfm, orifices_dfm]
                        
    #logger.info('Currently only simple weirs can be applied. From Hydamo the attributes \'laagstedoorstroomhoogte\' and \'kruinbreedte\' are used to define the weir dimensions.')

    

# def generate_orifices(orifices, afsluitmiddel=None, sturing=None):
#     """
#     Generate orifices from Hydamo input

#     Parameters
#     ----------
#     orifices : gpd.GeoDataFrame
#         GeoDataFrame with geometry and attributes for weirs.
    
#     Returns
#     -------
#     pd.DataFrame
#         DataFrame with orifice attributes suitable for dflowfm
#     """

#     orifices_dfm = orifices.copy().astype('object')
#     if 'maximaaldebiet' not in orifices_dfm:  
#         orifices_dfm['uselimitflow'] = 'false'
#         orifices_dfm['limitflow'] = 0.0
#     else:
#         orifices_dfm['uselimitflow'] = 'true'
#         orifices_dfm['limitflow'] = orifices_dfm['maximaaldebiet']
#     return orifices_dfm


def generate_uweirs(uweirs, yz_profiles=None):
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
            if 'stuwid' in yz_profiles:  
                prof = yz_profiles[yz_profiles['stuwid']==uweir.globalid]   
                if not prof.empty:                    
                    counts = len(prof.geometry.iloc[0].coords[:])
                    xyz = np.vstack(prof.geometry.iloc[0].coords[:])
                    length = np.r_[0, np.cumsum(np.hypot(np.diff(xyz[:, 0]), np.diff(xyz[:, 1])))]
                    yzvalues = np.c_[length, xyz[:, -1]-np.min(xyz[:,-1])]                            
     
           
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

def generate_bridges(bridges, yz_profiles=None):
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
    bridges_dfm['shift'] = [0.0 for _ in range(len(bridges_dfm))]
    
    for bridge in bridges.itertuples():        
        # first search in yz-profiles
        prof = yz_profiles[yz_profiles['brugid']==bridge.globalid]
        if len(prof) > 0:
            #bedlevel = np.min([c[2] for c in prof.geometry[0].coords[:]])  
            profile_id=prof.code.values[0]
        else:
            # return an error it is still not found
            raise ValueError(f'{bridge.code} is not found in any cross-section.')
        
        profile_id=prof.code.values[0]
        bridges_dfm.at[bridge.Index, 'crosssection'] = profile_id
        #bridges_dfm.at[bridge.Index, 'bedlevel'] = float(bedlevel)
             
    return bridges_dfm
          
def generate_culverts(culverts,management_device=None):

    culverts_dfm = culverts.copy()
    culverts_dfm['crosssection'] = [{} for _ in range(len(culverts_dfm))]
        
    for culvert in culverts.itertuples():

        # Generate cross section definition name
        if culvert.vormkoker == 'Rond' or culvert.vormkoker == 'Ellipsvormig':
            crosssection = {'shape': 'circle', 'diameter': culvert.hoogteopening}
            
        elif culvert.vormkoker == 'Rechthoekig' or culvert.vormkoker == 'Onbekend' or culvert.vormkoker=='Muilprofiel' or culvert.vormkoker=='Heulprofiel':
            crosssection = {'shape': 'rectangle', 'height': culvert.hoogteopening, 'width': culvert.breedteopening, 'closed': 1}
        
        else:
            crosssection = {'shape': 'circle', 'diameter': 0.40}
            print(f'Culvert {culvert.code} has an unknown shape: {culvert.vormkoker}. Applying a default profile (round - 40cm)')
        
        # Set cross section definition
        culverts_dfm.at[culvert.Index, 'allowedflowdir'] = 'both'
        culverts_dfm.at[culvert.Index, 'valveonoff'] = 0
        culverts_dfm.at[culvert.Index, 'numlosscoeff'] = 0
        culverts_dfm.at[culvert.Index, 'valveopeningheight'] = 0
        culverts_dfm.at[culvert.Index, 'relopening'] = 0
        culverts_dfm.at[culvert.Index, 'losscoeff'] = 0
        # check whether an afsluitmiddel is present and take action dependent on its settings
        if management_device is not None:
            
            if not management_device[management_device.duikersifonhevelid==culvert.globalid].empty:
                if len(management_device[management_device.duikersifonhevelid==culvert.globalid])==0:
                    raise IndexError(f'No instances of closing_device associated with culvert {culvert.code}')
                for _,i in management_device[management_device.duikersifonhevelid==culvert.globalid].iterrows():
                    if i['soortregelmiddel']=='terugslagklep':
                        culverts_dfm.at[culvert.Index, 'allowedflowdir'] = 'positive'
                    elif i['soortregelmiddel']=='schuif':
                        culverts_dfm.at[culvert.Index, 'valveonoff'] = 1
                        culverts_dfm.at[culvert.Index, 'valveopeningheight'] = float(i['hoogteopening'])
                        culverts_dfm.at[culvert.Index, 'numlosscoeff'] = 1
                        culverts_dfm.at[culvert.Index, 'relopening'] = float(i['hoogteopening'])/culvert.hoogteopening
                        culverts_dfm.at[culvert.Index, 'losscoeff'] = float(i['afvoercoefficient'])
                    else:
                        print(f'Type of closing device for culvert {culvert.code} is not implemented; only "schuif" and "terugslagklep" are allowed.')
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

def dwarsprofiel_to_yzprofiles(crosssections, roughness, branches, roughness_variant='Low'):
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
            if css.geometry.intersection(branche_geom[0]).geom_type=='MultiPoint':
                thalweg_xyz = css.geometry.intersection(branche_geom[0])[0].coords[:][0]                            
            else:
                thalweg_xyz = css.geometry.intersection(branche_geom[0]).coords[:][0]                
            # and the Y-coordinate of the thalweg
            thalweg = np.hypot( thalweg_xyz[0]-xyz[0,0], thalweg_xyz[1]-xyz[0,1])
        else: 
            thalweg = 0.0
        
        if roughness_variant == "High":
            ruwheid=roughness[roughness['profielpuntid']==css.globalid].ruwheidhoog
        if roughness_variant == "Low":
            ruwheid=roughness[roughness['profielpuntid']==css.globalid].ruwheidlaag
            
        # Add to dictionary
        cssdct[css.code] = {
            'branchid': css.branch_id,
            'chainage': css.branch_offset,
            'yz': yz,
            'thalweg':thalweg,
            'typeruwheid': roughness[roughness['profielpuntid']==css.globalid].typeruwheid.to_string(index=False),
            'ruwheid': float(ruwheid)
        }
    
    return cssdct

def parametrised_to_profiles(parametrised, parametrised_values,  branches, roughness_variant='Low'):
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
    
    cssdct = {}
    for param in parametrised.itertuples():
        branch = [branch for branch in branches if branch.globalid==param.hydroobjectid]
    
        values = parametrised_values[parametrised_values.normgeparamprofielid==param.normgeparamprofielid]
    
        #Drop profiles for which not enough data is available to write (as rectangle)
        # nulls = pd.isnull(parambranches[['bodembreedte', 'bodemhoogtebenedenstrooms', 'bodemhoogtebovenstrooms']]).any(axis=1).values
        # parambranches = parambranches.drop(ExtendedGeoDataFrame(geotype=LineString), parambranches.index[nulls], index_col='code',axis=0)
        # parambranches.drop(parambranches.index[nulls], inplace=True)
        
        if pd.isnull(values[values.soortparameter=='bodemhoogte benedenstrooms'].waarde).values[0]:
            print('bodemhoogte benedenstrooms not available for profile {}.'.format(param.globalid))
        if pd.isnull(values[values.soortparameter=='bodembreedte'].waarde).values[0]:
            print('bodembreedte not available for profile {}.'.format(param.globalid))
        if pd.isnull(values[values.soortparameter=='bodemhoogte bovenstrooms'].waarde).values[0]:
            print('bodemhoogte bovenstrooms not available for profile {}.'.format(param.globalid))
        
        # Determine characteristics
        botlev = (values[values.soortparameter=='bodemhoogte benedenstrooms'].waarde.values[0] + values[values.soortparameter=='bodemhoogte benedenstrooms'].waarde.values[0]) / 2.0       
        
        
        if pd.isnull(values[values.soortparameter=='taludhelling linkerzijde'].waarde).values[0]:
            csstype=='rectangle'
        else:
            css_type = 'trapezium'
            dh1 = values[values.soortparameter=='hoogte insteek linkerzijde'].waarde.values[0] - botlev
            dh2 = values[values.soortparameter=='hoogte insteek rechterzijde'].waarde.values[0] - botlev
            height = (dh1 + dh2) / 2.0
            # Determine maximum flow width and slope (both needed for output)
            maxflowwidth = values[values.soortparameter=='bodembreedte'].waarde.values[0] + values[values.soortparameter=='taludhelling linkerzijde'].waarde.values[0] * dh1 + values[values.soortparameter=='taludhelling rechterzijde'].waarde.values[0] * dh2
            slope = (values[values.soortparameter=='taludhelling linkerzijde'].waarde.values[0] + values[values.soortparameter=='taludhelling rechterzijde'].waarde.values[0]) / 2.0
           
        if roughness_variant=='Low':
            roughness =  values.ruwheidlaag.values[0]
        else:
            roughness = values.ruwheidhoog.values[0]
        # Determine name for cross section
        if css_type == 'trapezium':
            cssdct[branch[0].Index] = {
                'type': css_type,
                'slope': round(slope, 2),
                'maximumflowwidth': round(maxflowwidth, 1),
                'bottomwidth': round(values[values.soortparameter=='bodembreedte'].waarde.values[0], 3),
                'closed': 0,
                'thalweg': 0.0,
                'typeruwheid': values.typeruwheid.values[0],
                'ruwheid': roughness,
                'bottomlevel': botlev
            }
        elif css_type == 'rectangle':
            cssdct[branch[0].Index] = {
                'type': css_type,
                'height': 5.0,
                'width': round(values[values.soortparameter=='bodembreedte'].waarde.values[0], 3),
                'closed': 0,
                'thalweg': 0.0,
                'typeruwheid': values.typeruwheid.iloc[0],
                'ruwheid': roughness,        
                'bottomlevel': botlev
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

        if 'waterstand' in bndcnd.typerandvoorwaarde:
            bctype = 'waterlevel'
        elif 'debiet' in bndcnd.typerandvoorwaarde:
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
