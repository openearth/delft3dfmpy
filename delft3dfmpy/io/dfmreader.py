import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from delft3dfmpy.converters import hydamo_to_dflowfm
from delft3dfmpy.core import checks
from delft3dfmpy.datamodels.common import ExtendedGeoDataFrame
from shapely.geometry import LineString, MultiPolygon, Point, Polygon

logger = logging.getLogger(__name__)


class StructuresIO:

    def __init__(self, structures):
        self.structures = structures

    def pumps_from_hydamo(self, pompen, sturing, gemalen):
        """
        Method to generate dflowfm pumps from hydamo. Three dataframes are needed:
        pompen (pumps), sturing (control) and gemalen (pump houses).
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_pumps(pompen, sturing, gemalen)
        # Add to dict
        for pump in geconverteerd.itertuples():
            self.structures.add_pump(
                id=pump.code,
                branchid=pump.branch_id,
                chainage=pump.branch_offset,
                orientation='positive',
                numstages=1,
                controlside='suctionSide',
                capacity=pump.maximalecapaciteit,
                startlevelsuctionside=pump.startlevelsuctionside,
                stoplevelsuctionside=pump.stoplevelsuctionside                
            )

    def orifices_from_hydamo(self, orifices):
        """
        Method to generate dflowfm orifices.
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_orifices(orifices)
        # Add to dict
        for orifice in geconverteerd.itertuples():
            self.structures.add_orifice(
                id=orifice.code,
                branchid=orifice.branch_id,
                chainage=orifice.branch_offset,                
                crestlevel=orifice.laagstedoorstroomhoogte,
                crestwidth=orifice.laagstedoorstroombreedte,
                gateloweredgelevel=orifice.schuifhoogte,
                uselimitflowpos=orifice.uselimitflow,
                limitflowpos=orifice.limitflow,
                uselimitflowneg=orifice.uselimitflow,
                limitflowneg=orifice.limitflow,
                corrcoeff=orifice.afvoercoefficient                               
            )
            
    def weirs_from_hydamo(self, weirs, yz_profiles=None, parametrised_profiles=None, afsluitmiddel=None, sturing=None):
        """
        Method to convert dflowfm weirs from hydamo weirs.
        
        Split the HyDAMO weirs object into a subset with 'regular weirs' and 'universal weirs'. 
       
        For now the assumption is that weirs for which a profile is defined, are universal, other weirs are regular.
               
        """
        # Convert to dflowfm input
        
        regular = ExtendedGeoDataFrame(geotype=Point, required_columns=weirs.required_columns)
        index = np.zeros((len(weirs.code)))
        if yz_profiles is not None:
            if 'codegerelateerdobject' in yz_profiles:            
                index[np.isin(weirs.code , np.asarray(yz_profiles.codegerelateerdobject))]=1
        if parametrised_profiles is not None:
            if 'codegerelateerdobject' in parametrised_profiles:            
                index[np.isin(weirs.code , np.asarray(parametrised_profiles.codegerelateerdobject))]=1
           
        regular.set_data(weirs[index==0], index_col='code',check_columns=True)        
        regular_geconverteerd = hydamo_to_dflowfm.generate_weirs(regular, afsluitmiddel=afsluitmiddel)        
        regular.set_data(weirs, index_col='code',check_columns=True)        
        
        # Add to dict
        for weir in regular_geconverteerd.itertuples():
            self.structures.add_weir(
                id=weir.code,
                branchid=weir.branch_id,
                chainage=weir.branch_offset,                
                crestlevel=weir.laagstedoorstroomhoogte,
                crestwidth=weir.laagstedoorstroombreedte,
                corrcoeff=weir.afvoercoefficient
            )
    
        universal = ExtendedGeoDataFrame(geotype=Point, required_columns=weirs.required_columns)
        universal.set_data(weirs[index==1], index_col='code',check_columns=True)                                      
        universal_geconverteerd = hydamo_to_dflowfm.generate_uweirs(universal,yz_profiles=yz_profiles, parametrised_profiles=parametrised_profiles)
        
        if universal_geconverteerd.empty:
            return("No profile detected for universal weir)")
        # Add to dict
        for uweir in universal_geconverteerd.itertuples():
            self.structures.add_uweir(
                id=uweir.code,
                branchid=uweir.branch_id,
                chainage=uweir.branch_offset,                
                crestlevel=uweir.laagstedoorstroomhoogte,                
                numlevels=uweir.numlevels,
                yvalues=uweir.yvalues,
                zvalues=uweir.zvalues,
                allowedflowdir='both',
                dischargecoeff=uweir.afvoercoefficient                
            )    
        
    def bridges_from_hydamo(self, bridges, yz_profiles=None, parametrised_profiles=None):
        """
        Method to convert dflowfm bridges from hydamo bridges.
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_bridges(bridges, yz_profiles=yz_profiles, parametrised_profiles=parametrised_profiles)
        # Add to dict
        for bridge in geconverteerd.itertuples():
            self.structures.add_bridge(
                id=bridge.code,
                branchid=bridge.branch_id,
                chainage=bridge.branch_offset,
                length=bridge.lengte,
                bedlevel=bridge.bedlevel,
                upperheight=bridge.hoogtebovenzijde,
                lowerheight=bridge.hoogteonderzijde,                
                crosssection=bridge.crosssection,
                inletlosscoeff=bridge.intreeverlies,
                outletlosscoeff=bridge.uittreeverlies,
                frictiontype=hydamo_to_dflowfm.roughness_gml[bridge.ruwheidstypecode],
                frictionvalue=bridge.ruwheidswaarde
            )
            

    def culverts_from_hydamo(self, culverts, afsluitmiddel=None):
        """
        Method to convert dflowfm weirs from hydamo weirs.
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_culverts(culverts, afsluitmiddel)

        # Add to dict
        for culvert in geconverteerd.itertuples():            
                self.structures.add_culvert(
                    id=culvert.code,
        	        branchid=culvert.branch_id,
        	        chainage=culvert.branch_offset,
        	        leftlevel=culvert.hoogtebinnenonderkantbovenstrooms,
        	        rightlevel=culvert.hoogtebinnenonderkantbenedenstrooms,
        	        crosssection=culvert.crosssection,
        	        length=culvert.lengte,#geometry.length,
        	        inletlosscoeff=culvert.intreeverlies,
        	        outletlosscoeff=culvert.uittreeverlies,
                    allowedflowdir=culvert.allowedflowdir,
                    valveonoff=culvert.valveonoff,
                    numlosscoeff=culvert.numlosscoeff, 
                    valveopeningheight=culvert.valveopeningheight,
                    relopening=culvert.relopening,
                    losscoeff=culvert.losscoeff,                                    
                    frictiontype=hydamo_to_dflowfm.roughness_gml[culvert.ruwheidstypecode],
                    frictionvalue=culvert.ruwheidswaarde
                )
        
    def compound_structures(self, idlist, structurelist):
        """
        Method to add compound structures to the model.
        
        """
        geconverteerd = hydamo_to_dflowfm.generate_compounds(idlist, structurelist, self.structures)
        
         # Add to dict
        for compound in geconverteerd.itertuples():
            self.structures.add_compound(
                id=compound.code,
                numstructures=compound.numstructures,
    	        structurelist=compound.structurelist    	        
            )
        
                    

class CrossSectionsIO:

    def __init__(self, crosssections):
        self.crosssections = crosssections

    def from_hydamo(self, dwarsprofielen, parametrised=None, branches=None):
        """
        Method to add cross section from hydamo files. Two files
        can be handed to the function, the cross section file (dwarsprofiel) and the
        parametrised file (normgeparametriseerd). The
        hierarchical order is 1. dwarsprofiel, 2. normgeparametriseerd.
        Each branch will be assigned a profile following this order. If parametrised
        and standard are not given, branches can be without cross section. In that case
        a standard profile should be assigned
        """

        # first, make a selection as to use only the dwarsprofielen/parametrised that are related to branches, not structures
        if dwarsprofielen is not None and not dwarsprofielen.empty:
            if 'codegerelateerdobject' not in dwarsprofielen:
                dwarsprofielen['codegerelateerdobject'] = np.empty((len(dwarsprofielen)))*np.nan
            dp_branches = ExtendedGeoDataFrame(geotype=LineString, columns = dwarsprofielen.required_columns+['codegerelateerdobject'])
            dp_branches.set_data(dwarsprofielen[dwarsprofielen.codegerelateerdobject.isna()], index_col='code', check_columns=True)
        else:     
            dp_branches = ExtendedGeoDataFrame(geotype=LineString)
            
        if parametrised is not None and not parametrised.empty:
            if 'codegerelateerdobject' not in parametrised:
                parametrised['codegerelateerdobject'] = np.empty((len(parametrised)))*np.nan
            if len(parametrised)>0:        
                par_branches = ExtendedGeoDataFrame(geotype=LineString, columns = parametrised.required_columns+['codegerelateerdobject'])
                par_branches.set_data(parametrised[parametrised.codegerelateerdobject.isna()], index_col='code', check_columns=True)
            else:
                par_branches = ExtendedGeoDataFrame(geotype=LineString, columns = parametrised.required_columns+['codegerelateerdobject'])
            # Assign cross-sections to branches
            nnocross = len(self.crosssections.get_branches_without_crosssection())
            logger.info(f'Before adding the number of branches without cross section is: {nnocross}.')
        else:
            par_branches = ExtendedGeoDataFrame(geotype=LineString)
            
        if not dp_branches.empty:
            # 1. Collect cross sections from 'dwarsprofielen'
            crosssections = hydamo_to_dflowfm.dwarsprofiel_to_yzprofiles(dp_branches, branches)               
        
            for name, css in crosssections.items():
                # Add definition
                self.crosssections.add_yz_definition(yz=css['yz'], thalweg=css['thalweg'], name=name, roughnesstype=css['ruwheidstypecode'], roughnessvalue=css['ruwheidswaarde'])
                # Add location
                self.crosssections.add_crosssection_location(branchid=css['branchid'], chainage=css['chainage'], definition=name)
        
        # Check the number of branches with cross sections
        no_crosssection = self.crosssections.get_branches_without_crosssection()
        nnocross = len(no_crosssection)
        logger.info(f'After adding \'dwarsprofielen\' the number of branches without cross section is: {nnocross}.')
        if (nnocross == 0):
            print('No further branches without a profile.')
        elif par_branches.empty:
            print('No parametrised crossections available for branches.')
        else: 
            # Derive norm cross sections for norm parametrised
            crosssections = hydamo_to_dflowfm.parametrised_to_profiles(par_branches, no_crosssection)
            # Get branch information
            branchdata = self.crosssections.dflowfmmodel.network.branches.loc[list(crosssections.keys())]
            branchdata['chainage'] = branchdata.length / 2.
        
            # Add cross sections
            for branchid, css in crosssections.items():
                chainage = branchdata.at[branchid, 'chainage']
                    
                if css['type'] == 'rectangle':
                    name = self.crosssections.add_rectangle_definition(
                        height=css['height'],
                        width=css['width'],
                        closed=css['closed'],
                        roughnesstype=css['ruwheidstypecode'],
                        roughnessvalue=css['ruwheidswaarde']
                    )
    
                if css['type'] == 'trapezium':
                    name = self.crosssections.add_trapezium_definition(
                        slope=css['slope'],
                        maximumflowwidth=css['maximumflowwidth'],
                        bottomwidth=css['bottomwidth'],
                        closed=css['closed'],
                        roughnesstype=css['ruwheidstypecode'],
                        roughnessvalue=css['ruwheidswaarde']
                    ) 
    
                # Add location
                self.crosssections.add_crosssection_location(branchid=branchid, chainage=chainage, definition=name, shift=css['bottomlevel'])

            nnocross = len(self.crosssections.get_branches_without_crosssection())
            logger.info(f'After adding \'normgeparametriseerd\' the number of branches without cross section is: {nnocross}.')
        
        # Assign cross-sections to structures        
        nnocross = len(self.crosssections.get_structures_without_crosssection())        
        logger.info(f'Before adding the number of structures without cross section is: {nnocross}.')    
        
        if dwarsprofielen is not None:
            # and subsets for structures    
            dp_structures = ExtendedGeoDataFrame(geotype=LineString)
            dp_structures.set_data(dwarsprofielen[~dwarsprofielen.codegerelateerdobject.isna()], index_col='code', check_columns=True)
        
            # 1. Collect cross sections from 'dwarsprofielen'
            crosssections = hydamo_to_dflowfm.dwarsprofiel_to_yzprofiles(dp_structures, None)               
        
            for name, css in crosssections.items():        
                self.crosssections.add_yz_definition(yz=css['yz'], thalweg=css['thalweg'], name=name, roughnesstype=css['ruwheidstypecode'], roughnessvalue=css['ruwheidswaarde'])
        
        if parametrised is not None:
            par_structures = ExtendedGeoDataFrame(geotype=LineString, columns = parametrised.required_columns + ['codegerelateerdobject'])
            par_structures.set_data(parametrised[~parametrised.codegerelateerdobject.isna()], index_col='code', check_columns=True)            
        else:
            par_structures = ExtendedGeoDataFrame(geotype=LineString)
      
        no_crosssection = self.crosssections.get_structures_without_crosssection()        
        nnocross = len(no_crosssection)
        logger.info(f'After adding \'dwarsprofielen\' the number of branches without cross section is: {nnocross}.')
        if (nnocross == 0):
            print('No further structures without a profile.')
        elif par_structures.empty:
            print('No parametrised crosssections available for structures.')
        else:                 
            # Derive norm cross sections for norm parametrised
            crosssections = hydamo_to_dflowfm.parametrised_to_profiles(par_structures,[])                
            # Add cross section definitions
            for strucid, css in crosssections.items():             
                        
                    if css['type'] == 'rectangle':
                        name = self.crosssections.add_rectangle_definition(
                            height=css['height'],
                            width=css['width'],
                            closed=css['closed'],
                            roughnesstype=css['ruwheidstypecode'],
                            roughnessvalue=css['ruwheidswaarde'],
                            name=strucid,                        
                        )
        
                    if css['type'] == 'trapezium':
                        name = self.crosssections.add_trapezium_definition(
                            slope=css['slope'],
                            maximumflowwidth=css['maximumflowwidth'],
                            bottomwidth=css['bottomwidth'],
                            closed=css['closed'],
                            roughnesstype=css['ruwheidstypecode'],
                            roughnessvalue=css['ruwheidswaarde'],
                            name=strucid
                        ) 
                        
        nnocross = len(self.crosssections.get_structures_without_crosssection())
        logger.info(f'After adding \'normgeparametriseerd\' the number of structures without cross section is: {nnocross}.')
            
class ExternalForcingsIO:

    def __init__(self, external_forcings):
        self.external_forcings = external_forcings

    def from_hydamo(self, boundary_conditions):

        # Read from Hydamo
        bcdct = hydamo_to_dflowfm.generate_boundary_conditions(boundary_conditions, self.external_forcings.dflowfmmodel.network.schematised)

        # Add all items
        for key, item in bcdct.items():
            self.external_forcings.add_boundary_condition(key, item['geometry'], item['bctype'], item['value'], branchid=item['branchid'])
            # # Add to dataframe         
            # self.external_forcings.boundaries.loc[key] = item
            # Check if a 1d2d link should be removed
            #self.external_forcings.dflowfmmodel.network.links1d2d.check_boundary_link(self.external_forcings.boundaries.loc[key])

    def read_laterals(self, locations, lateral_discharges=None, rr_boundaries=None):
        """
        Process laterals

        Parameters
        ----------
        locations: gpd.GeoDataFrame
            GeoDataFrame with at least 'geometry' (Point) and the column 'code'
        lateral_discharges: pd.DataFrame
            DataFrame with lateral discharges. The index should be a time object (datetime or similar).
        rr_boundaries: pd.DataFrame
            DataFrame with RR-catchments that are coupled 
        """

        if rr_boundaries is None: rr_boundaries = []
        # Check argument
        checks.check_argument(locations, 'locations', gpd.GeoDataFrame, columns=['geometry'])
        if lateral_discharges is not None:
            checks.check_argument(lateral_discharges, 'lateral_discharges', pd.DataFrame)

        # Check if network has been loaded
        network1d = self.external_forcings.dflowfmmodel.network.mesh1d
        if not network1d.meshgeomdim.numnode:
            raise ValueError('1d network has not been generated or loaded. Do this before adding laterals.')

        # in case of 3d points, remove the 3rd dimension
        locations['geometry2'] = [Point([point.geometry.x, point.geometry.y]) for _,point in locations.iterrows()]    
        locations.drop('geometry', inplace=True, axis=1)
        locations.rename(columns={'geometry2':'geometry'}, inplace=True)
        
        # Find nearest 1d node per location and find the nodeid
        #lateral_crds = np.vstack([loc.geometry.coords[0] for loc in locations.itertuples()])             
        #nodes1d = network1d.get_nodes()
        #get_nearest = KDTree(nodes1d)
        #_, nearest_idx = get_nearest.query(lateral_crds[:,0:2])
                
        # Get time series and add to dictionary
        #for nidx, lateral in zip(nearest_idx, locations.itertuples()):
        for lateral in locations.itertuples():
            # crd = nodes1d[nearest_idx]
            #nid = f'{nodes1d[nidx][0]:g}_{nodes1d[nidx][1]:g}'
            
            # Check if a time is provided for the lateral
            if lateral.code in rr_boundaries:                
                # Add to dictionary
                self.external_forcings.laterals[lateral.code] = {
                    'branchid': lateral.branch_id,
                    'branch_offset': str(lateral.branch_offset)                    
                }
            else:
                if lateral_discharges is None:
                    logger.warning(f'No lateral_discharges provied. {lateral.code} expects them. Skipping.')
                    continue
                else:
                    if lateral.code not in lateral_discharges.columns:
                        logger.warning(f'No data found for {lateral.code}. Skipping.')
                        continue
                    
                # Get timeseries
                series = lateral_discharges.loc[:, lateral.code]
                
                # Add to dictionary
                self.external_forcings.laterals[lateral.code] = { 
                    'branchid': lateral.branch_id,
                    'branch_offset': str(lateral.branch_offset), 
                    'timeseries': series            
                }                                   
        
        
