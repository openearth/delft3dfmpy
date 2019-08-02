import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from scipy.spatial import KDTree

from delft3dfmpy.converters import hydamo_to_dflowfm
from delft3dfmpy.core import checks

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
                direction=1,
                nrstages=1,
                capacity=pump.maximalecapaciteit,
                startlevelsuctionside=pump.startlevelsuctionside,
                stoplevelsuctionside=pump.stoplevelsuctionside,
                locationfile=f'pump_{pump.code}.pli'
            )

    def weirs_from_hydamo(self, weirs):
        """
        Method to convert dflowfm weirs from hydamo weirs.
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_weirs(weirs)
        # Add to dict
        for weir in geconverteerd.itertuples():
            self.structures.add_weir(
                id=weir.code,
                branchid=weir.branch_id,
                chainage=weir.branch_offset,
                crestlevel=weir.laagstedoorstroomhoogte,
                crestwidth=weir.kruinbreedte,
                dischargecoeff=weir.afvoercoefficient
            )

    def culverts_from_hydamo(self, culverts):
        """
        Method to convert dflowfm weirs from hydamo weirs.
        """
        # Convert to dflowfm input
        geconverteerd = hydamo_to_dflowfm.generate_culverts(culverts)

        # Add to dict
        for culvert in geconverteerd.itertuples():
            self.structures.add_culvert(
                id=culvert.code,
    	        branchid=culvert.branch_id,
    	        chainage=culvert.branch_offset,
    	        leftlevel=culvert.hoogtebinnenonderkantbovenstrooms,
    	        rightlevel=culvert.hoogtebinnenonderkantbenedenstrooms,
    	        crosssection=culvert.crosssection,
    	        length=culvert.geometry.length,
    	        inletlosscoeff=culvert.intreeverlies,
    	        outletlosscoeff=culvert.uittreeverlies,
            )

class CrossSectionsIO:

    def __init__(self, crosssections):
        self.crosssections = crosssections

    def from_hydamo(self, dwarsprofielen, parameterised=None):
        """
        Method to add cross section from hydamo files. Two files and a default cross section
        can be handed to the function, the cross section file (dwarsprofiel), the
        parameterised file (normgeparametriseerd) and a standard cross section. The
        hierarchical order is 1. dwarsprofiel, 2. normgeparametriseerd, 3. standard.
        Each branch will be assigned a profile following this order. If parameterised
        and standard are not given, branches can be without cross section. In that case
        dflowfm will assign a standard profile.
        """
        nnocross = len(self.crosssections.get_branches_without_crosssection())
        logger.info(f'Before adding the number of branches without cross section is: {nnocross}.')

        # 1. Collect cross sections from 'dwarsprofielen'
        crosssections = hydamo_to_dflowfm.dwarsprofiel_to_yzprofiles(dwarsprofielen)
        for name, css in crosssections.items():
            # Add definition
            self.crosssections.add_yz_definition(yz=css['yz'], name=name, roughnesstype=css['ruwheidstypecode'], roughnessvalue=css['ruwheidswaarde'])
            # Add location
            self.crosssections.add_crosssection_location(branchid=css['branchid'], chainage=css['chainage'], definition=name)
        
        # Check the number of branches with cross sections
        no_crosssection = self.crosssections.get_branches_without_crosssection()
        nnocross = len(no_crosssection)
        logger.info(f'After adding \'dwarsprofielen\' the number of branches without cross section is: {nnocross}.')
        if (nnocross == 0) or parameterised is None:
            return None

        # Derive norm cross sections for norm parameterised
        crosssections = hydamo_to_dflowfm.parameterised_to_profiles(parameterised, no_crosssection)
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
                    roughnesstype=branchdata.at[branchid, 'ruwheidstypecode'],
                    roughnessvalue=branchdata.at[branchid, 'ruwheidswaarde']
                )

            if css['type'] == 'trapezium':
                name = self.crosssections.add_trapezium_definition(
                    slope=css['slope'],
                    maximumflowwidth=css['maximumflowwidth'],
                    bottomwidth=css['bottomwidth'],
                    closed=css['closed'],
                    roughnesstype=branchdata.at[branchid, 'ruwheidstypecode'],
                    roughnessvalue=branchdata.at[branchid, 'ruwheidswaarde']
                ) 

            # Add location
            self.crosssections.add_crosssection_location(branchid=branchid, chainage=chainage, definition=name)

        nnocross = len(self.crosssections.get_branches_without_crosssection())
        logger.info(f'After adding \'normgeparameteriseerd\' the number of branches without cross section is: {nnocross}.')

class ExternalForcingsIO:

    def __init__(self, external_forcings):
        self.external_forcings = external_forcings

    def from_hydamo(self, boundary_conditions):

        # Read from Hydamo
        bcdct = hydamo_to_dflowfm.generate_boundary_conditions(boundary_conditions, self.external_forcings.dflowfmmodel.network.schematised)

        # Add all items
        for key, item in bcdct.items():

            # Add to dataframe
            self.external_forcings.boundaries.loc[key] = item
            
            # Check if a 1d2d link should be removed
            self.external_forcings.dflowfmmodel.network.links1d2d.check_boundary_link(self.external_forcings.boundaries.loc[key])

    def read_laterals(self, locations, lateral_discharges):
        """
        Process laterals

        Parameters
        ----------
        locations: gpd.GeoDataFrame
            GeoDataFrame with at least 'geometry' (Point) and the column 'code'
        lateral_discharges: pd.DataFrame
            DataFrame with lateral discharges. The index should be a time object (datetime or similar).
        """

        # Check argument
        checks.check_argument(locations, 'locations', gpd.GeoDataFrame, columns=['geometry'])
        checks.check_argument(lateral_discharges, 'lateral_discharges', pd.DataFrame)

        # Check if network has been loaded
        network1d = self.external_forcings.dflowfmmodel.network.mesh1d
        if not network1d.meshgeomdim.numnode:
            raise ValueError('1d network has not been generated or loaded. Do this before adding laterals.')

        # Find nearest 1d node per location
        nodes1d = network1d.get_nodes()
        get_nearest = KDTree(nodes1d)
        lateral_crds = np.vstack([loc.geometry.coords[0] for loc in locations.itertuples()])
        
        _, nearest_idx = get_nearest.query(lateral_crds)
        
        # Get time series and add to dictionary
        for crd, lateral in zip(nodes1d[nearest_idx], locations.itertuples()):
            # Check if a time is provided for the lateral
            if lateral.code not in lateral_discharges.columns:
                logger.warning(f'No data found for {lateral.code}. Skipping.')
                continue

            # Get timeseries
            series = lateral_discharges.loc[:, lateral.code]

            # Add to dictionary
            self.external_forcings.laterals[lateral.code] = {
                'x': crd[0],
                'y': crd[1],
                'timeseries': series
            }
