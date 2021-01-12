import logging

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import LineString
from rasterstats import zonal_stats
from delft3dfmpy.core import checks, geometry
from delft3dfmpy.datamodels.common import ExtendedDataFrame
import rasterio
import warnings
from rasterio.transform import from_origin
import os
import imod
from tqdm.auto import tqdm
import logging

logger = logging.getLogger(__name__)
def generate_unpaved(catchments, landuse, surface_level, soiltype,  surface_storage, infiltration_capacity, initial_gwd, meteo_areas, zonalstats_alltouched=None):    
    """ 
    Combine all data to form a complete UNPAVED definition. ALso the coordinates for the networ topology are included.
    Zonal statistics are applied to land use, to get the areas per type. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed, and for soil type the dominant soil type in the catchment is used. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    
    """
    all_touched=False if zonalstats_alltouched is None else zonalstats_alltouched
               
    # required rasters
    warnings.filterwarnings('ignore')
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)   
    
    rast, affine = read_raster(soiltype, static=True)    
    soiltypes = zonal_stats(catchments, soiltype, affine = affine, stats='majority',all_touched=all_touched)
    
    rast, affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, rast, affine=affine, stats="median",all_touched=all_touched)    
        
    # optional rasters
    if isinstance(surface_storage, str):                
        rast,affine = read_raster(surface_storage, static=True)
        sstores = zonal_stats(catchments, rast, affine=affine, stats="mean",all_touched=True)            
    elif isinstance(surface_storage,int):
        surface_storage = float(surface_storage)
    if isinstance(infiltration_capacity, str):        
        rast,affine = read_raster(infiltration_capacity, static=True)
        infcaps = zonal_stats(catchments, rast, affine=affine, stats="mean",all_touched=True)    
    elif isinstance(infiltration_capacity,int):
        infiltration_capacity = float(infiltration_capacity)
    if isinstance(initial_gwd, str):                
        rast,affine = read_raster(initial_gwd, static=True)
        ini_gwds = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)       
    elif isinstance(initial_gwd,int):
        initial_gwd = float(initial_gwd)
        
    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    unpaved_drr = ExtendedDataFrame(required_columns=['code'])
    unpaved_drr.set_data( pd.DataFrame(np.zeros((len(catchments),12)), 
                                       columns=['code','total_area','lu_areas','mvlevel',
                                                'soiltype','surstor','infcap','initial_gwd','meteostat','px','py','boundary'
                                               ], dtype="str"), index_col='code')
    unpaved_drr.index = catchments.code            
    # HyDAMO Crop code; hydamo name, sobek index, sobek name: 
    # 1 aardappelen   3 potatoes
    # 2 graan         5 grain
    # 3 suikerbiet    4 sugarbeet
    # 4 mais          2 corn
    # 5 overige gew. 15 vegetables
    # 6 bloembollen  10 bulbous plants
    # 7 boomgaard     9 orchard
    # 8 gras          1 grass  
    # 9 loofbos      11 dediduous
    # 10 naaldbos    12 conferous
    # 11 natuuur     13 nature
    # 12 braak       14 fallow    
    sobek_indices = [3,5,4,2,15,10,9,1,11,12,13,14]
    for num, cat in enumerate(catchments.itertuples()):    
        # if no rasterdata could be obtained for this catchment, skip it.
        if mean_elev[num]['median'] is None: 
            logger.warning('No rasterdata available for catchment %s' % cat.code)
            continue        
        tm = [m for m in meteo_areas.itertuples() if m.geometry.contains(cat.geometry.centroid)]
        ms = meteo_areas.iloc[0,:][0] if tm==[] else tm[0].code   
        # find corresponding meteo-station
        #ms = [ms for ms in meteo_areas.itertuples() if ms.geometry.contains(cat.geometry.centroid)]
        #ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]
        mapping = np.zeros(16, dtype=int)
        for i in range(1,12):
            if i in lu_counts[num]: mapping[sobek_indices[i-1]-1] = lu_counts[num][i]*px_area            
        lu_map = ' '.join(map(str,mapping))        
        elev = mean_elev[num]['median']
        unpaved_drr.at[cat.code, 'code'] = str(cat.code)
        unpaved_drr.at[cat.code, 'total_area'] = f'{cat.geometry.area:.0f}'
        unpaved_drr.at[cat.code, 'lu_areas'] = lu_map
        unpaved_drr.at[cat.code, 'mvlevel'] = f'{elev:.2f}'
        unpaved_drr.at[cat.code, 'soiltype'] =f'{soiltypes[num]["majority"]+100.:.0f}'        
        if isinstance(surface_storage, float):    
            unpaved_drr.at[cat.code, 'surstor'] = f'{surface_storage:.3f}' 
        else:        
            unpaved_drr.at[cat.code, 'surstor'] = f'{sstores[num]["mean"]:.3f}'
        if isinstance(infiltration_capacity, float):    
            unpaved_drr.at[cat.code, 'infcap'] = f'{infiltration_capacity:.3f}' 
        else:        
            unpaved_drr.at[cat.code, 'infcap'] = f'{infcaps[num]["mean"]:.3f}'
        if isinstance(initial_gwd, float):    
            unpaved_drr.at[cat.code, 'initial_gwd'] = f'{initial_gwd:.2f}'  
        else:
            unpaved_drr.at[cat.code, 'initial_gwd'] = f'{ini_gwds[num]["mean"]:.2f}'
        unpaved_drr.at[cat.code, 'meteostat'] = ms
        unpaved_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]-10:.0f}'
        unpaved_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        unpaved_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                
    return unpaved_drr          
    
def generate_ernst(catchments, depths, resistance, infiltration_resistance, runoff_resistance):    
    """
    The lists with depths and resistances as well as the standard infiltration and runoff resistances, are converted, to a dataframe.
    """
    ernst_drr = ExtendedDataFrame(required_columns=['code'])
    ernst_drr.set_data( pd.DataFrame(np.zeros((len(catchments),5)), 
                                       columns=['code','reslist','lvs','cvi','cvs'], dtype="str"), index_col='code')
    ernst_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):            
        ernst_drr.at[cat.code, 'code'] = str(cat.code)
        ernst_drr.at[cat.code, 'reslist'] = ' '.join([str(res) for res in resistance])
        ernst_drr.at[cat.code, 'lvs'] = ' '.join([str(depth) for depth in depths])
        ernst_drr.at[cat.code, 'cvi'] = str(infiltration_resistance)
        ernst_drr.at[cat.code, 'cvs'] = str(runoff_resistance)
    return ernst_drr     

def generate_paved( catchments=None, 
                    overflows=None,
                    sewer_areas=None,                                   
                    landuse=None, 
                    surface_level=None,
                    street_storage=None,
                    sewer_storage=None,
                    pump_capacity=None, 
                    meteo_areas=None,
                    zonalstats_alltouched=None):
    """ 
        Combine all data to form a complete PAVED definition. ALso the coordinates for the networ topology are included.
    Zonal statistics are applied to land use, to get the paved area. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    """
    
    all_touched=False if zonalstats_alltouched is None else zonalstats_alltouched        
        
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)   
    sl_rast, sl_affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, sl_rast, affine=sl_affine, stats="median",all_touched=all_touched)         
    
    if isinstance( street_storage, str):                
        strs_rast, strs_affine = read_raster(street_storage, static=True)
        str_stors = zonal_stats(catchments, strs_rast, affine=strs_affine, stats="mean", all_touched=True)        
    elif isinstance( street_storage, int):
        street_storage = float(street_storage)
    if isinstance( sewer_storage, str):                
        sews_rast, sews_affine = read_raster(sewer_storage, static=True)
        sew_stors = zonal_stats(catchments, sews_rast, affine=sews_affine, stats="mean", all_touched=True)        
    elif isinstance( sewer_storage, int):
        sewer_storage = float(sewer_storage)    
    if isinstance(pump_capacity, str):     
        pump_rast, pump_affine = read_raster(pump_capacity, static=True)               
        pump_caps = zonal_stats(catchments, pump_rast, affine=pump_affine, stats="mean", all_touched=True)    
    elif isinstance( pump_capacity, int):
        pump_capacity = float(pump_capacity)    
        
    def update_dict(dict1, dict2):
        for i in dict2.keys():
            if i in dict1:
                dict1[i]+=dict2[i]
            else:
                dict1[i] = dict2[i]
        return dict1
    
    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    paved_drr = ExtendedDataFrame(required_columns=['code'])
    if sewer_areas is not None:        
        # if the parameters area rasters, do the zonal statistics per sewage area as well.
        if isinstance( street_storage, str):                
            str_stors_sa = zonal_stats(sewer_areas, strs_rast,affine=strs_affine,stats="mean", all_touched=True)
        if isinstance( sewer_storage, str):                
            sew_stors_sa = zonal_stats(sewer_areas, sews_rast,affine=sews_affine,stats="mean", all_touched=True)
        if isinstance(pump_capacity, str):     
            pump_caps_sa = zonal_stats(sewer_areas, pump_rast,affine=pump_affine,stats="mean", all_touched=True)
        mean_sa_elev = zonal_stats(sewer_areas, sl_rast, affine=sl_affine, stats="median",all_touched=True)         
        
        # initialize the array of paved nodes, which should contain a node for all catchments and all overflows
        paved_drr.set_data( pd.DataFrame(np.zeros((len(catchments)+len(overflows),10)), 
columns=['code','area','mvlevel', 'streetstor', 'sewstor', 'pumpcap','meteostat','px', 'py', 'boundary'], dtype="str"), index_col='code')
        paved_drr.index = catchments.code.append(overflows.code)
        
        # find the paved area in the sewer areas
        for isew, sew in enumerate(sewer_areas.itertuples()):            
            pav_area = 0
            for cat_ind, cat in enumerate(catchments.itertuples()):                
                # if no rasterdata could be obtained for this catchment, skip it.
                if mean_elev[cat_ind]['median'] is None:
                   logger.warning('No rasterdata available for catchment %s' % cat.code)     
                   continue
                if(cat.geometry.intersects(sew.geometry)):
                    test_intersect = cat.geometry.intersection(sew.geometry)
                    #print(cat.Index+' '+sew.Index+' '+test_intersect.type)
                    if test_intersect.type =='LineString':
                        logger.warning('Intersection in %s contains of LineStrings, not polygons. Skipping. '% cat.code)
                        continue
                    if test_intersect.type=='GeometryCollection':                                                
                        numpol = 0
                        logger.info('Intersection in %s contains a GeometryCollection - splitting into polygons.'% cat.code)
                        for int_ft in test_intersect:                            
                            if int_ft.type == 'Polygon':
                                if numpol==0:                                                                
                                    intersecting_pixels = zonal_stats(int_ft, lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)[0]
                                else:                                    
                                    temp_int = zonal_stats(int_ft, lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)[0]
                                    intersecting_pixels = update_dict(intersecting_pixels, temp_int)
                                numpol += 1                        
                    else:
                        # find the paved area within the intersection and add it to the sewer area sum
                        intersecting_pixels = zonal_stats(cat.geometry.intersection(sew.geometry), lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)[0]
                    if intersecting_pixels=={}:
                        continue
                    if 14.0 not in intersecting_pixels:
                        logger.warning('%s/%s: no paved area in sewer area intersection!' % (sew.code, cat.code))
                        continue
                        
                    pav_pixels = intersecting_pixels[14.0]
                    pav_area += pav_pixels*px_area
                    # subtract it fromthe total paved area in this catchment, make sure at least 0 remains
                    lu_counts[cat_ind][14.0] -=  pav_pixels
                    if lu_counts[cat_ind][14.0] < 0: lu_counts[cat_ind][14.0]  = 0
            
            elev = mean_sa_elev[isew]['median']
            # find overflows related to this sewer area
            ovf = overflows[overflows.codegerelateerdobject==sew.code]           
            for ov in ovf.itertuples():                
                # find corresponding meteo-station
                tm = [m for m in meteo_areas.itertuples() if m.geometry.contains(sew.geometry.centroid)]
                ms = meteo_areas.iloc[0,:][0] if tm==[] else tm[0].code                         
                #ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]                            
                # add prefix to the overflow id to create the paved-node id
                paved_drr.at[ov.code, 'code'] = str(ov.code)
                paved_drr.at[ov.code, 'area'] = str(pav_area * ov.fractie)
                paved_drr.at[ov.code, 'mvlevel'] = f'{elev:.2f}'         
                # if a float is given, a standard value is passed. If a string is given, a rastername is assumed to zonal statistics are applied.       
                if isinstance(street_storage, float):    
                    paved_drr.at[ov.code, 'streetstor'] = f'{street_storage:.2f}'
                else:
                    paved_drr.at[ov.code, 'streetstor'] = f'{str_stors_sa[isew]["mean"]:.2f}'                   
                if isinstance(sewer_storage, float):    
                    paved_drr.at[ov.code, 'sewstor'] = f'{sewer_storage:.2f}'
                else:
                    paved_drr.at[ov.code, 'sewstor'] = f'{sew_stors_sa[isew]["mean"]:.2f}'            
                if isinstance(pump_capacity, float):    
                    paved_drr.at[ov.code, 'pumpcap'] = f'{pump_capacity}'
                else:
                    paved_drr.at[ov.code, 'pumpcap'] = f'{pump_caps_sa[isew]["mean"]:.2f}'    
                paved_drr.at[ov.code,'meteostat'] = ms
                paved_drr.at[ov.code, 'px'] = f'{ov.geometry.coords[0][0]+10:.0f}'
                paved_drr.at[ov.code, 'py'] = f'{ov.geometry.coords[0][1]:.0f}'
                paved_drr.at[ov.code, 'boundary'] = ov.code              
    else:    
        # in this case only the catchments are taken into account. A node is created for every catchment nonetheless, but only nodes with a remaining area >0 are written.
        paved_drr.set_data( pd.DataFrame(np.zeros((len(catchments),10)), 
                                       columns=['code','area','mvlevel', 'streetstor', 'sewstor', 'pumpcap','meteostat', 'px', 'py', 'boundary'], dtype="str"), index_col='code')
        paved_drr.index = catchments.code
    
    for num, cat in enumerate(catchments.itertuples()):                    
        # if no rasterdata could be obtained for this catchment, skip it.
        if mean_elev[num]['median'] is None:
            logger.warning('No rasterdata available for catchment %s' % cat.code)     
            continue
        
        # find corresponding meteo-station
        tm = [m for m in meteo_areas.itertuples() if m.geometry.contains(cat.geometry.centroid)]
        ms = meteo_areas.iloc[0,:][0] if tm==[] else tm[0].code        
            
        elev = mean_elev[num]['median']
        paved_drr.at[cat.code, 'code'] = str(cat.code)
        paved_drr.at[cat.code, 'area'] = str(lu_counts[num][14]*px_area) if 14 in lu_counts[num] else '0'        
        paved_drr.at[cat.code, 'mvlevel'] = f'{elev:.2f}'         
        # if a float is given, a standard value is passed. If a string is given, a rastername is assumed to zonal statistics are applied.       
        if isinstance(street_storage, float):    
            paved_drr.at[cat.code, 'streetstor'] = f'{street_storage:.2f}'
        else:
            paved_drr.at[cat.code, 'streetstor'] = f'{str_stors[num]["mean"]:.2f}'                   
        if isinstance(sewer_storage, float):    
            paved_drr.at[cat.code, 'sewstor'] = f'{sewer_storage:.2f}'
        else:
            paved_drr.at[cat.code, 'sewstor'] = f'{sew_stors[num]["mean"]:.2f}'            
        if isinstance(pump_capacity, float):    
            paved_drr.at[cat.code, 'pumpcap'] = f'{pump_capacity}'
        else:
            paved_drr.at[cat.code, 'pumpcap'] = f'{pump_caps[num]["mean"]:.2f}'    
        paved_drr.at[cat.code,'meteostat'] = ms
        paved_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]+10:.0f}'
        paved_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        paved_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return paved_drr   
       
def generate_greenhouse(catchments, landuse, surface_level, roof_storage, meteo_areas, zonalstats_alltouched=None):    
    """ 
        Combine all data to form a complete GREENHSE defintion. ALso the coordinates for the network topology are included.
    Zonal statistics are applied to land use, to get the paved area. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    """
    all_touched=False if zonalstats_alltouched is None else zonalstats_alltouched
        
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=all_touched)   
    rast, affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, rast, affine=affine, stats="median", all_touched=all_touched) 
    # optional rasters
    if isinstance(roof_storage, str):                
        rast, affine = read_raster(roof_storage, static=True)
        roofstors = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)    
    elif isinstance(roof_storage, int):
        roof_storage = float(roof_storage)
        
    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    gh_drr = ExtendedDataFrame(required_columns=['code'])
    gh_drr.set_data( pd.DataFrame(np.zeros((len(catchments),8)), 
                                       columns=['code','area','mvlevel', 'roofstor', 'meteostat','px', 'py', 'boundary'], dtype="str"), index_col='code')
    gh_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):    
        # if no rasterdata could be obtained for this catchment, skip it.
        if mean_elev[num]['median'] is None:
            logger.warning('No rasterdata available for catchment %s' % cat.code)     
            continue
                
        # find corresponding meteo-station
        tm = [m for m in meteo_areas.itertuples() if m.geometry.contains(cat.geometry.centroid)]
        ms = meteo_areas.iloc[0,:][0] if tm==[] else tm[0].code   
            
        elev = mean_elev[num]['median']
        gh_drr.at[cat.code, 'code'] = str(cat.code)
        gh_drr.at[cat.code, 'area'] = str(lu_counts[num][15]*px_area) if 15 in lu_counts[num] else '0'        
        gh_drr.at[cat.code, 'mvlevel'] = f'{elev:.2f}' 
        if isinstance(roof_storage, float):    
            gh_drr.at[cat.code, 'roofstor'] = f'{roof_storage:.2f}'
        else:
            gh_drr.at[cat.code, 'roofstor'] = f'{roofstors[num]["mean"]:.2f}'    
        gh_drr.at[cat.code, 'meteostat'] = ms
        gh_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]+20:.0f}'
        gh_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        gh_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return gh_drr   

def generate_openwater(catchments, landuse, meteo_areas, zonalstats_alltouched=None):    
    """ 
        Combine all data to form a complete OPENWATE definotion. ALso the coordinates for the network topology are included.
    Zonal statistics are applied to land use, to get the open water area. The classificition described in the notebook is assumed. 
    """
    all_touched=False if zonalstats_alltouched is None else zonalstats_alltouched
    
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True, all_touched=all_touched)   

    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    ow_drr = ExtendedDataFrame(required_columns=['code'])
    ow_drr.set_data( pd.DataFrame(np.zeros((len(catchments),6)), 
                                       columns=['code','area','meteostat','px', 'py', 'boundary'], dtype="str"), index_col='code')
    ow_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):    
        # find corresponding meteo-station
        tm = [m for m in meteo_areas.itertuples() if m.geometry.contains(cat.geometry.centroid)]
        ms = meteo_areas.iloc[0,:][0] if tm==[] else tm[0].code   
        
        ow_drr.at[cat.code, 'code'] = str(cat.code)
        ow_drr.at[cat.code, 'area'] = str(lu_counts[num][13]*px_area) if 13 in lu_counts[num] else '0'        
        ow_drr.at[cat.code, 'meteostat'] = ms
        ow_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]-20:.0f}'
        ow_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        ow_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return ow_drr   

def generate_boundary(boundary_nodes, catchments, overflows=None):
    """
    Method to create boundary  nodes for RR.

    """
    if overflows is not None:
        numlats = len(catchments)+len(overflows)
    else:
        numlats = len(catchments)
    bnd_drr = ExtendedDataFrame(required_columns=['code'])
    bnd_drr.set_data( pd.DataFrame(np.zeros((numlats,3)), 
                                       columns=['code', 'px', 'py'], dtype="str"), index_col='code')                                        
    if overflows is not None:
        bnd_drr.index = catchments.code.append(overflows.code)
    else:
        bnd_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):    
        # print(num, cat.code)
        if boundary_nodes[boundary_nodes['code']==cat.lateraleknoopcode].empty:
            #raise IndexError(f'{cat.code} not connected to a boundary node. Skipping.')
            logger.warning('%s not connected to a boundary node. Skipping.' % cat.code)
            continue
        bnd_drr.at[cat.code, 'code'] = cat.lateraleknoopcode        
        bnd_drr.at[cat.code, 'px']  = str(boundary_nodes[boundary_nodes['code']==cat.lateraleknoopcode]['geometry'].x.iloc[0]).strip()                                                                 
        bnd_drr.at[cat.code, 'py']  = str(boundary_nodes[boundary_nodes['code']==cat.lateraleknoopcode]['geometry'].y.iloc[0]).strip()
    if overflows is not None:
        logger.info('Adding overflows to the boundary nodes.')
        for num, ovf in enumerate(overflows.itertuples()):
            bnd_drr.at[ovf.code, 'code'] = ovf.code
            bnd_drr.at[ovf.code, 'px'] = str(ovf.geometry.coords[0][0])
            bnd_drr.at[ovf.code, 'py'] = str(ovf.geometry.coords[0][1])       
    return bnd_drr    

def generate_seepage(catchments, seepage_folder):    
    """
    Method to obtain catchment-average seepage fluxes from rasters. The time step is deduced from the raster filenames. 
    
    We assume seepage is read from Metaswap (m3 per cell). It needs to be converted to mm/day.
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(seepage_folder)
    times = []    
    arr = np.zeros((len(file_list), len(catchments.code)))
    for ifile, file in tqdm(enumerate(file_list),total=len(file_list),desc='Reading seepage files'):
        array, affine, time = read_raster(os.path.join(seepage_folder, file))
        times.append(time)           
        stats = zonal_stats(catchments, array, affine=affine, stats="mean", all_touched=True)
        arr[ifile,:] = [s['mean'] for s in stats]        
    result = pd.DataFrame(arr,columns='sep_'+catchments.code)
    result.index = times 
    # convert units
    result_mmd = (result / (1e-3*(affine[0]*-affine[4])))/((times[2]-times[1]).total_seconds()/86400.)
    return result_mmd

def generate_precip(areas, precip_folder):
    """
    Method to obtain catchment-average seepage fluxes from rasters. The time step is deduced from the raster filenames.
    
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(precip_folder)
    times = []        
    arr = np.zeros((len(file_list), len(areas.code)))
    for ifile, file in tqdm(enumerate(file_list),total=len(file_list),desc='Reading precipitation files'):
        array, affine, time = read_raster(os.path.join(precip_folder, file))
        times.append(time)                   
        stats = zonal_stats(areas, array, affine=affine, stats="mean", all_touched=True)
        arr[ifile,:]= [s['mean'] for s in stats]        
    result = pd.DataFrame(arr, columns='ms_'+areas.code)
    result.index = times 
    return result

def generate_evap(areas, evap_folder):    
    """
    Method to obtain catchment-average evaporation fluxes from rasters. The time step is deduced from the raster filenames. Since only one timeeries is allowed, the meteo areas are dissolved to a user specifield field.
    
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(evap_folder)
    # aggregated evap
    areas['dissolve'] = 1
    agg_areas = areas.iloc[0:len(areas),:].dissolve(by='dissolve',aggfunc='mean')
    times = []    
    arr = np.zeros((len(file_list), 1))
    for ifile, file in tqdm(enumerate(file_list),total=len(file_list),desc='Reading evaporation files'):
        array, affine, time = read_raster(os.path.join(evap_folder, file))
        times.append(time)                  
        stats = zonal_stats(agg_areas, array, affine=affine, stats="mean",all_touched=True)
        arr[ifile,:] = [s['mean'] for s in stats]       
    result = pd.DataFrame(arr,columns=['ms_'+areas.iloc[0,0]])
    result.index = times 
    return result

def read_raster(file, static=False):
    """
    Method to read a raster. All rasterio types are accepted, plus IDF: in that case the iMod-package is used to read the IDF raster (IDF is cusomary for MODFLOW/SIMGRO models.)

    Parameters
    ----------
    file : raster
        
    static : BOOL, optional
        If static than no time information needs to be deduced.

    Returns
    -------
    rasterio grid and an affine object.        

    """
    if not static:         
        time = pd.Timestamp(os.path.split(file)[1].split('_')[1].split('.')[0])
    if file.lower().endswith('idf'):        
        dataset = imod.idf.open(file)
        header = imod.idf.header(file,pattern=None)        
        grid = dataset[0,0,:,:].values        
        affine =from_origin(header['xmin'], header['ymax'], header['dx'], header['dx'])                            
    else:        
        dataset = rasterio.open(file)
        affine = dataset.transform
        grid = dataset.read(1)
    if static:
        return grid, affine
    else:
        return grid, affine, time
        
