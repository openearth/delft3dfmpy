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

logger = logging.getLogger(__name__)
def generate_unpaved(catchments, landuse, surface_level, soiltype,  surface_storage, infiltration_capacity, initial_gwd, meteo_areas):    
    """ 
    Combine all data to form a complete UNPAVED definition. ALso the coordinates for the networ topology are included.
    Zonal statistics are applied to land use, to get the areas per type. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed, and for soil type the dominant soil type in the catchment is used. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    
    """
    # required rasters
    warnings.filterwarnings('ignore')
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True)   
    
    rast, affine = read_raster(soiltype, static=True)    
    soiltypes = zonal_stats(catchments, soiltype, affine = affine, stats='majority',all_touched=True)
    
    rast, affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, rast, affine=affine, stats="median",all_touched=True)    
        
    # optional rasters
    if not isinstance(surface_storage, float):                
        rast,affine = read_raster(surface_storage, static=True)
        sstores = zonal_stats(catchments, rast, affine=affine, stats="mean",all_touched=True)    
    if not isinstance(infiltration_capacity, float):        
        rast,affine = read_raster(infiltration_capacity, static=True)
        infcaps = zonal_stats(catchments, rast, affine=affine, stats="mean",all_touched=True)    
    if not isinstance(initial_gwd, float):                
        rast,affine = read_raster(initial_gwd, static=True)
        ini_gwds = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)       
    
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
        
        # find corresponding meteo-station
        ms = [ms for ms in meteo_areas.itertuples() if ms.geometry.contains(cat.geometry.centroid)]
        ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]
        mapping = np.zeros(16, dtype=int)
        for i in range(1,12):
            if i in lu_counts[num]: mapping[sobek_indices[i-1]-1] = lu_counts[num][i]*px_area            
        lu_map = ' '.join(map(str,mapping))        
        elev = mean_elev[num]['median']/100
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
        unpaved_drr.at[cat.code, 'meteostat'] = ms[0]
        unpaved_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]-10:.0f}'
        unpaved_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        unpaved_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                
    return unpaved_drr          
    
def generate_ernst(catchments, depths, resistance):    
    """
    The lists with depths and resistances are converted to a dataframe.
    """
    ernst_drr = ExtendedDataFrame(required_columns=['code'])
    ernst_drr.set_data( pd.DataFrame(np.zeros((len(catchments),3)), 
                                       columns=['code','reslist','lvs'], dtype="str"), index_col='code')
    ernst_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):            
        ernst_drr.at[cat.code, 'code'] = str(cat.code)
        ernst_drr.at[cat.code, 'reslist'] = ' '.join([str(res) for res in resistance])
        ernst_drr.at[cat.code, 'lvs'] = ' '.join([str(depth) for depth in depths])
    return ernst_drr     

def generate_paved(catchments, landuse, surface_level, street_storage, sewer_storage, pump_capacity, meteo_areas):
    """ 
        Combine all data to form a complete PAVED definition. ALso the coordinates for the networ topology are included.
    Zonal statistics are applied to land use, to get the paved area. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    """
    
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True)   
    rast, affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, rast, affine=affine, stats="median",all_touched=True)         
    
    if not isinstance( street_storage, float):                
        rast, affine = read_raster(street_storage, static=True)
        str_stors = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)        
    if not isinstance( sewer_storage, float):                
        rast, affine = read_raster(sewer_storage, static=True)
        sew_stors = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)        
    if not isinstance(pump_capacity, float):     
        rast, affine = read_raster(pump_capacity, static=True)               
        pump_caps = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)    
        
    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    paved_drr = ExtendedDataFrame(required_columns=['code'])
    paved_drr.set_data( pd.DataFrame(np.zeros((len(catchments),9)), 
                                       columns=['code','area','mvlevel', 'streetstor', 'sewstor', 'pumpcap', 'px', 'py', 'boundary'], dtype="str"), index_col='code')
    paved_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):             
        # find corresponding meteo-station
        ms = [ms for ms in meteo_areas.itertuples() if ms.geometry.contains(cat.geometry.centroid)]
        ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]
            
        elev = mean_elev[num]['median']/100
        paved_drr.at[cat.code, 'code'] = str(cat.code)
        paved_drr.at[cat.code, 'area'] = str(lu_counts[num][14]*int(px_area)) if 14 in lu_counts[num] else '0'        
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
        paved_drr.at[cat.code,'meteostat'] = ms[0]
        paved_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]+10:.0f}'
        paved_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        paved_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return paved_drr   
    
   
def generate_greenhouse(catchments, landuse, surface_level, roof_storage, meteo_areas):    
    """ 
        Combine all data to form a complete GREENHSE defintion. ALso the coordinates for the network topology are included.
    Zonal statistics are applied to land use, to get the paved area. The classificition described in the notebook is assumed. From the elevation grid, the median value per catchment is assumed. Other parameters can be prescribed as as float (spatially uniform) or as a raster name, in which case the mean value per catchment is used.
    """
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True)   
    rast, affine = read_raster(surface_level, static=True)
    mean_elev = zonal_stats(catchments, rast, affine=affine, stats="median", all_touched=True) 
    # optional rasters
    if not isinstance(roof_storage, float):                
        rast, affine = read_raster(roof_storage, static=True)
        roofstors = zonal_stats(catchments, rast, affine=affine, stats="mean", all_touched=True)    
    
    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    gh_drr = ExtendedDataFrame(required_columns=['code'])
    gh_drr.set_data( pd.DataFrame(np.zeros((len(catchments),7)), 
                                       columns=['code','area','mvlevel', 'roofstor', 'px', 'py', 'boundary'], dtype="str"), index_col='code')
    gh_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):             
        # find corresponding meteo-station
        ms = [ms for ms in meteo_areas.itertuples() if ms.geometry.contains(cat.geometry.centroid)]
        ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]
            
        elev = mean_elev[num]['median']/100
        gh_drr.at[cat.code, 'code'] = str(cat.code)
        gh_drr.at[cat.code, 'area'] = str(lu_counts[num][15]*int(px_area)) if 15 in lu_counts[num] else '0'        
        gh_drr.at[cat.code, 'mvlevel'] = f'{elev:.2f}' 
        if isinstance(roof_storage, float):    
            gh_drr.at[cat.code, 'roofstor'] = f'{roof_storage:.2f}'
        else:
            gh_drr.at[cat.code, 'roofstor'] = f'{roofstors[num]["mean"]:.2f}'    
        gh_drr.at[cat.code, 'meteostat'] = ms[0]
        gh_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]+20:.0f}'
        gh_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        gh_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return gh_drr   

def generate_openwater(catchments, landuse, meteo_areas):    
    """ 
        Combine all data to form a complete OPENWATE definotion. ALso the coordinates for the network topology are included.
    Zonal statistics are applied to land use, to get the open water area. The classificition described in the notebook is assumed. 
    """
    
    lu_rast, lu_affine = read_raster(landuse, static=True)
    lu_counts = zonal_stats(catchments, lu_rast, affine=lu_affine, categorical=True)   

    # get raster cellsize    
    px_area = lu_affine[0] * -lu_affine[4]
    
    ow_drr = ExtendedDataFrame(required_columns=['code'])
    ow_drr.set_data( pd.DataFrame(np.zeros((len(catchments),5)), 
                                       columns=['code','area','px', 'py', 'boundary'], dtype="str"), index_col='code')
    ow_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):    
        # find corresponding meteo-station
        ms = [ms for ms in meteo_areas.itertuples() if ms.geometry.contains(cat.geometry.centroid)]
        ms = ms[0] if ms != [] else meteo_areas.iloc[0,:][0]
        
        ow_drr.at[cat.code, 'code'] = str(cat.code)
        ow_drr.at[cat.code, 'area'] = str(lu_counts[num][13]*int(px_area)) if 13 in lu_counts[num] else '0'        
        ow_drr.at[cat.code, 'meteostat'] = ms[0]
        ow_drr.at[cat.code, 'px'] = f'{cat.geometry.centroid.coords[0][0]-20:.0f}'
        ow_drr.at[cat.code, 'py'] = f'{cat.geometry.centroid.coords[0][1]:.0f}'
        ow_drr.at[cat.code, 'boundary'] = cat.lateraleknoopcode                        
    return ow_drr   

def generate_boundary(boundary_nodes, catchments):
    """
    Method to create boundary  nodes for RR.

    """
    bnd_drr = ExtendedDataFrame(required_columns=['code'])
    bnd_drr.set_data( pd.DataFrame(np.zeros((len(catchments),3)), 
                                       columns=['code', 'px', 'py'], dtype="str"), index_col='code')                                        
    bnd_drr.index = catchments.code
    for num, cat in enumerate(catchments.itertuples()):    
        bnd_drr.at[cat.code, 'code'] = cat.lateraleknoopcode
        bnd_drr.at[cat.code, 'px']  = boundary_nodes[boundary_nodes['code']==cat.lateraleknoopcode]['X'].to_string(index=False).strip()
        bnd_drr.at[cat.code, 'py']  = boundary_nodes[boundary_nodes['code']==cat.lateraleknoopcode]['Y'].to_string(index=False).strip()
    return bnd_drr    

def generate_seepage(catchments, seepage_folder):    
    """
    Method to obtain catchment-average seepage fluxes from rasters. The time step is deduced from the raster filenames.
    
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(seepage_folder)
    times = []    
    for ifile, file in enumerate(file_list):
        array, affine, time = read_raster(os.path.join(seepage_folder, file))
        times.append(time)           
        stats = zonal_stats(catchments, array, affine=affine, stats="median", all_touched=True)
        if ifile==0:
            result = pd.DataFrame( [[s['median'] for s in stats]] , columns='sep_'+catchments.code)
        else:
            result = result.append(pd.DataFrame( [[s['median'] for s in stats]], dtype=float, columns='sep_'+catchments.code), ignore_index=True)    
    result.index = times 
    return result

def generate_precip(areas, precip_folder):
    """
    Method to obtain catchment-average seepage fluxes from rasters. The time step is deduced from the raster filenames.
    
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(precip_folder)
    times = []    
    for ifile, file in enumerate(file_list):
        array, affine, time = read_raster(os.path.join(precip_folder, file))
        times.append(time)           
        stats = zonal_stats(areas, array, affine=affine, stats="mean", all_touched=True)
        if ifile==0:
            result = pd.DataFrame( [[s['mean'] for s in stats]] , columns='ms_'+areas.code)
        else:
            result = result.append(pd.DataFrame( [[s['mean'] for s in stats]], dtype=float, columns='ms_'+areas.code), ignore_index=True)    
    result.index = times 
    return result

def generate_evap(areas, evap_folder, dissolve_field=None):    
    """
    Method to obtain catchment-average evaporation fluxes from rasters. The time step is deduced from the raster filenames. Since only one timeeries is allowed, the meteo areas are dissolved to a user specifield field.
    
    """
    warnings.filterwarnings('ignore')
    file_list = os.listdir(evap_folder)
    # aggregated evap
    agg_areas = areas.iloc[0:len(areas),:].dissolve(by=dissolve_field,aggfunc='mean')
    times = []    
    for ifile, file in enumerate(file_list):
        array, affine, time = read_raster(os.path.join(evap_folder, file))
        times.append(time)                  
        stats = zonal_stats(agg_areas, array, affine=affine, stats="mean",all_touched=True)
        if ifile==0:
            result = pd.DataFrame( [[s['mean'] for s in stats]] , columns=['ms_'+areas.iloc[0,0]])
        else:
            result = result.append(pd.DataFrame( [[s['mean'] for s in stats]], dtype=float, columns=['ms_'+areas.iloc[0,0]]), ignore_index=True)    
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
        