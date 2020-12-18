# coding: latin-1
import os
import re
import shutil
import datetime
import warnings
import numpy as np
from shapely.geometry import box, LineString, Point
import pandas as pd
import csv
import delft3dfmpy
from delft3dfmpy.core.geometry import orthogonal_line

import logging

logger = logging.getLogger(__name__)

class DFlowRRWriter:
    """Writer for RR files"""

    def __init__(self, rrmodel, output_dir, name):
        self.rrmodel = rrmodel
        # self.geometries = parent.geometries
        # self.boundary_conditions = parent.boundary_conditions
        self.output_dir = os.path.join(output_dir, 'rr')
        self.name = name

        self.version = delft3dfmpy.DFlowFMWriter.version
        
        if os.path.exists(self.output_dir):
            shutil.rmtree(self.output_dir)

        self.d3b_parameters = rrmodel.d3b_parameters
        
        self.precip_df  = {}
        self.run_dimrpad = rrmodel.dimr_path
        
    def write_all(self):  # write all RR files
        """
        Wrapper method to write all components
        """        
        # copy all 'standard' files from the resources folder
        self.copyRRFiles()

        self.write_topology()
        
        self.write_unpaved()
        
        self.write_paved()
        
        self.write_greenhouse()
        
        self.write_openwater()
        
        self.write_meteo()
        
        self.write_coupling()

        # Change mdu parameters
        for parameter, value in self.d3b_parameters.items():
            self.change_d3b_parameter(parameter, value)
        
    def copyRRFiles(self):
        """
        Method to copy the RR-files from resources.

        Returns
        -------
        bool
            Succesfull or not.

        """
        srcRR = os.path.join(os.path.dirname(__file__), 'resources', 'RRmodel')
        targetRR = os.path.join(self.output_dir)
        shutil.copytree(srcRR, targetRR)
        return True

    def write_topology(self):
        """
        Wrapper to write the topolgy files for RR. The following files are written:
            3B_NOD.TP: IDs and coordinates of all nodes: unpaved, paved, greenhouse, open water and boundary. Only nodes with a total area of >0 are written. Nodes are placed up to 20 m to the east and west of the catchment centroid.
            3B_LINK.TP: Links between all nodes and the correct boundaries.
            BOUND3B.3B: coupling information for every RR boundary.
            BoundaryConditions.bc: see bound3b.3b
        Returns
        -------
        None.

        """
        filepath = os.path.join(self.output_dir, '3B_NOD.TP')
        with open(filepath, 'w') as f:

            # Unpaved nodes
            if any(self.rrmodel.unpaved.unp_nodes):
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:                    
                        f.write('NODE id \''+dct['id']+'\' nm \''+dct['id']+'\' ri \'-1\' mt 1 \'2\' nt 44 ObID \'3B_UNPAVED\' px '+dct['px']+' py '+dct['py']+' node\n')
                              
            # Paved nodes
            if any(self.rrmodel.paved.pav_nodes):
                for _, dct in self.rrmodel.paved.pav_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        f.write('NODE id \''+dct['id']+'\' nm \''+dct['id']+'\' ri \'-1\' mt 1 \'1\' nt 43 ObID \'3B_PAVED\' px '+dct['px']+' py '+dct['py']+' node\n')
                  
            # Greenhouse nodes
            if any(self.rrmodel.greenhouse.gh_nodes):
                for _, dct in self.rrmodel.greenhouse.gh_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        f.write('NODE id \''+dct['id']+'\' nm \''+dct['id']+'\' ri \'-1\' mt 1 \'3\' nt 45 ObID \'3B_GREENHOUSE\' px '+dct['px']+' py '+dct['py']+' node\n')

            # Openwater nodes
            if any(self.rrmodel.openwater.ow_nodes):
                for _, dct in self.rrmodel.openwater.ow_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        f.write('NODE id \''+dct['id']+'\' nm \''+dct['id']+'\' ri \'-1\' mt 1 \'21\' nt 46 ObID \'OW_PRECIP\' px '+dct['px']+' py '+dct['py']+' node\n')
                      
            # Boundary_nodes
            if any(self.rrmodel.external_forcings.boundary_nodes):
                for _, dct in self.rrmodel.external_forcings.boundary_nodes.items():
                    f.write('NODE id \''+dct['id']+'\' nm \''+dct['id']+'\' ri \'-1\' mt 1 \'6\' nt 78 ObID \'SBK_SBK-3B-NODE\' px '+dct['px']+' py '+dct['py']+' node\n')
                    
        filepath = os.path.join(self.output_dir, '3B_LINK.TP')
        with open(filepath, 'w') as f:
          
            cnt = 0
            if any(self.rrmodel.unpaved.unp_nodes):                
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:                    
                        cnt += 1
                        f.write('BRCH id \''+str(cnt)+'\' ri \'-1\' mt 1 \'0\' bt 17 ObID \'3B_LINK\' bn \''+dct['id']+'\' en \''+dct['boundary_node']+'\' brch\n')
            if any(self.rrmodel.paved.pav_nodes):                
                for _, dct in self.rrmodel.paved.pav_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        cnt += 1
                        f.write('BRCH id \''+str(cnt)+'\' ri \'-1\' mt 1 \'0\' bt 17 ObID \'3B_LINK\' bn \''+dct['id']+'\' en \''+dct['boundary_node']+'\' brch\n')
            if any(self.rrmodel.greenhouse.gh_nodes):                
                for _, dct in self.rrmodel.greenhouse.gh_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        cnt += 1
                        f.write('BRCH id \''+str(cnt)+'\' ri \'-1\' mt 1 \'0\' bt 17 ObID \'3B_LINK\' bn \''+dct['id']+'\' en \''+dct['boundary_node']+'\' brch\n')
            if any(self.rrmodel.openwater.ow_nodes):                
                for _, dct in self.rrmodel.openwater.ow_nodes.items():
                    if float(dct['ar']) > 0.0:                    
                        cnt += 1
                        f.write('BRCH id \''+str(cnt)+'\' ri \'-1\' mt 1 \'0\' bt 17 ObID \'3B_LINK\' bn \''+dct['id']+'\' en \''+dct['boundary_node']+'\' brch\n')

        # bound3b.3b        
        filepath = os.path.join(self.output_dir, 'Bound3B.3B')
        with open(filepath, 'w') as f:
            for _, dct in self.rrmodel.external_forcings.boundary_nodes.items():                
                f.write('BOUN id \''+dct['id']+'\' bl 2 \'0\' is 0 boun\n')
                
        # BoundaryConditions.bc
        filepath = os.path.join(self.output_dir, 'BoundaryConditions.bc')
        header = {'majorVersion':'1', 'minorVersion':'0', 'fileType':'boundConds'}
        with open(filepath, 'w') as f:                        
            self._write_dict(f, header, 'General','\n')            
            for _, dct in self.rrmodel.external_forcings.boundary_nodes.items():                
                temp = {"name":''+dct['id'], 'function':'constant','quantity':'water_level','unit':'m'} 
                self._write_dict(f,temp,'Boundary','    0\n\n')
                                   
    def write_unpaved(self):
        """
        Method to write all files associated with unpaved nodes: UNPAVED.3B, UNPAVED.ALF, UNPAVED.STO, UNPAVED.INF and UNPAVED.SEP. All files contain a definition for every node  because they may or may not be spatially distributed.

        Returns
        -------
        None.

        """
        if any(self.rrmodel.unpaved.unp_nodes): 
            filepath = os.path.join(self.output_dir, 'UNPAVED.3B')
            with open(filepath, 'w') as f:
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:                    
                        f.write('UNPV id \''+dct['id']+'\' na '+dct['na']+' gw '+dct['ga']+' ar '+dct['ar']+' lv '+dct['lv']+' co '+dct['co']+' su '+dct['su']+' sd \'sto_'+dct['id']+'\' ic \'inf_'+dct['id']+'\' bt '+dct['bt']+' ed \''+dct['ed']+'\' sp \''+dct['sp']+'\' ig 0 '+dct['ig']+' mg '+dct['mg']+' gl '+dct['gl']+' is '+dct['is']+' ms \''+dct['ms']+'\' unpv\n')                    
                
            filepath = os.path.join(self.output_dir, 'UNPAVED.STO')
            with open(filepath, 'w') as f:            
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:
                        f.write('STDF id \'sto_'+dct['id']+'\' nm \'sto_'+dct['id']+'\' ml '+dct['sd']+' il 0 stdf\n')
                
            filepath = os.path.join(self.output_dir, 'UNPAVED.INF')
            with open(filepath, 'w') as f:            
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:
                        f.write('INFC id \'inf_'+dct['id']+'\' nm \'inf_'+dct['id']+'\' ic '+dct['ic']+' infc\n')
    
            filepath = os.path.join(self.output_dir, 'UNPAVED.ALF')
            with open(filepath, 'w') as f:            
                for _, dct in self.rrmodel.unpaved.ernst_defs.items():                    
                        f.write('ERNS id \''+dct['id']+'\' nm \''+dct['id']+'\' cvi '+dct['cvi']+' cvo 0 '+dct['cvo']+' cvs '+dct['cvs']+' lv '+dct['lv']+' erns\n')                    
                                                
            filepath = os.path.join(self.output_dir, 'UNPAVED.SEP')
            with open(filepath, 'w') as f:            
                for _, dct in self.rrmodel.unpaved.unp_nodes.items():                                            
                    if np.sum([float(d) for d in dct['ar'].split(' ')]) > 0.0:
                        f.write('SEEP id \''+dct['sp']+'\' nm \''+dct['sp']+'\' co 4 PDIN 1 1 \'365;00:00:00\' pdin ss 0\n')
                        f.write('TBLE\n')
                        for row in self.rrmodel.external_forcings.seepage[dct['sp']]['seepage'].iteritems():
                            f.write('\''+row[0].strftime('%Y/%m/%d;%H:%M:%S')+'\' '+f'{row[1]:.5f} <\n')
                        f.write('tble\nseep\n')              
            
    def write_paved(self):
        """
        Method to write all files associated with paved nodes: PAVED.3B, PAVED.STO and PAVED.DWA. The latter contains only one definition.

        Returns
        -------
        None.

        """
        if any(self.rrmodel.paved.pav_nodes):
            filepath = os.path.join(self.output_dir, 'PAVED.3B')
            with open(filepath, 'w') as f:          
                for _, dct in self.rrmodel.paved.pav_nodes.items():                    
                    if float(dct['ar']) > 0.0:                    
                        f.write('PAVE id \''+dct['id']+'\' ar '+dct['ar']+' lv '+dct['lv']+' sd \'sto_'+dct['id']+'\' ss 0 qc 0 '+dct['qc']+' 0 qo 1 0 ms \''+dct['ms']+'\' aaf 1 is '+dct['is']+' np '+dct['np']+' dw \'Def_DWA\' ro '+dct['ro']+' ru '+dct['ru']+' qh \'\' pave\n')
                                                                                                                                                                                                           
            filepath = os.path.join(self.output_dir, 'PAVED.STO')                       
            with open(filepath, 'w') as f:           
                for _, dct in self.rrmodel.paved.pav_nodes.items():                    
                    if float(dct['ar']) > 0.0:                                   
                        f.write('STDF id \'sto_'+dct['id']+'\' nm \'sto_'+dct['id']+'\' ms '+dct['strs']+' is 0 mr '+dct['sews']+' 0.0 ir 0.0 0.0 stdf\n')
                    
            filepath = os.path.join(self.output_dir, 'PAVED.DWA')             
            with open(filepath, 'w') as f:                           
                f.write('DWA id \'Def_DWA\' nm \'Def_DWA\' do 1 wc 0 wd 0 wh 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 dwa\n')
                    
    def write_greenhouse(self):
        """
        Method to write all files associated with greenhouse nodes: GREENHSE.3B, GREENHSE.RF and GREENHSE.SIL. The latter contains only one definition.

        Returns
        -------
        None.

        """
        if any(self.rrmodel.greenhouse.gh_nodes):
            filepath = os.path.join(self.output_dir, 'GREENHSE.3B')
            with open(filepath, 'w') as f:                    
                for _, dct in self.rrmodel.greenhouse.gh_nodes.items():    
                    if float(dct['ar']) > 0.0:                    
                        f.write('GRHS id \''+dct['id']+'\' na 10 ar 0 0 '+dct['ar']+' 0 0 0 0 0 0 0 sl '+dct['sl']+' as 0 si \'Def_silo\' sd \'sto_'+dct['id']+'\' ms \''+dct['ms']+'\' is '+dct['is']+' grhs\n')      
                       
            filepath = os.path.join(self.output_dir, 'GREENHSE.SIL')
            with open(filepath, 'w') as f:
                f.write('SILO id \'Def_silo\' nm \'Def_silo\' sc 0.0 pc 0.0 silo\n')
            
            filepath = os.path.join(self.output_dir, 'GREENHSE.RF')
            with open(filepath, 'w') as f:
                for _, dct in self.rrmodel.greenhouse.gh_nodes.items():    
                    if float(dct['ar']) > 0.0:                                    
                        f.write('STDF id \'sto_'+dct['id']+'\' nm \'sto_'+dct['id']+'\' mk '+dct['sd']+' ik 0 stdf\n')
    
    def write_openwater(self):
         """
         Method to write OPENWATE.3B file.

        Returns
        -------
        None.

        """
         # write openwater.3b
         if any(self.rrmodel.openwater.ow_nodes):
            filepath = os.path.join(self.output_dir, 'OPENWATE.3B')
            with open(filepath, 'w') as f:                    
                for _, dct in self.rrmodel.openwater.ow_nodes.items(): 
                    if float(dct['ar']) > 0.0:                    
                        f.write('OWRR id \''+dct['id']+'\' ar '+dct['ar']+' ms \''+dct['ms']+'\' owrr\n')                                 
         
    def write_meteo(self):      
        """
        Method to write meteofiles (DEFAULT.BUI and DEFAULT.EVP) based on values per catchment.

        Returns
        -------
        None.

        """            
         # precip         
        df = {}
        for ii,ms in enumerate( self.rrmodel.external_forcings.precip.items()):                 
            df[ms[0]] = pd.Series(ms[1]['precip'])
        filepath = os.path.join(self.output_dir, 'DEFAULT.BUI')
        f = open(filepath, 'w')
        f.write('*Name of this file: c:\Result\1058\DEFAULT.BUI\n')
        f.close()
        f = open(filepath,'a')
        f.write('*Date and time of construction: 00/00/2000 00:00:00.\n')
        f.write('1\n')
        f.write('*Aantal stations\n')
        f.write(str(len(self.rrmodel.external_forcings.precip))+'\n')
        f.write('*Namen van stations\n')
        for ms in self.rrmodel.external_forcings.precip.items():
            f.write('\''+ms[0]+'\''+'\n' )
        f.write('*Aantal gebeurtenissen (omdat het 1 bui betreft is dit altijd 1)\n')
        f.write('*en het aantal seconden per waarnemingstijdstap\n')
        f.write(f'1 {(ms[1]["precip"].index[1]-ms[1]["precip"].index[0]).total_seconds():.0f}\n')
        f.write('*Elke commentaarregel wordt begonnen met een * (asterisk).\n')
        f.write('*Eerste record bevat startdatum en -tijd, lengte van de gebeurtenis in dd hh mm ss\n')
        f.write('*Het format is: yyyymmdd:hhmmss:ddhhmmss')
        f.write('*Daarna voor elk station de neerslag in mm per tijdstap.\n')
        duration = (ms[1]['precip'].index[-1]-ms[1]['precip'].index[0]).total_seconds()
        f.write(ms[1]['precip'].index[0].strftime("%Y %#m %#d %#H %#M %#S ")+str(int(duration/86400.))+' '+str(int(np.remainder(duration,86400.)/3600.))+' 0 0\n')                        
        f.close()
        self._dict_to_df()
        self.precip_df.to_csv(filepath, header=False, index=False, sep=" ", float_format='%.3f', mode='a')
#%%
          
#%%%          
        # evaporation         
        filepath = os.path.join(self.output_dir, 'DEFAULT.EVP')
        f = open(filepath,'w')
        f.write('*Verdampingsfile\n')
        f.write('*Meteo data: evaporation intensity in mm/day\n')
        f.write('*First record: start date, data in mm/day\n')
        f.write('*Datum (year month day), verdamping (mm/dag) voor elk weerstation\n')
        f.write('*jaar maand dag verdamping[mm]\n')     
        f.close()
        table = list(self.rrmodel.external_forcings.evap.values())[0]['evap']          
        table.to_csv(filepath,float_format='%.3f', date_format='%#Y  %#m  %#d ', sep=" ", header=False ,quoting=csv.QUOTE_NONE, mode="a", escapechar=" ")         
               
    def _dict_to_df(self):        
        """
        Method to convert the precip dictionary to a dataframe
        Returns
        -------
        None.

        """
        df = {}
        for ii, stat in enumerate(self.rrmodel.external_forcings.precip.items()):              
            df[stat[0]] = pd.Series(stat[1]['precip'])        
        self.precip_df = pd.DataFrame(df)
    
    def _write_dict(self, f, dct, header, endline):
        """
        Write a dictionary to file, without transforming with mapping
        """
        # Write header
        f.write(f'[{header}]\n')
        for key, value in dct.items():
            f.write(f'    {key:20s} = {value:20s}\n')
        f.write(endline)
        
    def change_d3b_parameter(self, parameter, value):
        """
        Change parameter value in delft3b.ini file
        """
        with open(os.path.join(self.output_dir,'DELFT_3B.INI'), 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.lower().startswith(parameter.lower()):
                items = re.split('=|#|$\s+', line)
                key, oldvalue = items[0], items[1]
                if not key.strip().lower() == parameter.lower():
                    continue
                lines[i] = line.replace('='+oldvalue, f'={value}'.ljust(len(oldvalue)+1))
                break

        with open(os.path.join(self.output_dir, 'DELFT_3B.INI'), 'w') as f:
            f.write(''.join(lines))
      
    def write_coupling(self):        
        """
        Method to write the files required for coupling DFM and RR.
        """
        # run.bat
        with open(os.path.join(self.output_dir, '../run.bat'),'w') as f:
                 f.write('@ echo off\n')
                 f.write('set OMP_NUM_THREADS=2\n')
                 f.write('call \"'+self.run_dimrpad+'\"\n')
                 f.write('pause\n')

        # coupling XML
        DFM_comp_name = 'CoupledTest1_DFM'
        RR_comp_name = 'CoupledTest1_RR'
        with open(os.path.join(self.output_dir, '../dimr_config.xml'),'w') as f:       
                 f.write('<?xml version="1.0" encoding="utf-8" standalone="yes"?>\n')
                 f.write('<dimrConfig xmlns="http://schemas.deltares.nl/dimr" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://schemas.deltares.nl/dimr http://content.oss.deltares.nl/schemas/dimr-1.2.xsd">\n')
                 f.write('\t<documentation>\n')
                 f.write('\t\t<fileVersion>1.2</fileVersion>\n')
                 f.write('\t\t<createdBy>D-HyDAMO in delft3dfmpy v.'+self.version['number']+'</createdBy>\n')
                 f.write("\t\t<creationDate>"+self.version['date']+'</creationDate>\n')
                 f.write('\t</documentation>\n\n')        
                 f.write('\t<control>\n')
                 f.write('\t\t<parallel>\n')
                 f.write('\t\t\t<startGroup>\n')
                 tstr = f'0 {self.rrmodel.d3b_parameters["Timestepsize"]:.0f} {(pd.Timestamp(self.rrmodel.d3b_parameters["EndTime"])-pd.Timestamp(self.rrmodel.d3b_parameters["StartTime"])).total_seconds():.0f}'
                 f.write('\t\t\t\t<time>'+tstr+'</time>\n')
                 f.write('\t\t\t\t<coupler name="flow2rr" />\n')
                 f.write('\t\t\t\t<start name="'+RR_comp_name+'" />\n')
                 f.write('\t\t\t\t<coupler name="rr2flow" />\n')
                 f.write('\t\t\t</startGroup>\n')
                 f.write('\t\t\t<start name="'+DFM_comp_name+'" />\n')
                 f.write('\t\t</parallel>\n')
                 f.write('\t</control>\n\n')                                      

                 f.write('\t<component name="'+DFM_comp_name+'">\n')
                 f.write('\t\t<library>dflowfm</library>\n')
                 f.write('\t\t<workingDir>fm</workingDir>\n')
                 f.write('\t\t<inputFile>'+self.name+'.mdu</inputFile>\n')
                 f.write('\t</component>\n')                 
                 f.write('\t<component name="'+RR_comp_name+'">\n')
                 f.write('\t\t<library>rr_dll</library>\n')
                 f.write('\t\t<workingDir>rr</workingDir>\n')
                 f.write('\t\t<inputFile>Sobek_3b.fnm</inputFile>\n')
                 f.write('\t</component>\n\n')
        
                 f.write('\t<coupler name="flow2rr">\n')
                 f.write('\t\t<sourceComponent>'+DFM_comp_name+'</sourceComponent>\n')
                 f.write('\t\t<targetComponent>'+RR_comp_name+'</targetComponent>\n')		     
                 for i in self.rrmodel.external_forcings.boundary_nodes.items():
                     f.write('\t\t\t<item>\n')                     
                     # FOr now use one observation ponit- until water levels can be read from other nodes
                     #f.write('\t\t\t\t<sourceName>observations/'+str(i[0])+'/water_level</sourceName>\n')
                     f.write('\t\t\t\t<sourceName>laterals/'+str(i[0])+'/water_level</sourceName>\n')
                     f.write('\t\t\t\t<targetName>catchments/'+str(i[0])+'/water_level</targetName>\n')
                     f.write('\t\t\t</item>\n')                                  
                 f.write('\t\t<logger>\n')
                 f.write('\t\t\t<workingDir>.</workingDir>\n')
                 f.write('\t\t\t<outputFile>dflowfm_to_rr.nc</outputFile>\n')
                 f.write('\t\t</logger>\n')
                 f.write('\t</coupler>\n\n')
                 
                 f.write('\t<coupler name="rr2flow">\n')
                 f.write('\t\t<sourceComponent>'+RR_comp_name+'</sourceComponent>\n')
                 f.write('\t\t<targetComponent>'+DFM_comp_name+'</targetComponent>\n')		         
                 for i in self.rrmodel.external_forcings.boundary_nodes.items():
                     f.write('\t\t\t<item>\n')                     
                     f.write('\t\t\t\t<sourceName>catchments/'+str(i[0])+'/water_discharge</sourceName>\n')
                     f.write('\t\t\t\t<targetName>laterals/'+str(i[0])+'/water_discharge</targetName>\n')
                     f.write('\t\t\t</item>\n')                                 

                 f.write('\t\t<logger>\n')
                 f.write('\t\t\t<workingDir>.</workingDir>\n')
                 f.write('\t\t\t<outputFile>rr_to_dflowfm.nc</outputFile>\n')
                 f.write('\t\t</logger>\n')                 
                 f.write('\t</coupler>\n')
                 f.write('</dimrConfig>\n')
        
    