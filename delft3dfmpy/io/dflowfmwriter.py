# coding: latin-1
import os
import re
import shutil
import datetime

import numpy as np
from shapely.geometry import box, LineString, Point
import delft3dfmpy
from delft3dfmpy.io.UgridWriter import UgridWriter
from delft3dfmpy.core.geometry import orthogonal_line

import logging

logger = logging.getLogger(__name__)

class DFlowFMWriter:
    """Writer for FM files"""

    # versioning info
    version = { 'number'        : delft3dfmpy.__version__,
                'date'          : datetime.datetime.strftime(datetime.datetime.today(),'%Y-%m-%d %H:%M:%S'),
                'dfm_version'   : '1.2.105.67228M',
                'dimr_version'  : '2.13.02.67836',
                'suite_version' : '(Beta) 0.9.6.51435'} 
   
    def __init__(self, dflowfmmodel, output_dir, name):
        self.dflowfmmodel = dflowfmmodel
        self.output_dir = os.path.join(output_dir, 'fm')
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        self.name = name
       
        self.version = DFlowFMWriter.version
        
        self.mdu_parameters = dflowfmmodel.mdu_parameters
        self.mdu_parameters['NetFile'] = self.name + "_net.nc"

        self.mdu_parameters['ExtForceFileNew'] = self.name + "_new.ext"

        self.mdufile = os.path.join(self.output_dir, self.name +'.mdu')
        self.extfile_new = os.path.join(self.output_dir, self.mdu_parameters['ExtForceFileNew'])
        self.netfile = os.path.join(self.output_dir, self.mdu_parameters['NetFile'])
        
        if hasattr(self.dflowfmmodel, 'dimr_path'):
           self.run_dimrpad = self.dflowfmmodel.dimr_path
       
        with open(os.path.join(self.output_dir, 'boundaries.bc'), 'w') as f:
            self._write_header(f, 'boundConds', 1.01, extra_linebreak=False)

    def write_all(self):  # write all fm files from HyDAMO
        """
        Wrapper to write all components to DFM. Remove existing files and create new ones. Note that the old format ext-file is not used anymore,.
        """
        if os.path.exists(self.extfile_new):
            os.remove(self.extfile_new)

        with open(self.extfile_new,'w') as f:
            self._write_header(f, 'extForce', 2.01)

        # Write grid in Ugrid format to netcdf
        ugridWriter = UgridWriter()
        ugridWriter.write(dflowfmmodel=self.dflowfmmodel, path=self.netfile, version=self.version)

        # mdu file
        self.writeMdFiles()
        # Roughness files
        self.write_roughness()
        # Cross section definitions
        self.write_crosssection_definitions()
        # Cross section locations
        self.write_crosssection_locations()
        # Structures (weirs, pumps, culverts)
        self.write_structures()
        # Initial conditions
        self.write_initial_conditions()
        # Boundary conditions
        self.write_boundary_conditions()
        # Observation points
        self.write_observation_points()
        if hasattr(self, 'run_dimrpad'):
            # Files for running DIMR
            self.write_dimrfiles()
        
        # Change mdu parameters
        for parameter, value in self.mdu_parameters.items():
            self.change_mdu_parameter(parameter, value)

        # Write laterals
        if any(self.dflowfmmodel.external_forcings.laterals):
            self.write_laterals()

    def _write_dict(self, f, dct, header,extra_linebreak=True):
        """
        Write a dictionary to file, without transforming with mapping
        """
        # Write header
        f.write(f'[{header}]\n')
        for key, value in dct.items():
            f.write(f'\t{key} = {value}\n')
        if extra_linebreak: f.write('\n')

    def writeMdFiles(self):
        """
        Copy MDU file form resouces.
        """
        srcMdu = os.path.join(os.path.dirname(__file__), 'resources', 'FMmdu.txt')
        targetMdu = os.path.join(self.output_dir, f'{self.name}.mdu')

        shutil.copy2(srcMdu, targetMdu)

        return True

    def write_roughness(self):

        # Write roughness
        roughnessfiles = []

        # Create a direcroty for the roughness files
        roughnessdir = 'roughness'
        if not os.path.exists(os.path.join(self.output_dir, roughnessdir)):
            os.mkdir(os.path.join(self.output_dir, roughnessdir))

        # Write each roughness definition to a file
        for i, (_, roughness) in enumerate(self.dflowfmmodel.network.roughness_definitions.items()):
            # Write roughness file for channels
            #file = os.path.join(roughnessdir, f'r{i+1:03d}.ini')
            file = f'{roughnessdir}/r{i+1:03d}.ini'
            roughnessfiles.append(file)
            with open(os.path.join(self.output_dir, file), 'w') as f:
                # header roughness
                self._write_header(f, filetype='roughness', fileversion=3.00)
                dct = {
                    'frictionId': roughness['name'],
                    'frictionType': roughness['code'],
                    'frictionValue': roughness['value'],
                }
                self._write_dict(f, dct=dct, header='Global')

        self.mdu_parameters['frictFile'] = ';'.join(roughnessfiles)

    def write_crosssection_locations(self):
        # Write cross section locations
        with open(os.path.join(self.output_dir, 'cross_section_locations.ini'), 'w') as f:
            self._write_header(f, 'crossLoc', fileversion=1.01)
            for _, series in self.dflowfmmodel.crosssections.crosssection_loc.items():
                # Write the cross section item
                self._write_dict(f, dct=series, header='CrossSection')

            # Write the default profile (if specified) to branches that do not have a cross section
            if self.dflowfmmodel.crosssections.default_definition is not None:
                # Find branches without profile
                no_css = self.dflowfmmodel.crosssections.get_branches_without_crosssection()
                branches = (self.dflowfmmodel.network.branches.loc[no_css, 'geometry'].length / 2).to_dict()
                for branchid, chainage in branches.items():
                    dct = {
                        'id': f'{branchid}_{chainage:.1f}',
                        'branchid': branchid,
                        'chainage': chainage,
                        'shift': self.dflowfmmodel.crosssections.default_definition_shift,
                        'definitionId': self.dflowfmmodel.crosssections.default_definition
                    }
                    # Write the cross section item
                    self._write_dict(f, dct=dct, header='CrossSection')

    def _write_header(self, f, filetype, fileversion, extra_linebreak=True):
        f.write(f"# Generated with D-HyDAMO in delftd3fmpy v.{self.version['number']} on {self.version['date']}.\n")
        f.write(f"# Tested and compatible with D-Flow FM {self.version['dfm_version']}, DIMRset {self.version['dimr_version']} and D-HYDRO suite 1D2D {self.version['suite_version']}.\n")
        f.write('# \n')
        f.write('[General]\n')
        f.write(f'fileVersion = {fileversion:.2f}\n')
        f.write(f'fileType = {filetype}\n')
        if extra_linebreak:
            f.write('\n')

    def write_crosssection_definitions(self):

        with open(os.path.join(self.output_dir, 'cross_section_definitions.ini'), 'w') as f:
            # header
            self._write_header(f, filetype='crossDef', fileversion=3.00)
            # definitions
            for _, dct in self.dflowfmmodel.crosssections.crosssection_def.items():
                self._write_dict(f, dct=dct, header='definition')

    def write_structures(self):

        filepath = os.path.join(self.output_dir, 'structure.ini')
        with open(filepath, 'w') as f:
            # write header
            self._write_header(f, filetype='structure', fileversion=2.00)

            # Culverts
            if any(self.dflowfmmodel.structures.culverts):
                for _, dct in self.dflowfmmodel.structures.culverts.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Weirs
            if any(self.dflowfmmodel.structures.weirs):
                for _, dct in self.dflowfmmodel.structures.weirs.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Bridges
            if any(self.dflowfmmodel.structures.bridges):
                for _, dct in self.dflowfmmodel.structures.bridges.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Pumps
            if any(self.dflowfmmodel.structures.pumps):
                for pump_id, dct in self.dflowfmmodel.structures.pumps.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Universal weirs
            if any(self.dflowfmmodel.structures.uweirs):
                for _, dct in self.dflowfmmodel.structures.uweirs.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Orifices
            if any(self.dflowfmmodel.structures.orifices):
                for _, dct in self.dflowfmmodel.structures.orifices.items():
                    self._write_dict(f, dct=dct, header='Structure')

            # Compound structures
            if any(self.dflowfmmodel.structures.compounds):
                for _, dct in self.dflowfmmodel.structures.compounds.items():
                    self._write_dict(f, dct=dct, header='structure')

        # Write structure time series
        for bcnd in self.dflowfmmodel.external_forcings.structures.itertuples():
            # Write time series for structure
            filename = f"{bcnd.type}_{bcnd.id}.tim"
            # Get time series as data
            data = list(zip(bcnd.time, bcnd.value))
            # Write file
            write_fm_file(file=os.path.join(self.output_dir, filename), data=data)


    def write_laterals(self):
        """
        Write laterals to files. Note that the .pol files are deprececated; locations are now in the ext-file.
        """

        if not 'refdate' in self.mdu_parameters.keys():
            raise KeyError('Specify reference date before writing laterals.')

        for name, dct in self.dflowfmmodel.external_forcings.laterals.items():

            if 'timeseries' in dct:
                # Determine starttime from mdu parameters
                starttime = datetime.datetime.strptime(str(self.mdu_parameters['refdate']), '%Y%m%d')
                starttime += datetime.timedelta(seconds=int(self.mdu_parameters['tstart']))

                # Write time series - this should probably be written to a .bc file
                series = dct['timeseries'].copy()
                series.index -= starttime
                data = np.c_[series.index.total_seconds().values / 60, series.values]

                refd = str(self.dflowfmmodel.mdu_parameters["refdate"])

                 # and write data to the bc file
                with open(os.path.join(self.output_dir, 'boundaries.bc'), 'a') as f:
                    f.write(f'\n[Forcing]\n'
                        f'name       = {name}\n'
                        f'function   = timeseries\n'
                        f'timeInterpolation = linear\n'
                        f'quantity   = time\n'
                        f'unit       = minutes since {refd[0:4]}-{refd[4:6]}-{refd[6:9]} 00:00:00\n'
                        f'quantity   = lateral_discharge\n'
                        f'unit       = m3/s\n'
                       )
                write_fm_file(file=os.path.join(self.output_dir, 'boundaries.bc'), data=data, mode='a')
                discharge_kw = f'boundaries.bc'
            else:
                discharge_kw = 'REALTIME'

            # Add to external forcing file
            with open(self.extfile_new, 'a') as f:
                dct = {'id':f'{name}',
                       'name':f'{name}',
                       'type':'discharge',
                       'locationType':'1d',
                       'branchId': f'{dct["branchid"]}',
                       'chainage': f'{dct["branch_offset"]}',
                       'discharge':f'{discharge_kw}'
                      }
                self._write_dict(f,dct,'Lateral')

    def write_boundary_conditions(self):
        """
        Write boundary conditions to (new format) ext file and the time series to a .bc file.
        """

        for bc in self.dflowfmmodel.external_forcings.boundaries.values():

            # Add boundary to ext file
            with open(self.extfile_new, 'a') as f:
                # check whether 2D rainfall has been applied to the 2D model
                if any('rainfall_2D' in key for key in self.dflowfmmodel.external_forcings.boundaries):

                    if bc['bctype'] == 'rainfallbnd':
                        # get file name of rainfall netcdf
                        rainfall_netcdf = self.dflowfmmodel.external_forcings.boundaries['rainfall_2D']['file_name']

                        #new: add precipitation netcdf for 2D mesh [Meteo] header
                        dct_met = {'quantity': 'rainfall',
                                   'forcingFile': os.path.basename(rainfall_netcdf),
                                   'forcingFileType': 'netcdf'}
                        self._write_dict(f, dct_met, 'Meteo')

                        shutil.copy2(rainfall_netcdf, self.output_dir)

                    if bc['bctype'] == 'waterlevelbnd' or bc['bctype'] == 'dischargebnd':

                        dct = {'quantity':f"{bc['bctype']}",'nodeId': bc['nodeid'],'forcingFile':'boundaries.bc'}
                        self._write_dict(f, dct, 'Boundary')

                        # Get time series
                        if bc['time'] is None:
                            data = float(bc['value'])
                            headerType='constant'
                        else:
                            data = list(zip(bc['time'], bc['value']))
                            headerType='timeseries'

                else:
                    dct = {'quantity':f"{bc['bctype']}",'nodeId': bc['nodeid'],'forcingFile':'boundaries.bc'}
                    self._write_dict(f, dct, 'Boundary')

                    # Get time series
                    if bc['time'] is None:
                        data = float(bc['value'])
                        headerType='constant'
                    else:
                        data = list(zip(bc['time'], bc['value']))
                        headerType='timeseries'

                refd = str(self.dflowfmmodel.mdu_parameters["refdate"])
    
                # and write data to the bc file
                with open(os.path.join(self.output_dir, 'boundaries.bc'), 'a') as f:
                    if headerType=='timeseries':
                        f.write(f'\n[Forcing]\n'
                            f'name       = {bc["nodeid"]}\n'
                            f'function   = timeseries\n'
                            f'timeInterpolation = linear\n'
                            f'quantity   = time\n'
                            f'unit       = minutes since {refd[0:4]}-{refd[4:6]}-{refd[6:9]} 00:00:00\n'
                            f'quantity   = {bc["bctype"]}\n'
                            f'unit       = {bc["unit"]}\n'
                           )
                    if headerType=='constant':
                        f.write(f'\n[Forcing]\n'
                            f'name       = {bc["nodeid"]}\n'
                            f'function   = constant\n'
                            f'quantity   = {bc["bctype"]}\n'
                            f'unit       = {bc["unit"]}\n'
                            f'{data}\n'
                           )
                if headerType=='timeseries':
                    write_fm_file(file=os.path.join(self.output_dir, 'boundaries.bc'), data=data, mode='a')

    def write_initial_conditions(self):

        """
        Two options exist: waterlevels per polygon or water depths per point. Both are written in the new format - an iniFieldFie is created pointing to the data.
        """

        # Create folder for initial conditions if it does not exist yet
        initcondfolder = 'initialconditions'

        # initial files are not in the ext-file in the new format...
        self.dflowfmmodel.mdu_parameters['IniFieldFile'] = initcondfolder+'/initialFields.ini'

        initcondpath = os.path.join(self.output_dir, initcondfolder)
        if not os.path.exists(initcondpath):
            os.mkdir(initcondpath)

        with open(os.path.join(initcondpath,'initialFields.ini'),'w') as f:
            self._write_header(f, 'iniField',2.00)

        # we have a water level definition
        if len(self.dflowfmmodel.external_forcings.initial_waterlevel_polygons)>0:
            if self.dflowfmmodel.external_forcings.initial_waterlevel_polygons['geometry'][0] is None:
                # no polygons defined - global definition
                with open(os.path.join(initcondpath,'initialWaterLevel.ini'),'w') as f:
                    self._write_header(f, '1dField',2.00)
                    dct = {'quantity': 'waterLevel','unit': 'm+NAP','value': self.dflowfmmodel.external_forcings.initial_waterlevel_polygons.waterlevel[0]}
                    self._write_dict(f, dct, 'Global')

                with open(os.path.join(initcondpath,'initialFields.ini'),'a') as f:
                     dct  = {'quantity': 'waterLevel',
                             'dataFileType': '1dField',
                             'dataFile': initcondfolder+'/initialWaterLevel.ini'
                             }
                     self._write_dict(f,dct,'Initial')
            else:
                for row in self.dflowfmmodel.external_forcings.initial_waterlevel_polygons.itertuples():
                    # Write pol file
                    write_fm_file(file=os.path.join(initcondpath, row.Index+".pol"),
                                      data=[np.vstack(row.geometry.exterior.coords[:])],
                                      names=[row.Index]
                                  )
                    with open(os.path.join(initcondpath,'initialFields.ini'),'a') as f:
                        dct  = {'quantity': 'waterLevel',
                                'dataFileType': 'polygon',
                                'dataFile': f'{row.Index}.pol',
                                'interpolationMethod': 'constant',
                                'value': f'{row.waterlevel}'
                                }
                        self._write_dict(f,dct,'Initial')

        # we have a water depth definition
        if len(self.dflowfmmodel.external_forcings.initial_waterdepth_polygons)>0:
            if self.dflowfmmodel.external_forcings.initial_waterdepth_polygons['geometry'][0] is None:
                # no polygons defined - global definition
                with open(os.path.join(initcondpath,'initialWaterDepth.ini'),'w') as f:
                    self._write_header(f, '1dField',2.00)
                    dct = {'quantity': 'waterDepth','unit': 'm','value': self.dflowfmmodel.external_forcings.initial_waterdepth_polygons.waterdepth[0]}
                    self._write_dict(f, dct, 'Global')

                with open(os.path.join(initcondpath,'initialFields.ini'),'a') as f:
                     dct  = {'quantity': 'waterDepth',
                             'unit': 'm',
                             'dataFileType': '1dField',
                             'dataFile': initcondfolder+'/initialWaterDepth.ini'
                             }
                     self._write_dict(f,dct,'Initial')
            else:
                for row in self.dflowfmmodel.external_forcings.initial_waterdepth_polygons.itertuples():
                    # Write pol file
                    write_fm_file(file=os.path.join(initcondpath, row.Index+".pol"),
                                      data=[np.vstack(row.geometry.exterior.coords[:])],
                                      names=[row.Index]
                                  )
                    with open(os.path.join(initcondpath,'initialFields.ini'),'a') as f:
                        dct  = {'quantity': 'waterDepth',
                                'unit': 'm',
                                'dataFileType': 'polygon',
                                'dataFile': f'{row.Index}.pol',
                                'interpolationMethod': 'constant',
                                'value': f'{row.waterdepth}'
                                }
                        self._write_dict(f,dct,'Initial')

    def objects_to_ldb(self, scalefactor=1.0):
        """
        Method to write weirs, bridges, pumps, universal weires, and cross sections to lbd
        """
        # Write cross sections
        coords = []
        names = []

        network = self.dflowfmmodel.network
        crosssections = self.dflowfmmodel.crosssections
        structures = self.dflowfmmodel.structures.as_dataframe(weirs=True, pumps=True, bridges=True, culverts=True, uweirs=True, orifices=True, compounds=True)

        for name, css in crosssections.crosssection_loc.items():
            sch_branch = network.schematised.at[css['branchid'], 'geometry']
            geo_branch = network.branches.at[css['branchid'], 'geometry']
            coords.append(orthogonal_line(
                line=sch_branch,
                offset=sch_branch.project(geo_branch.interpolate(css['chainage'])),
                width=10.0
            ))
            names.append(name)

        # # observation points
        # if len(self.dflowfmmodel.observation_points)>0:
        #     obs1d = self.dflowfmmodel.observation_points[self.dflowfmmodel.observation_points.locationType=='1d']

        #     for obspt in obs1d.itertuples():
        #         sch_branch = network.schematised.at[obspt.branch_id, 'geometry']
        #         geo_branch = network.branches.at[obspt.branch_id, 'geometry']
        #         offset = sch_branch.project(geo_branch.interpolate(obspt.branch_offset))
        #         names.append(obspt.code)
        #         coords.append(sch_branch.interpolate(offset).buffer(2.0).simplify(0.5).exterior.coords[:])

        if structures is not None:
            for structure in structures.itertuples():
                # compounds have no coordinates themselves
                if structure.type == 'compound':
                    continue
                sch_branch = network.schematised.at[structure.branchid, 'geometry']
                geo_branch = network.branches.at[structure.branchid, 'geometry']
                # Get offset on schematised branch
                offset = sch_branch.project(geo_branch.interpolate(structure.chainage))

                if structure.structype == 'pump':
                    # Get circle coordinates
                    coords.append(sch_branch.interpolate(offset).buffer(20.0).simplify(0.5).exterior.coords[:])
                    names.append(structure.Index)

                elif structure.structype == 'weir':
                    # Create orthogonal line
                    line = orthogonal_line(line=sch_branch, offset=offset, width=10.0)
                    # get coordinates for weir
                    ringcoords = LineString(line).parallel_offset(-2.0).coords[:] + LineString(line).parallel_offset(2.0).coords[::-1]
                    ringcoords += [ringcoords[0]]
                    coords.append(ringcoords)
                    names.append(structure.Index)

                elif structure.structype == 'bridge':
                    line1 = orthogonal_line(line=sch_branch, offset=offset-5., width=10.0)
                    line2 = orthogonal_line(line=sch_branch, offset=offset+5., width=10.0)
                    ringcoords = LineString(line1).coords[:] + LineString(line2).coords[::-1]
                    ringcoords += [ringcoords[0]]
                    coords.append(ringcoords)
                    names.append(structure.Index)

                elif structure.structype == 'uweir':
                    line1 = orthogonal_line(line=sch_branch, offset=offset-5., width=15.0)
                    line2 = orthogonal_line(line=sch_branch, offset=offset+5., width=1.0)
                    ringcoords = LineString(line1).coords[:] + LineString(line2).coords[::-1]
                    ringcoords += [ringcoords[0]]
                    coords.append(ringcoords)
                    names.append(structure.Index)

                elif structure.structype == 'orifice':
                    line1 = orthogonal_line(line=sch_branch, offset=offset+5., width=15.0)
                    line2 = orthogonal_line(line=sch_branch, offset=offset-5., width=1.0)
                    ringcoords = LineString(line1).coords[:] + LineString(line2).coords[::-1]
                    ringcoords += [ringcoords[0]]
                    coords.append(ringcoords)
                    names.append(structure.Index)

                elif structure.structype == 'culvert':
                    # Get circle coordinates
                    coords.append(sch_branch.interpolate(offset).buffer(5.0).simplify(0.5).exterior.coords[:])
                    names.append(structure.Index)

                elif structure.structype == 'orifice':
                    line1 = orthogonal_line(line=sch_branch, offset=offset, width=5.0)
                    line2 = orthogonal_line(line=LineString(line1), offset=offset, width=5.0)
                    coords.append(LineString(line1).coords[:]+LineString(line2).coords[:])
                    names.append(structure.Index)

                else:
                    raise TypeError('Structure type not recognized')

        # Write to file
        filename = 'structures.ldb'
        write_fm_file(os.path.join(self.output_dir, filename), data=coords, names=names)

        # Modify in mdu file
        self.mdu_parameters['landboundaryfile'] = filename

    def write_observation_points(self):
        """
        Methode to write observation points to file, and set
        the observation points file name in the mdu file.
        """
        if not self.dflowfmmodel.observation_points.empty:
            # Define file name
            filename = 'obspts_obs.ini'
            with open(os.path.join(self.output_dir,filename),'w') as f:
                self._write_header(f, filetype='obsPoints', fileversion=2.00)
                self.dflowfmmodel.observation_points.drop('geometry', axis=1, inplace=True)
                if self.dflowfmmodel.observation_points['locationType'].str.count("1d").sum()>0:
                    obs1d = self.dflowfmmodel.observation_points[self.dflowfmmodel.observation_points['locationType']=='1d']
                    if self.dflowfmmodel.observation_points['locationType'].str.count("2d").sum()>0:
                        obs1d.drop(['x','y'], axis=1, inplace=True)
                    obs1ddct = obs1d.to_dict(orient='row')
                    for dct in obs1ddct:
                        self._write_dict(f, dct=dct, header='ObservationPoint')
                if self.dflowfmmodel.observation_points['locationType'].str.count("2d").sum()>0:
                    obs2d = self.dflowfmmodel.observation_points[self.dflowfmmodel.observation_points['locationType']=='2d']
                    obs2d.drop(['branchId','chainage'], axis=1, inplace=True)
                    obs2ddct = obs2d.to_dict(orient='row')
                    for dct in obs2ddct:
                        self._write_dict(f, dct=dct, header='ObservationPoint')

            # Write to fm format file
            #obspt_data = [[row.geometry.x, row.geometry.y, row.name] for row in self.dflowfmmodel.observation_points.itertuples()]
            #write_fm_file(os.path.join(self.output_dir, filename), obspt_data)

            # Change mdu
            self.change_mdu_parameter('ObsFile', filename)
        else:
            self.change_mdu_parameter('ObsFile', '')

    def change_mdu_parameter(self, parameter, value):
        """
        Change parameter value in mdu file
        """
        with open(os.path.join(self.output_dir, f'{self.name}.mdu'), 'r') as f:
            lines = f.readlines()

        for i, line in enumerate(lines):
            if line.lower().startswith(parameter.lower()):
                items = re.split('=|#|$\s+', line)
                key, oldvalue = items[0], items[1]
                if not key.strip().lower() == parameter.lower():
                    continue
                lines[i] = line.replace('='+oldvalue, f'= {value}'.ljust(len(oldvalue)+1))
                break

        with open(os.path.join(self.output_dir, f'{self.name}.mdu'), 'w') as f:
            f.write(''.join(lines))
    
    def write_dimrfiles(self):        
        """
        Method to write the files required for running in DIMR
        """
        # run.bat
        with open(os.path.join(self.output_dir, 'run.bat'),'w') as f:
                 f.write('@ echo off\n')
                 f.write('set OMP_NUM_THREADS=2\n')
                 f.write('call \"'+self.run_dimrpad+'\"\n')
                 f.write('pause\n')

        # coupling XML
        FM_comp_name = 'DFM'        
        with open(os.path.join(self.output_dir, 'dimr_config.xml'),'w') as f:       
                 f.write('<?xml version="1.0" encoding="utf-8" standalone="yes"?>\n')
                 f.write('<dimrConfig xmlns="http://schemas.deltares.nl/dimr" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://schemas.deltares.nl/dimr http://content.oss.deltares.nl/schemas/dimr-1.2.xsd">\n')
                 f.write('\t<documentation>\n')
                 f.write('\t\t<fileVersion>1.2</fileVersion>\n')
                 f.write('\t\t<createdBy>D-HyDAMO in delft3dfmpy v.'+self.version['number']+'</createdBy>\n')
                 f.write("\t\t<creationDate>"+self.version['date']+'</creationDate>\n')
                 f.write('\t</documentation>\n\n')        
                 f.write('\t<control>\n')
                 f.write('\t\t\t\t<start name="'+FM_comp_name+'" />\n')
                 f.write('\t</control>\n\n')                                      
                 f.write('\t<component name="'+FM_comp_name+'">\n')
                 f.write('\t\t<library>dflowfm</library>\n')
                 f.write('\t\t<workingDir>.</workingDir>\n')
                 f.write('\t\t<inputFile>'+self.name+'.mdu</inputFile>\n')
                 f.write('\t</component>\n')                                  
                 f.write('</dimrConfig>\n')
        
    

def _format_row(row):
    # Create format string
    string = ''
    for element in row:
        if isinstance(element, str):
            string += element + ' '
        elif isinstance(element, int):
            string += '{:6d} '.format(element)
        elif isinstance(element, float):
            string += '{:12.6f} '.format(element)
    return string

def write_fm_file(file, data, names=None, mode='w'):
    string = ''
    if names is not None:
        for geometry, name in zip(data, names):
            string += '{}\n'.format(name)
            string += '{:6d} {:6d}\n'.format(*np.shape(geometry))
            string += '\n'.join(list(map(_format_row, geometry))) + '\n'

    else:
        string = '\n'.join(list(map(_format_row, data))) + '\n'

    with open(file, mode) as f:
        f.write(string)

