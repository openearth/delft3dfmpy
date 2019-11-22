# coding: latin-1
import os
import re
import shutil
import datetime

import numpy as np
from shapely.geometry import box, LineString, Point

from delft3dfmpy.io.UgridWriter import UgridWriter
from delft3dfmpy.core.geometry import orthogonal_line

import logging

logger = logging.getLogger(__name__)

class DFlowFMWriter:
    """Writer for FM files"""

    def __init__(self, dflowfmmodel, output_dir, name):
        self.dflowfmmodel = dflowfmmodel
        # self.geometries = parent.geometries
        # self.boundary_conditions = parent.boundary_conditions
        self.output_dir = output_dir
        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)
        self.name = name


        self.mdu_parameters = dflowfmmodel.mdu_parameters
        self.mdu_parameters['NetFile'] = self.name + "_net.nc"
        self.mdu_parameters['ExtForceFile'] = self.name + ".ext"

        self.mdufile = os.path.join(output_dir, self.name +'.mdu')
        self.extfile = os.path.join(output_dir, self.mdu_parameters['ExtForceFile'])
        self.netfile = os.path.join(output_dir, self.mdu_parameters['NetFile'])

    def write_all(self):  # write all fm files from HyDAMO

        if os.path.exists(self.extfile):
            os.remove(self.extfile)

        # Write grid in Ugrid format to netcdf
        ugridWriter = UgridWriter()
        ugridWriter.write(dflowfmmodel=self.dflowfmmodel, path=self.netfile)

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

        # Change mdu parameters
        for parameter, value in self.mdu_parameters.items():
            self.change_mdu_parameter(parameter, value)

        # Write laterals
        if any(self.dflowfmmodel.external_forcings.laterals):
            self.write_laterals()

    def _write_dict(self, f, dct, header):
        """
        Write a dictionary to file, without transforming with mapping
        """
        # Write header
        f.write(f'[{header}]\n')
        for key, value in dct.items():
            f.write(f'{key} = {value}\n')
        f.write('\n')

    def writeMdFiles(self):
        srcMdu = os.path.join(os.path.dirname(__file__), '..', 'resources', 'FMmdu.txt')
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
            file = os.path.join(roughnessdir, f'r{i+1:03d}.ini')
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
            self._write_header(f, 'crossLoc', fileversion=2.01)
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
                        'definition': self.dflowfmmodel.crosssections.default_definition
                    }
                    # Write the cross section item
                    self._write_dict(f, dct=dct, header='crosssection')

    def _write_header(self, f, filetype, fileversion):
        f.write('[General]\n')
        f.write(f'fileVersion = {fileversion:.2f}\n')
        f.write(f'fileType = {filetype}\n')
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

            # Culverts
            if any(self.dflowfmmodel.structures.culverts):
                for _, dct in self.dflowfmmodel.structures.culverts.items():
                    self._write_dict(f, dct=dct, header='structure')

            # Weirs
            if any(self.dflowfmmodel.structures.weirs):
                for _, dct in self.dflowfmmodel.structures.weirs.items():
                    self._write_dict(f, dct=dct, header='structure')

            # Pumps
            if any(self.dflowfmmodel.structures.pumps):
                branches = self.dflowfmmodel.network.branches

                for pump_id, dct in self.dflowfmmodel.structures.pumps.items():
                    self._write_dict(f, dct=dct, header='structure')

                    # Write pli line. This line needs to cross the schematized branch, not the original geometry
                    nearest_pt = branches.at[dct['branchid'], 'geometry'].interpolate(dct['chainage'])
                    branch = self.dflowfmmodel.network.schematised.at[dct['branchid'], 'geometry']
                    line = orthogonal_line(line=branch, offset=branch.project(nearest_pt), width=0.1)
                    write_fm_file(os.path.join(self.output_dir, dct['locationfile']), [np.vstack(line)], [pump_id])

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
        Write laterals to files.
        """
        external_forcings = self.dflowfmmodel.external_forcings

        if not 'refdate' in self.mdu_parameters.keys():
            raise KeyError('Specify reference date before writing laterals.')

        # Create directory for laterals if it does not yet exist
        if not os.path.exists(os.path.join(self.output_dir, 'laterals')):
            os.mkdir(os.path.join(self.output_dir, 'laterals'))

        for name, dct in external_forcings.laterals.items():
            
            # Get coordinates of polygon around point
            point = Point(dct['x'], dct['y'])
            crds = np.vstack(box(*point.buffer(1.0).bounds).exterior.coords[:])

            # Write polygon to boundaries folder
            write_fm_file(os.path.join(self.output_dir, 'laterals', f'lateral_{name}.pol'), [crds], [name])

            # Determine starttime from mdu parameters
            starttime = datetime.datetime.strptime(str(self.mdu_parameters['refdate']), '%Y%m%d')
            starttime += datetime.timedelta(seconds=int(self.mdu_parameters['tstart']))

            # Write time series
            series = dct['timeseries'].copy()
            series.index -= starttime
            data = np.c_[series.index.total_seconds().values / 60, series.values]

            write_fm_file(
                file=os.path.join(self.output_dir, 'laterals', f"lateral_{name}_0001.tim"),
                data=data,
            )

            # Add to external forcing file
            with open(self.extfile, 'a') as dst:
                dst.write((
                    f"QUANTITY=lateraldischarge1d\n"
                    f"FILENAME=laterals/lateral_{name}.pol\n"
                    f"FILETYPE=10\n"
                    f"METHOD=1\n"
                    f"OPERAND=O\n\n"
                ))

    def write_boundary_conditions(self):
        """
        Write boundary conditions to ext file
        """

        # Copy external forcing file (*.ext)
        if not os.path.exists(self.extfile):
            shutil.copy2(os.path.join(os.path.dirname(__file__), 'resources', 'template.ext'), self.extfile)

        # For all boundary conditions
        for bc in self.dflowfmmodel.external_forcings.boundaries.itertuples():
            # Name
            name = f"{bc.bctype}_{bc.code}"

            # Get time series
            if bc.time is None:
                data = [[0, bc.value], [999999, bc.value]]
            else:
                data = list(zip(bc.time, bc.value))

            # Write pli
            if bc.filetype == 9:
                write_fm_file(
                    file=os.path.join(self.output_dir, name+".pli"),
                    data=[np.vstack(bc.geometry.coords[:])],
                    names=[name]
                )
                filename = f"FILENAME={name}.pli"
                # Write time series for pli (point 1)
                write_fm_file(file=os.path.join(self.output_dir, f"{name}_0001.tim"), data=data)

            elif bc.filetype == 1:
                # Write time series
                write_fm_file(file=os.path.join(self.output_dir, f"{bc.code}.tim"), data=data)
                filename = f"FILENAME={bc.code}.tim"

            else:
                raise NotImplementedError()

            # Add to ext file
            with open(self.extfile, 'a') as f:
                f.write((
                    f"QUANTITY={bc.bctype}\n"
                    f"{filename}\n"
                    f"FILETYPE={bc.filetype}\n"
                    f"METHOD={bc.method}\n"
                    f"OPERAND={bc.operand}\n\n"
                ))

    def write_initial_conditions(self):

        # Copy external forcing file (*.ext)
        if not os.path.exists(self.extfile):
            shutil.copy2(os.path.join(os.path.dirname(__file__), '..', 'resources', 'template.ext'), self.extfile)

        # Create folder for initial conditions if it does not exist yet
        initcondfolder = 'initialconditions'
        initcondpath = os.path.join(self.output_dir, initcondfolder)
        if not os.path.exists(initcondpath):
            os.mkdir(initcondpath)

        # Write water levels within polygons
        for row in self.dflowfmmodel.external_forcings.initial_waterlevel_polygons.itertuples():

            # Add to ext file
            with open(self.extfile, 'a') as f:
                f.write((
                    "QUANTITY=initialwaterlevel\n"
                    f"FILENAME={initcondfolder}\{row.Index}.pol\n"
                    "FILETYPE=10\n"
                    "METHOD=4\n"
                    "OPERAND=O\n"
                    f"VALUE={row.waterlevel}\n\n"
                ))

            # Write pol file
            write_fm_file(
                file=os.path.join(initcondpath, row.Index+".pol"),
                data=[np.vstack(row.geometry.exterior.coords[:])],
                names=[row.Index]
            )

        # Write water levels at xyz, if the water depth function is called (and xyz thus)
        if any(self.dflowfmmodel.external_forcings.initial_waterlevel_xyz):
            write_fm_file(
                file=os.path.join(self.output_dir, self.mdu_parameters['WaterLevIniFile']),
                data=self.dflowfmmodel.external_forcings.initial_waterlevel_xyz
            )

    def objects_to_ldb(self, scalefactor=1.0):
        """
        Method to write weirs, pumps and cross sections to lbd
        """
        # Write cross sections
        coords = []
        names = []

        network = self.dflowfmmodel.network
        crosssections = self.dflowfmmodel.crosssections
        structures = self.dflowfmmodel.structures.as_dataframe(weirs=True, pumps=True, culverts=True)

        for name, css in crosssections.crosssection_loc.items():
            sch_branch = network.schematised.at[css['branchid'], 'geometry']
            geo_branch = network.branches.at[css['branchid'], 'geometry']
            coords.append(orthogonal_line(
                line=sch_branch,
                offset=sch_branch.project(geo_branch.interpolate(css['chainage'])),
                width=10.0
            ))
            names.append(name)
        
        if structures is not None:
            for structure in structures.itertuples():

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

                elif structure.structype == 'culvert':
                    # Get circle coordinates
                    coords.append(sch_branch.interpolate(offset).buffer(5.0).simplify(0.5).exterior.coords[:])
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
            filename = 'obspts_obs.xyn'
            # Write to fm format file
            obspt_data = [[row.geometry.x, row.geometry.y, row.name] for row in self.dflowfmmodel.observation_points.itertuples()]
            write_fm_file(os.path.join(self.output_dir, filename), obspt_data)
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
