import logging
from itertools import product

import geopandas as gpd
import numpy as np
import pandas as pd
import PIL.Image
import PIL.ImageDraw
import rasterio
from matplotlib import path
from shapely import affinity
from shapely.geometry import MultiLineString, LineString, MultiPolygon, Polygon, MultiPoint, Point

from delft3dfmpy.core import checks
from delft3dfmpy.core.logging import ProgressLogger

logger = logging.getLogger(__name__)


def rotate_coordinates(origin, theta, xcrds, ycrds):
    """
    Rotate coordinates around origin (x0, y0) with a certain angle (radians)
    """
    x0, y0 = origin
    xcrds_rot = x0 + (xcrds - x0) * np.cos(theta) + (ycrds - y0) * np.sin(theta)
    ycrds_rot = y0 - (xcrds - x0) * np.sin(theta) + (ycrds - y0) * np.cos(theta)
    return xcrds_rot, ycrds_rot


def minimum_bounds_fixed_rotation(polygon, angle):
    """Get the minimum box for a polygon with a given axes rotation.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon
        Polygon that is rotated
    angle : int or float
        Rotation of the polygon in degrees

    Returns
    -------
    tuple
        Tuple with origin (x, y), xsize and ysize
    """
    # Determine spinning point
    spinpt = (polygon.envelope.bounds[0], polygon.envelope.bounds[1])

    # Rotate clip polygon with rotation, get envelope and rotate back.
    rotbox1 = affinity.rotate(polygon, angle=angle, origin=spinpt).envelope

    # Determine size of grid
    xsize = rotbox1.bounds[2] - rotbox1.bounds[0]
    ysize = rotbox1.bounds[3] - rotbox1.bounds[1]

    # Rotate again, and get origin
    rotbox2 = affinity.rotate(rotbox1, angle=-angle, origin=spinpt)
    origin = rotbox2.exterior.coords[0]

    return origin, xsize, ysize

def possibly_intersecting(dataframebounds, geometry, buffer=0):
    """
    Finding intersecting profiles for each branch is a slow process in case of large datasets
    To speed this up, we first determine which profile intersect a square box around the branch
    With the selection, the interseting profiles can be determines much faster.

    Parameters
    ----------
    dataframebounds : numpy.array
    geometry : shapely.geometry.Polygon
    """

    geobounds = geometry.bounds
    idx = (
        (dataframebounds[0] - buffer < geobounds[2]) &
        (dataframebounds[2] + buffer > geobounds[0]) &
        (dataframebounds[1] - buffer < geobounds[3]) &
        (dataframebounds[3] + buffer > geobounds[1])
    )
    # Get intersecting profiles
    return idx


def find_nearest_branch(branches, geometries, method='overal', maxdist=5):
    """
    Method to determine nearest branch for each geometry.
    The nearest branch can be found by finding t from both ends (ends) or the nearest branch from the geometry
    as a whole (overal), the centroid (centroid), or intersecting (intersect).

    Parameters
    ----------
    branches : geopandas.GeoDataFrame
        Geodataframe with branches
    geometries : geopandas.GeoDataFrame
        Geodataframe with geometries to snap
    method='overal' : str
        Method for determine branch
    maxdist=5 : int or float
        Maximum distance for finding nearest geometry
    minoffset : int or float
        Minimum offset from the end of the corresponding branch in case of method=equal
    """
    # Check if method is in allowed methods
    allowed_methods = ['intersecting', 'overal', 'centroid', 'ends']
    if method not in allowed_methods:
        raise NotImplementedError(f'Method "{method}" not implemented.')

    # Add columns if not present
    if 'branch_id' not in geometries.columns:
        geometries['branch_id'] = ''
    if 'branch_offset' not in geometries.columns:
        geometries['branch_offset'] = np.nan

    if method == 'intersecting':
        # Determine intersection geometries per branch
        geobounds = geometries.bounds.values.T
        for branch in branches.itertuples():
            selectie = geometries.loc[possibly_intersecting(geobounds, branch.geometry)].copy()
            intersecting = selectie.loc[selectie.intersects(branch.geometry).values]

            # For each geometrie, determine offset along branch
            for geometry in intersecting.itertuples():
                # Determine distance of profile line along branch
                geometries.at[geometry.Index, 'branch_id'] = branch.Index

                # Calculate offset
                branchgeo = branch.geometry
                mindist = min(0.1, branchgeo.length / 2.)
                offset = round(branchgeo.project(branchgeo.intersection(geometry.geometry).centroid), 3)
                offset = max(mindist, min(branchgeo.length - mindist, offset))
                geometries.at[geometry.Index, 'branch_offset'] = offset

    else:
        branch_bounds = branches.bounds.values.T
        # In case of looking for the nearest, it is easier to iteratie over the geometries instead of the branches
        for geometry in geometries.itertuples():
            # Find near branches
            nearidx = possibly_intersecting(branch_bounds, geometry.geometry, buffer=maxdist)
            selectie = branches.loc[nearidx]

            if method == 'overal':
                # Determine distances to branches
                dist = selectie.distance(geometry.geometry)
            elif method == 'centroid':
                # Determine distances to branches
                dist = selectie.distance(geometry.geometry.centroid)
            elif method == 'ends':
                # Since a culvert can cross a channel, it is
                crds = geometry.geometry.coords[:]
                dist = (selectie.distance(Point(*crds[0])) + selectie.distance(Point(*crds[-1]))) * 0.5

            # Determine nearest
            if dist.min() < maxdist:
                branchidxmin = dist.idxmin()
                geometries.at[geometry.Index, 'branch_id'] = dist.idxmin()
                if isinstance(geometry.geometry, Point):
                    geo = geometry.geometry
                else:
                    geo = geometry.geometry.centroid

                # Calculate offset
                branchgeo = branches.at[branchidxmin, 'geometry']
                mindist = min(0.1, branchgeo.length / 2.)
                offset = max(mindist, min(branchgeo.length - mindist, round(branchgeo.project(geo), 3)))
                geometries.at[geometry.Index, 'branch_offset'] = offset

def orthogonal_line(line, offset, width=1.0):
    """
    Parameters
    ----------
    line : shapely.geometry.LineString
        Line geometry object on which the orthogonal line is drawn
    offset : float
        Offset of the orthogonal line along line
    width : float
        Width of the orthogonal line

    Returns
    -------
    line : list
        List with coordinate tuples
    """

    # Determine angle at offset
    angle = np.angle(complex(*np.diff([
        line.interpolate(offset - 0.001).coords[0][:2],
        line.interpolate(offset + 0.001).coords[0][:2]
    ], axis=0)[0])) + 0.5 * np.pi

    # Create new line
    pt = line.interpolate(offset).coords[0]

    f = 0.5 * width
    line = [(pt[0] + np.cos(angle) * f, pt[1] + np.sin(angle) * f), (pt[0] - np.cos(angle) * f, pt[1] - np.sin(angle) * f)]

    return line

def extend_linestring(line, near_pt, length):

    # Get the nearest end
    nearest_end = (0, 1) if line.project(near_pt) < line.length / 2 else (-1, -2)

    # Extrapolate the end 1 meter, and create a perpendicular line
    coords = line.coords[:]
    x0, y0 = coords[nearest_end[0]]
    dx, dy = np.diff(np.vstack([coords[nearest_end[0]], coords[nearest_end[1]]]), axis=0)[0]
    segmentlength = (dx**2 + dy**2)**0.5
    dx /= (segmentlength * length)
    dy /= (segmentlength * length)

    # if nearest_end[0] == 0:
    #     coords = LineString[(x0 - dx, y0 - dy)] + coords
    # else:
    #     coords = coords + [(x0 - dx, y0 - dy)]

    return LineString([(x0, y0), (x0 - dx, y0 - dy)])

def points_in_polygon(points, polygon):
    """
    Determine points that are inside a polygon, taking
    holes into account.

    Parameters
    ----------
    points : numpy.array
        Nx2 - array
    polygon : shapely.geometry.Polygon
        Polygon (can have holes)
    """
    # First select points in square box around polygon
    ptx, pty = points.T
    mainindex = possibly_intersecting(
        dataframebounds=np.c_[[ptx, pty, ptx, pty]], geometry=polygon)
    boxpoints = points[mainindex]

    extp = path.Path(polygon.exterior)
    intps = [path.Path(interior) for interior in polygon.interiors]

    # create first index. Everything within exterior is True
    index = extp.contains_points(boxpoints)

    # set points in holes also to nan
    if intps:
        subset = boxpoints[index]
        # Start with all False
        subindex = np.zeros(len(subset), dtype=bool)

        for intp in intps:
            # update mask, set to True where point in interior
            subindex = subindex | intp.contains_points(subset)

        # Everything within interiors should be True
        # So, set everything within interiors (subindex == True), to True
        index[np.where(index)[0][subindex]] = False

    # Set index in main index to False
    mainindex[np.where(mainindex)[0][~index]] = False

    return mainindex

class RasterPart:

    def __init__(self, f, xmin, ymin, xmax, ymax):
        self.f = f
        # Indices, not coordinates
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        self.set_window()

        self.get_corners()

    @classmethod
    def from_bounds(cls, f, bnds):
        idxs = list(f.index(bnds[0], bnds[1]))[::-1] + list(f.index(bnds[2], bnds[3]))[::-1]
        return cls(f, min(idxs[0], idxs[2]), min(idxs[1], idxs[3]), max(idxs[0], idxs[2]), max(idxs[1], idxs[3]))

    def set_window(self):
        self.window = ((self.ymin, self.ymax), (self.xmin, self.xmax))
        self.shape = (self.ymax - self.ymin, self.xmax - self.xmin)

    def get_corners(self):
        x0 = self.f.xy(self.ymax, self.xmin)[0]
        x1 = self.f.xy(self.ymax, self.xmax)[0]
        y0 = self.f.xy(self.ymin, self.xmax)[1]
        y1 = self.f.xy(self.ymax, self.xmax)[1]

        self.lowerleft = (min(x0, x1), min(y0, y1))
        self.upperright = (max(x0, x1), max(y0, y1))

    def get_xy_range(self):
        x = np.linspace(self.lowerleft[0], self.upperright[0], (self.xmax - self.xmin), endpoint=False)
        y = np.linspace(self.lowerleft[1], self.upperright[1], (self.ymax - self.ymin), endpoint=False)[::-1]
        #TODO: FIX y-DIRECTION
        return x, y

    def read(self, layeridx):
        return self.f.read(layeridx, window=self.window)

    def get_pts_in_part(self, pts, buffer=0):
        self.get_corners()
        # Select points within part + buffer
        idx = (
            (pts[:, 0] > self.lowerleft[0] - buffer) &
            (pts[:, 0] < self.upperright[0] + buffer) &
            (pts[:, 1] > self.lowerleft[1] - buffer) &
            (pts[:, 1] < self.upperright[1] + buffer)
        )
        return idx

    def get_mask(self, polygon):

        valid = geometry_to_mask(polygon, self.lowerleft, abs(self.f.transform.a), self.shape)
        return valid

def geometry_to_mask(polygons, lowerleft, cellsize, shape):

    # Initialize mask
    mask = np.zeros(shape)

    for polygon in as_polygon_list(polygons):
        # Create from exterior
        mask += get_mask(polygon.exterior, lowerleft, cellsize, shape)
        # Subtract interiors
        for interior in polygon.interiors:
            mask -= get_mask(interior, lowerleft, cellsize, shape, outline=0)

    mask = (mask == 1)

    return mask


def get_mask(linestring, lowerleft, cellsize, shape, outline=1):

    # Create array from coordinate sequence
    path = np.vstack(linestring.coords[:])

    # Convert to (0,0) and step size 1
    path[:, 0] -= lowerleft[0]
    path[:, 1] -= lowerleft[1] + cellsize
    path /= cellsize
    # Convert from array to tuple list
    path = list(zip(*zip(*path)))

    # Create mask
    maskIm = PIL.Image.new('L', (shape[1], shape[0]), 0)
    PIL.ImageDraw.Draw(maskIm).polygon(path, outline=outline, fill=1)
    mask = np.array(maskIm)[::-1]

    return mask


def raster_in_parts(f, ncols, nrows, facedata=None):
    """
    Certain rasters are too big to read into memory at once.
    This function helps splitting them in equal parts of (+- ncols x nrows pixels)

    If facedata is given, each part is extended such that whole faces
    are covered by the parts
    """
    nx = max(1, f.shape[1] // ncols)
    ny = max(1, f.shape[0] // nrows)

    xparts = np.linspace(0, f.shape[1], nx+1).astype(int)
    yparts = np.linspace(0, f.shape[0], ny+1).astype(int)

    pts = facedata[['facex', 'facey']].values

    parts = []
    pl = ProgressLogger(logger, nx*ny, 10)

    for i, (ix, iy) in enumerate(product(range(nx), range(ny))):
        pl.set_step(i)

        part = RasterPart(f, xmin=xparts[ix], ymin=yparts[iy], xmax=xparts[ix+1], ymax=yparts[iy+1])

        if facedata is not None:
            # For each part, get points in part.
            # For narrow/diagonal shapes the points in part can be limited
            idx = part.get_pts_in_part(pts)
            if not idx.any():
                continue

            crds = facedata['crds']
            ll = list(zip(*[crds[i].min(axis=0) for i in np.where(idx)[0] + 1]))
            ur = list(zip(*[crds[i].max(axis=0) for i in np.where(idx)[0] + 1]))
            bounds = (min(ll[0]), min(ll[1]), max(ur[0]), max(ur[1]))

            # Get new part based on extended bounds
            part = RasterPart.from_bounds(f, bounds)

            # Add the cell centers within the window as index to the part
            part.idx = idx

        yield part

def rasterize_cells(facedata, prt):

    # Initialize mask
    # Create mask
    maskIm = PIL.Image.new('I', (prt.shape[1], prt.shape[0]), 0)
    todraw = PIL.ImageDraw.Draw(maskIm)

    cellnumber = np.zeros(prt.shape)
    cellsize = abs(prt.f.transform.a)

    for row in facedata.itertuples():

        # Create array from coordinate sequence
        path = row.crds.copy()
        # Convert to (0,0) and step size 1
        path[:, 0] -= (prt.lowerleft[0] - 0.5 * cellsize)
        path[:, 1] -= (prt.lowerleft[1] + 0.5 * cellsize)
        path /= cellsize
        # Convert from array to tuple list
        path = list(zip(*zip(*path)))

        # Create mask
        todraw.polygon(path, outline=row.Index, fill=row.Index)

    return np.array(maskIm, dtype=np.int32)[::-1]

def check_geodateframe_rasterstats(facedata):
    """
    Check for type, columns and coordinates
    """
    if not isinstance(facedata, gpd.GeoDataFrame):
        raise TypeError('facedata should be type GeoDataFrame')

    # Check if facedata has required columns
    if ('facex' not in facedata.columns) or ('facey' not in facedata.columns):
        xy = list(zip(*[pt.coords[0] for pt in facedata.geometry.centroid]))
        facedata['facex'] = xy[0]
        facedata['facey'] = xy[1]

    # Check if coordinates are present.
    if 'crds' not in facedata.columns:
        facedata['crds'] =[row.coords[:] for row in facedata.geometry]


def raster_stats_fine_cells(rasterpath, facedata, stats=['mean']):
    """
    Get raster stats

    Parameters
    ----------
    rasterpath : str
        Path to raster file
    facedata : geopandas.GeoDataFrame
        Dataframe with polygons in which the raster statistics are derived.
    stats : list
        List of statistics to retrieve. Should all be present as functions that
        can be applied to pandas group by
    """

    # Create empty array for stats
    stat_array = {stat: {} for stat in stats + ['count']}

    # Check geometries
    check_geodateframe_rasterstats(facedata)

    # Open raster file
    first = True
    i = 0
    with rasterio.open(rasterpath, 'r') as f:

        # Split file in parts based on shape
        parts = raster_in_parts(f, ncols=250, nrows=250, facedata=facedata)

        for prt in parts:

            # Get values from array
            arr = prt.read(1)
            valid = (arr != f.nodata)
            if not valid.any():
                continue

            # Rasterize the cells in the part
            cellidx_sel = rasterize_cells(facedata.loc[prt.idx], prt)
            cellidx_sel[~valid] = 0
            valid = (cellidx_sel != 0)

            for cell_idx in np.unique(cellidx_sel):
                if cell_idx == 0:
                    continue

                # Mask
                cellmask = (cellidx_sel == cell_idx)
                if not cellmask.any():
                    continue

                # Get bottom values
                bottom = arr[cellmask]

                # For each statistic, get the values and add
                for stat in stats:
                    stat_array[stat][cell_idx] = getattr(np, stat)(bottom)

                stat_array['count'][cell_idx] = len(bottom)

        # Cast to pandas dataframe
        df = pd.DataFrame.from_dict(stat_array).reindex(index=facedata.index)

    return df

def waterdepth_ahn(dempath, facedata, outpath, column):
    """
    Function that combines a dem and water levels to a water
    depth raster. No sub grid correction is done.

    Parameters
    ----------
    dempath : str
        Path to raster file with terrain level
    facedata : gpd.GeoDataFrame
        GeoDataFrame with at least the cell geometry and
        a column with water levels
    outpath : str
        Path to output raster file
    column : str
        Name of the column with the water level data
    """

    # Open raster file
    with rasterio.open(dempath, 'r') as f:
        first = True
        out_meta = f.meta.copy()

        cell_area = abs(f.transform.a * f.transform.e)

        # Split file in parts based on shape
        parts = raster_in_parts(f, ncols=250, nrows=250, facedata=facedata)

        for prt in parts:

            # Get values from array
            arr = prt.read(1)
            valid = (arr != f.nodata)
            if not valid.any():
                continue

            cellidx_sel = rasterize_cells(facedata.loc[prt.idx], prt)
            cellidx_sel[~valid] = 0
            valid = (cellidx_sel != 0)

            # Create array to assign water levels
            wlev_subgr = np.zeros(cellidx_sel.shape, dtype=out_meta['dtype'])

            for cell_idx in np.unique(cellidx_sel):
                if cell_idx == 0:
                    continue
                # Mask
                cellmask = (cellidx_sel == cell_idx)
                if not cellmask.any():
                    continue
                # Add values for cell to raster
                wlev_subgr[cellmask] = facedata.at[cell_idx, column]

            # Write to output raster
            with rasterio.open(outpath, 'w' if first else 'r+', **out_meta) as dst:
                # Determine water depth
                if not first:
                    wdep_subgr = dst.read(1, window=prt.window)
                else:
                    wdep_subgr = np.ones_like(wlev_subgr) * f.nodata
                    first = False
                wdep_subgr[valid] = np.maximum(wlev_subgr - arr, 0).astype(out_meta['dtype'])[valid]
                dst.write(wdep_subgr[None, :, :], window=prt.window)

    compress(outpath)


def compress(path):
    """
    Function re-save an existing raster file with compression.

    Parameters
    ----------
    path : str
        Path to raster file. File is overwritten with compress variant.
    """
    # Compress
    with rasterio.open(path, 'r') as f:
        arr = f.read()
        out_meta = f.meta.copy()
        out_meta['compress'] = 'deflate'
    with rasterio.open(path, 'w', **out_meta) as f:
        f.write(arr)

def as_geometry_list(geometry, singletype, multitype):
    """Convenience method to return a list with one or more

Polygons/LineString/Point from a given Polygon/LineString/Point
    or MultiPolygon/MultiLineString/MultiPoint.

    Parameters
    ----------
    polygon : list or Polygon or MultiPolygon
        Object to be converted

    Returns
    -------
    list
        list of Polygons
    """
    if isinstance(geometry, singletype):
        return [geometry]
    elif isinstance(geometry, multitype):
        return [p for p in geometry]
    elif isinstance(geometry, list):
        lst = []
        for item in geometry:
            lst.extend(as_geometry_list(item, singletype, multitype))
        return lst
    else:
        raise TypeError(f'Expected {singletype} or {multitype}. Got "{type(geometry)}"')
        
def as_linestring_list(linestring):
    return as_geometry_list(linestring, LineString, MultiLineString)

def as_polygon_list(polygon):
    return as_geometry_list(polygon, Polygon, MultiPolygon)

def as_point_list(point):
    return as_geometry_list(point, Point, MultiPoint)