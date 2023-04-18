from collections import namedtuple
from itertools import product
from scipy.optimize import brentq

import netCDF4
from matplotlib.collections import LineCollection, PatchCollection
from matplotlib.patches import Polygon as PolygonPatch
from shapely.geometry import Polygon, LineString
import numpy as np
import geopandas as gpd
import rasterio
import tqdm

import PIL.Image
import PIL.ImageDraw


class Results:
    def __init__(self, ncfile):
        self.file = ncfile
        self.open()
        self.variables = list(self.nc.variables.keys())
        self.close()
        self.network = Network(self)
        
    def open(self):
        self.nc = netCDF4.Dataset(self.file)
    
    def close(self):
        self.nc.close()
        
    def get_data(self, variable):
        self._check_variable(variable)
        self.open()
        values = self.nc.variables[variable][:]
        self.close()
        return values
        
    def get_Mesh2d_faces(self, as_polygon=False, as_patch=False):
        self.open()
        cells = np.dstack([self.nc.variables['Mesh2d_face_x_bnd'], self.nc.variables['Mesh2d_face_y_bnd']])
        nodata = (cells == self.nc.variables['Mesh2d_face_x_bnd']._FillValue)
        cells = [cell[idx] for idx, cell in zip(nodata.sum(axis=-1) != 2, cells.data)]
        self.close()
        if as_polygon:
            cells = [Polygon(cell) for cell in cells]
        elif as_patch:
            cells = [PolygonPatch(cell, closed=True) for cell in cells]
            
        return cells
    
    def get_mesh1d_edges(self, as_linestring=False):
        self.open()
        self.network.load_1d_network()

        edges = self.network.edges1d
        if as_linestring:
            edges = [LineString(edge) for edge in edges]
            
        return edges
    
    def get_Mesh2d_face_centers(self):
        self.open()
        coords = np.c_[self.nc.variables['Mesh2d_face_x'][:].data, self.nc.variables['Mesh2d_face_y'][:].data]
        self.close()
        return coords
    
    def get_max(self, variable, axis=0):
        return self.get_stat(variable, 'max', axis)

    def get_stat(self, variable, stat, axis=0):
        values = self.get_data(variable).data
        statvalue = getattr(np, stat)(values, axis=axis)
        return statvalue
        
    def get_network_count(self):
        self.open()
        counts = {
            'mesh1d': results.nc.dimensions['nmesh1d_node'].size,
            'Mesh2d': results.nc.dimensions['nMesh2d_face'].size,
            'links1d2d': results.nc.dimensions['nlink1d2d_connections'].size
        }
        self.close()
        return counts

    def _check_variable(self, variable):
        if variable not in self.variables:
            raise KeyError(f'Variable "{variable}" not in variables. Pick from: {", ".join(self.variables)}.')    

    def save_subgrid_results_fourrier(self, variable, statstep, ahnpath, outpath, shapeout=None, axis=0):
        """Method to save model results with subgrid correction. The variable should thus be
        the water depth, since this can be subgrid corrected.
        
        Parameters
        ----------
        variable : str
            Variable to rasterize, should be water depth or something similar.
        statstep : str or int
            Statistic (str, should be in numpy) or time step (int) to rasterize
        ahnpath : str
            Path to DEM from which the terrain level is read
        outpath : str
            Path to file to which the raster is written
        shapeout : path, optional
            If given (not None), the shapefile with 2D cells and the derived statistic/step is written to this file, by default None
        axis : int, optional
            Axis for which the statistic is derives, by default 0 since 'time' is on the first axis
        """

        # Check given variable
        self._check_variable(variable)
        
        xy = self.get_Mesh2d_face_centers()

        faces = self.get_Mesh2d_faces()
        facedata = gpd.GeoDataFrame(
            index=np.arange(len(xy), dtype=np.uint32) + 1,
            data=xy,
            columns=['x', 'y'],
            geometry=[Polygon(face) for face in faces],
            crs='epsg:28992'
        )
        facedata['crds'] = faces
        
        # Get water levels and water depths
#         if isinstance(statstep, str):
#             facedata['value'] = self.get_stat(variable, statstep, axis)
#         elif isinstance(statstep, int):
#             facedata['value'] = np.take(self.get_data(variable).data, statstep, axis=axis)
#         return(facedata['value'])
        facedata['value']=self.get_data(variable)
        facedata.value[facedata.value>30]=0
        #return(facedata)
        # Get waterdepth taking into account subgrid
        get_subgrid_waterdepth(ahnpath, facedata, outpath)
        
        # Save as shapefile (optional)
        if shapeout is not None:
            facedata.drop('crds', axis=1).to_file(shapeout)            
            
    def save_subgrid_results(self, variable, statstep, ahnpath, outpath, shapeout=None, axis=0):
        """Method to save model results with subgrid correction. The variable should thus be
        the water depth, since this can be subgrid corrected.
        
        Parameters
        ----------
        variable : str
            Variable to rasterize, should be water depth or something similar.
        statstep : str or int
            Statistic (str, should be in numpy) or time step (int) to rasterize
        ahnpath : str
            Path to DEM from which the terrain level is read
        outpath : str
            Path to file to which the raster is written
        shapeout : path, optional
            If given (not None), the shapefile with 2D cells and the derived statistic/step is written to this file, by default None
        axis : int, optional
            Axis for which the statistic is derives, by default 0 since 'time' is on the first axis
        """

        # Check given variable
        self._check_variable(variable)
        
        xy = self.get_Mesh2d_face_centers()

        faces = self.get_Mesh2d_faces()
        facedata = gpd.GeoDataFrame(
            index=np.arange(len(xy), dtype=np.uint32) + 1,
            data=xy,
            columns=['x', 'y'],
            geometry=[Polygon(face) for face in faces],
            crs='epsg:28992'
        )
        facedata['crds'] = faces
        
        # Get water levels and water depths
        if isinstance(statstep, str):
            facedata['value'] = self.get_stat(variable, statstep, axis)
        elif isinstance(statstep, int):
            facedata['value'] = np.take(self.get_data(variable).data, statstep, axis=axis)
        #return(facedata['value'])
        # Get waterdepth taking into account subgrid
        get_subgrid_waterdepth(ahnpath, facedata, outpath)
        
        # Save as shapefile (optional)
        if shapeout is not None:
            facedata.drop('crds', axis=1).to_file(shapeout)

    def save_rasterized_results(self, variable, statstep, rasterpath, outpath, axis=0):
        """Method to save model results with on raster. An example raster should be given to
        get the dimension and affine.
        
        Parameters
        ----------
        variable : str
            Variable to rasterize, should be water depth or something similar.
        statstep : str or int
            Statistic (str, should be in numpy) or time step (int) to rasterize
        rasterpath : str
            Path to raster from which the shape and affine is derived.
        outpath : str
            Path to file to which the raster is written
        axis : int, optional
            Axis for which the statistic is derives, by default 0 since 'time' is on the first axis
        """
        # Check given variable
        self._check_variable(variable)
        
        xy = self.get_Mesh2d_face_centers()

        faces = self.get_Mesh2d_faces()
        facedata = gpd.GeoDataFrame(
            index=np.arange(len(xy), dtype=np.uint32) + 1,
            data=xy,
            columns=['x', 'y'],
            geometry=[Polygon(face) for face in faces],
            crs='epsg:28992'
        )
        facedata['crds'] = faces
        
        # Get water levels and water depths
        if isinstance(statstep, str):
            facedata['value'] = self.get_stat(variable, statstep, axis)
        elif isinstance(statstep, int):
            facedata['value'] = np.take(self.get_data(variable).data, statstep, axis=axis)
        
        # Save the results without subgrid
        rasterize_results(rasterpath, facedata, outpath)
        
class Network:
    
    def __init__(self, parent):
        self.get_data = parent.get_data
        self.loaded_1d = False
        
    def load_1d_network(self):
        
        # Collect 1d network branches
        geomx = self.get_data('network_geom_x')
        geomy = self.get_data('network_geom_y')

        # Group by branches
        i = 0
        self.branches = []
        for nnodes in self.get_data('network_node_count'):
            self.branches.append(LineString([(geomx[j], geomy[j]) for j in range(i, i + nnodes)]))
            i += nnodes

        # Get id's and offsets
        branchids = self.get_data('mesh1d_nodes_branch_id')
        offsets = self.get_data('mesh1d_nodes_branch_offset')
        
        # Interpolate offsets on branches to get mesh1d nodes
        crds = []
        for i in range(len(branchids)):
            crds.append(self.branches[branchids[i]].interpolate(offsets[i]).coords[0])

        # Lees 1D links uit
        self.nodes1d = np.array(crds)
        self.edges1d = self.nodes1d[self.get_data('mesh1d_edge_nodes').data - 1]
        
        # Set flag loaded_1d to True
        self.loaded_1d = True
        
    def get_nearest(self, xy, type='node', one_based=False):
        
        # Load the network if it is not loaded yet
        if not self.loaded_1d:
            self.load_1d_network()
        
        # Get nearest node
        if type == 'node':
            imin = np.argmin(np.hypot(*(self.nodes1d - xy).T)) + int(one_based)
        
        elif type == 'edge':
            dist = np.hypot(*(self.edges1d[:, 0, :] - xy).T) + np.hypot(*(self.edges1d[:, 1, :] - xy).T)
            imin = np.argmin(dist) + int(one_based)
            
        else:
            raise ValueError(f'Type "{type}" not recognized. Specify node or edge.')
            
        return imin

def rasterize_results(rasterpath, facedata, outpath):
    """
    Get raster stats
    
    Parameters
    ----------
    dempath : str
        Path to raster file with terrain level
    outpath : str
        Path to output raster file
    pts : Nx3 array
        Array with x, y and z of results to interpolate
    buffer : float
        Max size between different points. This length is used to include points outside
        domain when splitting raster
    """

    # Check geometries
    check_geodateframe_rasterstats(facedata)

    # Open raster file
    with rasterio.open(rasterpath, 'r') as f:
        first = True
        out_meta = f.meta.copy()

        # Split file in parts based on shape
        parts = raster_in_parts(f, ncols=250, nrows=250, facedata=facedata)

        for prt in parts:

            # Get values from array
            cellidx_sel = rasterize_cells(facedata.loc[prt.idx], prt)

            # If no cells intersect the part, continue
            if (cellidx_sel == 0).all():
                continue
            
            # Create array to assign water levels
            results_subgr = np.full(cellidx_sel.shape, f.nodata, dtype=out_meta['dtype'])
            
            for cell_idx in np.unique(cellidx_sel):
                if cell_idx == 0:
                    continue
                # Mask
                cellmask = (cellidx_sel == cell_idx)
                if not cellmask.any():
                    continue
                # Add values for cell to raster
                results_subgr[cellmask] = facedata.at[cell_idx, 'value']
            
            # Write to output raster
            with rasterio.open(outpath, 'w' if first else 'r+', **out_meta) as dst:

                # Read already assigned results and add
                if not first:
                    old = dst.read(1, window=prt.window)
                    idx = old != f.nodata
                    results_subgr[idx] = old[idx]
                
                dst.write(results_subgr[None, :, :], window=prt.window)
                first = False
                
    compress(outpath)

def compress(path):
    # Compress
    with rasterio.open(path, 'r') as f:
        arr = f.read()
        out_meta = f.meta.copy()
        out_meta['compress'] = 'deflate'
    with rasterio.open(path, 'w', **out_meta) as f:
        f.write(arr)

def get_subgrid_waterdepth(rasterpath, facedata, outpath):
    """
    Get raster stats
    
    Parameters
    ----------
    rasterpath : str
        Path to raster file
    pts : Nx2 array
        Array with x and y location
    buffer : float
        Max size between different points. This length is used to include points outside
        domain when splitting raster
    df : pandas.DataFrame
        Dataframe with tables to be filled
    """

    # Check geometries
    check_geodateframe_rasterstats(facedata)

    # Open raster file
    with rasterio.open(rasterpath, 'r') as f:

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

            # # Get x and y range of part
            # x, y = prt.get_xy_range()
            
            # # Interpolate nearest cell
            # X, Y = np.meshgrid(x, y)
            # pts, values, mask = get_griddata_input(facedata.loc[prt.idx], prt)
            # valid = (mask & valid)
            # cellidx_sel = griddata(points=pts, values=values, xi=(X, Y), method='nearest')
            # cellidx_sel[~valid] = 0

            # Create array to assign water levels
            wlev_subgr = np.zeros(cellidx_sel.shape, dtype=out_meta['dtype'])
            
            # Create function to determine volume given the level
            def vol(level, bottom, cell_area, fm_wdep):
                volume = fm_wdep * len(bottom) * cell_area
                return np.maximum(0, level - bottom).sum() * cell_area - volume

            for cell_idx in np.unique(cellidx_sel):
                if cell_idx == 0:
                    continue
                
                # Mask
                cellmask = (cellidx_sel == cell_idx)
                if not cellmask.any():
                    continue
                
                # Get bottom values
                bottom = arr[cellmask]
                
                # Iterate to depth
                fm_wdep = facedata.at[cell_idx, 'value']
                bottom_mean = bottom.mean()
                facedata.at[cell_idx, 'avg_bedlev'] = bottom_mean
                facedata.at[cell_idx, 'std_bedlev'] = bottom.std()
                # Determine subgrid waterlevle
                if fm_wdep <= 0.0:
                    subgr_wlev = bottom.min()

                else:
                    # Determine subgrid waterlevle
                    subgr_wlev = brentq(
                        f=vol,
                        a=bottom.min() + fm_wdep - 1e-3,
                        b=bottom.mean() + fm_wdep + 1e-3,
                        args=(bottom, cell_area, fm_wdep),
                        xtol=2e-4
                    )
                
                # Add waterdepth to dataframe
                facedata.at[cell_idx, 'subgr_wlev'] = subgr_wlev

                # Add values for cell to raster
                wlev_subgr[cellmask] = subgr_wlev

            # Write to output raster
            with rasterio.open(outpath, 'w' if first else 'r+', **out_meta) as dst:
            
                # Determine water depth
                if not first:
                    wdep_subgr = dst.read(1, window=prt.window)
                else:
                    wdep_subgr = np.ones_like(wlev_subgr) * getattr(np, out_meta['dtype'])(f.nodata)
                    first = False
                
                wdep_subgr[valid] = np.maximum(wlev_subgr - arr, 0).astype(out_meta['dtype'])[valid]
                # wdep_subgr[arr == f.nodata] = f.nodata
            
                dst.write(wdep_subgr[None, :, :], window=prt.window)

    compress(outpath)

def check_geodateframe_rasterstats(facedata):
    """
    Check for type, columns and coordinates
    """
    if not isinstance(facedata, gpd.GeoDataFrame):
        raise TypeError('facedata should be type GeoDataFrame')

    if 0 in facedata.index.tolist():
        raise ValueError('0 occurs in index. This value if reserved.')

    # Check if facedata has required columns
    if ('facex' not in facedata.columns) or ('facey' not in facedata.columns):
        xy = list(zip(*[pt.coords[0] for pt in facedata.geometry.centroid]))
        facedata['facex'] = xy[0]
        facedata['facey'] = xy[1]

    # Check if coordinates are present.
    if 'crds' not in facedata.columns:
        facedata['crds'] =[row.coords[:] for row in facedata.geometry]


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
    for ix, iy in tqdm.tqdm(product(range(nx), range(ny)), total=nx*ny):

        part = RasterPart(f, xparts[ix], yparts[iy], xparts[ix+1], yparts[iy+1])
        
        if facedata is not None:
            # For each part, get points in part
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

        # print(ix, iy)
        yield part

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
        corners = self.get_corners()
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