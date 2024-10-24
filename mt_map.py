import matplotlib.pyplot as plt
import matplotlib.colors
import rasterio
from rasterio.mask import mask

import numpy as np
import pandas as pd
import rasterio
import geopandas
import warnings
import shapely

matplotlib.rcParams['font.family'] = 'Arial'

earth_radius = 6377.6
deg = earth_radius * np.pi / 180  # 111.3 km
arcmin = deg / 60  # 1.8 km
arcsec = deg / 3600.  # 31 m
cell_size = 9.  # km

def _names_match(name, mtn_name):
    ref = name.lower()
    query = mtn_name.lower()
    if ' ' in ref:
        if (ref.split(' ')[-1][:4] not in 
            ['moun', 'mts.', 'plat', 'isle', 'isla', 'peni', 'shie'] or
            ref.split(' ')[0] == 'new'):
            ref, query = ref.split(' ')[-1], query.split(' ')[-1]
    return ref[:4] == query[:4]

def _bin_map(matrix, factor=5, value='relief'):
    mat = np.where(matrix > 0., matrix, np.nan)
    height, width = mat.shape
    res_h, res_w = int(height / factor), int(width / factor)
    tail_h, tail_w = height % factor, width % factor
    in_mat = mat[:height - tail_h, :width - tail_w].reshape(res_h, factor, res_w, factor)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        if value == 'mean':
            return np.nanmean(np.nanmean(in_mat, axis=1), axis=2)
        elif value == 'std':
            return np.sqrt(np.nanmean(np.nanmean(in_mat ** 2, axis=1), axis=2) -
                           np.nanmean(np.nanmean(in_mat, axis=1), axis=2) ** 2)
        elif value == 'min':
            return np.nanmin(np.nanmin(in_mat, axis=1), axis=2)
        elif value == 'max':
            return np.nanmax(np.nanmax(in_mat, axis=1), axis=2)
        else:  # value == 'relief'
            return (np.nanmean(np.nanmean(in_mat, axis=1), axis=2) -
                     np.nanmin(np.nanmin(in_mat, axis=1), axis=2))


def _get_map(mt_name, topo):
    if type(mt_name) == int:
        name = shp.name[mt_name]
        geom = shp.geometry[mt_name]
    else:
        for name, geom in zip(shp.name, shp.geometry):
            if _names_match(name, mt_name):
                break
    
    geom_map = [shapely.geometry.mapping(geom)]
    maps_topo = mask(topo, geom_map, crop=True)[0][0]

    return name, np.where((maps_topo > 0), maps_topo, np.nan), geom

def _solve(p, A):
    Delta = -16*A + p**2
    if Delta < 0:
        return np.nan, np.nan
    else:
        return p/4 + np.sqrt(Delta)/4, p/4 - np.sqrt(Delta)/4

def _cal_eros(maps):  # From Vance et al. (2003)
    factor = int(cell_size / arcmin) + 1
    rlf = _bin_map(maps, factor=factor, value='relief')
    eros_map = 1e-2 * 10 ** (rlf / 500.)
    return np.nanmean(eros_map), eros_map

def plot_mt(mt_name, savfig=0, show=1):
    fig, ax = plt.subplots(1, 1, figsize=(6.5, 5.5))

    name, maps_topo, geom = _get_map(mt_name, topo)
    bounds = geom.bounds
    extent = [bounds[0], bounds[2], bounds[1], bounds[3]]
    
    colors_land = plt.cm.terrain(np.linspace(0, 1, 200))
    # combine them and build a new colormap
    cut_terrain_map = matplotlib.colors.LinearSegmentedColormap.from_list('cut_terrain', colors_land)


    plt_map = np.where(maps_topo > 0., maps_topo, np.nan)
    fig = plt.imshow(plt_map, extent=extent, cmap=cut_terrain_map, vmin=0, vmax=8000)
    # cbar = plt.colorbar(fig, orientation='horizontal')

    ax.set_xlabel('Longitude (deg)')
    ax.set_ylabel('Latitude (deg)')

    if savfig:
        plt.savefig(fig_path + mt_name + '.png', dpi=600)
    if show:
        plt.show()
    
shp = geopandas.read_file(data_path + 'ne_50m_geography_regions_polys/ne_50m_geography_regions_polys.shp')
topo = rasterio.open(data_path + 'ETOPO1_Bed_g_geotiff.tif')

# plot_mt("Tibet", savfig=1)
# plot_mt("Himalayas")

# plot_mt("Alps", savfig=1)
plot_mt("Zagros", savfig=0)
