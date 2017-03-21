import matplotlib as mpl
mpl.use('Agg')
import tracpy
import tracpy.calcs
from tracpy.tracpy_class import Tracpy
import os
from datetime import datetime, timedelta
import numpy as np
import netCDF4 as netCDF
import matplotlib.pyplot as plt
from glob import glob
from scipy import ndimage
from matplotlib.path import Path
import cartopy
import cartopy.crs as ccrs
import xarray as xr
import pandas as pd
import paper

loc = glob('/rho/raid/dongyu/blended*.nc')

# from time import time
# from shapely.geometry import Point
# from shapely.geometry.polygon import Polygon
# start = time()
# inds = ~baypathll.contains_points(np.vstack([lonp.flatten(), latp.flatten()]).T).reshape(lonp.shape)
# end = time()
#
# polygon = Polygon(zip(lonpath, latpath))
# points = [Point(lon, lat) for (lon, lat) in zip(lonp.flatten(), latp.flatten())]
# start = time()
# end = time()



def make_baypath():

    # getting path encompassing NOT bay only
    baypathfile = 'calcs/baypathfile.npz'
    # from init()
    if not os.path.exists(baypathfile):
        grid_filename = 'blended_grid.nc'
        proj = tracpy.tools.make_proj(setup='nwgom')

        # Read in grid
        grid = tracpy.inout.readgrid(grid_filename, proj,
                                     vert_filename=loc[0],
                                     usespherical=True)
        # use pts in projected space OUTSIDE of the coastpath I clicked out before
        outerpathxy = np.load('../shelf_transport/calcs/coastpaths.npz')['outerpathxy'].item()

        # remove points outside the path
        ov = outerpathxy.vertices
        xpath = np.hstack((ov[0:300,0], ov[300,0], ov[0,0], ov[0,0]))
        ypath = np.hstack((ov[0:300,1], 390, 390, ov[0,1]))
        newv = np.vstack((xpath, ypath)).T
        baypathxy = Path(newv)  # make path encompassing shelf not bay for differentiating
        lonpath, latpath = grid.proj(xpath, ypath, inverse=True)
        newvll = np.vstack((lonpath, latpath)).T
        baypathll = Path(newvll)
        np.savez(baypathfile, baypathxy=baypathxy, baypathll=baypathll)
    else:
        baypathll = np.load(baypathfile)['baypathll'].item()

    return baypathll


def io():

    baypathll = make_baypath()

    basename = '_14days_dx300'
    refdate = datetime(2011, 2, 15, 0, 0)  # datetime(2010, 7, 15, 0, 0) datetime(2010, 2, 1, 0, 0) datetime(2010, 7, 1, 0, 0)

    if refdate.day == 15:
        basename = '_backward' + basename
    elif refdate.day == 1:
        basename = '_forward' + basename

    if 'forward' in basename:
        start = refdate
        end = start + timedelta(days=14*2)
    elif 'backward' in basename:
        start = refdate - timedelta(days=14)
        end = start + timedelta(days=2*14)
    Files = glob('tracks/' + start.isoformat()[:7] + '*' + basename + '*.nc')
    dfdates = pd.date_range(start=start.isoformat()[:10] + ' 00:00:01', end=end.isoformat()[:10] + ' 00:00:01', freq='900S')
    df = pd.DataFrame(index=dfdates)
    # import pdb; pdb.set_trace()
    for File in Files:
        # print(File)
        d = netCDF.Dataset(File)
        lonp = d['lonp'][:]; latp = d['latp'][:]; tp = d['tp'][0,:]
        d.close()

        # points that are outside of bay
        inds = baypathll.contains_points(np.vstack([lonp.flatten(), latp.flatten()]).T).reshape(lonp.shape)

        numexit = inds.sum(axis=0)  # number of drifter outside of bay by time
        # add column to dataframe with this information
        simstartdate = File.split('/')[-1].split('14days')[0][:-1]

        dftemp = pd.DataFrame(index=netCDF.num2date(tp, 'seconds since 1970-01-01  00:00:00'), data={simstartdate: numexit})
        df = df.join(dftemp)  # add column to dataframe

        df.to_csv('calcs/enterexit_sim3_' + start.isoformat()[:7] + basename + '.csv')  # save every time step


def plot_io():

    # df = pd.read_csv('calcs/enterexit_sim3_2010-02_forward_14days_dx300.csv', parse_dates=True, index_col=0)  # analysis of drifter tracks
    # cols = df.filter(like='forward').columns  # list of columns
    # number of files being used per time, as proxy for # available drifters
    # start = '2009-05-14' # '2009-04-14'
    # stop = '2009-06-10'  # '2009-05-10'
    # nfiles = (~np.isnan(df[start:stop][cols])).astype(int).sum(axis='columns')
    df = pd.read_csv('calcs/enterexit_sim3_2010-02_forward_14days_dx300.csv', parse_dates=True, index_col=0)  # analysis of drifter tracks
    start = df.index[0].isoformat()[:10]
    stop = df.index[-1].isoformat()[:10]
    nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
    loc = sorted(glob('/rho/raid/dongyu/blended*.nc'))
    dm = xr.open_mfdataset(loc)  # blended model product

    # Plot forward, high flow
    fig, axes = plt.subplots(4, 1, figsize=(16,8), sharex=True)

    # tides
    dm['zeta'].sel(ocean_time=slice(start, stop)).isel(yr=235, xr=435).plot(ax=axes[0])

    # river discharge
    paper.get_from_suntans(start, stop)['boundary_Q'].plot(ax=axes[1])

    # winds


    # forward: drifter exiting bay?
    df.sum(axis='columns').divide(nfiles).plot(ax=axes[3])  # sum across columns


def tracks():

    Files = glob('tracks/*_15days_*.nc')
    for File in Files:
        figname = 'figures/' + File.split('/')[-1].split('.')[0] + '.png'
        if not os.path.exists(figname):
            d = netCDF.Dataset(File)
            lonp = d['lonp'][:]; latp = d['latp'][:];
            d.close()
            fig = plt.figure(figsize=(12,10))
            ax = plt.axes(projection=ccrs.Mercator(central_longitude=-85.0))
            gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
            ax.set_extent([-95.4, -94.2, 28.8, 29.9], ccrs.PlateCarree())
            ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
            ax.plot(lonp[:,0], latp[:,0], 'r.', transform=ccrs.PlateCarree());
            ax.plot(lonp[:].T, latp[:].T, 'k', lw=0.3, alpha=0.7, transform=ccrs.PlateCarree());
            fig.savefig(figname, bbox_inches='tight')
            plt.close(fig)


if __name__ == "__main__":
    io()
