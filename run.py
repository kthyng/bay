import matplotlib as mpl
# mpl.use('Agg')
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

# loc = sorted(glob('/rho/raid/dongyu/blended*.nc'))  # original bay? (Rayson, on blended grid)
# loc = sorted(glob('/rho/raid/dongyu/superposition/blended*.nc'))
# loc = '/rho/raid/dongyu/201007_new/blended201007.nc'  # new bay model output
loc = '/rho/raid/dongyu/superposition/blended201007.nc'  # new superposition

def check_tides(date, ff):
    '''Check tidal cycle for starting date and following 2 weeks.'''

    d = xr.open_mfdataset(loc)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    # data from entrance to bay
    # if ff == 1:  # forward
    d['zeta'].sel(ocean_time=slice(date, date+timedelta(days=14))).isel(yr=235, xr=435).plot()
    # elif ff == -1:  # backward
    #     d['zeta'].sel(ocean_time=slice(date-timedelta(days=14), date)).isel(yr=235, xr=435).plot()


def init(name, ff):
    time_units = 'seconds since 1970-01-01  00:00:00'

    nsteps = 10

    # Number of steps to divide model output for outputting drifter location
    N = 4

    # Number of days
    ndays = 14

    # # This is a forward-moving simulation
    # ff = 1

    # Time between outputs
    tseas = 3600.  # time between output in seconds
    ah = 0.
    av = 0.  # m^2/s

    # surface drifters
    z0 = 's'

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0
    doturb = 0

    # for periodic boundary conditions in the x direction
    doperiodic = 0

    # Flag for streamlines. All the extra steps right after this are for streamlines.
    dostream = 0

    # for 3d flag, do3d=0 makes the run 2d and do3d=1 makes the run 3d
    do3d = 0

    ## Choose method for vertical placement of drifters
    z0 = 's' # I know the size from checking #'s after eliminating those outside domain ' #'z' #'salt' #'s'
    num_layers = 1
    zpar = 0

    grid_filename = 'blended_grid.nc'

    proj = tracpy.tools.make_proj(setup='nwgom-pyproj')

    # Read in grid
    grid = tracpy.inout.readgrid(grid_filename, proj,
                                 vert_filename=loc[0],
                                 usespherical=True)

    # Initialize Tracpy class
    tp = Tracpy(loc, grid, name=name, tseas=tseas,
                ndays=ndays, nsteps=nsteps, N=N, ff=ff, ah=ah, av=av,
                doturb=doturb, do3d=do3d, z0=z0, zpar=zpar,
                time_units=time_units,
                usespherical=True)

    dx = 300  # distance (meters) between drifters in x and y
    ll0fname = 'calcs/ll0_inbay_dx' + str(dx) + '.npz'
    if not os.path.exists(ll0fname):
        llcrnrlon = -95.2; urcrnrlon = -94.5;
        llcrnrlat = 29; urcrnrlat = 29.78;
        xcrnrs, ycrnrs = grid.proj([llcrnrlon, urcrnrlon], [llcrnrlat, urcrnrlat])
        X, Y = np.meshgrid(np.arange(xcrnrs[0], xcrnrs[1], dx), np.arange(ycrnrs[0], ycrnrs[1], dx))

        lon0, lat0 = grid.proj(X, Y, inverse=True)

        # lon0, lat0 in projected coordinates
        x0, y0 = grid.proj(lon0, lat0)

        baypathxy = calcs.baypath(which='xy')
        # points that are in the bay (outside of shelf)
        inds = ~baypathxy.contains_points(np.vstack([x0.flatten(), y0.flatten()]).T).reshape(x0.shape)

        # remove shelf points
        x0 = x0[inds]; y0 = y0[inds]

        # back to latlon
        lon0, lat0 = grid.proj(x0, y0, inverse=True)

        # Eliminate points that are outside domain or in masked areas
        lon0, lat0 = tracpy.tools.check_points(lon0, lat0, tp.grid)

        # save starting locations for future use
        np.savez(ll0fname, lon0=lon0, lat0=lat0)

        # plot
        fig = plt.figure()
        figname = 'figures/starting_points_bay_dx' + str(dx) + '.png'
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-85.0))
        gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
        ax.set_extent([-95.4, -94.2, 28.8, 29.9], ccrs.PlateCarree())
        ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'
        ax.plot(lon0, lat0, 'r.', transform=ccrs.PlateCarree());
        # import pdb; pdb.set_trace()
        fig.savefig(figname, bbox_inches='tight')
    else:
        dtemp = np.load(ll0fname)
        lon0 = dtemp['lon0']; lat0 = dtemp['lat0']

    return tp, lon0, lat0


def run():

    ffs = [1]  # [1, -1]  # forward and backward moving simulations
    # to keep consistent between sims
    refdates = [datetime(2010, 7, 1, 0, 0)]
    # refdates = [datetime(2010, 2, 1, 0, 0), datetime(2010, 7, 1, 0, 0),
    #             datetime(2011, 2, 1, 0, 0), datetime(2011, 7, 1, 0, 0)]

    for refdate in refdates:
        for ff in ffs:
            if ff == 1:
                overallstartdate = refdate
                overallstopdate = overallstartdate + timedelta(days=14)
                basename = '_forward'
            elif ff == -1:
                overallstartdate = refdate + timedelta(days=14)
                overallstopdate = overallstartdate + timedelta(days=14)
                basename = '_backward'

            date = overallstartdate

            # Start from the beginning and add days on for loop
            # keep running until we hit the next month
            while date < overallstopdate:

                name = 'newsuperposition/' + date.isoformat()[0:13] + basename

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc') and \
                    not os.path.exists('tracks/' + name + 'gc.nc'):

                    tp, lon0, lat0 = init(name, ff)
                    lonp, latp, zp, t, T0, U, V = tracpy.run.run(tp, date, lon0, lat0)

                # Increment by 4 hours for next loop
                date = date + timedelta(hours=4)



if __name__ == "__main__":
    run()
