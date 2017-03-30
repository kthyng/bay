'''
Plot functions for bay work.
'''

# import matplotlib as mpl
# mpl.use('Agg')
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
from glob import glob
import xarray as xr
import cartopy.crs as ccrs
import cmocean.cm as cmo
# from tracpy import op
import seaborn as sns
import matplotlib as mpl
import pandas as pd
import numpy as np
from tracpy import op
import os
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.ticker as mticker


mpl.rcParams.update({'font.size': 12})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


def conditions(season, direction='forward'):
    '''Plot wind, river, tide for selected time periods.'''

    plotdrifters = True  # whether or not to plot summary drifter results for same time period

    basename = '_14days_dx300'
    if season == 'winter':
        month = 2
    elif season == 'summer':
        month = 7
    if direction == 'forward':
        day = 1
        basename = '_forward' + basename
        dtitle = 'Number drifters\noutside bay'
    elif direction == 'backward':
        day = 15
        basename = '_backward' + basename
        dtitle = 'Number drifters\noutside bay'
    refdates = [datetime(2010, month, day, 0, 0), datetime(2011, month, day, 0, 0)]

    dys = [-10, 10]
    dy2s = [-0.4, 0.4]
    colors = ['crimson', 'darkcyan']
    xtexts = [0.85, 0.925]

    if plotdrifters == True:
        nplots = 4
    else:
        nplots = 3

    if len(refdates) == 1:  # don't shift if only one dataset plotted
        dys = [0, 0]
        dy2s = [0, 0]

    fig, axes = plt.subplots(nplots, 1, sharex=True, figsize=(14,8))
    for refdate, dy, dy2, color, xtext in zip(refdates, dys, dy2s, colors, xtexts):

        if direction == 'forward':
            # datetime format
            start = refdate
            stop = start + timedelta(days=14*2)
        elif direction == 'backward':
            start = refdate - timedelta(days=14)
            stop = start + timedelta(days=2*14)

        # string format
        start = start.isoformat()[:10]
        stop = stop.isoformat()[:10]
        name = start + '-' + stop + '_'

        # read in dataframe
        File = 'calcs/df_' + refdate.isoformat()[:7] + '_' + direction + '.csv'
        df = pd.read_csv(File, parse_dates=True, index_col=0)

        # tide: from blended model output
        doy = df.index.dayofyear + df.index.hour/24. + df.index.minute/3600.
        axes[0].plot(doy, df['zeta'].values + dy2, color=color, lw=1)
        axes[0].set_ylabel('Sea surface\nheight [m]')

        # wind: from shelf model
        axes[1].quiver(doy, np.zeros(len(doy))+dy, df['uwind'].data, df['vwind'].data, color=color,
                       headaxislength=0, headlength=0, width=0.2, units='y', scale_units='y', scale=1)
        axes[1].set_ylabel('Wind [m/s]')
        axes[1].set_ylim(-22, 22)

        # river: from suntans model
        axes[2].plot(doy, df['river'], color=color, lw=2)
        axes[2].set_ylabel('River discharge\n[m$^3$/s]')
        axes[2].text(xtext, 0.9, start[:4], color=color, transform=axes[2].transAxes)

        # drifters
        if plotdrifters:
            axes[3].plot(doy, df['drifters'], color=color)
            axes[3].set_ylabel(dtitle)
            if 'drifters' not in name:
                name += 'drifters'

    axes[-1].xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %d'))
    axes[-1].axis('tight')
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, hspace=0.17)
    fig.autofmt_xdate(bottom=0.07)
    fig.savefig('figures/io/' + name + '.png', bbox_inches='tight')
    fig.savefig('figures/io/' + name + '.pdf', bbox_inches='tight')
    plt.show()


# def io():
#
#     # df = pd.read_csv('calcs/enterexit_sim3_2010-02_forward_14days_dx300.csv', parse_dates=True, index_col=0)  # analysis of drifter tracks
#     # cols = df.filter(like='forward').columns  # list of columns
#     # number of files being used per time, as proxy for # available drifters
#     # start = '2009-05-14' # '2009-04-14'
#     # stop = '2009-06-10'  # '2009-05-10'
#     # nfiles = (~np.isnan(df[start:stop][cols])).astype(int).sum(axis='columns')
#     df = pd.read_csv('calcs/enterexit_sim3_2010-02_forward_14days_dx300.csv', parse_dates=True, index_col=0)  # analysis of drifter tracks
#     start = df.index[0].isoformat()[:10]
#     stop = df.index[-1].isoformat()[:10]
#     nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
#     loc = sorted(glob('/rho/raid/dongyu/blended*.nc'))
#     dm = xr.open_mfdataset(loc)  # blended model product
#
#     # Plot forward, high flow
#     fig, axes = plt.subplots(4, 1, figsize=(16,8), sharex=True)
#
#     # tides
#     dm['zeta'].sel(ocean_time=slice(start, stop)).isel(yr=235, xr=435).plot(ax=axes[0])
#
#     # river discharge
#     paper.get_from_suntans(start, stop)['boundary_Q'].plot(ax=axes[1])
#
#     # winds
#
#
#     # forward: drifter exiting bay?
#     df.sum(axis='columns').divide(nfiles).plot(ax=axes[3])  # sum across columns


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


def drifters():
    '''Plot drifters in time from multiple simulations.'''

    year = '2010'
    month = '02'

    # model output
    m = xr.open_dataset('/rho/raid/dongyu/blended' + year + month + '.nc')
    dates = m['ocean_time'].to_pandas().dt.to_pydatetime()
    lon = m['lon'].isel(xr=slice(1,-1), yr=slice(1,-1)).data
    lat = m['lat'].isel(xr=slice(1,-1), yr=slice(1,-1)).data
    dd = 8  # quiver
    datestr0 = dates[0].isoformat()[:13]  # starting date
    datestr1 = dates[-1].isoformat()[:13]  # ending date

    # load previously-saved mechanism time series
    # hourly, for this time frame
    zeta = pd.read_csv('calcs/zeta2009-2011.csv', parse_dates=True, index_col=0)[datestr0:datestr1].resample('60T').interpolate()
    uwind = pd.read_csv('calcs/uwind2009-2011.csv', parse_dates=True, index_col=0)['-6.20427322388']
    uwind = uwind[datestr0:datestr1].resample('60T').interpolate()
    vwind = pd.read_csv('calcs/vwind2009-2011.csv', parse_dates=True, index_col=0)['-2.68748617172']
    vwind = vwind[datestr0:datestr1].resample('60T').interpolate()
    wmax = abs(vwind).max()  # vertical only
    wmin = -wmax
    river = pd.read_csv('calcs/river2009-2011.csv', parse_dates=True, index_col=0)['108.394273505']
    rmin = river.min()
    rmax = river.max()
    river = river[datestr0:datestr1].resample('60T').interpolate()
    # surface speed max
    smax = 1.2

    dss = []  # forward moving drifter sims
    dsbase = '_forward_14days_dx300'
    dsfiles = glob('tracks/' + year + '-' + month + '*' + dsbase + '.nc')
    # # name of drifter file that has starting date of date
    # fname = dffiles[np.where([datestr in dffiles[i] for i in range(len(dffiles))])[0][0]]
    for dsfile in dsfiles:
        dstemp = xr.open_dataset(dsfile)
        dstemp['tp'] = (('nt'), dstemp['tp'].isel(ntrac=0))  # change tp to 1d
        dstemp = dstemp.swap_dims({'nt': 'tp'})  # so that can index off tp
        # choose only drifters that enter/exit domain with df['lonp'].sel(tp=datestr).isel(ntrac=[1,100])
        File = 'calcs/enterexit/enterexit_sim3_' + dsfile.split('/')[1][:13] + '_forward_14days_dx300.npz'
        idrifters = list(np.load(File)['idrifters'].item())
        dss.append(dstemp.isel(ntrac=idrifters))

    import cartopy.feature as cfeature
    land_10m = cfeature.NaturalEarthFeature('physical', 'land', '10m',
                                            edgecolor='face',
                                            facecolor=cfeature.COLORS['land'])

    os.makedirs('figures/' + datestr0)
    for date in dates:
        datestr = date.isoformat()[:13] # e.g. '2010-02-01T00'
        datestrlong = date.isoformat()  # since tracks are every 15 min, to be more specific
        datestrlong2 = (date+timedelta(seconds=1)).isoformat()  # since tracks are every 15 min, to be more specific
        datenice = date.strftime('%b %d, %Y %H:%M')
        figname = 'figures/' + datestr0 + '/' + datestr + '.png'

        # start plot
        fig = plt.figure(figsize=(9.5, 10))
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-85.0))
        gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
        # the following two make the labels look like lat/lon format
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        # gl.xlocator = mticker.FixedLocator([-105, -95, -85, -75, -65])  # control where the ticks are
        # gl.xlabel_style = {'size': 15, 'color': 'gray'}  # control how the tick labels look
        # gl.ylabel_style = {'color': 'red', 'weight': 'bold'}
        gl.xlabels_top = False  # turn off labels where you don't want them
        gl.ylabels_right = False

        ax.set_extent([-95.42, -94.4, 28.95, 29.8], ccrs.PlateCarree())
        # ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'

        # plot model output
        u = op.resize(m['u'].sel(ocean_time=datestr).isel(yr=slice(1,-1)).data, 1)
        v = op.resize(m['v'].sel(ocean_time=datestr).isel(xr=slice(1,-1)).data, 0)
        # salt = m['salt'].sel(ocean_time=datestr).isel(yr=slice(1,-1), xr=slice(1,-1), s_rho=-1).data
        s = np.sqrt(u**2 + v**2)
        mappable = ax.pcolormesh(lon, lat, s, cmap=cmo.speed, transform=ccrs.PlateCarree(), vmin=0, vmax=smax)
        # lower resolution part
        ax.quiver(lon[:145:dd/4,::dd], lat[:145:dd/4,::dd],
                  u[:145:dd/4,::dd], v[:145:dd/4,::dd], scale=15,
                  color='k', transform=ccrs.PlateCarree(), pivot='middle')
        # higher resolution part in bay
        Q = ax.quiver(lon[144::dd,::dd], lat[144::dd,::dd],
                  u[144::dd,::dd], v[144::dd,::dd], scale=15,
                  color='k', transform=ccrs.PlateCarree(), pivot='middle')
        qk = ax.quiverkey(Q, 0.15, 0.25, 0.5, r'0.5 m$\cdot$s$^{-1}$ current', labelcolor='k', fontproperties={'size': '10'})

        ax.add_feature(land_10m, facecolor='0.9')

        # colorbar
        cax = fig.add_axes([0.15, 0.8, .25, 0.02])
        cb = fig.colorbar(mappable, cax=cax, orientation='horizontal')#, pad=0.1)
        cb.ax.tick_params(labelsize=12, length=2, color='k', labelcolor='k')
        cb.set_ticks(np.arange(0, 1.4, 0.2))
        cb.set_label(r'Surface speed [m$\cdot$s$^{-1}$]', fontsize=12, color='k')

        # tide signal
        axtide = fig.add_axes([0.14, 0.62, .27, 0.08], frameon=False)#, transform=ax.transAxes)
        zeta.plot(ax=axtide, color='k', legend=False, linewidth=1.5)
        zeta[datestr:datestr].plot(ax=axtide, marker='o', color='r', legend=False)
        axtide.get_yaxis().set_visible(False)
        axtide.get_xaxis().set_visible(False)
        axtide.text(0.03, 0.95, 'sea surface', transform=axtide.transAxes, fontsize=12)

        # river discharge
        axriver = fig.add_axes([0.14, 0.52, .27, 0.08], frameon=False)#, transform=ax.transAxes)
        river.plot(ax=axriver, color='k', legend=False, linewidth=1.5)
        river[datestr:datestr].plot(ax=axriver, marker='o', color='r', legend=False)
        axriver.get_yaxis().set_visible(False)
        axriver.get_xaxis().set_visible(False)
        axriver.text(0.03, 0.95, 'river', transform=axriver.transAxes, fontsize=12)
        axriver.set_ylim(rmin, rmax)  # to compare with other months
        # axriver.axis('tight')

        # wind time series
        doy = uwind.index.dayofyear + uwind.index.hour/24. + uwind.index.minute/3600.
        uwindtemp = pd.DataFrame(uwind)[datestr:datestr]  # need it to be a dataframe not series for next line indexing
        doynow = uwindtemp.index.dayofyear + uwindtemp.index.hour/24. + uwindtemp.index.minute/3600.
        axwind = fig.add_axes([0.14, 0.41, .27, 0.1], frameon=False)#, transform=ax.transAxes)
        axwind.quiver(doy, np.zeros(len(doy)), uwind, vwind, color='k',
                       headaxislength=0, headlength=0, width=0.1, units='y', scale_units='y', scale=1)
        axwind.quiver(doynow, 0, uwind[datestr], vwind[datestr], color='r',
                       headaxislength=0, headlength=0, width=1.0, units='y', scale_units='y', scale=1)
        axwind.get_yaxis().set_visible(False)
        axwind.get_xaxis().set_visible(False)
        axwind.text(0.03, 0.95, 'wind', transform=axwind.transAxes, fontsize=12)
        axwind.set_ylim(wmin, wmax)  # to compare with other months
        axwind.set_xlim([doy.min(), doy.max()])

        # plot time
        fig.text(0.15, 0.375, datenice, transform=fig.transFigure, fontsize=14)

        # plot drifters
        for ds in dss:
            try:  # drifter files are available every 4 hours not hourly
                ds = ds.sel(tp=slice(datestrlong,datestrlong2))
                ax.plot(ds['lonp'], ds['latp'], 'o', color='#782277', markersize=5, alpha=0.7, transform=ccrs.PlateCarree());
            except:
                # print(ds.tp[0])
                continue
        fig.savefig(figname, bbox_inches='tight')
        plt.close(fig)


def fit(model, df, cols):
    '''Plot statistical fit from calcs.stats().'''

    plt.figure(figsize=(14,6))
    df['drifters_subtidal'].plot(color='k', linewidth=3)
    model.fittedvalues.plot(color='r')
    model.resid.plot(color='g')
    names = ['drifters', 'fit', 'residual']
    for col, ls in zip(cols, ['-','--','-.']):
        df[col].plot(color='darkorange', linewidth=1, linestyle=ls)
        names.append(col)
    plt.legend(names, fontsize=12, loc='best')


if __name__ == "__main__":
    drifters()
