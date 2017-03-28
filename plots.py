'''
Plot functions for bay work.
'''

# import matplotlib as mpl
# mpl.use('Agg')
from datetime import datetime
import matplotlib.pyplot as plt
from glob import glob
import xarray as xr
import cartopy.crs as ccrs
import cmocean.cm as cmo
from tracpy import op
import seaborn as sns


mpl.rcParams.update({'font.size': 14})
mpl.rcParams['font.sans-serif'] = 'Arev Sans, Bitstream Vera Sans, Lucida Grande, Verdana, Geneva, Lucid, Helvetica, Avant Garde, sans-serif'
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.cal'] = 'cursive'
mpl.rcParams['mathtext.rm'] = 'sans'
mpl.rcParams['mathtext.tt'] = 'monospace'
mpl.rcParams['mathtext.it'] = 'sans:italic'
mpl.rcParams['mathtext.bf'] = 'sans:bold'
mpl.rcParams['mathtext.sf'] = 'sans'
mpl.rcParams['mathtext.fallback_to_cm'] = 'True'


def conditions():
    '''Plot wind, river, tide for selected time periods.'''

    refdates = [datetime(2010, 2, 1, 0, 0), datetime(2011, 2, 1, 0, 0)]
    plotdrifters = False  # whether or not to plot summary drifter results for same time period


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

    name = ''

    fig, axes = plt.subplots(nplots, 1, sharex=True, figsize=(14,8))
    for refdate, dy, dy2, color, xtext in zip(refdates, dys, dy2s, colors, xtexts):

        basename = '_14days_dx300'
        if refdate.day == 15:
            basename = '_backward' + basename
        elif refdate.day == 1:
            basename = '_forward' + basename
            dtitle = 'Number drifters\noutside bay'

        if 'forward' in basename:
            # datetime format
            start = refdate
            stop = start + timedelta(days=14*2)
        elif 'backward' in basename:
            start = refdate - timedelta(days=14)
            stop = start + timedelta(days=2*14)
            # print(start, stop)
        # string format
        start = start.isoformat()[:10]
        stop = stop.isoformat()[:10]

        name += start + '-' + stop + '_'

        # tide: from blended model output
        zeta = get_from_blended('zeta', start, stop)
        doy = zeta['ocean_time'].to_index().dayofyear + zeta['ocean_time'].to_index().hour/24.
        axes[0].plot(doy, zeta.values + dy2, color=color, lw=1)
        axes[0].set_ylabel('Sea surface\nheight [m]')

        # wind: from shelf model
        uwind = get_from_shelf('Uwind', start, stop)
        vwind = get_from_shelf('Vwind', start, stop)
        doy = uwind['ocean_time'].to_index().dayofyear + uwind['ocean_time'].to_index().hour/24.
        # import pdb; pdb.set_trace()
        axes[1].quiver(doy, np.zeros(len(doy))+dy, uwind.data, vwind.data, color=color,
                       headaxislength=0, headlength=0, width=0.2, units='y', scale_units='y', scale=1)
        axes[1].set_ylabel('Wind [m/s]')
        axes[1].set_ylim(-22, 22)

        # river: from suntans model
        river = get_from_suntans(start, stop)
        doy = river.index.dayofyear + river.index.hour/24.
        axes[2].plot(doy, river.boundary_Q, color=color, lw=2)
        # axes[2].xaxis_date()
        axes[2].set_ylabel('River discharge\n[m$^3$/s]')

        axes[2].text(xtext, 0.01, start[:4], color=color, transform=axes[2].transAxes)
        # drifters
        if plotdrifters:
            # import pdb; pdb.set_trace()
            df = pd.read_csv('calcs/enterexit_sim3_' + start[:7] + basename + '.csv', parse_dates=True, index_col=0)  # analysis of drifter tracks
            nfiles = len(df.columns)
            # nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
            # nperfile = np.load('calcs/ll0_inbay_dx300.npz')['lon0'].item().shape
            doy = df.index.dayofyear + df.index.hour/24.
            axes[3].plot(doy, df.sum(axis='columns').divide(nfiles), color=color)
            axes[3].set_ylabel(dtitle)
            # df.sum(axis='columns').divide(nfiles).plot(ax=axes[3])  # sum across columns
            if 'drifters' not in name:
                name += 'drifters'

    # import pdb; pdb.set_trace()
    axes[-1].xaxis.set_major_formatter(mpl.dates.DateFormatter('%b %d'))
    axes[-1].axis('tight')
    fig.subplots_adjust(left=0.08, right=0.99, top=0.97, hspace=0.17)
    fig.autofmt_xdate(bottom=0.07)
    fig.savefig('figures/io/' + name + '.png', bbox_inches='tight')
    fig.savefig('figures/io/' + name + '.pdf', bbox_inches='tight')


def io():

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


def drifters():
    '''Plot drifters in time from multiple simulations.'''

    year = '2010'
    month = '02'

    # model output
    m = xr.open_dataset('/rho/raid/dongyu/blended' + year + month + '.nc')
    dates = m['ocean_time'].to_pandas().dt.to_pydatetime()
    lon = m['lon'].isel(xr=slice(1,-1), yr=slice(1,-1)).data
    lat = m['lat'].isel(xr=slice(1,-1), yr=slice(1,-1)).data
    dd = 5  # quiver

    df = []  # forward moving drifter sims
    dfbase = '_forward_14days_dx300'
    dffiles = glob('tracks/' + year + '-' + month + '*' + dfbase + '.nc')
    for date in dates:
        datestr = date.isoformat()[:13] # e.g. '2010-02-01T00'
        # datestr2 = date.strftime('%Y-%m-%d %H:%M:%S')  # e.g. '2010-02-01 00:00:00'
        # name of drifter file that has starting date of date
        fname = dffiles[np.where([datestr in dffiles[i] for i in range(len(dffiles))])[0][0]]
        # add on next forward drifters sim
        df.append(xr.open_dataset(fname))
        df['tp'] = (('nt'), df['tp'].isel(ntrac=0))  # change tp to 1d
        df = df.swap_dims({'nt': 'tp'})  # so that can index off tp


        # # backward drifters sim 1
        # db1 = xr.open_dataset('tracks/')

        # start plot
        fig = plt.figure(figsize=(12,10))
        ax = plt.axes(projection=ccrs.Mercator(central_longitude=-85.0))
        gl = ax.gridlines(linewidth=0.2, color='gray', alpha=0.5, linestyle='-', draw_labels=True)
        ax.set_extent([-95.0, -94.5, 29.2, 29.6], ccrs.PlateCarree())
        ax.coastlines(resolution='10m')  # coastline resolution options are '110m', '50m', '10m'

        # plot model output
        # u = op.resize(m['u'].sel(ocean_time=datestr).isel(yr=slice(1,-1)).data, 1)
        # v = op.resize(m['v'].sel(ocean_time=datestr).isel(xr=slice(1,-1)).data, 0)
        # s = np.sqrt(u**2 + v**2)
        ax.pcolormesh(lon, lat, s, cmap=cmo.speed, transform=ccrs.PlateCarree())
        ax.quiver(lon[::5,::5], lat[::5,::5], u[::5,::5], v[::5,::5], transform=ccrs.PlateCarree())

        # plot drifters
        # choose only drifters that enter/exit domain with df['lonp'].sel(tp=datestr).isel(ntrac=[1,100])
        ax.plot(df['lonp'].sel(tp=datestr), df['latp'].sel(tp=datestr), 'r.', transform=ccrs.PlateCarree());
        # ax.plot(lonp[:].T, latp[:].T, 'k', lw=0.3, alpha=0.7, transform=ccrs.PlateCarree());
        fig.savefig(figname, bbox_inches='tight')
        plt.close(fig)


def fit():
    '''Plot statistical fit from calcs.stats().'''

    figure()
    model.fittedvalues.plot(color='r')
    model.resid.plot(color='g')
    df['drifters_tidal'].plot(color='k', linewidth=3)
    (df['zeta_shifted']/(df['zeta_shifted'].max()/df['drifters_tidal'].max())).plot(color='m')
