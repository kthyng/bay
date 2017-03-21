'''
Supporting code for paper
'''

import matplotlib.pyplot as plt
from glob import glob
import xarray as xr
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.dates import date2num
from datetime import datetime, timedelta
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm
import statsmodels.api as sm
# import statsmodels.api as sm
# tsa = sm.tsa
from scipy.signal import argrelextrema


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


def get_from_blended(var, start, stop, yrho=230, xrho=445):
    '''Get some information from blended model product.

    Default xr/yr location is in entrance to Galveston Bay.
    '''

    # assuming time periods are consecutive (two months in a row)
    loc = glob('/rho/raid/dongyu/blended' + start[:4] + start[5:7] + '*.nc')
    loc.extend(glob('/rho/raid/dongyu/blended' + stop[:4] + stop[5:7] + '*.nc'))
    # remove nonunique values (in case all within one month)
    d = xr.open_mfdataset(sorted(list(set(loc))), concat_dim='ocean_time')

    return d[var].sel(ocean_time=slice(start, stop)).isel(yr=yrho, xr=xrho)


def get_from_shelf(var, start, stop, xrho=274, yrho=165):
    '''Get something from shelf model.

    Default xr/yr location is in entrance to Galveston Bay.
    '''

    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_hindcast_agg'
    d = xr.open_dataset(loc)

    if var == 'wind':
        var1 = d['Uwind'].sel(ocean_time=slice(start, stop)).isel(eta_rho=yrho, xi_rho=xrho)
        var2 = d['Vwind'].sel(ocean_time=slice(start, stop)).isel(eta_rho=yrho, xi_rho=xrho)
        var = np.sqrt(var1**2 + var2**2)
    elif var == 'Uwind' or var == 'Vwind':
        var = d[var].sel(ocean_time=slice(start, stop)).isel(eta_rho=yrho, xi_rho=xrho)
    # print(var)
    return var


def get_from_TWDB(start, stop):
    '''Get river data.

    from https://waterdatafortexas.org/coastal/sites/TRIN
    Not using this anymore
    '''

    df = pd.read_csv('TRIN_seawater_salinity_daily.csv', comment='#', parse_dates=True, index_col=0)
    return df[start:stop]


def get_from_suntans(start, stop):
    '''Get river discharge from suntans boundary file.'''

    ds = xr.open_dataset('/rho/raid/dongyu/GalvCoarse_BC_' + start[:4] + '.nc')
    ds.swap_dims({'Nt': 'time'}, inplace=True)  # so can index based on time
    return ds['boundary_Q'].sel(time=slice(start, stop)).sum(axis=1).to_dataframe()


def plot_conditions():
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


def stats():
    '''Stats between drifters and forcing mechanisms.'''

    File = 'calcs/enterexit_sim3_2010-07_forward_14days_dx300.csv'
    df = pd.read_csv(File, parse_dates=True, index_col=0)
    nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
    y = df.sum(axis='columns').divide(nfiles).resample('15min', base=0).interpolate()

    start = df.index[0].isoformat()[:10]
    stop = df.index[-1].isoformat()[:10]

    river = get_from_suntans(start, stop).resample('15min').interpolate()['boundary_Q'][:y.index[-1].isoformat()]

    # scipy.stats.pearsonr(y.values, river[:y.index[-1].isoformat()]['boundary_Q'].values)
    # scipy.stats.pearsonr(y, river)


    zeta = get_from_blended('zeta', start, stop).to_dataframe().resample('15min').interpolate()['zeta'][:y.index[-1].isoformat()]
    # scipy.stats.pearsonr(y, zeta)



    uwind = get_from_shelf('Uwind', start, stop).to_dataframe().resample('15min').interpolate()['Uwind'][:y.index[-1].isoformat()]
    vwind = get_from_shelf('Vwind', start, stop).to_dataframe().resample('15min').interpolate()['Vwind'][:y.index[-1].isoformat()]
    s = np.sqrt(uwind**2 + vwind**2)
    theta = np.unwrap(np.arctan2(vwind, uwind))
    # scipy.stats.pearsonr(y, s)
    # scipy.stats.pearsonr(y, theta)


    df = pd.DataFrame({'drifters':y})
    df.index.name = 'datetime'
    df['drifternomean'] = (df['drifters']-df['drifters'].rolling(160).mean())
    df['driftermean'] = (df['drifters'].rolling(160).mean())
    df['river'] = river
    df['zeta'] = zeta
    # imax = argrelextrema(zeta.values, np.greater)[0]  # local maxima for zeta
    # imin = argrelextrema(zeta.values, np.less)[0]  # local minima for zeta
    # df['flood'] = zeta[]
    # df['zetap'] = zeta[zeta>0]
    # df['zetam'] = zeta[zeta<0]
    df['u'] = uwind
    df['v'] = vwind
    df['s'] = s
    df['theta'] = theta

    name = 'calcs/df_' + start[:7]
    if 'forward' in File:
        name += '_forward'
    elif 'backward' in File:
        name += '_backward'
    df.to_csv(name + '.csv')  # save every time step
    # df = pd.read_csv(name, parse_dates=True, index_col=0)

    # # model = ols("drifters ~ theta", df).fit()
    # # print(model.rsquared)
    # # model = ols("drifters ~ s", df).fit()
    # # print(model.rsquared)
    # # model = ols("drifters ~ u", df).fit()
    # # print(model.rsquared)
    # # model = ols("drifters ~ v", df).fit()
    # # print(model.rsquared)
    # # model = ols("drifters ~ zeta", df).fit()
    # # print(model.rsquared)
    # model = ols("drifters ~ river + zeta + s + theta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + zeta + u + theta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + zeta + s", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + zeta + u", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + zeta + theta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + theta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river + zeta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifternomean ~ zeta", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("drifters ~ river", df).fit()
    # print(model.rsquared,model.bic)
    # model = ols("driftermean ~ river", df).fit()
    # print(model.rsquared,model.bic)
    #
    # # Peform analysis of variance on fitted linear model
    # anova_results = anova_lm(model)


if __name__ == "__main__":
    stats()
