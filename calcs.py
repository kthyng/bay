'''
Calculation functions for bay work.
'''

from glob import glob
import pandas as pd
import read
import itertools
import numpy as np
from statsmodels.formula.api import ols
import statsmodels.api as sm
tsa = sm.tsa
from matplotlib.path import Path
from datetime import datetime, timedelta
import os
import netCDF4 as netCDF


baypathll = np.load('calcs/pathbayandjetty.npz', encoding='latin1')['pathll'].item()

def io(ind=-1):
    '''Calculates - in time - the number of drifters inside the bay.

    This is determined using baypath.'''

    basename = 'newsuperposition/'
    if not os.path.exists('calcs/enterexit/' + basename):
        os.makedirs('calcs/enterexit/' + basename)
    refdate = datetime(2010, 7, 1, 0, 0)  # datetime(2010, 7, 15, 0, 0) datetime(2010, 2, 1, 0, 0) datetime(2010, 7, 1, 0, 0)

    if refdate.day == 15:
        direct = '_backward'
    elif refdate.day == 1:
        direct = '_forward'

    if 'forward' in direct:
        start = refdate
        end = start + timedelta(days=14*2)
    elif 'backward' in direct:
        start = refdate - timedelta(days=14)
        end = start + timedelta(days=2*14)
    Files = glob('tracks/' + basename + start.isoformat()[:7] + '*' + direct + '*.nc')
    dfdates = pd.date_range(start=start.isoformat()[:10] + ' 00:00:01', end=end.isoformat()[:10] + ' 00:00:01', freq='900S')
    df = pd.DataFrame(index=dfdates)
    for i, File in enumerate(Files):

        # short-circuit if ind is not -1 and just use that index instead
        if ind != -1:
            if i != ind:
                continue
        # print(File)
        d = netCDF.Dataset(File)
        lonp = d['lonp'][:]; latp = d['latp'][:]; tp = d['tp'][0,:]
        d.close()
        # points that are inside of bay
        inds = baypathll.contains_points(np.vstack([lonp.flatten(), latp.flatten()]).T).reshape(lonp.shape)
        # import pdb; pdb.set_trace()

        # save indices of drifters that are outside bay
        outsideinfo = np.where(~inds)  # tuple: [drifter index, time] for drifter outside bay
        idrifters = set(outsideinfo[0])  # indices of drifters that are outside bay at some point (overall)
        numoutside = (~inds).sum(axis=0)  # number of drifter outside of bay by time
        # import pdb; pdb.set_trace()
        # add column to dataframe with this information
        simstartdate = File.split('/')[-1].split('14days')[0][:-1]

        dftemp = pd.DataFrame(index=netCDF.num2date(tp, 'seconds since 1970-01-01  00:00:00'), data={simstartdate: numoutside})
        df = df.join(dftemp)  # add column to dataframe

        if ind == -1:
            df.to_csv('calcs/enterexit/' + basename + start.isoformat()[:7] + '.csv')  # save every time step to update
        else:
            df.to_csv('calcs/enterexit/' + basename + File.split('/')[-1][:13] + '.csv')  # save every time step to update
        # the following is for plotting drifters that exit domain
        np.savez('calcs/enterexit/' + basename + simstartdate[:13] + '.npz', idrifters=idrifters, outsideinfo=outsideinfo)


def make_dfs(which='2010-07.csv'):
    '''Make dataframes between drifters and forcing mechanisms for running stats.'''

    basename = 'newsuperposition/'
    Files = sorted(glob('calcs/enterexit/' + basename + which))
    for File in Files:
        # File = 'calcs/enterexit_sim3_2010-07_backward_14days_dx300.csv'
        df = pd.read_csv(File, parse_dates=True, index_col=0)
        # nfiles = (~np.isnan(df)).astype(int).sum(axis='columns')
        # y = df.sum(axis='columns').divide(nfiles).resample('15min', base=0).interpolate()
        # delta drifters
        # import pdb; pdb.set_trace()
        # delta drifters, averaged over file and smoothed
        y = df.diff(axis=0).mean(axis=1).rolling(window=16).mean()
        ind = np.isnan(y)
        y = y.resample('15T', base=0).interpolate()  # moving away from 1 minute base
        # insert nan's back
        y[np.where(ind)[0]] = np.nan

        start = df.index[0].isoformat()[:10]
        stop = df.index[-1].isoformat()[:10]

        river = read.from_suntans(start, stop).resample('15min').interpolate()['boundary_Q'][:y.index[-1].isoformat()]
        # import pdb; pdb.set_trace()
        zeta = read.from_blended('zeta', start, stop).to_dataframe().resample('15min').interpolate()['zeta'][:y.index[-1].isoformat()]
        uwind = read.from_shelf('Uwind', start, stop).to_dataframe().resample('15min').interpolate()['Uwind'][:y.index[-1].isoformat()]
        vwind = read.from_shelf('Vwind', start, stop).to_dataframe().resample('15min').interpolate()['Vwind'][:y.index[-1].isoformat()]
        s = np.sqrt(uwind**2 + vwind**2)
        theta = np.unwrap(np.arctan2(vwind, uwind))
        sustr = read.from_shelf('sustr', start, stop).to_dataframe().resample('15min').interpolate()['sustr'][:y.index[-1].isoformat()]
        svstr = read.from_shelf('svstr', start, stop).to_dataframe().resample('15min').interpolate()['svstr'][:y.index[-1].isoformat()]


        df = pd.DataFrame({'drifters':y})
        df.index.name = 'datetime'
        df['drifters_tidal'] = (df['drifters']-df['drifters'].rolling(192, center=True).mean())
        df['drifters_subtidal'] = (df['drifters'].rolling(192, center=True).mean())
        df['river'] = river
        df['zeta'] = zeta
        df['dzeta'] = df['zeta'].diff()  # derivative of zeta - proxy for current

        df['uwind'] = uwind
        df['vwind'] = vwind
        df['s'] = s
        df['theta'] = theta
        df['dtheta'] = df['theta'].diff()  # change in wind direction
        df['sustr'] = sustr
        df['svstr'] = svstr

        # add some columns on for other analysis

        # What if river discharge is related to drifters but with a time lag?
        # istart, iend are indices in time series where subtidal signal starts and
        # ends being available due to filtering making nan's
        istart = np.where(np.diff(df['drifters_subtidal'].fillna(-999) == -999))[0][0] + 1
        iend = np.where(np.diff(df['drifters_subtidal'].fillna(-999) == -999))[0][1] + 1
        ccf_river = tsa.ccf(df['drifters_subtidal'][istart:iend], df['river'][istart:iend])
        # this gives the shift forward necessary for the drifters to correlate best
        # with the river discharge
        imax = ccf_river.argmax()
        # this column has accounted for shift in river discharge by pushing river time forward
        df['river_shifted'] = df['river'].shift(imax)

        # found no correlation shift for feb 2010
        # ccf_uwind = tsa.ccf(df['drifters_subtidal'][istart:iend], df['uwind'][istart:iend])
        # imax = ccf_uwind.argmax()
        # df['uwind_shifted'] = df['uwind'].shift(imax)

        ccf_theta = tsa.ccf(df['drifters_subtidal'][istart:iend], df['theta'][istart:iend])
        imax = ccf_theta.argmax()
        df['theta_shifted'] = df['theta'].shift(imax)

        ipos = df['dzeta'] >= 0
        df['dzeta_floor0'] = df['dzeta'][ipos]

        ccf_zeta = tsa.ccf(df['drifters_tidal'][istart:iend], df['zeta'][istart:iend])
        imax = ccf_zeta.argmax()
        df['zeta_shifted'] = df['zeta'].shift(imax)

        ipos = df['zeta'] >= 0
        df['zeta_floor0'] = df['zeta'][ipos]

        name = 'calcs/enterexit/' + basename + 'df_' + which
        # name = 'calcs/enterexit/' + basename + 'df_' + File.split('/')[-1].split('.')[0] + start[:7]
        if 'forward' in File:
            name += '_forward'
        elif 'backward' in File:
            name += '_backward'
        df.to_csv(name)  # save every time step


def powerset(inlist):
    '''Returns list of all unique combinations of items in input list.

    Example:
    powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    '''

    s = list(inlist)
    return list(itertools.chain.from_iterable(itertools.combinations(s, r) for r in range(len(s)+1)))


def filter_powerset(df, which):
    '''Find appropriate combinations of mechanisms.'''

    if which == 'subtidal':
        # names of forcing mechanism columns
        mechanisms = df.columns[['drifters' not in column and 'zeta' not in column for column in df.columns]]
        combos = powerset(mechanisms)[1:]  # skip empty

        # filter allowed combinations
        iremove = list(np.where([('river' in combo and 'river_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('theta' in combo and 'theta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('theta' in combo and 'dtheta' in combo) for combo in combos])[0])
        iremove.extend(np.where([('dtheta' in combo and 'theta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('uwind' in combo and 'sustr' in combo) or ('vwind' in combo and 'svstr' in combo) for combo in combos])[0])
        iremove.extend(np.where([('uwind' in combo and 'svstr' in combo) or ('vwind' in combo and 'sustr' in combo) for combo in combos])[0])
        # don't want uwind/vwind and sustr+svstr together
        iremove.extend(np.where([(('uwind' in combo or 'vwind' in combo) and ('sustr' in combo and 'svstr' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo and 'vwind' in combo) and ('sustr' in combo or 'svstr' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo or 'vwind' in combo) and ('s' in combo and 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo and 'vwind' in combo) and ('s' in combo or 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('sustr' in combo or 'svstr' in combo) and ('s' in combo and 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('sustr' in combo and 'svstr' in combo) and ('s' in combo or 'theta' in combo)) for combo in combos])[0])

        # remove these instances
        combos = np.delete(np.asarray(combos), list(set(iremove)))

    elif which == 'tidal':

        mechanisms = df.columns[['drifters' not in column and 'zeta' in column for column in df.columns]]
        combos = powerset(mechanisms)[1:]  # skip empty

        iremove = list(np.where([('dzeta' in combo and 'dzeta_floor0' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta' in combo and 'zeta_floor0' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta' in combo and 'zeta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta_floor0' in combo and 'zeta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta_shifted' in combo and 'zeta_shifted_floor0' in combo) for combo in combos])[0])
        # remove these instances
        combos = np.delete(np.asarray(combos), list(set(iremove)))

    elif which == 'drifters':

        mechanisms = df.columns[['drifters' not in column and 'tidal' not in column and 'subtidal' not in column for column in df.columns]]
        combos = powerset(mechanisms)[1:]  # skip empty
        # import pdb; pdb.set_trace()
        # filter allowed combinations
        iremove = list(np.where([('river' in combo and 'river_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('theta' in combo and 'theta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('theta' in combo and 'dtheta' in combo) for combo in combos])[0])
        iremove.extend(np.where([('dtheta' in combo and 'theta_shifted' in combo) for combo in combos])[0])
        iremove.extend(np.where([('uwind' in combo and 'sustr' in combo) or ('vwind' in combo and 'svstr' in combo) for combo in combos])[0])
        iremove.extend(np.where([('uwind' in combo and 'svstr' in combo) or ('vwind' in combo and 'sustr' in combo) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo or 'vwind' in combo) and ('sustr' in combo and 'svstr' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo and 'vwind' in combo) and ('sustr' in combo or 'svstr' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo or 'vwind' in combo) and ('s' in combo and 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('uwind' in combo and 'vwind' in combo) and ('s' in combo or 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('sustr' in combo or 'svstr' in combo) and ('s' in combo and 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([(('sustr' in combo and 'svstr' in combo) and ('s' in combo or 'theta' in combo)) for combo in combos])[0])
        iremove.extend(np.where([('dzeta' in combo and 'dzeta_floor0' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta' in combo and 'zeta_floor0' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta_shifted' in combo and 'zeta_floor0' in combo) for combo in combos])[0])
        # iremove.extend(np.where([('zeta' in combo and 'dzeta' in combo) for combo in combos])[0])
        iremove.extend(np.where([('zeta' in combo and 'zeta_shifted' in combo) for combo in combos])[0])
        # iremove.extend(np.where([('dzeta' in combo) for combo in combos])[0])
        # iremove.extend(np.where([('dzeta_floor0' in combo) for combo in combos])[0])
        # iremove.extend(np.where([('zeta_floor0' in combo) for combo in combos])[0])
        # remove these instances
        combos = np.delete(np.asarray(combos), list(set(iremove)))
    # import pdb; pdb.set_trace()

    return combos


def transform(df, which):
    '''Transform forcing mechanism using outside information beyond just signal.'''

    if 'river' in which:
        varfull = pd.read_csv('calcs/river2009-2011.csv', parse_dates=True, index_col=0)['108.394273505']
    elif 'uwind' in which:
        varfull = pd.read_csv('calcs/uwind2009-2011.csv', parse_dates=True, index_col=0)['-6.20427322388']
    elif 'zeta' in which:
        varfull = pd.read_csv('calcs/zeta2009-2011.csv', parse_dates=True, index_col=0)['0.0549901995635']
    elif 'vwind' in which:
        varfull = pd.read_csv('calcs/vwind2009-2011.csv', parse_dates=True, index_col=0)['-2.68748617172']
    elif 'sustr' in which:
        varfull = pd.read_csv('calcs/sustr2009-2011.csv', parse_dates=True, index_col=0)['-0.0590621456504']
    elif 'svstr' in which:
        varfull = pd.read_csv('calcs/svstr2009-2011.csv', parse_dates=True, index_col=0)['-0.0251256525517']
    elif which == 's':
        varfull1 = pd.read_csv('calcs/uwind2009-2011.csv', parse_dates=True, index_col=0)['-6.20427322388']
        varfull2 = pd.read_csv('calcs/vwind2009-2011.csv', parse_dates=True, index_col=0)['-2.68748617172']
        varfull = np.sqrt(varfull1**2 + varfull2**2)
    elif 'theta' in which:
        varfull1 = pd.read_csv('calcs/uwind2009-2011.csv', parse_dates=True, index_col=0)['-6.20427322388']
        varfull2 = pd.read_csv('calcs/vwind2009-2011.csv', parse_dates=True, index_col=0)['-2.68748617172']
        varfull = np.unwrap(np.arctan2(varfull1, varfull2))
    # the following are specific to one case only
    elif which == 'drifters_subtidal':
        Files = glob('calcs/enterexit/newsuperposition/df_*.csv')
        varfull = []
        for File in Files:
            varfull.extend(pd.read_csv(File, parse_dates=True, index_col=0)['drifters_subtidal'])
        varfull = np.asarray(varfull)
    elif which == 'drifters_tidal':
        Files = glob('calcs/enterexit/newsuperposition/df_*.csv')
        varfull = []
        for File in Files:
            varfull.extend(pd.read_csv(File, parse_dates=True, index_col=0)['drifters_tidal'])
        varfull = np.asarray(varfull)
    elif which == 'drifters':
        Files = glob('calcs/enterexit/newsuperposition/df_*.csv')
        varfull = []
        for File in Files:
            varfull.extend(pd.read_csv(File, parse_dates=True, index_col=0)['drifters'])
        varfull = np.asarray(varfull)

    # import pdb; pdb.set_trace()
    if np.isnan(varfull).sum()>0:
        return (df[which] - np.nanmean(varfull))/np.nanstd(varfull)
    else:
        return (df[which] - varfull.mean())/varfull.std()


def scaled(df, combo, which='subtidal'):
    '''Wrapper to create scaled version of dataframe with desired columns.'''

    dfscaled = pd.DataFrame()
    for item in list(combo):
        dfscaled[item] = transform(df, which=item)  # scale using years-long data
    if 'tidal' in which:
        item = 'drifters_' + which
    else:
        item = 'drifters'
    dfscaled[item] = (df[item] - df[item].mean())/df[item].std()

    return dfscaled


def stats(which='subtidal', direction='forward', simname='newsuperposition'):
    '''Calculate stats.

    Maybe only do subtidal since tidal appears to have no impact on how many
    drifters stay outside.'''

    # Files = glob('calcs/df_????-??_*ward.csv')
    Files = glob('calcs/enterexit/' + simname + '/df_????-??.csv')
    # Files = glob('calcs/enterexit/' + simname + '/df_????-??_' + direction + '.csv')
    # Files = glob('calcs/df_2010-02_forward.csv')
    # import pdb; pdb.set_trace()
    for File in Files:
        # print(File)
        df = pd.read_csv(File, parse_dates=True, index_col=0)

        # which = 'subtidal'
        combos = filter_powerset(df, which=which)

        models = []; dfscaleds = []
        for combo in combos:

            # run ordinary least squares analysis on drifter column vs. combination of mechanisms
            # create temporary new dataframe for this loop's analysis where all
            # variables are scaled
            dfscaled = scaled(df, combo, which=which)
            # import pdb; pdb.set_trace()
            if 'tidal' in which:
                models.append(ols('drifters_' + which + ' ~ ' + " + ".join(list(combo)), dfscaled).fit())
            else:
                models.append(ols('drifters' + ' ~ ' + " + ".join(list(combo)), dfscaled).fit())
            dfscaleds.append(dfscaled)

        # find top N r^2, lowest BIC values
        N = 50
        ir2 = np.argsort(-np.asarray([model.rsquared_adj for model in models]))[:N]
        ibic = np.argsort(np.asarray([model.bic for model in models]))[:N]
        # indstemp = list(set(np.concatenate((ir2, ibic))))  # has to be both highest r^2 and lowest BIC
        inds = []
        for ind in ir2:
            if (models[ind].pvalues[1:]>0.1).sum() > 0 or models[ind].bic < 0:
                pass
            else:
                inds.append(ind)
        if len(inds) > 5:
            inds = inds[:5]

        # fig, axes = plt.subplots(N, 1, sharex=True, figsize=(16,10))
        # Summary
        print('')
        print('Simulation set: %s' % File)
        print('Number of combinations checked: %i' % len(combos))
        print('Top adjusted r^2, lowest BIC performers, no p>0.1: %i' % len(inds))
        for i, ind in enumerate(inds):
            model = models[ind]
            dfscaled = dfscaleds[ind]
            print('Adjusted r^2: %0.2f' % model.rsquared_adj)
            print('BIC: %d' % model.bic)
            print('coefficients: ')
            print(model.params.to_string())
            print('pvalues: ')
            print(model.pvalues[1:].to_string())
            print('')

        # plot
        # plt.figure()
        # model.fittedvalues.plot(color='r', ax=axes[i], label='fit')
        # model.resid.plot(color='g', ax=axes[i], label='residual')
        # dfscaled['drifters_subtidal'].plot(ax=axes[i])
        #     # dfscaled['drifters_' + which].plot(color='k', linewidth=3, ax=axes[i])
        #     # for key in model.params.keys()[1:]  # skips intercept
        #     #     dfscaled[key].plot(ax=axes[i])
        # fig.tight_layout()



    #
    # figure()
    # plot(df['drifters_subtidal'][istart:iend][imax:-imax], df['theta'][istart:iend][imax:-imax], 'g.')
    # plot(df['drifters_subtidal'][istart:iend][imax:-imax], df['theta'][istart:iend][imax+imax:], 'r.')
    # plot(df['drifters_subtidal'][istart:iend][imax+imax:], df['theta'][istart:iend][imax:-imax], 'b.')


# if __name__ == "__main__":
    # io()
    # make_dfs()
