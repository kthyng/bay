'''
Read functions for bay work.
'''

from glob import glob
import xarray as xr
import pandas as pd
import numpy as np


# loc = glob('/rho/raid/dongyu/blended*.nc')


def from_blended(var, start, stop, yrho=230, xrho=445):
    '''Get some information from blended model product.

    Default xr/yr location is in entrance to Galveston Bay.
    '''

    # assuming time periods are consecutive (two months in a row)
    loc = glob('/rho/raid/dongyu/blended' + start[:4] + start[5:7] + '*.nc')
    loc.extend(glob('/rho/raid/dongyu/blended' + stop[:4] + stop[5:7] + '*.nc'))
    # remove nonunique values (in case all within one month)
    d = xr.open_mfdataset(sorted(list(set(loc))), concat_dim='ocean_time')

    return d[var].sel(ocean_time=slice(start, stop)).isel(yr=yrho, xr=xrho)


def from_shelf(var, start, stop, xrho=274, yrho=165):
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
    elif var == 'sustr':
        var = d[var].sel(ocean_time=slice(start, stop)).isel(eta_u=yrho, xi_u=xrho)
    elif var == 'svstr':
        var = d[var].sel(ocean_time=slice(start, stop)).isel(eta_v=yrho, xi_v=xrho)
    # print(var)
    return var


def from_TWDB(start, stop):
    '''Get river data.

    from https://waterdatafortexas.org/coastal/sites/TRIN
    Not using this anymore
    '''

    df = pd.read_csv('TRIN_seawater_salinity_daily.csv', comment='#', parse_dates=True, index_col=0)
    return df[start:stop]


def from_suntans(start, stop):
    '''Get river discharge from suntans boundary file.'''

    # Always read in all and then remove from there. Files overlap so had to
    # deal with specially.
    years = np.array([2009, 2010, 2011])
    base = '/rho/raid/dongyu/GalvCoarse_BC_'
    river = pd.DataFrame()
    for year in years:
        ds = xr.open_dataset(base + str(year) + '.nc')
        ds.swap_dims({'Nt': 'time'}, inplace=True)  # so can index based on time
        river = river.append(ds['boundary_Q'].sel(time=slice(str(year) + '-01-01T00:00:00', str(year) + '-12-31T23:00:00')).sum(axis=1).to_dataframe())
    return river[start:stop]
