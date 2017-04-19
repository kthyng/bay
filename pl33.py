'''
Good tidal signal filter from Octant.
'''

import numpy as np


class plfilt(object):
    """
    pl33 filter class, to remove tides and intertial motions from timeseries

    USAGE:
    ------

    >>> pl33 = plfilt(dt=4.0)   # 33 hr filter
    >>> pl33d = plfilt(dt=4.0, cutoff_period=72.0)  # 3 day filter

    dt is the time resolution of the timeseries to be filtered in hours.  Default dt=1
    cutoff_period defines the timescales to low pass filter. Default cutoff_period=33.0

    Calling the class instance can have two forms:

    >>> uf = pl33(u)   # returns a filtered timeseries, uf.  Half the filter length is
                       # removed from each end of the timeseries

    >>> uf, tf = pl33(u, t)  # returns a filtered timeseries, uf, plus a new time
                             # variable over the valid range of the filtered timeseries.

    """

    _pl33 = np.array([-0.00027, -0.00114, -0.00211, -0.00317, -0.00427, -0.00537,
                      -0.00641, -0.00735, -0.00811, -0.00864, -0.00887, -0.00872,
                      -0.00816, -0.00714, -0.0056 , -0.00355, -0.00097,  0.00213,
                       0.00574,  0.0098 ,  0.01425,  0.01902,  0.024  ,  0.02911,
                       0.03423,  0.03923,  0.04399,  0.04842,  0.05237,  0.05576,
                       0.0585 ,  0.06051,  0.06174,  0.06215,  0.06174,  0.06051,
                       0.0585 ,  0.05576,  0.05237,  0.04842,  0.04399,  0.03923,
                       0.03423,  0.02911,  0.024  ,  0.01902,  0.01425,  0.0098 ,
                       0.00574,  0.00213, -0.00097, -0.00355, -0.0056 , -0.00714,
                      -0.00816, -0.00872, -0.00887, -0.00864, -0.00811, -0.00735,
                      -0.00641, -0.00537, -0.00427, -0.00317, -0.00211, -0.00114,
                      -0.00027], dtype='d')

    _dt = np.linspace(-33, 33, 67)

    def __init__(self, dt=1.0, cutoff_period=33.0):

        if np.isscalar(dt):
            self.dt = float(dt) * (33.0 / cutoff_period)
        else:
            self.dt = np.diff(dt).mean() * (33.0 / cutoff_period)

        filter_time = np.arange(0.0, 33.0, self.dt, dtype='d')
        self.Nt = len(filter_time)
        self.filter_time = np.hstack((-filter_time[-1:0:-1], filter_time))

        self.pl33 = np.interp(self.filter_time, self._dt, self._pl33)
        self.pl33 /= self.pl33.sum()

    def __call__(self, u, t=None, mode='valid'):
        uf = np.convolve(u, self.pl33, mode=mode)
        if t is None:
            return uf
        else:
            tf = t[self.Nt-1:-self.Nt+1]
            return uf, tf
