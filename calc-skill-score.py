#!/usr/bin/env python

'''
This script is intended for calculating an energy-aware skill score for a given climate model.

Separate scores are calculated, one for each property of interest, e.g. temperature, precipitation, humidity.

The script will therefore need to take a number of arguments; climate model identifier, energy consumed by model run,
property list, start and end dates, time step and ranges for the atmospheric levels, longitude and latitude.

The last three arguments define a spatial volume over which the properties are averaged.
And so, for each property in the list, the script calculates a time series of property values extracted from the model output.

The start and end of this series is defined by the start and end date arguments.
The interval between each element in the series is set by the time step argument.

The same property time series are calculated but this time using data from the ERA5 archive.
This dataset provides a ground truth. Fortunately, the ERA5 archive can be accessed using the get_jasmin_era5 Python package.

The calculation of the skill score for each property requires us to compare the model outputs with the ERA5 ground truth.
Specifically, the script calculates from the ERA5 data, two sets of values per property: the first set represents the MASE scaling
and second is the ground truth autocorrelation. Each element in these sets corresponds to a different lag value.
The lag values are given as a command line argument.


This script is under development: so far, I've implemented the functions to calculate the MASE scalings and autocorrelations
from the ground truth ERA5 data accessed via the get_jasmin_era5 Python package. This package reads the ERA5 data from the
netCDF files stored off the "/badc/ecmwf-era5/data/oper/an_ml" path on the JASMIN server.

Michael Bareford
m.bareford@epcc.ed.ac.uk
Feb 2026
'''



import argparse
import math
import numpy as np

from get_jasmin_era5 import Find_era5



# Calculate the Mean Absolute Scaled Error (MASE) for a ground truth ('gt')
# with respect to the supplied autocorrelation lag values ('h').
#
# Assume that 'gt' is a numpy array containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of MASE scalings based on "gt" where the lag value is the key.
# The "lags" parameter is used to decide which scalings are stored in the returned dictionary.
def calc_mase_scalings(gt, lags):

    scalings = {}
    gt_len = len(gt)

    for h in lags:

        if h in range(1,gt_len//2):
            gt_subset = gt[0:gt_len:h]

            n = len(gt_subset)

            d_sum = 0.0
            for i in range(1,n):
                d_sum += abs(gt_subset[i]-gt_subset[i-1])

            scalings[h] = d_sum / (n-1)

        else:
            print('calc_mase_scalings: Lag value ' + str(h) + ' out of range for ground truth (1-'+str(gt_len)+').')

    return scalings



# Calculate the autocorrelations for a data series ('ds') with respect to the supplied lag values.
# The data series represents a ground truth.
#
# Assume that 'ds' is a numpy array containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of autocorrelations based on "ds" where the lag value is the key.
# The "lags" parameter is used to decide which autocorrelations are stored in the returned dictionary.
def calc_auto_correlations(ds, lags):

    corrs = np.correlate(ds, ds, mode='full')
    corrs = corrs[corrs.size//2:] 
    corrs = corrs[:len(ds)//2 - 1]
    corrs = corrs / corrs.max()

    auto_corrs = {}

    for h in lags:
    
        if  h-1 in range(corrs.size):
            auto_corrs[h] = corrs[h-1]
        else:
            print('calc_auto_correlations: Lag value ' + str(h) + ' out of range for data series (1-'+str(len(ds))+').')
            continue

    return auto_corrs



'''
An example call of calc_skill_score.py for the current version.

./calc_skill_score.py --properties 't' \
                      --start-date '2020-06-01' \
                      --end-date '2020-06-02' \
                      --time-step '1h' \
                      --levels  100 137 \
                      --longitudes 90 140 \
                      --latitudes -20 20 \
                      --lags 1 6 10
'''

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# TODO: add arguments that identify the climate model and how to access its output on the JASMIN server.
# TODO: add argument that gives the energy consumed during the model run

parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-p', '--properties', nargs='+', default='t', help='Property list: format "<p1> <p2> ... <pn>"')
parser.add_argument('-s', '--start-date', type=str, default='2020-06-01', help='The start date: format is "YYYY-MM-DD"')
parser.add_argument('-e', '--end-date', type=str, default='2020-06-02', help='The end date: format is "YYYY-MM-DD"')
parser.add_argument('-t', '--time-step', type=str, default='1h', help='Time step')
parser.add_argument('-l', '--levels', nargs='+', type=int, default='1 137', help='ECMWF L137 model level range')
parser.add_argument('-o', '--longitudes', nargs='+', type=int, default='0 360', help='The longitudinal range')
parser.add_argument('-a', '--latitudes', nargs='+', type=int, default='-90 90', help='The latitudinal range')
parser.add_argument('-g', '--lags', nargs='+', type=int, default='1', help='The list of autocorrelation lag values')


args = parser.parse_args()


# Locate the ERA5 archive and assign the result to ground truth archive
gt_arc = Find_era5()


# Example code showing how to extract a subset of data.
#ds = gt_arc[('t'), '2020-06-01':'2020-06-02':'1h', 100:137, 90:140, -20:20]

# Extract from the ERA5 archive the required property values over the required time frame and spatial region.
gt_ds = gt_arc[args.properties,
               args.start_date:args.end_date:args.time_step,
               args.levels[0]:args.levels[1],
               args.longitudes[0]:args.longitudes[1],
               args.latitudes[0]:args.latitudes[1]]


# Iterate over the properties and for each property average the values over the atmospheric level, latitude and longitude ranges.
for prop_name, prop_values in gt_ds.items():

    gt_avg = prop_values.mean(dim=['level', 'longitude', 'latitude'])

    mase_scalings = calc_mase_scalings(gt_avg.values, args.lags)
    print('MASE Scalings\n' + str(mase_scalings) + '\n\n')

    auto_corrs = calc_auto_correlations(gt_avg.values, args.lags)
    print('Autocorrelations\n' + str(auto_corrs) + '\n\n')



'''
TODO

1. Locate a model output archive (e.g. from within "/badc/cmip6/data/CMIP6/HighResMIP/NERC/HadGEM3-GC31-HH") on the JASMIN server.
2. Extract the property values according to the start and end dates, time step and atmospheric level, latitude and longitude ranges.
3. Average the property values over the indicated spatial volume.
4. Calculate the skill score for each property for each lag value passed via the command line.
5. Make these scores energy aware by dividing each score by the energy consumed by the model run.
6. Plot the scores by property and by lag value.
'''


