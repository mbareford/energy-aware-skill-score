#!/usr/bin/env python

'''
This script is intended for calculating an energy-aware skill score for a given climate model.

Separate scores are calculated, one for each property of interest, e.g. temperature, wind speed, humidity.

The script will therefore need to take a number of arguments; climate model identifier, energy consumed by model run,
property list, start and end dates, time step and ranges for the longitude and latitude (and possibly the atmospheric levels too).

The longlat ranges define a spatial volume over which the properties are averaged.
And so, for each property in the list, the script calculates a time series of property values extracted from the model output.

The start and end of this series is defined by the start and end date arguments.
The interval between each element in the series is set by the time step argument.

The same property time series are calculated but using data from the ERA5 archive.
This dataset provides a ground truth. Fortunately, the ERA5 archive can be accessed using the get_jasmin_era5 Python package.

The calculation of the skill score for each property requires us to compare the model outputs with the ERA5 ground truth.
Specifically, the script calculates from the ERA5 data, two sets of values per property: the first set represents the MASE scaling
and second is the ground truth autocorrelation. Each element in these sets corresponds to a different lag value.
The lag values are given as a command line argument.


This script is under development: so far, I've implemented the functions to calculate the skill scores.

The ground truth is provided by the ERA5 data accessed via the get_jasmin_era5 Python package.
This package reads the ERA5 data from the netCDF files stored off the "/badc/ecmwf-era5/data/oper/an_ml" path on the JASMIN server.

Example model output comes from the CMIP6 HadGEM3-GC31 model, specifically, the netcdf file containing 3-hourly air surface temperatures.

Michael Bareford
m.bareford@epcc.ed.ac.uk
Mar 2026
'''



import argparse
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from get_jasmin_era5 import Find_era5
from cdo import *



class MappingFunctions:

    MAP_FUNC_DEFAULT=1
    def default(x):
        return x / (x + 1)

    MAP_FUNC_EXPONENT=2 
    def exponent(x):
        return 1 - math.exp(-x)

    MAP_FUNC_INVERSE_TAN=3
    def inverse_tan(x):
        return 2*math.pi*math.atan(x)

    MAP_FUNC_HYPERBOLIC_TAN=4
    def hyperbolic_tan(x):
        return math.tanh(x)

    MAP_FUNC_SQUARE_ROOT=5
    def square_root(x):
        return math.sqrt(x) / (math.sqrt(x) + 1)

    MAP_FUNC_SQUARE=6
    def square(x):
        return x^2 / (x^2 + 1)
    
    def get_function(this, func_type):
        mfunc = MappingFunctions.default

        match func_type:
            case MappingFunctions.MAP_FUNC_EXPONENT:
                mfunc = MappingFunctions.exponent
            case MappingFunctions.MAP_FUNC_INVERSE_TAN:
                mfunc = MappingFunctions.inverse_tan
            case MappingFunctions.MAP_FUNC_HYPERBOLIC_TAN:
                mfunc = MappingFunctions.hyperbolic_tan
            case MappingFunctions.MAP_FUNC_SQUARE_ROOT:
                mfunc = MappingFunctions.square_root
            case MappingFunctions.MAP_FUNC_SQUARE:
                mfunc = MappingFunctions.square
            case _:
                mfunc = MappingFunctions.default

        return mfunc




# Calculate the Mean Absolute Scaled Error (MASE) scalings for a ground truth ('gt')
# with respect to the supplied autocorrelation lag values ('h').
#
# Assume 'gt' is a numpy array containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of MASE scalings based on 'gt' where the lag value is the key.
# The 'lags' parameter is used to decide which scalings are stored in the returned dictionary.
def calc_mase_scalings(gt, gt_interval, lags):

    scalings = {}
    gt_len = len(gt)

    for h in lags:

        if h % gt_interval != 0:
            # lag value is not a multiple of the ground truth interval
            continue
        
        if h in range(1,(gt_len*gt_interval)//2):
            gt_subset = gt[0:gt_len:h//gt_interval]

            n = len(gt_subset)

            d_sum = 0.0
            for i in range(1,n):
                d_sum += abs(gt_subset[i]-gt_subset[i-1])

            scalings[h] = d_sum / (n-1)

        else:
            print('calc_mase_scalings: Lag value ' + str(h) + ' out of range for ground truth (1-'+str(gt_len*gt_interval)+').')

    return scalings


# Calculate the autocorrelations for a ground truth ('gt') with respect to the supplied lag values.
#
# Assume 'gt' is a numpy array containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of autocorrelations based on 'gt' where the lag value is the key.
# The 'lags' parameter is used to decide which autocorrelations are stored in the returned dictionary.
def calc_auto_correlations(ds, ds_interval, lags):

    auto_corrs = {}
    
    ds_len = len(ds)

    corrs = np.correlate(ds, ds, mode='full')
    corrs = corrs[corrs.size//2:] 
    corrs = corrs[:ds_len//2 - 1]
    corrs = corrs / corrs.max()

    for h in lags:
    
        if h % ds_interval != 0:
            # lag value is not a multiple of the dataset interval
            continue

        if  (h/ds_interval)-1 in range(corrs.size):
            auto_corrs[h] = corrs[(h//ds_interval)-1]
        else:
            print('calc_auto_correlations: Lag value ' + str(h) + ' out of range for data series (1-'+str(ds_len*ds_interval)+').')
            continue

    return auto_corrs


# Calculate the Mean Absolute Scaled Error (MASE) for a model truth with respect to a
# ground truth ('gt') and the supplied autocorrelation lag values ('h').
#
# Assume 'mt' and 'gt' are numpy arrays containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of MASE values derived from comparing 'mt' with 'gt'; the lag value is the key.
# The 'lags' parameter is used to decide which error values are stored in the returned dictionary.
def calc_mase_errors(mt, gt, interval, mase_scalings, lags):

    errors = {}

    mt_len = len(mt)
    gt_len = len(gt)

    if mt_len == gt_len:
        # assume model truth and ground truth series have same interval

        for h in lags:

            if h in range(1,(mt_len*interval)//2):

                mt_subset = mt[0:mt_len:h//interval]
                gt_subset = gt[0:gt_len:h//interval]

                n = len(mt_subset)

                e_sum = 0.0
                for i in range(1,n):
                    e_sum += abs(gt_subset[i]-mt_subset[i]) / mase_scalings[h]

                errors[h] = e_sum / n

            else:
                print('calc_mase_errors: Lag value ' + str(h) + ' out of range for model/ground truth (1-'+str(mt_len*interval)+').')
    else:
        print('calc_mase_errors: The number of elements in model truth (' + str(mt_len) + ') does not match ground truth (' + str(gt_len) + ').')

    return errors


# Calculate the skill scores by combining MASE values and ground truth autocorrelations that
# correspond to a range of lag values.
def calc_skill_scores(mase_errors, auto_corrs, map_func, lags):

    scores = {}

    for h in lags:

        #scores[h] = (1 - (mase_errors[h] / (mase_errors[h] + 1))) * (1 - abs(auto_corrs[h]))
        scores[h] = (1 - map_func(mase_errors[h])) * (1 - abs(auto_corrs[h]))

    return scores



'''
An example call of calc-skill-score.py for the current version.

The longitudinal and latitudinal ranges describe a region that covers part of continental Europe,
that runs from the south east corner of France up to north Netherlands, across to the eastern Belarus
and down to east Romania.

./calc-skill-score.py --properties '2t' --start-date '2020-06-01' --end-date '2020-07-01' --time-step '3h' --longitudes 6 28 --latitudes 46 53 --lags 6 12 24 168

./calc-skill-score.py --properties '2t' \
                      --start-date '2020-06-01' \
                      --end-date '2020-07-01' \
                      --time-step '3h' \
                      --longitudes 6 28 \
                      --latitudes 46 53 \
                      --lags 6 12 24 168 \
                      --map-func-type 1
'''

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# TODO: add arguments that identify the climate model and how to access its output on the JASMIN server.
# TODO: add argument that gives the energy consumed during the model run

parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-p', '--properties', nargs='+', default='t', help='ERA5 property list: format "<p1> <p2> ... <pn>"')
parser.add_argument('-d', '--model', type=str, default='CMIP-HadGEM3', help='The climate model to be evaluated.')
parser.add_argument('-s', '--start-date', type=str, default='2020-06-01', help='The start date: format is "YYYY-MM-DD"')
parser.add_argument('-e', '--end-date', type=str, default='2020-06-02', help='The end date: format is "YYYY-MM-DD"')
parser.add_argument('-t', '--time-step', type=str, default='1h', help='Time step')
parser.add_argument('-o', '--longitudes', nargs='+', type=int, default='0 360', help='The longitudinal range')
parser.add_argument('-a', '--latitudes', nargs='+', type=int, default='-90 90', help='The latitudinal range')
parser.add_argument('-g', '--lags', nargs='+', type=int, default='1', help='The list of autocorrelation lag values')
parser.add_argument('-m', '--map-func-type', type=int, default=1, help='The type of mapping function used to ensure that the skill score is in the range 0-1.')


# Parse command line arguments
args = parser.parse_args()

time_interval = int(args.time_step[:-1])

long_min = args.longitudes[0]
long_max = args.longitudes[1]

lat_min = args.latitudes[0]
lat_max = args.latitudes[1]

lags = args.lags

map_func = MappingFunctions().get_function(args.map_func_type)
    

fig = plt.gcf()


# Prepare CDO tmp directory
tempPath = '~/tmp/'
cdo = Cdo(tempdir=tempPath)
cdo.cleanTempDir()


# Locate the ERA5 archive and assign the result to ground truth archive
gt_arc = Find_era5()


# Extract from the ERA5 archive the required property values over the required time frame and spatial region.
gt_ds = gt_arc[args.properties,
               args.start_date:args.end_date:args.time_step,
               None, # ignore pressure levels
               long_min:long_max,
               lat_min:lat_max]


# Iterate over the properties, averaging the property values over the longitude and latitude ranges.
for prop_name, prop_values in gt_ds.items():

    print('Calculating skill scores for property ' + str(prop_name) + '...')

    gt_avg = prop_values.mean(dim=['longitude', 'latitude'])

    gt_vals = np.float32(gt_avg.values)

    mase_scalings = calc_mase_scalings(gt_vals, time_interval, lags)
    print('MASE Scalings\n' + str(mase_scalings) + '\n\n')

    auto_corrs = calc_auto_correlations(gt_vals, time_interval, lags)
    print('Autocorrelations\n' + str(auto_corrs) + '\n\n')


    # hardcode path to HadGEM model output file
    model_out_fn = '/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-MM/historical/r1i1p1f3/3hr/tas/gn/latest/tas_3hr_HadGEM3-GC31-MM_historical_r1i1p1f3_gn_200001010300-200101010000.nc'
    # hardcode CMIP variable selection - ERA5 variable '2t' maps to CMIP variable 'tas'
    var_sel = 'tas'
    # hardcode time series index selection - 153/392 covers the range 2020-06-01 to 2020-06-30 inclusive from the CMIP HadGEM file for the year 2000, i.e. '200001010300-200101010000'.
    tmi_sel = '153/392'


    lonlat_extent = str(long_min) + ',' + str(long_max) + ',' + str(lat_min) + ',' + str(lat_max)
    cdo_in_sel = '-selvar,' + var_sel + ' -seltimestep,' + tmi_sel + ' ' + model_out_fn

    mt_ds = cdo.sellonlatbox(lonlat_extent, input=cdo_in_sel, returnXArray='tas')

    mt_avg = mt_ds.mean(dim=['lat','lon'])

    # assume model truth is in 32-bit float precision
    mt_vals = mt_avg.values


    mase_errors = calc_mase_errors(mt_vals, gt_vals, time_interval, mase_scalings, lags)
    print('MASE Errors\n' + str(mase_errors) + '\n\n')

    skill_scores = calc_skill_scores(mase_errors, auto_corrs, map_func, lags)
    print('Skill Scores\n' + str(skill_scores) + '\n\n')

    lag_labels = [ '6 hrs', '12 hrs', '1 day', '1 wk' ]
    plt.bar(lag_labels, skill_scores.values(), color='#66c2a5')
    
    plt.xlabel('Autocorrelation Lag')
    plt.ylabel('Skill Score')
    plt.title(r'CMIP6 HadGEM3: Air Surface Temperature Skill Score')

    fig.savefig('cmip6-hadgem3-tas-skill-score.jpg', format='jpg', dpi=1000)

    plt.clf()

'''
TODO

1. Locate a model output archive (e.g. from within "/badc/cmip6/data/CMIP6/HighResMIP/NERC/HadGEM3-GC31-HH") on the JASMIN server.
2. Extract the property values according to the start and end dates, time step and atmospheric level, latitude and longitude ranges.
3. Average the property values over the indicated spatial volume.
4. Calculate the skill score for each property for each lag value passed via the command line.
5. Make these scores energy aware by dividing each score by the energy consumed by the model run.
6. Plot the scores by property and by lag value.

/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-MM/historical/r1i1p1f3/3hr/tas/gn/latest/tas_3hr_HadGEM3-GC31-MM_historical_r1i1p1f3_gn_200001010300-200101010000.nc -> ../files/d20200720/tas_3hr_HadGEM3-GC31-MM_historical_r1i1p1f3_gn_200001010300-200101010000.nc
'''


