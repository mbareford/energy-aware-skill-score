#!/usr/bin/env python

'''
This script is intended for calculating an energy-aware skill score for the various properties computed
by a specified weather forecast code.

Separate scores are calculated, one for each property of interest, e.g. temperature, wind speed, humidity.

The script will therefore need to take a number of arguments; weather forecast identifier, energy consumed by code run,
property list, start and end dates, time step and ranges for the longitude and latitude (and possibly the atmospheric levels too).

The longlat ranges define a spatial volume over which the properties are averaged.
And so, for each property in the list, the script calculates a time series of property values extracted from the model output.

The start and end of this series is defined by the start and end date arguments.
The interval between each element in the series is set by the time step argument.

The same property time series are extracted from the ERA5 archive.
This dataset provides a ground truth.
Fortunately, the ERA5 archive can be accessed using the get_jasmin_era5 Python package.

The calculation of the skill score for each property requires us to compare the model outputs with the ERA5 ground truth.
The script calculates from the ERA5 data, two sets of values per property: the first set represents the MASE scaling
and second is the ground truth autocorrelation.
Each element in these sets corresponds to a different lag value.
The lag values are given as a command line argument.

The ground truth is provided by the ERA5 data accessed via the get_jasmin_era5 Python package.
This package reads the ERA5 data from the netCDF files stored off the "/badc/ecmwf-era5/data/oper/an_ml" path on the JASMIN server.

Example model output comes from the CMIP6 HadGEM3-GC31 model, specifically, the netcdf files containing 6-hourly ("6hrPlev") data, see the
"/badc/cmip6/data/CMIP6/CMIP/MOHC/HadGEM3-GC31-MM/historical/r1i1p1f3/6hrPlev/<variable_name>/gn/latest/<var_name>_6hrPlev_HadGEM3-GC31-MM_historical_r1i1p1f3_gn_200001010300-200012302100.nc"
path accessible from the JASMIN server.

TODO: add argument that gives the energy consumed during the model run

Michael Bareford
m.bareford@epcc.ed.ac.uk
Apr 2026
'''



import os
import argparse
import math
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from get_jasmin_era5 import Find_era5
from cdo import *


# A selection of mapping functions for ensuring that skill scores are in the range 0-1.
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
# The 'lag_values' parameter is used to decide which scalings are stored in the returned dictionary.
def calc_mase_scalings(gt, gt_interval, lag_values):

    scalings = {}
    gt_len = len(gt)

    for h in lag_values:

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
# The 'lag_values' parameter is used to decide which autocorrelations are stored in the returned dictionary.
def calc_auto_correlations(gt, gt_interval, lag_values):

    auto_corrs = {}
    
    gt_len = len(gt)
    corrs = np.correlate(gt, gt, mode='full')
    corrs = corrs[corrs.size//2:] 
    corrs = corrs / corrs.max()

    for h in lag_values:
    
        if h % gt_interval != 0:
            # lag value is not a multiple of the ground truth interval
            continue

        i = h//gt_interval
        if i in range(corrs.size):
            auto_corrs[h] = corrs[i]
        else:
            print('calc_auto_correlations: Lag value ' + str(h) + ' out of range for data series (1-'+str(gt_len*gt_interval)+').')
            continue

    return auto_corrs


# Calculate the Mean Absolute Scaled Error (MASE) for a model truth with respect to a
# ground truth ('gt') and the supplied autocorrelation lag values ('h').
#
# Assume 'mt' and 'gt' are numpy arrays containing property values separated by a fixed interval,
# which can be temporal or spatial.
#
# Returns a dictionary of MASE values derived from comparing 'mt' with 'gt'; the lag value is the key.
# The 'lag_values' parameter is used to decide which error values are stored in the returned dictionary.
def calc_mase_errors(mt, gt, interval, mase_scalings, lag_values):

    errors = {}

    mt_len = len(mt)
    gt_len = len(gt)

    if mt_len == gt_len:
        # assume model truth and ground truth series have same interval
        for h in lag_values:

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
def calc_skill_scores(mase_errors, auto_corrs, map_func, lag_values):

    scores = {}

    for h in lag_values:

        scores[h] = (1 - map_func(mase_errors[h])) * (1 - abs(auto_corrs[h]))

    return scores


# Calculate the near-surface relative humidity ground truth from the
# 2m temperature and 2m dewpoint temperature extracted from the ERA5 archive
#
# 'hurs' is the CMIP variable name for this form of humidity.
def calc_hurs_gt(era5_values_dict):
    t2m = era5_values_dict['t2m']  # 2m temperatue
    d2m = era5_values_dict['d2m']  # 2m dew point temperatue

    if len(t2m) == len(d2m):

        # define Magnus coefficents
        mc_beta = 17.625
        mc_lambda = 243.04

        hurs = []
        for i in range(len(t2m)):
            t = t2m[i]
            dpt = d2m[i]

            hurs.append(100.0*(math.exp((mc_beta*dpt)/(mc_lambda+dpt))
                               /
                               math.exp((mc_beta*t)/(mc_lambda+t))))
            
        return np.float32(hurs)
    
    else:
        print('calc_hurs_gt: Cannot calculate "hurs" ground truth: mismatch between ERA5 temperature arrays.')
        return []


# Calculate the 10m wind speed magnitude using the 10m eastward (U) and
# 10m northward (V) wind speed components extracted from the ERA5 archive.
#
# 'sfcWind' is the CMIP variable name for this form of wind speed.
def calc_sfcwind_gt(era5_values_dict):
    u10 = era5_values_dict['10u']  # 10 metre eastward (U) wind component
    v10 = era5_values_dict['10v']  # 10 metre northward (V) wind component

    if len(u10) == len(v10):

        sfcwind = []
        for i in range(len(u10)):
            u = u10[i]
            v = v10[i]

            sfcwind.append(np.sqrt(u**2 + v**2))
            
        return np.float32(sfcwind)
    
    else:
        print('calc_sfcwind_gt: Cannot calculate "sfcWind" ground truth: mismatch between ERA5 wind speed arrays.')
        return []
    

# The 'list_fn' file contains lines that follow a "<model output file>,<time index selection>" format.
# It is assumed that the files are listed in ascending date order.
# The model output files are NetCDF files.
# The name of each file contains the temporal resolution and CMIP variable.
def parse_model_output_list(list_fn, temporal_res, var_list):

    model_output_fns = {}
    model_output_tis = {}

    print('\nParsing model output list: ' + str(list_fn) + '...\n')
    output_lines = os.popen('cat ' + list_fn).read().split('\n')
    
    for out_line in output_lines:
        out_line_split = out_line.split(',')
        out_fn = out_line_split[0]
        out_ti = out_line_split[-1]

        if temporal_res in out_fn:
            # data in output file has correct temporal resolution
            for var in var_list:
                if var in out_fn:
                    # data in output file concerns variable var
                    if var in model_output_fns:
                        model_output_fns[var].append(out_fn)
                        model_output_tis[var].append(out_ti)
                    else:
                        model_output_fns[var] = [out_fn]
                        model_output_tis[var] = [out_ti]
        
    for var in var_list:
        if var not in model_output_fns:
            print('parse_model_output_list: CMIP variable ' + str(var) + ' not captured at ' + str(temporal_res) + ' intervals by any of the model output files listed in ' + list_fn + '.')

    return model_output_fns, model_output_tis


'''
An example call of calc-skill-score-wf.py for the current version.

The longitudinal and latitudinal ranges describe a region that covers part of continental Europe,
that runs from the south east corner of France up to north Netherlands, across to the eastern Belarus
and down to east Romania.

./calc-skill-score-wf.py -vn 'tas' -vc '#1b9e77' -s '2000-06-01T03:00' -e '2000-06-30T22:00' -t '6h' -o 6 28 -a 46 53 -d 'CMIP6 HadGEM3' -u './model-output-v1.lst' -m 1 -lv 6 12 24 168 -ll '6 hrs' '12 hrs' '1 day' '1 wk'
./calc-skill-score-wf.py -vn 'tas' -vc '#1b9e77' -s '2000-01-01T03:00' -e '2000-12-30T22:00' -t '6h' -o 6 28 -a 46 53 -d 'CMIP6 HadGEM3' -u './model-output-v2.lst' -m 1 -lv 6 12 24 168 720 -ll '6 hrs' '12 hrs' '1 day' '1 wk' '30 days'
./calc-skill-score-wf.py -vn 'tas' -vc '#1b9e77' -s '2000-01-01T03:00' -e '2010-12-30T22:00' -t '6h' -o 6 28 -a 46 53 -d 'CMIP6 HadGEM3' -u './model-output-v3.lst' -m 1 -lv 6 12 24 168 720 8766 -ll '6 hrs' '12 hrs' '1 day' '1 wk' '30 days' '1 yr'

./calc-skill-score-wf.py -vn 'tas' 'psl' 'hurs' -vc '#1b9e77' '#d95f02' '#7570b3' -s '2000-06-01T03:00' -e '2000-06-30T22:00' -t '6h' -o 6 28 -a 46 53 -d 'CMIP6 HadGEM3' -u './model-output.lst' -m 1 -lv 6 12 24 168 -ll '6 hrs' '12 hrs' '1 day' '1 wk'

./calc-skill-score-wf.py --variable-names 'tas' 'psl' 'hurs' \
                         --variable-colours '#1b9e77' '#d95f02' '#7570b3' \
                         --start-date '2000-06-01T03:00' \
                         --end-date '2000-06-30T21:00' \
                         --time-step '6h' \
                         --longitudes 6 28 \
                         --latitudes 46 53 \
                         --model-name 'CMIP6 HadGEM3' \
                         --model-output-list './model-output-v1.lst'
                         --map-func-type 1 \
                         --lag-values 6 12 24 168 \
                         --lag-labels '6 hrs' '12 hrs' '1 day' '1 wk' 
'''

cmip_property_dict = { 
        'tas': { 'desc': 'near-surface air temperature',   'units': 'K',   'era5': ['2t'] },
        'psl': { 'desc': 'sea level pressure',             'units': 'Pa',  'era5': ['msl'] },
       'hurs': { 'desc': 'near-surface relative humidity', 'units': '%',   'era5': ['2t', '2d'] },
    'sfcWind': { 'desc': 'near-surface wind speed',        'units': 'm/s', 'era5': ['10u', '10v'], 'cmip': ['uas', 'vas'] }
}



parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
parser.add_argument('-vn', '--variable-names', nargs='+', type=str, default='tas', help='CMIP variable name list: format "<v1> <v2> ... <vn>"')
parser.add_argument('-vc', '--variable-colours', nargs='+', type=str, default='red', help='CMIP variable colour list: format "<c1> <c2> ... <cn>"')
parser.add_argument('-s', '--start-date', type=str, default='2020-06-01', help='The start date: format is "YYYY-MM-DDThh:mm:ss"')
parser.add_argument('-e', '--end-date', type=str, default='2020-06-02', help='The end date: format is "YYYY-MM-DDThh:mm:ss"')
parser.add_argument('-t', '--time-step', type=str, default='1h', help='Time step')
parser.add_argument('-o', '--longitudes', nargs='+', type=int, default='0 360', help='The longitudinal range')
parser.add_argument('-a', '--latitudes', nargs='+', type=int, default='-90 90', help='The latitudinal range')
parser.add_argument('-d', '--model-name', type=str, default='Unknown Model', help='The name of the weather forecast model to be evaluated')
parser.add_argument('-u', '--model-output-list', type=str, default='', help='The name of a file that contains a list of model output files and time index selections.' )
parser.add_argument('-m', '--map-func-type', type=int, default=1, help='The type of mapping function used to ensure that the skill score is in the range 0-1')
parser.add_argument('-lv', '--lag-values', nargs='+', type=int, default='1', help='The list of autocorrelation lag values in hours')
parser.add_argument('-ll', '--lag-labels', nargs='+', type=str, default='1 hr', help='The list of autocorrelation lag labels')


# Parse command line arguments
args = parser.parse_args()

variable_names = args.variable_names
variable_colours = args.variable_colours

time_interval = int(args.time_step[:-1])

long_min = args.longitudes[0]
long_max = args.longitudes[1]

lat_min = args.latitudes[0]
lat_max = args.latitudes[1]

lag_values = args.lag_values
lag_labels = args.lag_labels


# Retreive the model output filenames and corresponding time indexes 
model_output_fns, model_output_tis = parse_model_output_list(args.model_output_list, args.time_step, variable_names)

# Retrieve mapping function used to convert skill score to 0-1 range
map_func = MappingFunctions().get_function(args.map_func_type)
    

# Prepare CDO tmp directory
tempPath = '~/tmp/'
cdo = Cdo(tempdir=tempPath)
cdo.cleanTempDir()

# The skill scores dictionary holds a list of skill scores per variable, where each score
# corresponds to an autocorrelation lag value
skill_scores = {}

# Locate the ERA5 archive and assign the result to the ground truth archive
gt_arc = Find_era5()


# Iterate over the CMIP variables
for cmip_var in variable_names:

    cmip_prop = cmip_property_dict[cmip_var]
    print('Calculating skill scores for ' + cmip_prop['desc'] + '...')

    print('Obtaining ground truth...')
    # Extract from the ERA5 archive the required ground truth values for the specified time frame and spatial region
    gt_ds = gt_arc[cmip_prop['era5'],
                   args.start_date:args.end_date:args.time_step,
                   None, # ignore pressure levels
                   long_min:long_max,
                   lat_min:lat_max]
    
    print('Averaging ground truth...')
    gt_values_dict = {}
    for gt_var, gt_values in gt_ds.items():
        print('gt_var='+str(gt_var))
        # Average the ground truth over the specified spatial region
        gt_avg = gt_values.mean(dim=['longitude', 'latitude'])
        gt_values_dict[gt_var] = np.float32(gt_avg.values)

    # Obtain the ground truth corresponding to the CMIP variable
    # In some cases, the ground truth may need to be calculated from the extracted ERA5 data
    match cmip_var:
        case 'tas':
            gt_values = gt_values_dict['t2m']
        case 'psl':
            gt_values = gt_values_dict['msl']
        case 'hurs':
            gt_values = calc_hurs_gt(gt_values_dict)
        case 'sfcWind':
            gt_values = calc_sfcwind_gt(gt_values_dict)
        case _:
            print('Error, CMIP variable ' + str(cmip_var) + ' cannot be derived from available ERA5 ground truth data.')
            gt_values = []

    if 0 == len(gt_values):
        continue
            
    mase_scalings = calc_mase_scalings(gt_values, time_interval, lag_values)
    print('MASE Scalings\n' + str(mase_scalings) + '\n\n')

    auto_corrs = calc_auto_correlations(gt_values, time_interval, lag_values)
    print('Autocorrelations\n' + str(auto_corrs) + '\n\n')


    print('Obtaining model truth...')
    if cmip_var in model_output_fns:
        mt_ds = None
        for i, out_path in enumerate(model_output_fns[cmip_var]):
            ti_sel = model_output_tis[cmip_var][i]
    
            lonlat_extent = str(long_min) + ',' + str(long_max) + ',' + str(lat_min) + ',' + str(lat_max)
            cdo_in_sel = '-selvar,' + cmip_var + ' -seltimestep,' + ti_sel + ' ' + out_path

            mt_ds_i = cdo.sellonlatbox(lonlat_extent, input=cdo_in_sel, returnXArray=cmip_var)
            if mt_ds == None:
                mt_ds = mt_ds_i
            else:
                mt_ds = xr.combine_by_coords([mt_ds_i, mt_ds])
    else:
        continue

    print('Averaging ground truth...')    
    mt_avg = mt_ds.mean(dim=['lat','lon'])

    # assume model truth is in 32-bit float precision
    mt_values = mt_avg.values

    mase_errors = calc_mase_errors(mt_values, gt_values, time_interval, mase_scalings, lag_values)
    print('MASE Errors\n' + str(mase_errors) + '\n\n')

    skill_scores[cmip_var] = calc_skill_scores(mase_errors, auto_corrs, map_func, lag_values)
    print('Skill Scores\n' + str(skill_scores[cmip_var]) + '\n\n')


# Plot skill scores as a grouped bar chart
bar_width = 0.25

br = np.arange(len(lag_values))

fig, ax = plt.subplots(dpi=1000)
for i, cmip_var in enumerate(skill_scores):    
    ax.bar(br+(i*bar_width), skill_scores[cmip_var].values(), color=variable_colours[i], width=bar_width, edgecolor='white', label=cmip_property_dict[cmip_var]['desc'])

ax.set_ylabel('Skill Score', fontweight='bold')
ax.set_xlabel('Autocorrelation Lag', fontweight='bold')
ax.set_xticks(br + bar_width)
ax.set_xticklabels(lag_labels)

temporal_range = args.start_date.split('T')[0] + ':' + args.end_date.split('T')[0] + ':' + args.time_step
spatial_range = 'long: ' + str(long_min) + '-' + str(long_max) + ', lat: ' + str(lat_min) + '-' + str(lat_max)
ax.set_title(args.model_name + ' Skill Scores\n(' + temporal_range + ', ' + spatial_range + ')')

ax.legend()

fig_name = args.model_name
fig_name = fig_name.replace(' ', '-')
fig_name = fig_name.lower() + '-skill-scores.jpg'
fig.savefig(fig_name, format='jpg', dpi=1000)
