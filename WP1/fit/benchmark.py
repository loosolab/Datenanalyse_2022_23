import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit
from scipy.stats import entropy

from fragment import utils
from fragment import calc
from fit import calc
from fit import benchmark as bc
from fit import function as func
from fit import score


def load_data(path: str):
    print('loading fragments...')
    fragments = utils.read_fragment_file(path) 
    df = utils.create_dataframe(fragments)
    
    df['Fragment-Count'] = [len(x) for x in df['Fragments']]
    
    print('calculate distribution...')
    df['Distribution'] = get_distribution(df)
        
    print('calculate maxima...')
    df['Maxima'] = get_maxima(df)
    df['Maxima-Count'] = [len(x) for x in df['Maxima']]
     
    print('calculate score...')
    df['Score'] = score.compute_score_df(df)
    return df

def get_distribution(df, bins=30, min_value = 50, max_value = 700):
    min_value = min(df['Fragments'].apply(min))
    max_value = max(df['Fragments'].apply(max))
    value_range = max_value - min_value
    bin_scale = int(value_range/bins)
    bin_index = [x for x in range(min_value,max_value,bin_scale)]
    distribution = []
    for i in df['Fragments']:
        inds = np.digitize(i, bin_index)
        y_values = []
        for x in range(len(bin_index)):
            y_values.append(np.count_nonzero(inds == x+1))
        distribution.append(y_values)       
    return distribution

def get_maxima(df, distribution = 'Distribution'):
    maxima = []
    for i in df[distribution]:
        if i is np.nan:
            maxima.append(np.nan)
        else:
            maxima.append(score.get_real_maxima(i))
    return maxima



##### STILL IN DEVELOPMENT #####

def get_fit(df, function = func.fitBaseline, bins=30, log_plot = False, max_val_num = 0.1):

    fit_error = 0
    x_scale = calc.get_x_scale(df)
    fits = []
    for distribution in df['Distribution']:    
             
        # normalize data
        normalized_distribution = calc.normalize(distribution)
        # calculate best fit_parameters
        if function == func.fitBaseline:
            try:
                parameters, covariance = curve_fit(function,x_scale,normalized_distribution, bounds = (0,[2.,2.,100.]),maxfev=200000)
            
                fit_A = parameters[0]
                fit_B = parameters[1]    
                fit_C = parameters[2]
                fit_y = function(x_scale, fit_A, fit_B, fit_C)
                fits.append(fit_y)
            except:
                fit_error += 1
                fits.append(np.nan)
            continue
        else:
            try:
                parameters, covariance = curve_fit(function,x_scale,normalized_distribution, bounds = (0,[2.,2.]),maxfev=200000)
                fit_A = parameters[0]
                fit_B = parameters[1] 
                fit_y = function(x_scale, fit_A, fit_B)
                fits.append(fit_y)
            except:
                fit_error += 1
                fits.append(np.nan)
                continue
    print(f'Count of not fittet cells: {fit_error}')
    return fits

def get_kl_divergence(df):
    kl_divergence = []
    data = zip(np.asarray(df['Distribution']), np.asarray(df['Fit']))
    for d, f in data:
        if f is np.nan:
            kl_divergence.append(np.nan)
        else:
            kl_divergence.append(entropy(d,f))
    return kl_divergence


def get_maxima_position(df):
    scores = []
    score_fails = 0
    baseline = [x for x in df['Basemaxima']]
    fit = [x for x in df['Maxima']]
    for b, f in zip(baseline, fit):
        try:
            scores.append(calc.compute_score(b, f))
        except:
            score_fails += 1
            scores.append(np.nan)
    print(f'failed scores: {score_fails}')
    return scores
