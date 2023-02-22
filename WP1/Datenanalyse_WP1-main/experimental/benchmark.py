import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import entropy

import src.plot
import src.calc as calc
import experimental.function as func


def analyse_data(df, list_of_column, output = None):
    # plots a histogram for each and a comparisson between each variables, defined by a list, in a dataframe
    for i in list_of_column:
        src.plot.histplt(df, output, i)
    for i in list_of_column:
        for z in list_of_column:
            src.plot.compplt(df, output, i, z)
    return


def get_fit(df, function=func.fitBaseline, bins=30, log_plot=False, max_val_num=0.1):
    fit_error = 0
    x_scale = calc.calculate_xscale(df)
    fits = []
    for distribution in df['Distribution']:

        # normalize data
        normalized_distribution = calc.normalize(distribution)
        # calculate best fit_parameters
        if function == func.fitBaseline:
            try:
                parameters, covariance = curve_fit(function, x_scale, normalized_distribution,
                                                   bounds=(0, [2., 2., 100.]), maxfev=200000)

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
                parameters, covariance = curve_fit(function, x_scale, normalized_distribution, bounds=(0, [2., 2.]),
                                                   maxfev=200000)
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
            kl_divergence.append(entropy(d, f))
    return kl_divergence
