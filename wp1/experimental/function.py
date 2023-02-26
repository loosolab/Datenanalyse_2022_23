import numpy as np


def getBenchmarkPlotValuesForRangeX(x_start=50, x_end=700, x_step=1):
    # Use this function to get the benchmark( y-values) for a given range of x_values.
    # This function returns a list x_values, y_values which can be used for a plot.
    # x_start : the start value of x. Set to 50 by default.
    # x_end : the last value of x. Set to 700 by default.
    # x_step : the difference between two x values. This is set to 1 by default.

    x_values = np.arange(x_start, x_end, x_step)
    y_values = []
    for x in x_values:
        y_values.append(benchmark(x))
    return [x_values, y_values]


def benchmark(x):
    # The original benchmark function with all the parameters preset.
    # Edit the parameters in this function and also in the 'fitBaselineBenchmark' to reset the parameters.

    first_peak = 50
    return fitBaselineBenchmark(x, first_peak)


def fitBaselineBenchmark(x, first_peak):
    # A wrapper function to preset values for peak and bottom for the'fitBaseline' function.
    # This can also be used to calculate the fit.
    # Peak set to 10**-2 by default.
    # Set bottom to 10**-5, when the 'end' in the 'fitBaseline' is set to 1400.
    # Set bottom to 2*10**-4, when the 'end' in the 'fitBaseline' is set to 700.

    peak = 10 ** -2
    bottom = 2 * 10 ** -4
    return fitBaseline(x, peak, bottom, first_peak)


def baseline_function(x, bottom, peak):
    y = (bottom * 10 ** (np.log10(peak / bottom) * (1 - (x - 50) / (700 - 50)))) * (
                1 - (0.6) ** (np.log10(10 + ((x - 50) / 150) ** 4)) * np.sin(2 * np.pi * (x - 50) / 200))
    return y


def fitBaseline(x, peak, bottom, first_peak):
    """
    IMORTANT: While using the function to  calculate a fit, set bounds for the respective parameters.

    This function tries to replicate the structure of the read density vs fragment length curve.
    The function was constructed in regards to 3 structural properties:
    1) The exponential decay of read density wrt increase in fragment length.
    2) Periodical peaks and bottoms with a peiod length of approx. 200 bp.
    3) The flattening of these peaks and bottoms wrt increase in fragment length.

    To construct this function, the following variables have been used:

    1) x : fragment length.
    2) y : read density.
    3) period : the period with which these peaks/bottoms occur.
    4) end : the end x value for the function.
    5) peak : the max y value at the first peak.
    6) bottom : the min y value observed.
    7) first_peak : the x value of the first peak. It is the point where the function starts giving the outut.

    The following variables are taken as parameters for the function:
    1) peak
    2) bottom
    3) first_peak

    """

    # The following varaibles are defined in the function itself:
    period = 160
    end = 700

    # Defining the 3 structural properties of the function:

    # 1) The exponential decay:
    decay = bottom * 10 ** (np.log10(peak / bottom) * (1 - (x - first_peak) / (end - first_peak)))

    # 2) The period function:
    # Here a sin function is used to simulate the period.
    period_function = np.sin(2 * np.pi * (x - first_peak) / period)

    # 3) The flattening or the amlitude function for the period
    amplitude = (0.6) ** (np.log10(10 + ((x - first_peak) / period) ** 4))

    # To create the complete period function the amplitude is first multpilied to the period_function,
    # and then the result is subtracted from 1.
    # Theoretically this should prevent the complete period function from reaching values below zero.
    # Then we multiply this complete period function with the decay function so that this period function
    # swings around this exonential decay function as x increases.
    # The final result:
    y = decay * (1 - amplitude * period_function)

    return y