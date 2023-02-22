# This file contains methods used for
# different calculations.
import math

import numpy as np


def calculate_mean(value_list: list, decimal_places=2):
    """
    This method computes the mean value out of a list
    of values. The result is rounded to two decimal places.
    :param value_list: list of values.
    :param decimal_places: rounded to these decimal places.
    :return: mean of the given list of values.
    """

    # Create variable to be returned
    mean = 0

    # Iterate over value list
    for length in value_list:
        mean += length

    mean = round(mean / len(value_list), decimal_places)

    return mean


def calculate_median(value_list: list):
    """
    This method computes the median value out of a list
    of values.
    :param value_list: list of values.
    :return: median of the given list of values.
    """

    # Create sorted copy of list
    sorted_list = sorted(value_list)

    # Check if list length is even or odd
    if is_even(sorted_list) is False:
        # return middle value
        middle_index = int(math.ceil(len(sorted_list)/2))
        return sorted_list[middle_index-1]
    else:
        # Get the two values above/below middle index
        middle_index = int(len(sorted_list)/2)

        below = sorted_list[middle_index-1]
        above = sorted_list[middle_index]

        # return median value
        return (above+below)/2


def is_even(value_list: list):
    """
    This method checks if the length of a list is even.
    :param value_list: list of values.
    :return: True, if the length of the list is even, odd otherwise.
    """
    if len(value_list) % 2 == 0:
        return True
    return False


def calculate_fragment_length(start: int, stop: int):
    """
    This method simply calculates the fragment length of the
    two specified start/stop base pairs.
    :param start: start position of the fragment (bp).
    :param stop: stop position of the fragment (bp).
    :return: length of the fragment.
    """
    return abs(start-stop)


def normalize(array):
    """
    This method normalizes an array with the formula (a-min(a)/(max(a)-min(a)).
    :param array: array to normalize
    :return: normalized array.
    """
    array = np.asarray(array)
    array_min = array.min()
    array_max = array.max()
    normalized_array = (array - array_min) / (array_max - array_min)
    return normalized_array


def calculate_xscale(df, bins=30):
    """
    This method normalizes an array with the formula (a-min(a)/(max(a)-min(a)).
    :param df: dataframe with 'Fragments' column.
    :return: x_scale with bin points over the range of fragment lengths.
    """
    # calculate bin scale over the range of fragment lengths
    max_value = max(df['Fragments'].apply(max))
    min_value = min(df['Fragments'].apply(min))
    value_range = max_value - min_value
    bin_scale = value_range/bins
    
    # return a numpy array with the beginning value of each bin
    x_scale = [x for x in np.arange(min_value,max_value,bin_scale)]
    x_scale = np.asarray(x_scale)
    return x_scale


def calculate_maxima(value_list):
    """
    This method computes all local maxima in a given list of numerical
    values using a "sliding window" to filter false-positive local maxima
    resulting from noisy data. This window has a fixed size of 5.
    A list of all indices of such local maxima found in the provided
    value list is then returned.
    :param value_list: list of numerical values
    :return: list of indices of all local maxima
    """

    # Create local maxima list
    local_maxima = []

    # Store length of distribution parameter
    distr_len = len(value_list)

    # Iterate over distribution
    for index, count in enumerate(value_list):

        # Create "empty" sliding window
        window = np.array([None, None, None, None, None])

        # Set window
        if index - 2 >= 0:
            window[0] = value_list[index - 2]
        if index - 1 >= 0:
            window[1] = value_list[index - 1]
        if index >= 0:
            window[2] = value_list[index]
        if index + 1 <= distr_len - 1:
            window[3] = value_list[index + 1]
        if index + 2 <= distr_len - 1:
            window[4] = value_list[index + 2]

        # Check if value in window is local maxima
        poss_maxima = window[2]
        is_maxima = False
        # If window[1] is None, then window[0] is also None
        if window[1] is None:
            if poss_maxima > window[3] and poss_maxima > window[4]:
                is_maxima = True
        # Check if only window[0] is None
        elif window[0] is None:
            if window[1] < poss_maxima and poss_maxima > window[3] and poss_maxima > window[4]:
                is_maxima = True
        # Check if window[4] is not None, then window[3] must also be not None
        # This is the usual case
        elif window[4] is not None:
            if window[0] < poss_maxima and window[1] < poss_maxima and poss_maxima > window[3] and poss_maxima > window[4]:
                is_maxima = True

        # Check if possible maxima is a real local maxima
        # If yes, then append it to local_maxima list
        if is_maxima:
            local_maxima.append(index)

    # Return list of all local maxima
    return local_maxima


def calculate_score(peaks, min_frag, bin_size, bins=30, period=160):
    """
    Computes the average difference between the provided peak indices and
    the specified period. The calculation is divided into three cases:

    1. The peak list is empty:
        -> Return (positive) infinity
    2. The peak list contains just one peak:
        -> Return absolute difference between peak and period+min_frag
    3. The peak list contains more than one peak:
        -> Return average difference between all peaks w.r.t. specified period

    Hence, the closer this value conforms to 0 the better, making 0 the best possible
    value.
    :param peaks: list of peaks indices
    :param period: period in which the peaks should occur
    :param min_frag: length of the smallest fragment
    :param bin_size: size of each bin
    :param bins: number of bins the peak indices will be mapped to
    """

    # Create list for later computation
    temp = range(0, bins + 1)

    # Create empty list for computed locations
    peak_bins = []


    for i in temp:
        if temp.index(i) in peaks:
            peak_bins.append(i * bin_size + min_frag)

    # Compare Difference
    diff = 0

    # Check for 0 peaks
    if len(peaks) == 0:
        return np.inf

    # Check for 1 peak
    if len(peaks) == 1:
        peak_diff = abs(peak_bins[0] - (period + min_frag))
        return round(peak_diff, 2)

    # More than 1 peak
    for index in range(len(peak_bins) - 1):
        bin_diff = abs(peak_bins[index] - peak_bins[index + 1])
        diff += abs(bin_diff - period)

    diff /= len(peaks) - 1

    return round(diff, 2)
