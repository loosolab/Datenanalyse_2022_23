import math

import numpy as np
from scipy.stats import entropy
from scipy.integrate import simpson
import similaritymeasures as similaritymeasures
from fit.calc import findLocalMaximaMinima, compute_local_maxima_diffMean, get_x_scale, get_sublist, closest


# TODO: In load_data einbinden !
def compute_score_df(dataframe):
    """
    Computes a score for each row of the dataframe and adds
    a corresponding column "Score" containing each individual score.
    The dataframe needs to have a column "Distribution" with
    lists of numerical values to enable a reliable scoring.
    If a score of a row could not be calculated, the value in
    the new cell will be NaN.
    :param dataframe: the dataframe to be extended
    """

    # Check if column "Distribution" exists
    if "Distribution" not in dataframe.columns:
        print("Dataframe does not have a column Distribution!")
        return

    # Check if column "Score" exists, if not add it
    if "Score" not in dataframe.columns:
        dataframe["Score"] = np.NaN

    # Create empty score list
    score_list = []

    # Iterate over dataframe rows
    for index, frame in dataframe.iterrows():

        # Compute Score for this row
        score = peak_avg_loc_diff(get_real_maxima(dataframe["Distribution"][index]))

        # Append score list with computed value
        score_list.append(score)

    # Set the values of the column "Score" in the dataframe
    dataframe["Score"] = score_list


def get_real_maxima(value_list):
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


def peak_avg_loc_diff(peaks, period=160, min_frag=50, max_frag=700, bins=30):
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
    :param max_frag: length of the largest fragment
    :param bins: number of bins the peak indices will be mapped to
    """

    # Create list for later computation
    temp = range(0, bins + 1)

    # Create empty list for computed locations
    peak_bins = []

    # Compute bin size
    bin_size = (max_frag - min_frag) / bins
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
