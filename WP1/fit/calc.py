import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.optimize import curve_fit
from scipy.stats import entropy

from fragment import utils
from fragment import calc


def calculate_maxima(arr, stdout_print=True):
    """
    This function computes all local maxima of a given function. The data of
    the function has to be provided as an array. This function returns all indices
    of the found local maxima.
    Note: If one has a function f(x) = x (or any polynomial containing x), then all values
          of the "y-dimension" have to be given to this function.
    :param arr: data of the function
    :param stdout_print: if set to false, the indices will not be printed
    :return: array containing all local maxima
    """
    # Empty lists to store points of
    # local maxima and minima
    mx = []
    mi = []
    # Checking whether the first point is
    # local maxima
    if arr[0] > arr[1]:
        mx.append(0)
    if arr[0] < arr[1]:
        mi.append(0)

    # Iterating over all points to check
    # local maxima
    for i in range(1, len(arr) - 1):

        # Condition for local maxima
        if arr[i - 1] < arr[i] > arr[i + 1]:
            mx.append(i)
        if arr[i - 1] > arr[i] < arr[i + 1]:
            mi.append(i)
        
    # Checking whether the last point is
    # local maxima
    if arr[-1] > arr[-2]:
        mx.append(len(arr) - 1)
    if arr[-1] < arr[-2]:
        mi.append(len(arr) - 1)
        
#     max_val_num = 0.1
#     difference = []
#     validated_mx = []
#     if min(mx) < min(mi):
#         if arr[mx[0]]-arr[mi[0]] > max_val_num:
#             validated_mx.append(mx[0])
        
#         for i in range(1,len(mx)-1):
#             if min(arr[mx[i]]-arr[mi[i]],arr[mx[i]]-arr[mi[i-1]]) > max_val_num:
#                 validated_mx.append(mx[i])
         
#     if min(mx) > min(mi):
#         for i in range(len(mx)-1):
#             if min(arr[mx[i]]-arr[mi[i]],arr[mx[i]]-arr[mi[i+1]]) > max_val_num:
#                 validated_mx.append(mx[i])
        
#     if arr[mx[-1]]-arr[mi[-1]] > max_val_num:
#         validated_mx.append(mx[-1])
        
    mx = [x*30+50 for x in mx]

    return mx


def normalize(distribution):
    distribution = np.asarray(distribution)
    distribution_min = distribution.min()
    distribution_max = distribution.max()
    normalized_distribution = (distribution - distribution_min) / (distribution_max - distribution_min)
    return normalized_distribution

def get_kl_divergence(cells):
    return entropy(cells['Distribution']/ cells['Fit'])


def get_x_scale (cells, bins=30):
    max_value = max(cells['Fragments'].apply(max))
    min_value = min(cells['Fragments'].apply(min))
    value_range = max_value - min_value
    if bins > max_value:
        print(f'Fragment range to small for {bins} bins. Bins get adjusted to {int(max_value)}.')
        bins = max_value
    bin_scale = int(value_range/bins)
    x_scale = [x for x in range(min_value,max_value,bin_scale)]
    x_scale = np.asarray(x_scale)
    return x_scale



def compute_local_maxima_diffMean(benchmark, fit):

    if len(benchmark) == 0 or len(fit) == 0:
        return 0

    # Create differences list
    diff_list = []

    # Iterate over benchmark function
    for maxima in benchmark:
        # First check if benchmark_maxima overextends the length of fit_maxima
        index = benchmark.index(maxima)
        if index <= len(fit)-1:
            diff_list.append(abs(maxima-fit[index]))

    # Compute Score and return it
    return statistics.mean(diff_list)


def get_sublist(function, index_list):

    # Create sublist
    sublist = []

    # Iterate over index_list and append sublist
    for index in index_list:
        sublist.append(function[index])

    # Return sublist
    return sublist


def findLocalMaximaMinima(n, arr):
    # Empty lists to store points of
    # local maxima and minima
    mx = []
    mn = []

    # Checking whether the first point is
    # local maxima or minima or neither
    if (arr[0] > arr[1]):
        mx.append(0)
    elif (arr[0] < arr[1]):
        mn.append(0)

    # Iterating over all points to check
    # local maxima and local minima
    for i in range(1, n - 1):

        # Condition for local minima
        if (arr[i - 1] > arr[i] < arr[i + 1]):
            mn.append(i)

        # Condition for local maxima
        elif (arr[i - 1] < arr[i] > arr[i + 1]):
            mx.append(i)

    # Checking whether the last point is
    # local maxima or minima or neither
    if (arr[-1] > arr[-2]):
        mx.append(n - 1)
    elif (arr[-1] < arr[-2]):
        mn.append(n - 1)

    return mx


def closest(lst, K):
    return lst[min(range(len(lst)), key=lambda i: abs(lst[i] - K))]