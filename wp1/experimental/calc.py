import statistics


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