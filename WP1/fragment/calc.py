# This file contains methods used for
# different calculations.
import math



def compute_mean(value_list: list, decimal_places=2):
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


def compute_median(value_list: list):
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


# This method checks if the length of a list is even.
def is_even(value_list: list):
    """
    This method checks if the length of a list is even.
    :param value_list: list of values.
    :return: True, if the length of the list is even, odd otherwise.
    """
    if len(value_list) % 2 == 0:
        return True
    return False


# This method simply calculates the fragment length of the
# two specified start/stop base pairs.
def calculate_fragment_length(start: int, stop: int):
    """
    This method simply calculates the fragment length of the
    two specified start/stop base pairs.
    :param start: start position of the fragment (bp).
    :param stop: stop position of the fragment (bp).
    :return: length of the fragment.
    """
    return abs(start-stop)
