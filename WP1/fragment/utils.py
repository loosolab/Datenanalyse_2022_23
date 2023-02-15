# This file contains wrapper and utility functions for
# working with fragment files.

import logging
import os
import pandas as pd
import anndata as ad
import scanpy as sc
import fragment.calc as calc


def read_fragment_file(abs_path: str):
    """
    This method reads a fragment file (.bed) and returns
    a dictionary with the cellbarcodes as keys and the
    computed fragment lengths as the corresponding values.
    The fragment lengths are stored as a list.
    Returns None when file extension is not .bed.
    :param abs_path: (absolute) path to the .bed file.
    :return: dictionary (keys=cellbarcode/values=list of fragment lengths)
    """

    # Check if the file extension is .bed
    if abs_path[-4:] != ".bed":
        logging.warning("provided file does not have '.bed' extension:\t"+abs_path)
        return None

    # Get file from disk
    fragment_file = open(abs_path, "r")

    # Create dictionary
    frag_dictionary = {}

    # Read file
    for line in fragment_file:

        # Split current line
        line_values = line.split()

        # Get cell barcode
        cellbarcode = line_values[3]

        # Get both fragments
        start = int(line_values[1])
        stop = int(line_values[2])

        # Check if cellbarcode key is in dictionary
        if cellbarcode in frag_dictionary:
            # Append list of fragment lengths
            frag_dictionary[cellbarcode].append(calc.calculate_fragment_length(start, stop))
        else:
            # Create new Key/Value (Cellbarcode/Fragment length list)
            frag_dictionary[cellbarcode] = [calc.calculate_fragment_length(start, stop)]

    # Close open access to file
    fragment_file.close()

    # Return dictionary
    return frag_dictionary


def read_fragment_directory(abs_path: str):
    """
    This method reads all .bed files from the
    specified path (directory).
    :param abs_path: path to the directory.
    :return: list of dictionaries. One for each
             .bed file found.
    """

    # Create list of dictionaries to be returned
    dictionary_list = []

    # Iterate over files in directory
    for file in os.listdir(abs_path):
        frag_dict = read_fragment_file(abs_path + "/" + file)
        if frag_dict is not None:
            dictionary_list.append(frag_dict)

    return dictionary_list


def combine_fragment_dictionaries(dictionary_list: list):
    """
    This method takes a list of dictionaries (e.g. through usage of
    read_fragment_directory) and combines them into a single dictionary.
    The returned dictionary contains all keys out of the dictionary list and
    their corresponding fragment length lists, merged into a single list.
    Note: Duplicate fragment length values are still added to the
          fragment length list of a cellbarcode.
    :param dictionary_list: list of dictionaries.
    :return: single dictionary with merged keys/values.
    """

    # Create dictionary to be returned
    frag_dictionary = {}

    # Iterate over dictionary list
    for dictionary in dictionary_list:

        # Iterate through dictionary (cellbarcode = key)
        for cellbarcode in dictionary:

            # Check if cellbarcode is not in dictionary
            if cellbarcode not in frag_dictionary:
                # Create new entry
                frag_dictionary[cellbarcode] = []

            # Extend values of the cellbarcode
            frag_dictionary[cellbarcode].extend(dictionary[cellbarcode])

    return frag_dictionary


def create_dataframe(frag_dictionary: dict, tissue=""):
    """
    This method creates a new dataframe object with the values
    taken from the input dictionary. The resulting object only
    contains the cellbarcode from the dictionary as an index and
    the corresponding mean/median values as well as the fragment lists
    from the provided frag_dictionary.
    :param frag_dictionary: dictionary from which the object is created from.
    :param tissue: string of the belonging tissue.
    :return: dataframe object with cellbarcodes, means, medians and fragment length lists.
    """

    # Create empty dictionary
    mean_median_dictionary = {}

    # Iterate over keys in frag_dictionary
    if tissue == "":
        for cellbarcode in frag_dictionary:
            mean_median_dictionary[cellbarcode] = {}
            mean_median_dictionary[cellbarcode]["Mean"] = calc.compute_mean(frag_dictionary[cellbarcode])
            mean_median_dictionary[cellbarcode]["Median"] = calc.compute_median(frag_dictionary[cellbarcode])
            mean_median_dictionary[cellbarcode]["Fragments"] = frag_dictionary[cellbarcode]
    else:
        for cellbarcode in frag_dictionary:
            mean_median_dictionary[tissue+"+"+cellbarcode] = {}
            mean_median_dictionary[tissue+"+"+cellbarcode]["Mean"] = calc.compute_mean(frag_dictionary[cellbarcode])
            mean_median_dictionary[tissue+"+"+cellbarcode]["Median"] = calc.compute_median(frag_dictionary[cellbarcode])
            mean_median_dictionary[cellbarcode]["Fragments"] = frag_dictionary[cellbarcode]

    # Transform dictionary into dataframe
    data_frame = pd.DataFrame(mean_median_dictionary).T

    # Return dataframe
    return data_frame


def add_mean_to_h5ad(ann_data: ad.AnnData, fragment_dictionary: dict, tissue: str):
    """
    This method adds the mean value of a given list of fragment lengths
    to an existing annData object. Only those rows (obs) with matching
    tissue+cellbarcode (from fragment_dictionary) as index are extended.
    :param ann_data: annData object to be extended.
    :param fragment_dictionary: dictionary containing cellbarcodes/list of fragment lengths.
    :param tissue: string of the belonging tissue.
    """
    # Iterate over fragment_dictionary
    for cellbarcode in fragment_dictionary:
        if tissue+"+"+cellbarcode in ann_data.obs_names:
            ann_data.to_df().at[tissue+"+"+cellbarcode, "Mean"] = calc.compute_mean(fragment_dictionary[cellbarcode])


def add_median_to_h5ad(ann_data: ad.AnnData, fragment_dictionary: dict, tissue: str):
    """
    This method adds the median value of a given list of fragment lengths
    to an existing annData object. Only those rows (obs) with matching
    tissue+cellbarcode (from fragment_dictionary) as index are extended.
    :param ann_data: annData object to be extended.
    :param fragment_dictionary: dictionary containing cellbarcodes/list of fragment lengths.
    :param tissue: string of the belonging tissue.
    """

    # Iterate over fragment_dictionary
    for cellbarcode in fragment_dictionary:
        if tissue+"+"+cellbarcode in ann_data.obs_names:
            ann_data.to_df().at[tissue+"+"+cellbarcode, "Mean"] = calc.compute_median(fragment_dictionary[cellbarcode])
