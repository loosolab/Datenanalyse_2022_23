import logging
import os
import itertools

import numpy as np
import pandas as pd
import anndata as ad
import HTSeq
import src.calc as calc


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
            mean_median_dictionary[cellbarcode]["Mean"] = calc.calculate_mean(frag_dictionary[cellbarcode])
            mean_median_dictionary[cellbarcode]["Median"] = calc.calculate_median(frag_dictionary[cellbarcode])
            mean_median_dictionary[cellbarcode]["Fragments"] = frag_dictionary[cellbarcode]
    else:
        for cellbarcode in frag_dictionary:
            mean_median_dictionary[tissue+"+"+cellbarcode] = {}
            mean_median_dictionary[tissue+"+"+cellbarcode]["Mean"] = calc.calculate_mean(frag_dictionary[cellbarcode])
            mean_median_dictionary[tissue+"+"+cellbarcode]["Median"] = calc.calculate_median(frag_dictionary[cellbarcode])
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
            ann_data.to_df().at[tissue+"+"+cellbarcode, "Mean"] = calc.calculate_mean(fragment_dictionary[cellbarcode])


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
            ann_data.to_df().at[tissue+"+"+cellbarcode, "Mean"] = calc.calculate_median(fragment_dictionary[cellbarcode])


def load_data(path: str, bins = 30, penalty = 200):
    """
    This method creates an dataframe with cell barcode as index and colums
    for fragment lengths ('Fragments'), fragment count ('Fragment-Count'),
    fragment length distribution ('Distribution'), maxima position ('Maxima'),
    maxima count ('Maxima-Count') und einen quality score ('Score').
    :param path: path to a .bed file.
    :param bins: resolution of the calculated distribution further calculations based on.
    :return: dataframe with the colums specified above.
    """
    
    # loding fragment_file and creating dataframe with fragments, fragment_count, mean and median
    print('loading fragments...')
    fragments = read_fragment_file(path)
    df = create_dataframe(fragments)
    df['Fragment-Count'] = [len(x) for x in df['Fragments']]

    # calculate distribution in predifined bins and add it to the dataframe
    print('calculate distribution...')
    df['Distribution'] = get_distribution(df, bins = bins)
    
    # calculate maxima and their count and add them to the dataframe
    print('calculate maxima...')
    df['Maxima'] = get_maxima(df)
    df['Maxima-Count'] = [len(x) for x in df['Maxima']]
    
    # calculate score and and add it to the dataframe
    print('calculate score...')
    get_score(df, bins = bins, penalty = penalty)
    return df


def get_distribution(df, bins=30):
    """
    This method scales a given dataframe with a column of fragment lengths and binnes them
    to get the y-values of the distribution.
    :param df: dataframe with a column for fragment lengths.
    :param bins: resolution of the calculated distribution further calculations based on.
    """
    
    # calculate list of bin_indexes over the range of fragment lengths
    min_value = min(df['Fragments'].apply(min))
    max_value = max(df['Fragments'].apply(max))
    value_range = max_value - min_value
    bin_scale = value_range / bins
    bin_index = [x for x in np.arange(min_value, max_value, bin_scale)]
    
    # digitize fragment length vor every dataframe index into predefined bins
    distribution = []
    for i in df['Fragments']:
        inds = np.digitize(i, bin_index)
        
        # count the elements in every bin to get the y-value every point in the distribution
        y_values = []
        for x in range(len(bin_index)):
            y_values.append(np.count_nonzero(inds == x + 1))
        distribution.append(y_values)
    return distribution


def get_maxima(df, distribution='Distribution'):
    """
    This method calculates the local maxima for every index in a dataframe.
    :param df: dataframe with a kind of distribution column
    :param distribution: column in a dataframe with arrays of values to calculate the maximas from
    """
    
    # itterate over whole dataframe
    maxima = []
    for i in df[distribution]:
        
        # filter empty cells in the dataframe
        if i is np.nan:
            maxima.append(np.nan)
        else:
            # calculate local maxima
            maxima.append(calc.calculate_maxima(i))
            
    return maxima


def get_score(df, bins = 30, penalty = 200):
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
    if "Distribution" not in df.columns:
        print("Dataframe does not have a column Distribution!")
        return

    # Check if column "Score" exists, if not add it
    if "Score" not in df.columns:
        df["Score"] = np.NaN
        
    # Calculate bin_size
    min_frag = min(df['Fragments'].apply(min))
    max_frag = max(df['Fragments'].apply(max))
    bin_size = (max_frag - min_frag) / bins

    # Create empty score list
    score_list = []

    # Iterate over dataframe rows
    for index, frame in df.iterrows():

        # Compute Score for this row
        score = calc.calculate_score(calc.calculate_maxima(df["Distribution"][index]), min_frag, bin_size, bins = bins)
        
        # add penalty for fragment count
        if df['Fragment-Count'][index] < penalty:
            score = score + (penalty - df['Fragment-Count'][index])
        score = score*(1/np.log(df['Fragment-Count'][index]))

        # Append score list with computed value
        score_list.append(score)

    # Set the values of the column "Score" in the dataframe
    df["Score"] = score_list
    return score_list


def get_dataframe_from_BED_file(bed_file):
    
    #read bed file
    bedfile = open(bed_file, "r")
    
    #get data from file
    frag_list = []
    for line in bedfile:
        chrom, start, end, cell, value, etc = line.strip().split()
        frag_list.append([chrom, start, end, cell, value, etc])
        
    #create a dataframe
    frag_df = pd.DataFrame(data = frag_list, columns = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand'])

    #Typecasting
    frag_df.Start = pd.to_numeric(frag_df.Start, downcast='integer')
    frag_df.End = pd.to_numeric(frag_df.End, downcast='integer')
    frag_df.Score = pd.to_numeric(frag_df.Score, downcast='float')
    
    return frag_df


def get_genomic_array_from_BED_dataframe(data_frame):
    
    #Create Genomic Array (coverage)
    coverage = HTSeq.GenomicArray("auto", stranded=False, typecode="i")
    
    #add entries from dataframe to genomic array
    for row in data_frame.itertuples():
        coverage[HTSeq.GenomicInterval(row.Chromosome, int(row.Start), int(row.End), row.Strand)] += 1
        
    return coverage


def get_TSS_set_from_GTF_file(gtf_file):
    
    #read gtf file
    gtffile = HTSeq.GFF_Reader(gtf_file)
    
    #Add TSS data in a set.
    tsspos = set()
    for feature in gtffile:
        if feature.type == "exon" and feature.attr["exon_number"] == "1":
            tsspos.add(feature.iv.start_d_as_pos)
    
    return tsspos


def get_profile_for_plots(genomic_array, tss_positions, half_window_width = 500):
    
    #profile is an array that stores y-values for the plot
    #initialize profile
    profile = np.zeros(2*half_window_width, dtype='i')
    
    #update profile at all TSS positions.
    for p in tss_positions:
        start = p.pos - half_window_width
        if(start<0): continue
        window = HTSeq.GenomicInterval(p.chrom, p.pos - half_window_width, p.pos + half_window_width, ".")
        wincvg = np.fromiter(genomic_array[window], dtype='i', count=2*half_window_width)
        if p.strand == "+":
            profile += wincvg
        else:
            profile += wincvg[::-1]
    
    return profile
