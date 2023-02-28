import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import src.calc as calc
import src.utils as utils


def bindistplt(df, data='Fragments', column_name='Mean', bins=1, mode='equal', plot_bins=50, show=True, output_path=None):
    """
    This method bins data from a dataframe by an variable and plot every bin as an Histogram.
    :param df: dataframe with a column for data and the variable
    :param data: column of lists of values
    :param variable: column of values to sort by
    :param bins: count of bins the lists get sorted in
    :param mode: how the lists get sorted; equal: every bin same size, linear: every bin same range
    :output_path: if defined, path to save images as .png
    """

    # initialise a list containing an empty list for every bin
    grouped_data = []
    for t in range(bins):
        grouped_data.append([])
    bin_index = 0
    # sort dataframe by variable and define max and min of the variable sorted by
    sorted_df = df.loc[df[column_name] != float("inf")].sort_values(by=[column_name])
    min_variable = sorted_df[column_name][0]
    max_variable = sorted_df[column_name][-1]

    # equal mode
    if mode == 'equal':
        count = 1
        df_length = len(df.index)
        bin_size = int(df_length / bins)
        bin_start = [min_variable]
        bin_end = []

        # iterate over an sorted fragment_dictionary and fill the previously defined list with the
        # fragment lengths. Every bin contains the same number of cells.
        for i, c in enumerate(sorted_df[data]):

            # fill fragments in bins and count cells finished
            if count != bin_size + 1:
                grouped_data[bin_index] += c
                count += 1

            # count == bin_size => change bin and reset count
            else:
                count = 1
                bin_index += 1
                grouped_data[bin_index] += c
                bin_start.append(sorted_df[column_name][i])
                bin_end.append(sorted_df[column_name][i])

        bin_end.append(max_variable)

    # linear mode
    elif mode == 'linear':
        variable_range = max_variable - min_variable
        bin_size = variable_range / bins
        bin_start = [min_variable]
        bin_end = []

        # iterate over an sorted fragment_dictionary and fill the previously defined list with the
        # fragment lengths. Every bin contains the same variable range.
        for i, c in enumerate(sorted_df[data]):

            # fill fragments in bins till the variable limit is reached
            if sorted_df[column_name][i] <= min_variable + (bin_index + 1) * (bin_size):
                grouped_data[bin_index] += c

            # variable is over the limit => change bin
            else:
                bin_index += 1
                grouped_data[bin_index] += c
                bin_start.append(sorted_df[column_name][i])
                bin_end.append(sorted_df[column_name][i])
        bin_end.append(max_variable)

    # error-message, if mode not correctly defined
    else:
        print('Please only use "equal" or "linear" as variables for "mode"!')
        return

    # visualise data for every bin
    for i, s in enumerate(grouped_data):
        sns.histplot(s, bins=plot_bins, kde=True)
        plt.title(f'Distribution of {data}\n({column_name}: {bin_start[i]}-{bin_end[i]})')
        plt.xlabel(data)
        plt.ylabel('Count')

        # check output settings
        if output_path != None:
            try:
                title = f'distribution_{data}_{column_name}_{bin_start[i]}_{bin_end[i]}'
                plt.savefig(f'{output_path}{title}.png', bbox_inches='tight')
            except:
                print('Output_path not found!')
        if show:
            plt.show()
        else:
            plt.clf()
    return


def multiplt(df, column_name='Fragment-Count', distribution='Distribution', bins=1, mode='base', lower_limit=None,
             upper_limit=None, output_path=None):
    """
    This method sort data from a dataframe by an variable and plot the distribution
    of every data entry together in one plot for every bin.
    :param df: dataframe with a column for distribution and the variable
    :param variable: column of values to sort by
    :param distribution: column of y-values to plot
    :param bins: count of bins the data get sorted in
    :param mode: how the data get displayed; base: absolute value, normalize: data normalized, percent: relative frequency
    :param lower_limit: lowest value of the variable that get visualized
    :param upper_limit: highest value of the variable that get visualized
    :output_path: if defined, path to save images as .png
    """
    # define missing variables
    if lower_limit is None:
        lower_limit = min(df[column_name])
    if upper_limit is None:
        upper_limit = max(df[column_name])
    if upper_limit == float("inf"):
        upper_limit = max(df.loc[df[column_name] != float("inf"), column_name])
    if lower_limit == upper_limit:
        print("lower limit can't be the same value as upper limit.")
        return
    bin_scale = (upper_limit - lower_limit) / bins
    x_bins = len(df[distribution][0])
    x_scale = calc.calculate_xscale(df, bins=x_bins)
    count = 0

    # two loops, first for every bin, second for every datapoint in the defined bin space
    for i in np.arange(lower_limit, upper_limit, bin_scale):
        if i == upper_limit - bin_scale:
            last_bin = True
        else:
            last_bin = False
        for d in df.loc[(df[column_name] >= i) & (df[column_name] < i + bin_scale), distribution]:

            #  adjust y-values for the specific mode
            ydata = np.asarray(d)
            if mode == 'normalize':
                ydata_min = ydata.min()
                ydata_max = ydata.max()
                ydata = (ydata - ydata_min) / (ydata_max - ydata_min)
            elif mode == 'percent':
                ydata = ydata / sum(ydata)

            # save the subplot
            plt.plot(x_scale, ydata)
            count += 1
            
        if last_bin:
            for d in df.loc[(df[column_name] == i + bin_scale), distribution]:
                ydata = np.asarray(d)
                if mode == 'normalize':
                    ydata_min = ydata.min()
                    ydata_max = ydata.max()
                    ydata = (ydata - ydata_min) / (ydata_max - ydata_min)

                elif mode == 'percent':
                    ydata = ydata / sum(ydata)
                plt.plot(x_scale, ydata)
                count += 1

        # if there are no subplots, continue with new bin without plotting an empy plot
        if count == 0:
            continue

        # plot all subplots
        plt.title(f'Fragment length distribution of {count} cells\n({column_name}: {i}-{i + bin_scale})')
        plt.xlabel('fragment_length')
        if mode == 'normalize':
            plt.ylabel('normalized count per cell')
        elif mode == 'percent':
            plt.ylabel('Frequency (%)')
        else:
            plt.ylabel('count per cell')

        # output settings
        if output_path != None:
            try:
                title = f'{mode}_distribution_{count}_cells_{column_name}_{i}_{i + bin_scale}'
                plt.savefig(f'{output_path}{title}.png', bbox_inches='tight')
            except:
                print('Output_path not found!')
        plt.show()
        count = 0
    return


def split_and_multiplot_profile(df, tss_positions, half_window_width, split_point_1 = None, split_point_2 = None):
    """
    This function splits fragment into categories based on the
    fragment lengths and then plots their profile around TSS.
    
    :param df: dataframe containing chrom, start and end as columns.
    :param tss_positions: a set of TSS positions.
    :split_point_1: an integer value or an array of values, 
        used to determine the splitting point (length) of small and long fragments.
    :split_point_2: an integer value or an array of values, 
        used to determine the splitting point (length) of long and very long fragments.
    """
    
    if split_point_1 == None:
        if split_point_2 == None:
            genomic_array = utils.get_genomic_array_from_bed_dataframe(df)
            profile  = utils.get_profile_array_for_plots([genomic_array], tss_positions, half_window_width)[0]
            label = "All fragment Lengths"
            profplt([profile], [label], half_window_width)
            return
        else: 
            print('FAILED: Split point 2 cannot be passed without a value for split point 1.')
            return
        
    elif type(split_point_1) == int:
        if split_point_2 == None:
            array_of_genomic_arrays = utils.split_df_and_get_genomic_arrays(df, [split_point_1])
            [profile_1, profile_2] = utils.get_profile_array_for_plots(array_of_genomic_arrays, tss_positions, half_window_width)
            label_1 = "Short fragments : <= "+ str(split_point_1) +"bp"
            label_2 = "Long fragments : > "+ str(split_point_1) +"bp"
            profplt([profile_1, profile_2], [label_1, label_2], half_window_width)
            return
            
        elif type(split_point_2) == int:
            array_of_genomic_arrays = utils.split_df_and_get_genomic_arrays(df, [split_point_1, split_point_2])
            [profile_1, profile_2, profile_3] = utils.get_profile_array_for_plots(array_of_genomic_arrays, tss_positions, half_window_width)
            label_1 = "Short fragments : <= "+ str(split_point_1) +"bp"
            label_2 = "Long fragments : "+ str(split_point_1) +" - "+ str(split_point_2) +"bp"
            label_3 = "Very long fragments : > "+ str(split_point_2) +"bp"
            profplt([profile_1, profile_2, profile_3], [label_1, label_2, label_3], half_window_width)
            return
        
        elif type(split_point_2) == list:
            
            genomic_array_1 = utils.get_genomic_array_from_bed_dataframe(df.loc[df['End']-df['Start'] <= split_point_1])
            array_of_genomic_arrays = utils.split_df_and_get_genomic_arrays(df.loc[df['End']-df['Start'] > split_point_1], split_point_2)
            
            profile_1  = utils.get_profile_array_for_plots([genomic_array_1], tss_positions, half_window_width)[0]
            label_1 = "Short fragments : <= "+ str(split_point_1) +"bp"
            
            profile_array = utils.get_profile_array_for_plots(array_of_genomic_arrays, tss_positions, half_window_width)
            for index, split_point in enumerate(split_point_2):
                label_2 = "Long fragments : "+ str(split_point_1) +" - "+ str(split_point) +"bp"
                label_3 = "Very long fragments : > "+ str(split_point) +"bp"
                profile_2 = np.zeros(2*half_window_width, dtype='i')
                profile_3 = np.zeros(2*half_window_width, dtype='i')
                for profile in profile_array[0:index+1]:
                    profile_2 += profile
                for profile in profile_array[index+1:]:
                    profile_3 += profile
                profplt([profile_1, profile_2, profile_3], [label_1, label_2, label_3], half_window_width)
            return
    
    elif type(split_point_1) == list:
        if split_point_2 == None:
            array_of_genomic_arrays = utils.split_df_and_get_genomic_arrays(df, split_point_1)
            profile_array = utils.get_profile_array_for_plots(array_of_genomic_arrays, tss_positions, half_window_width)
            for index, split_point in enumerate(split_point_1):
                label_1 = "Short fragments : <= "+ str(split_point) +"bp"
                label_2 = "Long fragments : > "+ str(split_point) +"bp"
                profile_1 = np.zeros(2*half_window_width, dtype='i')
                profile_2 = np.zeros(2*half_window_width, dtype='i')
                for profile in profile_array[0:index+1]:
                    profile_1 += profile
                for profile in profile_array[index+1:]:
                    profile_2 += profile
                profplt([profile_1, profile_2], [label_1, label_2], half_window_width)
            return
        
        else: 
            print('FAILED: Split point 2 cannot be passed when an array is passed for split point 1.')
            return
    return


def histplt(df, column_name, output=None, bins=50):
    # plots a histogram for a specific variable in a dataframe
    # Uses Seaborn library to plot.
    # the output parameter is the path for saving the plot.
    
    sns.histplot(x=df[column_name], bins=bins)
    
    plt.xlabel(column_name)
    plt.ylabel('Count')
    
    if output != None:
        plt.savefig(f'{output}{column_name}_hist.png', bbox_inches='tight')
        
    plt.show()
    return


def compplt(df, column_name_1, column_name_2, output=None):
    # plots a comparison between two variables in a dataframe
    # Uses Seaborn library to plot.
    # the output parameter is the path for saving the plot.
    
    x = [x for x in df[column_name_1]]
    y = [x for x in df[column_name_2]]
    
    sns.regplot(x=x, y=y, line_kws={"color": "red"})
         
    plt.xlabel(column_name_1)
    plt.ylabel(column_name_2)
         
    if output != None:
        plt.savefig(f'{output}{column_name_1}_{column_name_2}.png', bbox_inches='tight')
         
    plt.show()
    return


def vioplt(df, column_name, output=None):
    # Plots violin plots for specific columns in the DataFrame.
    # Uses Seaborn library to plot.
    # the output parameter is the path for saving the plot.

    sns.violinplot(df.loc[(df[column_name]!=float(np.inf)),column_name])

    plt.ylabel(column_name)
         
    if output != None:
        plt.savefig(f'{output}{column_name}_vio.png', bbox_inches='tight')
         
    plt.show()
    return


def profplt(profile_array, label_array, half_window_width, name = None, output = None):
    """
    Plots multiple profiles from the array over a single plot.
    
    :param profile_array: array of profiles to plot.
    :param label_array: array of labels corresonding to profiles.
    :param half_window_width: half-size of the window used to create profile.
    :param name: name of the plot.
    :param output: path for the output file. 
    """
    
    fig, ax = plt.subplots()             
    x = np.arange(-half_window_width, half_window_width) 
    
    #Plots for individual profiles
    for index, profile in enumerate(profile_array):
        ax.plot(x, profile, label = label_array[index])
        
    #Set Labels
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc = 'upper right')
    plt.xlabel('Distance from TSS')
    plt.ylabel('Fragment Count')  
    
    #Set Grid
    plt.grid(which='major', axis = 'x',color='#DDDDDD', linewidth=0.8)
    plt.grid(which='minor', axis = 'x', color='#EEEEEE', linestyle=':', linewidth=0.5)
    plt.minorticks_on()
    
    #Set Name
    if name != None:
        ax.set_title(name)
    
    #Save Plot
    if output != None:
        plt.savefig(f'{output}_TSS_profile.png', bbox_inches='tight')
        
    plt.show()
    return
