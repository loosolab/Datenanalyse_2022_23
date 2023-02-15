import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from fragment import utils
from fragment import calc
from fit import calc
    
def bin_distribution(df, data = 'Fragments', variable = 'Mean', bins = 20, mode = 'equal', plot_bins = 50, show = True, output_path = None):
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
    sorted_df = df.sort_values(by=[variable])
    min_variable = sorted_df[variable][0]
    max_variable = sorted_df[variable][-1]
    
    # equal mode
    if mode == 'equal':
        count = 1
        df_length = len(df.index)
        bin_size = int(df_length / bins)
        bin_start = [min_variable]
        bin_end = []
        
        # iterate over an sorted fragment_dictionary and fill the previously defined list with the 
        # fragment lengths. Every bin contains the same number of cells.
        for i,c in enumerate(sorted_df[data]):
            
            # fill fragments in bins and count cells finished
            if count != bin_size + 1:  
                grouped_data[bin_index] += c
                count += 1
            
            # count == bin_size => change bin and reset count
            else:
                count = 1
                bin_index += 1
                grouped_data[bin_index] += c
                bin_start.append(sorted_df[variable][i]) 
                bin_end.append(sorted_df[variable][i])
                
        bin_end.append(max_variable)

    # linear mode
    elif mode == 'linear':
        variable_range = max_variable - min_variable
        bin_size = variable_range / bins
        bin_start = [min_variable]
        bin_end = []
        
        # iterate over an sorted fragment_dictionary and fill the previously defined list with the 
        # fragment lengths. Every bin contains the same variable range.
        for i,c in enumerate(sorted_df[data]):
            
            # fill fragments in bins till the variable limit is reached
            if sorted_df[variable][i] <= min_variable + (bin_index+1)*(bin_size):
                grouped_data[bin_index] += c
            
            # variable is over the limit => change bin
            else:
                bin_index += 1
                grouped_data[bin_index] += c 
                bin_start.append(sorted_df[variable][i])
                bin_end.append(sorted_df[variable][i])
        bin_end.append(max_variable)
    
    # error-message, if mode not correctly defined
    else:
        print('Please only use "equal" or "linear" as variables for "mode"!')
        return
    
    # visualise data for every bin
    for i,s in enumerate(grouped_data):
        sns.histplot(s, bins = plot_bins, kde=True)
        plt.title(f'Distribution of {data}\n({variable}: {bin_start[i]}-{bin_end[i]}')
        plt.xlabel(data)
        plt.ylabel('Count')
        
        # check output settings
        if output_path != None:
            try:
                title = f'distribution_{data}_{variable}_{bin_start[i]}_{bin_end[i]}'
                plt.savefig(f'{output_path}{title}.png', bbox_inches='tight')  
            except:
                print('Output_path not found!')
        if show:
            plt.show()
        else:
            plt.clf()
    return


        
def multi_multiplot_distribution(df, variable = 'Fragment-Count', distribution = 'Distribution', bins = 1, mode = 'base',  lower_limit = None,  upper_limit= None, output_path = None):
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
        lower_limit = min(df[variable])
    if upper_limit is None:
        upper_limit = max(df[variable])
    if upper_limit == float("inf"):
        upper_limit = max(df.loc[df[variable] != float("inf"), variable])
    if lower_limit == upper_limit:
        print("lower limit can't be the same value as upper limit.")
        return
    bin_scale = (upper_limit - lower_limit) / bins
    x_bins =len(df[distribution][0])
    x_scale = calc.get_x_scale(df, x_bins = x_bins)
    count = 0
    
    # two loops, first for every bin, second for every datapoint in the defined bin space
    for i in np.arange(lower_limit, upper_limit, bin_scale):
        for d in df.loc[(df[variable] >= i) & (df[variable] < i+bin_scale), distribution]:
            
            #  adjust y-values for the specific mode
            ydata = np.asarray(d)            
            if mode == 'normalize':
                ydata_min = ydata.min()
                ydata_max = ydata.max()
                ydata = (ydata - ydata_min) / (ydata_max - ydata_min)
            elif mode == 'percent':
                ydata = ydata/sum(ydata)
         
            # save the subplot   
            plt.plot(x_scale, ydata)
            count += 1
            
        # if there are no subplots, continue with new bin without plotting an empy plot
        if count == 0:
            continue
            
        # plot all subplots
        plt.title(f'Fragment length distribution of {count} cells\n({variable}: {i}-{i+bin_scale})')
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
                title = f'distribution_{count}_cells_{variable}_{i}_{i+bin_scale}'
                plt.savefig(f'{output_path}{title}.png', bbox_inches='tight')  
            except:
                print('Output_path not found!')
        plt.show()
        count = 0
    return
       
    
def plot_hist(df, variable,output = None, bins = 50):
    # plots a histogram for a specific variable in a dataframe
    sns.histplot(x = df[variable], bins = bins)
    plt.xlabel(variable)
    plt.ylabel('Count')
    if output != None:
        plt.savefig(f'{output}{variable}.png', bbox_inches='tight')
    plt.show()
    return
    
def plt_comp(df,x_data,y_data,output = None):
    # plots a comparison between two variables in a dataframe
    
    y = [x for x in df[y_data]]
    x = [x for x in df[x_data]]
    
    sns.regplot(x = x, y = y,line_kws={"color": "red"} )
    plt.xlabel(x_data)
    plt.ylabel(y_data)
    if output != None:
        plt.savefig(f'{output}{x_data}_{y_data}.png', bbox_inches='tight')
    plt.show()
    return

def plot_violin(df,column_name,output = None):
    # Plots violin plots for specific columns in the DataFrame.
    # Uses Seaborn library to plot.
    # the output parameter is the path for saving the plot.
    
    sns.violinplot(dataframe[column_name])

    plt.ylabel(column_name)
    if output != None:
        plt.savefig(f'{output}{column_name}.png', bbox_inches='tight')
    plt.show()
    return
    
def analyse_data(df, list_of_column, output = None):
    # plots a histogram for each and a comparisson between each variables, defined by a list, in a dataframe
    for i in list_of_column:
        plot_hist(df,output,i)
    for i in list_of_column:
        for z in list_of_column:
            plt_comp(df,output,i,z)
    return
