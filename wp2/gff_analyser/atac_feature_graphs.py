import os
import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt


def filter_by_cell_count(adata, threshold, key = "SCSA_pred_celltype"):
    """Filter the adata object by cells that have a higher count then threshold in the column described by key"""
    key_values = set(adata.obs[key])
    
    counts = {}
    for key_value in key_values:
        counts[key_value] = adata.obs[adata.obs[key] == key_value].shape[0]
    
    for k, count in counts.items():
        if count < threshold:
            adata = adata.obs[adata.obs[key] != k]

    return adata

def violin_plots(adata, output, group=None, filtered=None, multi_panel=True, log=False):
    
    if filtered == None:
        item_list = []

        for item in adata.obs:
            if adata.obs[item].dtype == 'int64':
                item_list.append(item)
            elif adata.obs[item].dtype == 'float64' and not item.startswith("n_"):
                item_list.append(item)  
    else:
        item_list = filtered
    
    create_output_directory(f"{output}violins/")
            
    for item in item_list:
        sc.pl.violin(adata, keys=item, groupby=group, multi_panel=multi_panel, rotation=90, log=log)
        
        try:
            tmp_label = group.split(' ')
            label = ''.join(tmp_label)
        except:
            label = group
            
        plt.savefig(f"{output}violins/{item}_{label}_violin.png")
        
def scatter_plots(adata, output, feature_filtered=None, celltype_filtered=None, figsize=[60,3]):
    
    if celltype_filtered != None:
        celltype_of_interest = celltype_filtered
        adata_tmp = adata[adata.obs['cell type'].isin(celltype_of_interest)]
    else:
        adata_tmp = adata
    
    if feature_filtered == None:
        item_list = []
        for item in adata_tmp.obs:
            if adata_tmp.obs[item].dtype == 'float64' and not item.startswith("n_"):
                item_list.append(item)
                
    else:
        item_list = feature_filtered
        
    create_output_directory(f"{output}images/")

    for itemx in item_list:

        fig, ax= plt.subplots(1, len(item_list), figsize=figsize)
        x = adata_tmp.obs[itemx]
        for n, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata_tmp.obs[itemy]
                ax = plt.subplot(1, len(item_list), n + 1)
                ax.scatter(x, y)
                plt.title(celltype_of_interest)
                plt.xlabel(itemx, fontsize=10)
                plt.ylabel(itemy, fontsize=15)
                fig.tight_layout()
        plt.savefig(f"{output}images/{itemx}.{itemy}.png")
        plt.show()
        plt.close()
        
def simple_scatter(adata, output, feature=['n_fragments_in_CDS', 'n_total_fragments'], display_all=False):
    item_list = []
    for item in adata.obs:
        if adata.obs[item].dtype == 'int64':
                item_list.append(item)
        elif adata.obs[item].dtype == 'float64' and not item.startswith("n_"):
            item_list.append(item)
    if display_all:
        feature = item_list
        
    create_output_directory(f"{output}scatter")

    for itemx in feature:

        x = adata.obs[itemx]
        for _, itemy in enumerate(item_list):
            if not itemy == itemx:
                y = adata.obs[itemy]
                plt.scatter(x, y)
                plt.xlabel(itemx, fontsize=10)
                plt.ylabel(itemy, fontsize=10)
                plt.show()
                plt.savefig(f"{output}scatter/{itemx}_{itemy}.png")
                plt.close()
                
                
def multi_plot(adata, feature1, feature2, group, out, celltype_filtered=None, multi_panel=False):
    output= f"{out}{feature1}_{feature2}/"
    try:
        if not out.endswith('/'):
            print("Please specify valid Path ('/')")
            return
        create_output_directory(out)
    except:
        print('First out exists')
        
    create_output_directory(output)
    
    if celltype_filtered != None:
        celltype_of_interest = celltype_filtered
        adata_tmp = adata[adata.obs['cell type'].isin(celltype_of_interest)]
    else:
        adata_tmp = adata
        
    sc.settings.figdir = "./"
    sc.pl.violin(adata_tmp, keys=[feature1, feature2], groupby=group, multi_panel=multi_panel, rotation=90)
    plt.savefig(f"{output}violin.png")
    
    for cell_type_raw in adata_tmp.obs['cell type'].unique():
        
        tmp = cell_type_raw.split(' ')
        cell_type = ''.join(tmp)
        
        # for getting rid of spaces, this object only includes filtered cell types since its filtered before
        adata_tmp = adata_tmp[adata_tmp.obs['cell type'].isin([cell_type_raw])]
        
        
        plt.scatter(adata_tmp.obs[feature1], adata_tmp.obs[feature2])
        plt.title(cell_type_raw)
        plt.xlabel(feature1)
        plt.ylabel(feature2)
        plt.savefig(f"{output}{feature1}_{feature2}_{cell_type}.png")
        plt.show()
        

def get_list(element) -> list:
    """ Return a list with element (list or any)"""
    list = []
    if type(element) is list:
        list.extend(element)
    else:
        list.append(element)
    return list

def compare_dimensionreductions(adata, key: str or list, comparator: str or list, out=None, legend_loc="right margin"):
    """Create a Dimension reduction for each key and comparator, save to out"""
    # TODO: Test this
    keys = get_list(key)
    keys.append(get_list(comparator))
    sc.pl.umap(adata, color=keys, legend_loc=legend_loc)
    if out:
        create_output_directory(out)
        if type(key) is list:
            filename = keys.join('_')
        else:
            filename = key
        plt.savefig(f'{out}umaps/{filename}.png')
    plt.close()

def compare_feature_to_celltypes(adata, feature: list or str, comparator: str, out=None):
    """
    Wrapper function for `render_compare_feature_to_celltypes. Draws a violin plot for each element in feature and one violin plot per feature that is grouped by the comparator

    Arguments:
    adata - Object that should be plotted
    feature - can be either a list of keys or a key. Is being compared to the comparator
    comparator - the element to groupby making a comparison easier
    out - location to save the plots to, can be None
    """
    if type(feature) is list:
        features = feature
        for feature in features:
            print(f'{feature}:')
            render_compare_feature_to_celltypes(adata, feature, comparator, out=out)
    else:
        render_compare_feature_to_celltypes(adata, feature, comparator, out=out)

def render_compare_feature_to_celltypes(adata, feature, comparator = 'cell type', out = None):
    """Make a violinplot of feature and feature grouped by comparator and optionally save to out"""
    label = 'Percent' if feature.startswith('pct') else 'Count'
    groups = sort_groups(adata, group_celltypes(adata, 5, comparator))
    fig, axs = plt.subplots(nrows=1, ncols=len(groups)+1, gridspec_kw={'wspace':0.4})
    sc.pl.violin(adata, feature, ax = axs[0], show=False, ylabel=label)
    for index, group in enumerate(groups):
        sc.pl.violin(adata, feature, groupby=group, ax = axs[index+1], rotation=90, show=False, ylabel="")
    if out:
        create_output_directory(out)
        fig.savefig(f"{out}violins/{feature}_{comparator}_violin.png")

def sort_groups(adata, groups, feature, key = 'cell type'):
    adata_list = []
    for group in groups:
        adata_tmp = adata[adata.obs[key] == group]
        max = adata_tmp.obs[feature].max()
        min = adata_tmp.obs[feature].min()
        adata_list.append({'data': adata_tmp, 'max': max, 'min': min})
    adata_list.sort(lambda g: g['max'])



def group_celltypes(adata, max_size = 5, key = 'cell type'):
    groups = {}
    index = 0
    for item in adata.obs[key].unique():
        if index not in groups:
            groups[index] = []
        if len(groups[index]) <= max_size:
            groups[index].append(item)
        else:
            index += 1
    return groups

def create_output_directory(out):
    """Creates the declared directory."""
    try:
        os.mkdir(f"{out}")
    except:
        print(f"Folder {out} exists.")
