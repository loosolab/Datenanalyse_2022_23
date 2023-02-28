import os
import numpy as np
import scanpy as sc
import pandas as pd
from matplotlib import pyplot as plt


def filter_by_cell_count(adata, threshold, key = "cell type"):
    """Filter the adata object by cells that have a higher count then threshold in the column described by key"""
    key_values = set(adata.obs[key])
    
    counts = {}
    for key_value in key_values:
        counts[key_value] = adata.obs[adata.obs[key] == key_value].shape[0]
    
    for k, count in counts.items():
        if count < threshold:
            adata = adata[adata.obs[key] != k]

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
    """ Return a list with element """
    element_list = []
    if type(element) is not list:
        element_list.append(element)
    else:
        element_list = element
    return element_list

def compare_dimensionreductions(adata, key: str or list, comparator: str or list, out=None, legend_loc="right margin"):
    """Create a Dimension reduction for each key and comparator, save to out"""
    keys = [*get_list(key), *get_list(comparator)]
    sc.pl.umap(adata, color=keys, legend_loc=legend_loc)
    if out:
        create_output_directory(out)
        if type(key) is list:
            filename = keys.join('_')
        else:
            filename = key
        plt.savefig(f'{out}umaps/{filename}.png')
    plt.close()

def compare_feature_to_celltype(adata, feature, celltype, name, sharey=True, out=None):
    if type(feature) is list:
        features = feature
        for feature in features:
            print(f'{feature}:')
            render_feature_to_celltypes(adata, feature, celltype, name, out, sharey)
    else:
        render_feature_to_celltype(adata, feature, celltype, name, out, sharey)

def render_feature_to_celltype(adata, feature, celltype, name, out=None, sharey=True):
    adata_tmp = adata.copy()
    adata_tmp = adata_tmp[adata_tmp.obs['cell type'] == celltype]
    fig, axs = plt.subplots(nrows=1, ncols=2, gridspec_kw={'wspace':0.4, 'hspace':0.5}, figsize=(10,5), sharey=sharey)
    sc.pl.violin(adata, feature, ax = axs[0], ylabel=name, show=False)
    sc.pl.violin(adata_tmp, feature, ax = axs[1], ylabel=name, show=False)
    if out:
        create_output_directory(out)
        create_output_directory(out+'violins/')
        fig.savefig(f"{out}violins/{feature}_{celltype}_violin.png")

def compare_feature_to_celltypes(adata, feature: list or str, comparator: str, max_size, groups: list = None, name=None, out=None, sharey=True):
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
            render_compare_feature_to_celltypes(adata=adata, feature=feature, groups=groups, comparator=comparator, name=name, max_size=max_size, out=out, sharey=sharey)
    else:
        render_compare_feature_to_celltypes(adata=adata, feature=feature, groups=groups, comparator=comparator, name=name, max_size=max_size, out=out, sharey=sharey)

def render_compare_feature_to_celltypes(adata, feature, groups, comparator = 'cell type', name=None, max_size=5, out=None, sharey=True):
    """Make a violinplot of feature and feature grouped by comparator and optionally save to out"""
    label = 'Percent' if feature.startswith('pct') else 'Count'
    adata = filter_by_cell_count(adata, 5, comparator)

    groups = get_sorted_groups_with_size(adata=adata, feature=feature, max_size=max_size, key=comparator, groups=groups)
    fig, axs = plt.subplots(nrows=1, ncols=len(groups)+1, gridspec_kw={'wspace':0.4, 'hspace':0.5}, figsize=(40,20), sharey=sharey)

    sc.pl.violin(adata, feature, ax = axs[0], show=False, ylabel=label)
    for index, group in groups.items():
        adata_tmp = adata.copy()
        adata_tmp = adata_tmp[adata_tmp.obs[comparator].isin(group)]
        sc.pl.violin(adata_tmp, feature, groupby=comparator, ax = axs[index+1], rotation=90, show=False, ylabel=name)
    if out:
        create_output_directory(out+'violins/')
        fig.savefig(f"{out}violins/{feature}_{comparator}_violin.png")

def get_sorted_groups_with_size(adata, feature: str, max_size: int, key: str = 'cell type', groups=None):
    sorting = {}
    if 'cell groups' not in adata.uns:
        adata.uns['cell groups'] = {}
    if not groups:
        groups = adata.obs[key].unique()
    for group in groups:
        adata_tmp = adata.copy()
        adata_tmp = adata[adata.obs[key] == group]

        sorting[group] =  {
            'min':  adata_tmp.obs[feature].min(),
            'max':  adata_tmp.obs[feature].max(),
            'mean': adata_tmp.obs[feature].mean()
        }
            
    sorted_groups = sorted(sorting.items(), reverse=True, key=lambda v: v[1]['max'] - v[1]['min'])
    
    bins = {}
    index = 0
    bins[index] = []
    for group in sorted_groups:
        if len(bins[index]) == max_size:
            index += 1
            bins[index] = []
        bins[index].append(group[0])
    return bins
                         
def create_output_directory(out):
    """Creates the declared directory."""
    try:
        os.mkdir(f"{out}")
    except:
        print(f"Folder {out} exists.")
